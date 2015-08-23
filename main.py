#!/home/kochbj/bin/bin/python2.7

##import subprocess, operator, os
##from math import log, exp

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from sys import stderr, exit, stdin, stdout
from os import chdir, getcwd, mkdir, listdir
from shutil import rmtree
from time import localtime
from Queue import Queue

from model_index import model_index
from rax_runner import rax_runner
from sample_size import sample_size_obj
from tree_maker import tree_maker


from os import path
my_path=path.dirname( path.realpath( __file__ ) )

RAXMLPATH=my_path+'/raxmlHPC-HYBRID'

def parse_args():
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('input',default=stdin, type=str, help="Input file.")
    parser.add_argument('--alnformat','-a',default='phylip',choices=['nexus','phylip','fasta','clustal'], type=str, help="Input Format")
    
    parser.add_argument('--threads', '-t', type=int, default=2, help="Number of threads for RAxML.")
    parser.add_argument('--models','-m',nargs='*',
                        default=['BLOSUM62', 'CPREV', 'DAYHOFF', 'DCMUT', 'FLU', 'HIVB', 'JTT', 'JTTDCMUT', 'LG', 'MTART', 'MTMAM', 'MTREV', 'MTZOA', 'PMB', 'RTREV', 'VT', 'WAG'],
                        help="Models to use.")
    ##ADD TEST
    parser.add_argument('--framework','-f',type=str,default='AIC',choices=['AIC','BIC','DTC'],help="Sort by criterion: AIC (Akaike), BIC (Bayesian Inference), DTC (Decision Theory)")
    parser.add_argument('--correction','-c',type=str,default='raw',choices=['length','shannon','c3'],
                        help="Sort by correction: length (alignment length), shannon (summed shannon entropies), c3 (number of taxa * alignment length * normalized shannon entropy X)")
    parser.add_argument('--dtc', action='store_true',default=False,help="Include the decision theory criterion.  Note this could increase run time for large alignments")
    parser.add_argument('--cong', action='store_true',default=False,help="Calculate symmetric distances for the best 5 models.  Note this could increase run time for large alignments")
    parser.add_argument('--opath','-o', default='.',help="Output path for RAxML_files folder")
    parser.add_argument('--startingtree','-s', default=None, help="Starting tree name")
    parser.add_argument('--treeformat','-r',default='newick',choices=['newick','nexus'],type=str,help="If an input tree is specified, what format is it?")
    parser.add_argument('--seed','-b', default=None,type=str, help="int seed for starting tree")

    args=parser.parse_args()
    if args.framework in ['BIC','DTC'] and args.correction=='raw': args.correction='length'
    return args

    

def print_banner():
    stdout.write('''       _     _        __  __        _     _         
  _ _ /_\   /_\ __ __|  \/  |___ __| |___| |___ _ _ 
 | '_/ _ \ / _ \\\ \ /| |\/| / _ | _` / -_) / -_) '_|
 |_|/_/ \_|_/ \_|_\_\|_|  |_\___|__,_\___|_\___|_| v0.9\n\n''')
                                               
    stdout.flush()
    stdout.write("This is free software, so there is no guarantee that it works properly,\nhas any error checks, or does anything other than print the cool banner above.\n")
    stdout.write("Stolen from Alexandros Stamatakis (RAxML) and Federico Abascal, Rafael Zardoya,\nand David Posada (ProtTest) by Bernie.\n\n")
    stdout.flush()    
    

def check_already_ran():  
    if "RAxML_files" in listdir('.'):
	stderr.write("ERROR:  rAAxModeler found a 'RAxML_files' folder in cwd.\n\trAAxModeler needs to create this directory to run.\n")
    	stderr.write("\tIs it okay to delete it? (y or n): ")
	stderr.flush()
    	decision = stdin.read(1)
        if decision=='y':
            rmtree("RAxML_files")
            return
        else:
            stderr.write("Exiting...\n")
            stderr.flush()
            exit(0)


def print_rax_cmd(ifile,opath,threads,best_model):
        stdout.write('\n'+'='*127+'\n')
        stdout.write("Suggested RAxML command:\nraxmlHPC-PTHREADS -T {threads} -f d -N 100 -s {ifile} -w {opath} -m {modelstr} -n [OUT NAME]\n".format(threads=threads,
                                                                                                                                                       opath=opath,
                                                                                                                                                       ifile=ifile,
                                                                                                                                                       modelstr=best_model.rax_name
                                                                                                                                                       )) 
        stdout.write('='*127+'\n')
        stdout.flush()
    

def main():
    print_banner()

    args=parse_args()  ### PASSS ARGS TO EVERYTHING
    chdir(args.opath)
    check_already_ran()
    time=localtime()
    timestr='_'.join([str(time.tm_sec),
                      str(time.tm_min),
                      str(time.tm_hour),
                      str(time.tm_mday),
                      str(time.tm_mon),
                      str(time.tm_year)])
    cwd=getcwd()
    args.opath=cwd+'/RAxML_files_'+timestr
    mkdir(args.opath)
    args.clist=['length','shannon','c3']
    

    
    s_size=sample_size_obj(args)
    s_size.print_stats(args.models)
    
    args.input=s_size.input_to_phy(iformat=args.alnformat,path=args.opath,ifile=args.input)
    
    
    m_index=model_index(args)


    t_maker=tree_maker(raxml_path=RAXMLPATH,args=args)

    mkdir(args.opath+"/STARTING_TREE")

    if args.startingtree!=None:
        args.startingtree=t_maker.nexus_to_newick()
    else:
        args.startingtree=t_maker.make_starting_tree()

    m_queue=Queue()
    model_list=[]

    for matrix in args.models:
        mkdir(args.opath+'/'+matrix)
        for params in ['G','GI','GF','GIF']:
            model_tuple=(matrix,params)
            m_queue.put(model_tuple)
    
    for thread in range(args.threads):
        runner=rax_runner(queue=m_queue,
                            model_index=m_index,
                            raxml_path=RAXMLPATH,
                            args=args,
                            s_size=s_size)
        runner.setDaemon(True)
        runner.start()

    


    m_queue.join()
            
            
    for criterion in ['AIC','BIC']:
        m_index.calc_deltas(criterion)
        m_index.calc_weights(criterion)
        m_index.rank_criteria_by_weight(criterion)

    if args.dtc==True:
        m_index.calc_dtcs()
        m_index.calc_deltas('DTC')
        m_index.calc_weights('DTC')
        m_index.rank_criteria_by_weight('DTC')

    m_index.sort_models(metric=args.framework,crit=args.correction)
    
    if args.cong==True:
        m_index.calc_congruence(num_best_models=4)
        
        
    
    m_index.print_delta_table(metric=args.framework,crit=args.correction)
    m_index.print_big_table(metric=args.framework,crit=args.correction)
    m_index.print_optional_table(metric=args.framework,crit=args.correction, num_best_models=4)


    print_rax_cmd(threads=args.threads,ifile=args.input,opath=cwd,best_model=m_index.ranked_models[0])
        
if __name__ == "__main__":        
	main()
        

    
