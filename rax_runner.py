from sys import stdout, exit
from subprocess import Popen, PIPE
from model import model
from dendropy import Tree
from random import randint
from threading import Thread



class rax_runner(Thread):

    def __init__(self,queue,model_index,raxml_path,args,s_size):
        Thread.__init__(self)
        self.args=args
        self.s_size=s_size
        self.stree=args.startingtree
        self.queue=queue
        self.raxml_path=raxml_path
        self.output=args.opath
        self.ifile_name=args.input
        self.raxstdout=''
        self.raxstderr=''
        self.model_index=model_index

    def run(self):
        while True:
            model_tuple=self.queue.get()
            model_name=model_tuple[0]
            param_names=model_tuple[1]

            model=self.__raxml_to_model(model_name=model_name,param_names=param_names,args=self.args,s_size=self.s_size)
                
            printable_lnL="{0:.6f}".format(model.lnL)
            printable_time="{0:.6f}".format(model.time)

            stdout.write("{0:<22} ".format(model_name+" +"+param_names+'...'))        
            stdout.write("{0:>15}\n".format("-lnL= "+printable_lnL))
            stdout.write("{0:>40}\n\n".format(printable_time+" secs"))
            stdout.flush()
            
            self.model_index.add_model(model)
            
            self.queue.task_done()
        
    def __model_to_raxmodelstr(self,model_name,param_names):
        if param_names=='G':
            modelstr="PROTGAMMA"+model_name
        elif param_names=='GI':
            modelstr='PROTGAMMAI'+model_name
        elif param_names=='GF':
            modelstr="PROTGAMMA"+model_name+'F'
        elif param_names=='GIF':
            modelstr='PROTGAMMAI'+model_name+'F'
        return modelstr

    def __raxml_to_model(self, model_name, param_names,args,s_size):
        raxmlmodelstr=self.__model_to_raxmodelstr(model_name,param_names)
        ##think about changing this so it catches crashes and continues
        rax_cmd=Popen('{raxpath} -f e -m {modelstr} -s {phy_file} -t {starting_tree} -n {output_names} -w {out_path}/{model_name}'.format(
                raxpath=self.raxml_path,
                modelstr=raxmlmodelstr,
                output_names=param_names,
                out_path=self.output,
                model_name=model_name,
                phy_file=self.ifile_name,
                starting_tree=self.stree),
                      
            shell=True,
            stdout=PIPE,
            stderr=PIPE)
				                                      
        rax_out = rax_cmd.communicate()
        full_name=model_name+param_names
        self.raxstdout+='='*80+'\n'+full_name.center(80)+'\n'+'='*80+'\n'
        self.raxstdout+=rax_out[0]

        self.raxstderr+='='*80+'\n'+full_name.center(80)+'\n'+'='*80+'\n'
        self.raxstderr+=rax_out[0]
        
        new_model=model(model_name=model_name,
                        param_names=param_names,
                        rax_name=raxmlmodelstr,
                        args=self.args,
                        s_size=s_size
                        )
        
        return new_model


    
        
                   
        

