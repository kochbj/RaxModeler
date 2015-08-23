from sys import stdout,stderr, exit
from subprocess import Popen, PIPE
from Bio import Phylo
from random import randint



class tree_maker:

    def __init__(self,raxml_path,args):
        self.raxml_path=raxml_path
        self.output=args.opath
        self.ifile_name=args.input
        self.tree_name=args.startingtree
        self.tree_format=args.treeformat
        self.seed=args.seed
        self.raxstdout=open(self.output+'/RAxML_stdout','w')
        self.raxstderr=open(self.output+'/RAxML_stderr', 'w')

    def nexus_to_newick(self):
        try:
            tree=Phylo.read(trees=self.tree_name,format=self.tree_format)
        except (ValueError):
            stderr.write("ERROR: Tree file "+self.tree_name+"is not in "+self.tree_format+'\n')
            exit(0)
        if self.tree_format=="newick":
            return self.tree_name  ##OPEN AND CLOSE TO CHECK IF NEWICK
        else:
            self.tree_name=self.output+'/STARTING_TREE/TREE_CONVERTED.tre'
            self.tree_format="newick"
            Phylo.write(trees=self.tree_name,format=self.tree_format)
            return self.tree_name


    def make_starting_tree(self):
        if self.seed==None:
            self.seed=''
            for i in range(10):
                self.seed+=str(randint(0,9))
            
            stdout.write("Created random seed: "+self.seed+'\n')
        else:
            stdout.write("Using given seed: "+self.seed+'\n')
        stdout.flush()
        
        stdout.write("Making starting tree... ")
        stdout.flush()
        STproc=Popen('{path} -p {seed} -y -m PROTCATJTT -s {ifile_name} -n ST.tre -w {output}/STARTING_TREE'.format(
                path=self.raxml_path,
                seed=self.seed,
                ifile_name=self.ifile_name,
                output=self.output),
            shell=True,
            stdout=PIPE,
            stderr=PIPE)
        
        STout= STproc.communicate()
        raxstdout=open(self.output+'/RAxML_stdout','w')
        raxstderr=open(self.output+'/RAxML_stderr', 'w')

        
        raxstdout.write('='*80+'\n'+'STARTING TREE'.center(80)+'\n'+'='*80+'\n')
        for line in STout[0].split('\n'):
            raxstdout.write(line+'\n')

        raxstderr.write('='*80+'\n'+'STARTING TREE'.center(80)+'\n'+'='*80+'\n')
        for line in STout[1].split('\n'):
            raxstderr.write(line+'\n')
        
        stdout.write("done.\n\n")
        stdout.flush()
        raxstdout.close()
        raxstderr.close()
        return self.output+'/STARTING_TREE/RAxML_parsimonyTree.ST.tre'
