from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Alphabet import generic_protein

from math import log

from sys import stderr, stdout, exit

from dendropy import TaxonSet

class sample_size_obj:


    def __init__(self,args):
        self.aln=self.__get_aln(ifile_name=args.input,iformat=args.alnformat)

        self.branches=2*len(self.aln)-3
        self.alength=self.aln.get_alignment_length()
        self.n_taxa=len(self.aln)

        if args.alnformat=='phylip':
            self.__phyheader_equals_aln(args.input)
            
            
        
        self.pssm=self.__get_pssm()
        
        self.sum_shannon=self.__calc_sum_shannon()
        self.norm_shannon=self.__calc_norm_shannon()
        self.c3=self.__calc_c3()
        self.correction_dict={}
        
        for cfactor in args.clist:
            if cfactor=='length': 
                self.correction_dict['length']=self.alength
            elif cfactor=='shannon':
                self.correction_dict['shannon']=self.sum_shannon
            elif cfactor=='c3':
                self.correction_dict['c3']=self.c3

        self.taxon_set=None
        if args.dtc==True or args.cong==True:
            self.taxon_set=self.__create_taxon_set()
        
    def __get_aln(self,ifile_name,iformat):
        try:
            alnfile=open(ifile_name,'r')
        except(IOError):
            stderr.write("ERROR: "+ifile_name+" does not exist.\nExiting...\n")
            stderr.flush()
            exit(0)
        try:
            aln=AlignIO.read(handle=alnfile,format=iformat,alphabet=generic_protein)
        except(ValueError):
            stderr.write("ERROR: "+ifile_name+" is not in "+iformat+" format.\nExiting...\n")
            stderr.flush()
            exit(0)
        alnfile.close()
        return aln

    def input_to_phy(self,iformat,path,ifile):
        if iformat=='phylip':
            return ifile
        else:
            newalnfile=open(path+'/ALN_CONVERTED.phy','w')
            newalnfile.write(self.aln.format('phylip'))
            newalnfile.close()
            return path+'/ALN_CONVERTED.phy'
            
            


    # I am using this check because for some reason Biopython does not check
    # that the alignment length is correct.
    def __phyheader_equals_aln(self,ifile_name):
        alnfile=open(ifile_name,'r')
        test_n_taxa, test_len =alnfile.next().split()
        alnfile.close()
        if int(test_n_taxa)!=self.n_taxa or int(test_len)!=self.alength:
            stderr.write("ERROR: Dimensions in "+ifile_name+" header are incorrect.\nExiting...\n")
            stderr.flush()
            exit(0)
        for record in self.aln:
            for char in '!@#$%^&*()\'\"+=':
                    if record.description.find(char)!=-1:
                         stderr.write("ERROR: Sequence headers in "+ifile_name+" ("+record.description+") contain non-alphanumerical characters other than dashes and underscores.\nExiting...\n")
                         stderr.flush()
                         exit(0)
        



    def __get_pssm(self):
        summary_aln = AlignInfo.SummaryInfo(self.aln)
        consensus_pssm=summary_aln.pos_specific_score_matrix()
        pssm=consensus_pssm.pssm
        return pssm
    

    def __calc_sum_shannon(self):
        protein_alphabet=['A', 'C', 'E', 'D', 'G', 'F', 'I', 'H', 'K', 'M', 'L', 'N', 'Q', 'P', 'S', 'R', 'T', 'W', 'V', 'Y']
        total_entropy=0.0
        for col in range(self.alength):
            col_entropy=0.0
            for aa in protein_alphabet:
                try: freqs=self.pssm[col][1][aa]/sum(self.pssm[col][1].values())
                except KeyError: freqs=0.0  
                if freqs > 0.0:
                    col_entropy+=(freqs*log(freqs,2))
                print "AA",aa, freqs
            print "COL ENTROPY", col_entropy
            total_entropy-=col_entropy
        print "TOTAL_ENTROPY", total_entropy
        return total_entropy

    def __calc_norm_shannon(self):
        norm_shannon=self.sum_shannon/self.alength/4.32
        return norm_shannon

    def __calc_c3(self):
        c3=self.n_taxa*self.alength*self.norm_shannon
        return c3

    def __create_taxon_set(self):
        

        taxonset=TaxonSet()
        for seq in self.aln:
            taxonset.new_taxon(label=seq.id)
        return taxonset

    def print_stats(self,model_list):
        stdout.write('='*127+'\n')

        stdout.write("Analyzing models {0}\n".format((', ').join(model_list)))
        stdout.flush()
        stdout.write('='*127+'\n')
        stdout.write("Number of parameters:\n")
        str_branches=str(self.branches).zfill(4) 
        stdout.write("\tBranches:\t{0:>59}\n".format(str_branches))
        stdout.write("\tGamma:                                              "+' '*13+"1+{branches}= {total}\n".format(branches=str_branches,total=self.branches+1))
        stdout.write("\tGamma + Invariant Sites:                            "+' '*11+"1+1+{branches}= {total}\n".format(branches=str_branches,total=self.branches+2))
        stdout.write("\tGamma + Empirical Base Frequencies:                 "+' '*10+"1+19+{branches}= {total}\n".format(branches=str_branches,total=self.branches+20))
        stdout.write("\tGamma + Invariant Sites+Empirical Base Frequencies: "+' '*8+"1+1+19+{branches}= {total}\n".format(branches=str_branches,total=self.branches+21))
        stdout.flush()
        stdout.write('='*127+'\n')
        stdout.flush()

        stdout.write("Correction factors:\n")
        str_length=str(self.correction_dict['length']).zfill(4)
        stdout.write("\tlength (alignment length): {0:>48}\n".format(str_length))
        str_shannon="{0:6.6f}".format(self.correction_dict['shannon'])
        stdout.write("\tshannon (summed shannon entropies): {0:>39}\n".format(str_shannon))
        str_c3="{0:6.6f}".format(self.correction_dict['c3'])
        stdout.write("\tc3 (normalized_shannon / length * max_entropy): {0:>27}\n".format(str_c3))
        stdout.write('='*127+'\n')
        stdout.flush()
