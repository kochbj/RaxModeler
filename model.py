from math import log, exp
from dendropy import Tree

class model:

    def __init__(self,model_name,param_names,rax_name,args,s_size):
        self.model_name=model_name
        self.param_names=param_names
        self.rax_name=rax_name
        self.name=model_name+param_names
        self.path=args.opath+'/'+model_name
        
        self.time, self.lnL=self.__get_lnL()
        
        self.n_params =self.__calc_n_params(branches=s_size.branches)
        self.tree=None
        self.cong_list=[]
        if args.dtc==True or args.cong==True:
            self.tree=self.__get_tree(taxon_set=s_size.taxon_set)
 

        if args.cong==True: self.cong_list=[] 

        self.criterion=self.__create_crit_dict(clist=args.clist)
        self.delta=self.__create_crit_dict(clist=args.clist)
        self.preweights=self.__create_crit_dict(clist=args.clist)
        self.weights=self.__create_crit_dict(clist=args.clist)
        self.rank=self.__create_crit_dict(clist=args.clist)

        
        if args.dtc==True:
            self.criterion['DTC']={}
            self.delta['DTC']={}
            self.preweights['DTC']={}
            self.weights['DTC']={}
            self.rank['DTC']={}
            
            for c in args.clist:
                self.criterion['DTC'][c]=0.0 
                self.delta['DTC'][c]=0.0 
                self.preweights['DTC'][c]=0.0
                self.weights['DTC'][c]=0.0
                self.rank['DTC'][c]=0.0
        
            
        self.criterion['AIC']['raw']=self.__calc_aic()
        for c in s_size.correction_dict:
            print c
            correction=self.__calc_correction(sample_size=s_size.correction_dict[c])
            self.criterion['AIC'][c]= self.criterion['AIC']['raw']+correction
            self.criterion['BIC'][c]=self.__calc_bic(sample_size=s_size.correction_dict[c])

        
            
    def __get_lnL(self):
        lnL_file_name=self.path+"/RAxML_log."+self.param_names
        lnLfile=open(lnL_file_name,'r')
        time, lnL = lnLfile.next().split()
        lnLfile.close()
        time=float(time)
        lnL=-1.0*float(lnL)
        return time, lnL

    def __get_tree(self,taxon_set):
        tree_file_name=self.path+"/RAxML_result."+self.param_names
        tree=Tree.get_from_path(tree_file_name,'newick',encode_splits=True,taxon_set=taxon_set)
        return tree

    def __calc_n_params(self, branches):
        if self.param_names=='G':
            params=1
        elif self.param_names=='GI':
            params=2
        elif self.param_names=='GF':
            params=20
        elif self.param_names=='GIF':
            params=21
        n_params=params+branches
        return n_params

    def __create_crit_dict(self, clist):
        dict={}
        dict['AIC']={}
        dict['AIC']['raw']=0.0
        dict['BIC']={}
        
        
        for c in clist:
            dict['AIC'][c]= 0.0
            dict['BIC'][c]= 0.0
        
        return dict


    def __calc_aic(self):
        return 2.0*self.lnL+2.0*self.n_params


    def __calc_bic(self,sample_size):
        return 2.0*self.lnL+self.n_params*log(sample_size)

    def __calc_correction(self,sample_size):
        print "sample_size", sample_size, self.n_params
        #A DIVIDE BY ZERO ERROR IS CAUSED BY THE ALN LENGTH CORRECTION, FOR NOW JUST SET CORRECTION REALLY HIGH
        if sample_size-self.n_params-1.0<=0.0: correction=1000000000
        else: correction=(2.0*self.n_params*(self.n_params+1.0))/(sample_size-self.n_params-1.0)
        return correction
