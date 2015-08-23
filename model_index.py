from math import exp,log
from sys import stdout
from dendropy.treecalc import false_positives_and_negatives, euclidean_distance

##PUT ALL PUBLICS UP TOP
class model_index:

    def __init__(self,args):
        self.models={}
        self.ranked_models=[]
        self.use_dtc=args.dtc
        self.use_cong=args.cong
        self.clist=args.clist

        
        self.min_criterion=self.__create_crit_dict()
        self.sumweights=self.__create_crit_dict()


        if self.use_dtc==True:
            self.min_criterion['DTC']={}
            self.sumweights['DTC']={}
            for cfactor in self.clist:
                self.min_criterion['DTC'][cfactor]=0.0
                self.sumweights['DTC'][cfactor]=0.0
            self.dtc_sum_dict={}            

    def __create_crit_dict(self):
        dict={}
        dict['AIC']={}
        dict['AIC']['raw']=0.0
        dict['BIC']={}
        
        
        for cfactor in self.clist:
            dict['AIC'][cfactor]= 0.0
            dict['BIC'][cfactor]= 0.0
        
        return dict


    def add_model(self,model):
        self.models[model.name]=model


    ##helper for calc dtcs
    def calc_dtcs(self): 
        for cfactor in self.clist:  
            min_bic=self.min_criterion['BIC'][cfactor]
            for model in self.models:
                dtc_sum=0.0
                model_tree=self.models[model].tree
                for other_model in self.models:
                    if other_model==model: continue
                    
                    other_tree=self.models[other_model].tree

                    distance=euclidean_distance(model_tree,other_tree)
                    if distance <= 0.0: continue
                    other_bic=self.models[other_model].criterion['BIC'][cfactor]
                    power = log(distance) - other_bic + min_bic
                    if power > -30.0:
                        dtc_sum+=exp(power)
                self.models[model].criterion['DTC'][cfactor]=dtc_sum    
   

    #helper for fill deltas
    def __get_minima(self,metric):
        for crit in self.min_criterion[metric]:
            best_model = min(self.models.itervalues(),key=lambda x: x.criterion[metric][crit])
            self.min_criterion[metric][crit]=best_model.criterion[metric][crit]


    def calc_deltas(self,metric):
        self.__get_minima(metric)
        for crit in self.min_criterion[metric]:
            min=self.min_criterion[metric][crit]
            for model in self.models:
                crit_of_model=self.models[model].criterion[metric][crit]
                self.models[model].delta[metric][crit]=crit_of_model-min
       
     ###Helper for calc weights           
    def __fill_pre_weights_get_sum_weights(self,metric):
        for crit in self.sumweights[metric]:
            for model in self.models:

                delta=self.models[model].delta[metric][crit]
                #print model, "DELTA", delta
                preweight=exp(-1.0*delta/2.0)
                #print model, "PREWIGHT",preweight
                self.models[model].preweights[metric][crit]=preweight
                self.sumweights[metric][crit]+=preweight

    def calc_weights(self,metric):
        self.__fill_pre_weights_get_sum_weights(metric)
        for crit in self.sumweights[metric]:
            sumweight=self.sumweights[metric][crit]
            #print "SUM WEIGHT", sumweight
            for model in self.models:
                preweight=self.models[model].preweights[metric][crit]
                self.models[model].weights[metric][crit]=preweight/sumweight


    def rank_criteria_by_weight(self,metric):
        for crit in self.sumweights[metric]:
            self.ranked_models=self.models.values()
            self.ranked_models.sort(key=lambda model: model.weights[metric][crit],reverse=True)
            rank_count=1
            for model in self.ranked_models:
                self.models[model.name].rank[metric][crit]=rank_count
                rank_count+=1
                
    def sort_models(self,metric,crit):
        self.ranked_models=self.models.values()
        self.ranked_models.sort(key=lambda x: x.weights[metric][crit],reverse=True)

    def calc_congruence(self,num_best_models):
        for model in self.ranked_models:
            for best_index in range(num_best_models):
                if self.ranked_models[best_index]==model:
                    model.cong_list.append(' * ')
                else:
                    congruence=false_positives_and_negatives(reference_tree=self.ranked_models[best_index].tree,test_tree=model.tree)
                    model.cong_list.append(congruence)           
                
                

    def print_delta_table(self,metric,crit):
        stdout.write('\nSelected statistics for '+metric+'-'+crit+':\n')
        stdout.flush()

        stdout.write('='*127+'\n')
        stdout.write(" "*10+'\t')
        stdout.write("{0:<20}\t".format('delta '+metric+'-'+crit))
        stdout.write("{0:<20}\t".format(metric+'-'+crit))
        stdout.write("{0:<20}\t".format(metric+'-'+crit+' weight'))
        stdout.write('\n'+'='*127+'\n')
        stdout.flush()
        
        for model in self.ranked_models:
            stdout.write("{0:<10}\t".format(model.model_name+' +'+model.param_names))
            delta="{0:.6f}".format(model.delta[metric][crit])
            criterion="{0:.6f}".format(model.criterion[metric][crit])
            weight="{0:.6f}".format(model.weights[metric][crit])
            stdout.write("{0:<20}\t{1:<20}\t{2:<20}\t\n".format(delta,criterion,weight))
            stdout.flush()


    #This probably needs reworking
    def print_big_table(self,metric,crit):
        stdout.write('\nWeight table sorted by '+metric+'-'+crit+':\n')
        stdout.flush()
        
        stdout.write('='*127+'\n')
        stdout.write(" "*10+'\t')
        for criterion in ['AIC-raw','AIC-length','AIC-shannon','AIC-c3','BIC-length','BIC-shannon','BIC-c3']:
            stdout.write("{0:<15}\t".format(criterion))
        stdout.write('\n'+'='*127+'\n')
        stdout.flush()

        for model in self.ranked_models:
            stdout.write("{0:<8}\t".format(model.model_name+' +'+model.param_names))

            #print AIC-raw first
            weight="{0:.6f}".format(model.weights['AIC']['raw'])
            rank=str(model.rank['AIC']['raw'])
            rank="({0})".format(rank.zfill(2))
            stdout.write("{weight:>8} {rank}\t".format(weight=weight,rank=rank))
            

            for framework in ['AIC','BIC']:
                for criterion in self.clist:
                    weight="{0:.6f}".format(model.weights[framework][criterion])
                    #jumping through hoops because I cant figure out how to get str.format
                    #to both truncate floats and justify

                    rank=str(model.rank[framework][criterion])
                    rank="({0})".format(rank.zfill(2))
                    stdout.write("{weight:<8} {rank}\t".format(weight=weight,rank=rank))

            stdout.write('\n')
            stdout.flush()

        
        
    ##helper for print big table
    def print_optional_table(self,metric,crit,num_best_models):
        if self.use_dtc==False and self.use_cong==False:
            return
        stdout.write('='*127+'\n')
        stdout.write(" "*10+'\t')
        if self.use_dtc==True:
            for cfactor in self.clist:
                criterion='DTC-'+cfactor
                stdout.write("{0:<15}\t".format(criterion))            
        if self.use_cong==True:
            stdout.write('  ')
            for model in self.ranked_models[:num_best_models]:
                stdout.write("{0:<10}\t".format(model.model_name+' +'+model.param_names+' SD'))
        stdout.write('\n'+'='*127+'\n')
        stdout.flush()

        for model in self.ranked_models:
            stdout.write("{0:<10}\t".format(model.model_name+' +'+model.param_names))
            if self.use_dtc==True:
                for criterion in self.clist:

                    weight="{0:.6f}".format(model.weights['DTC'][criterion])
                    #jumping through hoops because I cant figure out how to get str.format
                    #to both truncate floats and justify

                    rank=str(model.rank['DTC'][criterion])
                    rank="({0})".format(rank.zfill(2))
                    stdout.write("{weight:<8} {rank}\t".format(weight=weight,rank=rank))
               
            if self.use_cong==True:
                if self.use_dtc==True:
                    stdout.write('| ')
                for i in range(num_best_models):
                    stdout.write("{cong:<10}\t".format(cong=model.cong_list[i]))

            stdout.write('\n')
            stdout.flush()

