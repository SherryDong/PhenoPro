import cPickle as pickle
from scipy.stats.mstats import ks_2samp
from numpy import array
import sys


class PhenoBayes_Data:
            def __init__(self,HPOs,Diseases,Genes,Diseases_all_Ps,Diseases_Genes):
                  self.HPOs=HPOs
                  self.Diseases=Diseases
                  self.Genes=Genes
                  self.Diseases_all_Ps=Diseases_all_Ps
                  self.Diseases_Genes=Diseases_Genes


class HPO_Class:
            def __init__(self, _id=[], _name=[], _alt_id=[],  _def=[], _comment=[], _synonym=[], _xref=[], _is_a=[],_alt_Hs={}):
                  self._id = _id
                  self._name = _name
                  self._alt_id = _alt_id
                  self._def = _def
                  self._comment = _comment
                  self._synonym = _synonym
                  self._xref = _xref
                  self._is_a = _is_a
                  self._father=set()
                  self._disease=set()
                  self._child_self=set()
                  self._alt_Hs= _alt_Hs


def loading(obo_file,HPO_disease_gene_file,save_file_dir='./'):
      
      HPOs={}
      _id=[];_name=[];_alt_id=[];_def=[];_comment=[];_synonym=[];_xref=[];_is_a=[];
      fi=open(obo_file)
      obo_terms=fi.read().split('[Term]')
      fi.close()
      alt_Hs={}
      for term in obo_terms:
            if 'id: ' in term:
                  seq=term.split('\n')
                  for one in seq:
                        if ': ' in one:    
                              if 'id: ' in one and 'alt_id: '  not in one:
                                    _id.append(one.split(': ')[1])
                              if 'name: ' in one:
                                    _name.append(one.split(': ')[1])
                              if 'alt_id: ' in one:
                                    alt_Hs[one.split(': ')[1]] = _id[-1]
                                    _alt_id.append(one.split(': ')[1])
                              if 'def: ' in one:
                                    _def.append(one.split(': ')[1])
                              if 'comment: ' in one:
                                    _comment.append(one.split(': ')[1])
                              if 'synonym: ' in one:
                                    if '"' in one:
                                        _synonym.append(one.split(': ')[1].split('"')[1])
                                        #print one
                                    else:
                                        _synonym.append(one.split(': ')[1])
                              if 'xref: ' in one:
                                    _xref.append(one.split(': ')[1])
                              if 'is_a: ' in one:
                                    _is_a.append(one.split(': ')[1].split('!')[0].replace(' ',''))
      
                  HPOs[_id[0]]=HPO_Class(_id,_name,_alt_id,_def,_comment,_synonym,_xref,_is_a)
                  _id=[];_name=[];_alt_id=[];_def=[];_comment=[];_synonym=[];_xref=[];_is_a=[];
      
      HPOs['HP:0000118']._alt_Hs=alt_Hs
      def find_father(_ori_id,_id):
             
            if  HPOs[_id]._name=='All':
                  pass
                  
            else:
                  for one in HPOs[_id]._is_a:
                              HPOs[_ori_id]._father.add(one)
                              find_father(_ori_id,one)  
      
      for _id in HPOs:
            find_father(_id,_id)
      
      
      fi=open(HPO_disease_gene_file)
      Diseases={}
      Genes={}
      Diseases_Genes={}
      for line in fi:
            if line[0] != '#':
                  seq=line.replace('\n','').split('\t')
                  for one in seq:
                        if 'OMIM:' in one:
                              D=one
                        elif 'HP:' in one:
                              H=one
                  G=seq[1]
                  try:
                              Diseases[D].add(H)
                  except Exception, e:
                              Diseases[D]=set()
                              Diseases[D].add(H)
                  try:
                              Diseases_Genes[D].add(G)
                  except Exception, e:
                              Diseases_Genes[D]=set()
                              Diseases_Genes[D].add(G)
                  try:
                              Genes[G].add(D)
                  except Exception, e:
                              Genes[G]=set()
                              Genes[G].add(D)
                  
                  HPOs[H]._disease.add(D)
                  D='';H='';

      #Remove ans
      print "remove ans ..."
      fo=open('Disease_spHPO.txt','w')
      for D in Diseases:
            tmp=set()
            for H in Diseases[D]:
                if H in alt_Hs:
                    tmp.add(alt_Hs[H])
                else:
                    tmp.add(H)
            Diseases[D]=tmp
            tmp_ans=set() 
            tmp=set()
            for H in Diseases[D]:
                  tmp_ans = tmp_ans | HPOs[H]._father
            for H in Diseases[D]:
                  if H not in tmp_ans:
                        tmp.add(H)
            Diseases[D]=tmp
            fo.write(D)
            for H in Diseases[D]:
                  fo.write('\t'+H)
            fo.write('\n')



      j=0 
      for asHPO in HPOs:
            #if asHPO in HPOs['HP:0000118']._child_self:    
                  for HPO in HPOs:
                        if asHPO in HPOs[HPO]._father:
                              HPOs[asHPO]._child_self.add(HPO)
                  HPOs[asHPO]._child_self.add(asHPO) 
                  j=j+1;print " Finding HPOs' child nodes: "+str(j)+' HPOs \r',



      def P_asHPO_on_spHPO(asHPO,spHPO):
            if asHPO in HPOs[spHPO]._father or asHPO == spHPO:
                  return 1/float(len(HPOs[spHPO]._father)+1)
            else:
                  return 0.0

      def P_Disease_on_spHPO(Disease,spHPO):
            if spHPO in Diseases[Disease]:
                  return 1/float(len(HPOs[spHPO]._disease))
            else:
                  return 0.0


      SUM_HPO=0
      for one in HPOs:
            SUM_HPO=SUM_HPO+len(HPOs[one]._disease)
      def P_spHPO(spHPO):
            return float(len(HPOs[spHPO]._disease))/SUM_HPO



      def P_Disease_on_asHPO(Disease,asHPO):
            UP=0.0
            DOWN=0.0
            for spHPO in HPOs[asHPO]._child_self:
                  if spHPO in Diseases[Disease]:
                        UP=UP + P_Disease_on_spHPO(Disease,spHPO) * P_asHPO_on_spHPO(asHPO,spHPO) * P_spHPO(spHPO)
                  DOWN=DOWN + P_asHPO_on_spHPO(asHPO,spHPO) * P_spHPO(spHPO)
            if DOWN==0:
                  return 0

            P=UP/DOWN
            return  P

      P_Disease_on_asHPOs={}
      print ''
      t=0
      for asHPO in HPOs:
            
            if asHPO in HPOs['HP:0000118']._child_self:
                  t=t+1 ; print " Calculating P(Disease|asHPO): "+str(t)+' HPOs \r',
                  for Disease in Diseases:
                        P_Disease_on_asHPOs[Disease+':'+asHPO]=P_Disease_on_asHPO(Disease,asHPO)
                   #     print len(P_Disease_on_asHPOs)
      print ''


      Diseases_all_Ps={}
      q=0
      for Disease in Diseases:
            Diseases_all_Ps[Disease]={}
            q=q+1;print " Assigning P(Disease|asHPO): "+str(q)+' Diseases \r',
            HPOs_Ps=[]
            for HPO in HPOs:
                  if HPO in HPOs['HP:0000118']._child_self:
                        HPOs_Ps.append([P_Disease_on_asHPO(Disease,HPO),HPO])
            HPOs_Ps.sort(reverse=True)
           
            i=0
            tmp=1.1
            while i<len(HPOs_Ps):
                 if HPOs_Ps[i][0] != 0:
                       Diseases_all_Ps[Disease][HPOs_Ps[i][1]]=HPOs_Ps[i][0]
                 i=i+1
          
            Diseases_all_Ps[Disease]['LENGTH']=i
            

      fo=open(save_file_dir+'PhenoBayes_Data.pk','w')
      data = PhenoBayes_Data(HPOs,Diseases,Genes,Diseases_all_Ps,Diseases_Genes)
      pickle.dump(data,fo)
      fo.close()
      return data
      


def Get_record(given_HPOs,given_Ps,all_Ps):
    tmp_Ps=all_Ps
    tmp_Ps.sort()
    tmp_Ps=tmp_Ps[::-1]
    N=float(len(all_Ps))
    #print tmp_Ps[0]
    n=len(given_Ps)
    tmp_record=['']*n
    
    j=0
    while j < n:
        gvp=given_Ps[j]
        i=0
        
        while i < N:
            
            if tmp_Ps[i] <= gvp and gvp==0:
                tmp_record[j]=''#given_HPOs[j]+'_PR_'+str(given_Ps[j])+'_rank_'+str(int((float(i)+float(N))/2.0))+'('+str(int(N))+')'
                break
            if tmp_Ps[i] <= gvp:
                t=0
                tt=0
                for one in tmp_Ps[i+1:]:
                    if tmp_Ps[i]==one:
                        t=t+1
                        tt=tt+1
                    elif one != 0:
                        tt=tt+1
                    else:
                        break
                this_rank=max([i,1])#int(float(i+i+t)/(2))
                none_zero=int(tt+i)
                tmp_record[j]=given_HPOs[j]+'_Pr_'+str(round(given_Ps[j]*100,2))+'%_Rank_'+str(this_rank)  + '('+str(int(none_zero))+')'
                break
            i += 1 
        j += 1
    tmp=tmp_record
    tmp_record=[]
    for one in tmp:
        if one!='':
            tmp_record.append(one) 
    if tmp_record==[]:
        tmp_record=['None']
    return ';'.join(tmp_record)

def Ranked_Score_Disease_Pheno(given_HPOs,data):

#      fdata = open('/home/zhangfeng/server/applications/phenobayes/scripts/PhenoBayes_Data.pk')
#      data = pickle.load(fdata)
#      fdata.close()
      HPOs=data.HPOs
      Diseases=data.Diseases
      Genes=data.Genes
      Diseases_all_Ps=data.Diseases_all_Ps
      Diseases_Genes=data.Diseases_Genes
      tmp=[]
      for one in given_HPOs:
          tmp.append(one)
      given_HPOs=tmp
      def  P_value(Disease):
            all_Ps=[]
            given_Ps=[]
            for HPO in Diseases_all_Ps[Disease]:
                  if  HPO!='LENGTH':
                        all_Ps.append(Diseases_all_Ps[Disease][HPO])
            all_Ps=all_Ps+[0]*(Diseases_all_Ps[Disease]['LENGTH']-len(all_Ps))
            for given_HPO in given_HPOs:
                  if given_HPO in HPOs['HP:0000118']._alt_Hs:
                        given_HPO = HPOs['HP:0000118']._alt_Hs[given_HPO]
                  if given_HPO in Diseases_all_Ps[Disease]:
                        given_Ps.append(Diseases_all_Ps[Disease][given_HPO])
                  else:
                        given_Ps.append(0)

            all_Ps=array(all_Ps)
            given_Ps=array(given_Ps)
            record=Get_record(given_HPOs,given_Ps,all_Ps)
            return [ks_2samp(all_Ps,given_Ps,alternative='greater')[1],record]
      #fi.close()
      P_values=[]
      d=0
      for Disease in Diseases:
      #      d += 1;print d
            if len(Diseases[Disease])>0:
                 tmp_result=P_value(Disease)
                 P_values.append([tmp_result[0],Disease,tmp_result[1]])
    #             print " Calculating p-values: "+str(len(P_values))+' Dieases \r',
      P_values.sort()
      output=[one+[','.join([d for d in Diseases_Genes[one[1]]])] for one in P_values]
      return output


#loading(sys.argv[1],sys.argv[2],'./')









def wordscore(term1,term2):
        term1=term1.replace('','').replace('\r','').replace('\n','').lower()
        term2=term2.replace('','').replace('\r','').replace('\n','').lower()
        overlap=[]
        j_cutoff=min(max(len(term1)/2,1),5)
        i=0
        tmp_end=0
        while i < len(term1):
                j=j_cutoff
                tmp=[0]
                tmp_j=[0]
                tmp_end_lst={}
                tmp_end_lst[0]=[0]
                while j < len(term1)-i+1:
                        if term1[i:i+j] in term2[tmp_end:]:
                                flag=0
                                if flag!=1:
                                        tmp_j.append(j)
                                        tmp.append(j/float(i+1))
                                        try:
                                                tmp_end_lst[j].append(term2.find(term1[i:i+j]))
                                        except Exception, e:
                                                tmp_end_lst[j]=[term2.find(term1[i:i+j])]
                        j=j+1
                overlap.append(max(tmp))
                i=i+max(tmp_j)+1
                tmp_end=tmp_end_lst[max(tmp_j)][0]+max(tmp_j)+1
                wscore=sum(overlap)/float(len(term1)+len(term2)-sum(overlap))
                if wscore<0.2:
                        wscore=0
        return wscore


def compareterm(term1,term2):
        words1=term1.lower().replace('"','').split(' ')
        words2=term2.lower().replace('"','').split(' ')
        score={}
        for word1 in words1:
                tmp_score=[]
                tmp_word={}
                for word2 in words2:
                        wscore=wordscore(word1,word2)
                        tmp_score.append(wscore)
                        tmp_word[wscore]=word2
                if max(tmp_score) >= 0:
                        score[word1]=[max(tmp_score),tmp_word[max(tmp_score)]]
        SCORE=sum( [ score[w][0] for w in score])
        return SCORE


def pure_words(words):
        words=words.rstrip()
        words=words.replace('-',' ')
        words=words.replace(',','').replace(':','').replace('/','').replace('(','').replace(')','')
        words=":".join(words.split())
        seq=words.split(':')
        remove_lst=['of','the','an','a','with','in']
        i=0
        while i< len(seq):
                if seq[i] in remove_lst:
                        del seq[i]
                else:
                        i=i+1
        words=""
        for one in seq:
                words=words+one+' '
        words=words[:-1]
        #words=words.replace(':',' ')
        return words




def interpreting(keywords,data):
      keywords=pure_words(keywords.lower())
      if keywords=='':
          return ['None']
      score=[]
      HPOs=data.HPOs
      for HPO in HPOs:
            tmp=[0]
            check_list=HPOs[HPO]._name+HPOs[HPO]._synonym
            for term in check_list:
                  tmp.append(compareterm(keywords,term.lower()))
                  
            score.append([max(tmp),HPO])
#            if HPO=='HP:0001250':
#	               print HPOs[HPO]._synonym
 #                      print tmp
      score.sort(reverse=True)
      output=[]
      final_output=[] 
      if score[0][0]<=0.5:
            final_output=['None']
      else:
            tmp=score[0][0]
            i=0
            ori_output=[]
            while score[i][0]==tmp:
                  ori_output.append(score[i][1])
                  i=i+1
            child_group=set()
            #print HPOs['HP:0012758']._child_self
            #print HPOs['HP:0001263']._child_self
            for HPO in ori_output:
                this_HPO=set()
                this_HPO.add(HPO)
                tmp = HPOs[HPO]._child_self - this_HPO
                child_group=child_group | tmp
                #print child_group
            #print child_group
            for one in ori_output:
                if one not in child_group:
                    output.append(one)
            # deal with local opt
            limit_length=1
            if len(output)>limit_length:
                tmp=[]
                for one in output:
                    father_len=len(HPOs[one]._father)
                    tmp.append([father_len/float(father_len+len(HPOs[one]._child_self)),one])
                tmp.sort()
                for one in tmp[0:limit_length]:
                    final_output.append(one[1])
                combined_father=set()
                this_HPO=set()
                this_HPO.add(tmp[limit_length][1])
                combined_father=HPOs[tmp[limit_length][1]]._father | this_HPO
                for one in tmp[limit_length:]:
                    this_HPO=set()
                    this_HPO.add(one[1])
                #    print this_HPO
                    combined_father=combined_father & (HPOs[one[1]]._father | this_HPO)
                tmptmp=[]
                for one in combined_father:
                    tmptmp.append([len(HPOs[one]._father),one])
                #print tmptmp
                tmptmp.sort()
                final_output.append(tmptmp[-1][1])
                tmp=[]
                for one in  final_output:
                   if 'HP:0000118' not in HPOs[one]._child_self:
                       tmp.append(one)
                final_output=tmp
                #print final_output
            else:
                final_output=output

      return final_output


            

