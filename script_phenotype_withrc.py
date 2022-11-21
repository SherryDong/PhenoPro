import cPickle as pickle
import sys,nltk,re
import phenobayes_withrc as phenobayes

class PhenoBayes_Data:
            def __init__(self,HPOs,Diseases,Genes,Diseases_all_Ps):
                  self.HPOs=HPOs
                  self.Diseases=Diseases
                  self.Genes=Genes
                  self.Diseases_all_Ps=Diseases_all_Ps


class HPO_Class:
            def __init__(self, _id=[], _name=[], _alt_id=[],  _def=[], _comment=[], _synonym=[], _xref=[], _is_a=[]):
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


fi=open(sys.argv[1])
given_terms=set()
for line in fi:
    if line[0]!='#':
        seq=re.split(';|,|\*|/|\.',line.replace('\n','').rstrip().lstrip())
        for one in seq:
              if 'HP:' in one:
                    one=one.replace(' ','')
              given_terms.add(one)

fdata = open('PhenoBayes_Data.pk')
data = pickle.load(fdata)
fdata.close()

trim_lst=[',','related','activity','deficiency','susceptibility','included','autosomal','abnormality','abnormal','deficiencies','type','due','disorder','features','year','years','malformation','malformations','failure','syndrome','disease']
def trim_term(term):
    seq=term.lower().split(' ')
    output=[]
    for one in seq:
        if one not in trim_lst:
            output.append(one)

    return ' '.join(output)


def get_tag(sen):
    output=[]
    seq=sen.split(' ')
    for one in seq:
        output.append(nltk.pos_tag(nltk.word_tokenize(one))[0])
    return output




def select_words(line):
    line=line.lower()
    tags = get_tag(line.strip())
    tmp=[[]]
    tmp_jj=""
    flag=0
    tmp_jj_jj=""
#    print tags
    for one in tags:
        if  one[1] =='JJ' or 'JJ' in one[1]:
          if len(one[0])  <4 or len(one[0])>3 and one[0][-3:] != "ing":
            tmp_jj_jj=tmp_jj
            tmp_jj=one[0]
          
        elif (len(one[0])>3 and one[0][-3:]== "ing") or one[1] =='NN' or one[1] == 'NNS' or 'NN' in one[1]:
            if tmp_jj!='':
                tmp[-1].append(tmp_jj)
            tmp[-1].append(one[0])
            tmp_jj_jj=tmp_jj
            tmp_jj=''
        elif one[1] =='CC' or one[1] =='IN' or one[1] =='TO':
            if tmp_jj!='':
                if tmp_jj_jj !='':
                    tmp[-1].append(tmp_jj_jj)
                    tmp[-1].append(tmp_jj)
                else:
                    tmp[-1].append(tmp_jj)

            tmp_jj=""
            tmp_jj_jj=""
            tmp.append([])
    if tmp_jj!='':
         if tmp_jj_jj !='': 
             tmp[-1].append(tmp_jj_jj)
             tmp[-1].append(tmp_jj)
         else: 
             tmp[-1].append(tmp_jj)

    output=[]
    for one in tmp:
        if len(one)>0:
            output.append(' '.join(one))
    return output






import re


given_HPOs=[]
interpret_record=[]
for one in given_terms:
  if len(one)>0:  
    if 'HP:' not in one :
        one = trim_term(one)
        a1 = re.compile('\(.*\)' )
        one = a1.sub('',one)
        a2 = re.compile('\{.*\}' )
        one = a2.sub('',one)
        a3 = re.compile('\[.*\]' )
        one = a3.sub('',one)

        
        for words in select_words(one):
            intd=phenobayes.interpreting(words,data)
            if intd[0]!='None':
                given_HPOs=given_HPOs+intd
            interpret_record.append([words,','.join(intd)])



    elif "|" in one:
        given_HPOs.append(one.split('|')[1].split('(')[0])
    else:
        given_HPOs.append(one)


old=set()
tmp=[]
for one in given_HPOs:
      if one not in old and one in data.HPOs:
            tmp.append(one)
            old.add(one)

given_HPOs=tmp


if len(given_HPOs)==0:
      result=[[1]]
else:
      result=phenobayes.Ranked_Score_Disease_Pheno(given_HPOs,data)


OMIM_name={}
fa=open('OMIMID_Name.txt')
for line in fa:
      seq=line.strip().split('\t')
      OMIM_name[seq[0]]=seq[1]
fa.close()

def OMIM_NAME(OMIM):
      name=''
      if OMIM in OMIM_name:
            name=OMIM_name[OMIM].split(';')[0].lower().rstrip().lstrip()
            if OMIM.split(':')[1] in name:
                  name=name.replace(OMIM,'')
                  name=name.replace(OMIM.split(':')[1],'')
                  name=name[1:]
                  name=name.rstrip().lstrip()
            tmp=""
            for i in range(len(name)):
                  if i==0:
                        tmp += name[0].upper()
                  elif name[i-1]==' ':
                        tmp += name[i].upper()
                  else:
                        tmp += name[i]
            name=tmp
      else:
            name='None'
      return name





fo=open(sys.argv[2],'w')
fo.write('#Given_term\n')
for one in given_terms:
    fo.write(one+'\n')
fo.write('#Interpreted_term\tHPOs\n')    
for one in interpret_record:
    fo.write(one[0]+'\t'+one[1]+'\n')
fo.write('#Given_HPO\tHPO_name\n')
for one in given_HPOs:
    fo.write(one+'\t'+data.HPOs[one]._name[0]+'\n')
fo.write('#Disease\tDisease_name\tP_value\tCausal_genes\tHPO_P_R\n')




if given_HPOs==['HP:0000118']:
    for one in result:
        fo.write(one[1]+'\t'+OMIM_NAME(one[1])+'\t'+'0.0'+'\t'+one[3]+'\t'+one[2]+'\n')
else:
    for one in result:
    #if one[0]!=1.0:
        fo.write(one[1]+'\t'+OMIM_NAME(one[1])+'\t'+str(one[0])+'\t'+one[3]+'\t'+one[2]+'\n')



fi.close()
fo.close()
