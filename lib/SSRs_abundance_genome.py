import numpy as np
import pandas as pd
import os
from io import StringIO
import string
from collections import Counter
import sys

folder_output = sys.argv[1]

if folder_output[-1:] == "/":
    list_SEQ = os.listdir(folder_output)
else:
    list_SEQ = os.listdir(folder_output + "/")


Sequencias = {}
for i in list_SEQ:
    if i[-3:] == '.fa':
        df_name = i[:-3]
    if i[-6:] == '.fasta':
        df_name = i[:-6]
    if i[-4:] == '.fna':
        df_name = i[:-4]
    if i[-3:] == '.gb':
        df_name = i[:-3]
    if i[-4:] == '.gbk':
        df_name = i[:-4]
    
    df_data = pd.read_csv(folder_output + "/" + str(i) +  "/" + "output.satin", delimiter = '	', skiprows=3)
    
    Sequencias[df_name]=df_data

#------------------------------------
tupla_SSRs = []
df_Seq =[]
for k in list_SEQ:
    if k[-3:] == '.fa':
        df_name = k[:-3]
    if k[-6:] == '.fasta':
        df_name = k[:-6]
    if k[-4:] == '.fna':
        df_name = k[:-4]
    if k[-3:] == '.gb':
        df_name = k[:-3]
    if k[-4:] == '.gbk':
        df_name = k[:-4]
    df_Seq.append(df_name)
    df_name1 = df_name + "_SSRs"
    SSR_list = Sequencias[df_name]["SSR"]
    Sequencias[df_name1] = []
    for j in SSR_list:
        #motif = ''.join([t for t in j if not t.isdigit()])
        Sequencias[df_name1].append(j)
    tupla_SSRs.append(df_name1)
for a in tupla_SSRs:
    seq= "Sequencias['" + str(a) + "']"
    index = tupla_SSRs.index(a)
    seq =  eval(seq)
    tupla_SSRs[index] = seq
tupla_SSRs = tuple(tupla_SSRs)
freq = list(map(Counter, tupla_SSRs))
res = {motif: [count[motif] for count in freq] for motif in {motif for count in freq for motif in count}} 

#--------------------------

abundance = []
chaves = list(res.keys()) #motivos
list_index = []
abundance_unique =[]
for i in chaves:
    abundance.append(res[i])
abundance_sort = sorted(abundance, reverse= True)
for i in abundance_sort:
    if i not in abundance_unique:
        abundance_unique.append(i)

for t in abundance_unique:
    index = [ i for i in range(len(abundance_sort)) if abundance[i] == t ]
    list_index.append(index)

###--------------------------------

Index_Genome=[]
for i in range(len(df_Seq)):
    Seq = df_Seq[i]
    name = str(Seq) + "_Index_Genome"
    Index_Genome.append(name)
    Sequencias[name] = []
    for u in list_index:
        for v in u:
            li=[]
            index_Seq = np.where(np.array(tupla_SSRs[i]) == chaves[v])
            index_Seq = list(index_Seq)
            li.append(index_Seq)
        Sequencias[name].append(li)
        
##------------------------------------------
SSR ={}
li =[]
for t in df_Seq:
    SSR_seq = "SSR_set_" + str(t)
    SSR_set = set(list(Sequencias[t]["SSR"]))
    li.append(SSR_seq)
    SSR[SSR_seq] = SSR_set
Gen =[]

for k in li:
    Gen.append(list(SSR[k]))
freq = list(map(Counter, Gen))
ssr_freq = {gen: [cnt[gen] for cnt in freq] for gen in {gen for cnt in freq for gen in cnt}}

abun = []
chav = list(ssr_freq.keys())
li_index = []

for i in chav:
    abun.append(ssr_freq[i])
for c in range(len(abun)):
    percent = ((sum(abun[c]))/(len(abun[c])))*100
    ssr_freq[chav[c]] = percent
    



####--------------------------------------
df_SSR_Seq = []
for i in df_Seq:
    l = i  + "_SSRs"
    df_SSR_Seq.append(l)


######----------------------------------------------------------------------------------   
lista_arquivos =[]
for i in list_SEQ:
    if i[-3:] == '.fa':
        filename = i[:-3]
    if i[-6:] == '.fasta':
        filename = i[:-6]
    if i[-4:] == '.fna':
        filename = i[:-4]
    if i[-3:] == '.gb':
        filename = i[:-3]
    if i[-4:] == '.gbk':
        filename = i[:-4]
    
    filename_saida1 = open(folder_output + "/" + i + "/" + str(filename) + "_SSRs_in_genome_abundance.txt", 'w')
    lista_arquivos.append(folder_output + "/" + i + "/" + str(filename) + "_SSRs_in_genome_abundance.txt")



###--------------------------------------------



cabecalho = ["Start", "End","SSR", "SSR%freqin_genome/allgenomes", "SSR-presence-ausence", "SSR%_-presence-ausence-genomes"]
for i in range(len(lista_arquivos)):
    filename_saida1 = open(lista_arquivos[i], 'w')
    linha = '\t'.join(cabecalho)
    filename_saida1.write(linha + '\n')
    for x in range(len(list_index)):
        vec_lido = []
        
        for t in np.array(Sequencias[Index_Genome[i]][x][0][0]):
            posi_conc = str(Sequencias[df_Seq[i]]["Start"][t]) + str(Sequencias[df_Seq[i]]["End"][t])
            if(posi_conc not in vec_lido):
                vec_lido.append(posi_conc)
            soma_rel = 0
            for m in res[Sequencias[df_SSR_Seq[i]][t]]:
                if  m != 0:
                    soma_rel += 1
                #print(soma_rel)
                corpo = [str(Sequencias[df_Seq[i]]["Start"][t]), str(Sequencias[df_Seq[i]]["End"][t]), str(Sequencias[df_Seq[i]]["SSR"][t]),
                        str((res[Sequencias[df_SSR_Seq[i]][t]][i]/sum(res[Sequencias[df_SSR_Seq[i]][t]]))*100) , str(ssr_freq[Sequencias[df_Seq[i]]["SSR"][t]]), str(soma_rel/len(res[Sequencias[df_SSR_Seq[i]][t]]))]
                #print(str(soma_rel/len(res[Sequencias[df_SSR_Seq[i]][t]])))
                linha = '\t'.join(corpo)
                filename_saida1.write(linha + '\n')
                

##------------------------------------------


lista_arquivos_SSR =[]
for i in list_SEQ:
    if i[-3:] == '.fa':
        filename = i[:-3]
    if i[-6:] == '.fasta':
        filename = i[:-6]
    if i[-4:] == '.fna':
        filename = i[:-4]
    if i[-3:] == '.gb':
        filename = i[:-3]
    if i[-4:] == '.gbk':
        filename = i[:-4]
    
    filename_saida2 = open(folder_output + "/" + i + "/" + str(filename) + "_SSRs_frequence_per_genome.txt", 'w')
    lista_arquivos_SSR.append(folder_output + "/" + i + "/" + str(filename) + "_SSRs_frequence_per_genome.txt")           


###--------------------------------------------

cabecalho2 = ["SSR"]
for i in df_Seq:
    cabecalho2.append(i)


for i in range(len(lista_arquivos_SSR)):
    filename_saida2 = open(lista_arquivos_SSR[i], 'w')
    linha = '\t'.join(cabecalho2)
    filename_saida2.write(linha + '\n')
        
    for t in set(Sequencias[df_SSR_Seq[i]]):
        corpo = [str(t)]
        for u in res[t]:
            corpo.append(str(u))
        
        linha = '\t'.join(corpo)
        filename_saida2.write(linha + '\n')
        
