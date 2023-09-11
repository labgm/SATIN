import numpy as np
import pandas as pd
import os
from io import StringIO
import string
from collections import Counter
import sys

folder_output = sys.argv[1]


list_SEQ = os.listdir(folder_output+"/")
Sequencias = {}
for i in list_SEQ:
    df_name = i[:-3]
    df_data = pd.read_csv(folder_output + "/" + str(i) +  "/" + str(i) + ".coding", delimiter = '	')
    Sequencias[df_name]=df_data


tupla_SSRs = []
df_Seq =[]
for k in list_SEQ:
    df_name = k[:-3]
    df_Seq.append(df_name)
    df_name1 = k[:-3] + "_SSRs"
    SSR_list = Sequencias[df_name]["SSR"]
    Sequencias[df_name1] = []
    for j in SSR_list:
        motif = ''.join([t for t in j if not t.isdigit()])
        Sequencias[df_name1].append(motif)
    tupla_SSRs.append(df_name1)
for a in tupla_SSRs:
    seq= "Sequencias['" + str(a) + "']"
    index = tupla_SSRs.index(a)
    seq =  eval(seq)
    tupla_SSRs[index] = seq
tupla_SSRs = tuple(tupla_SSRs)

freq = list(map(Counter, tupla_SSRs))
res = {motif: [count[motif] for count in freq] for motif in {motif for count in freq for motif in count}} 

#----------

abundance = []
chaves = list(res.keys())
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

###

Index_Genes=[]
for i in range(len(tupla_SSRs)):
    Seq = df_Seq[i]
    name = str(Seq) + "_Index_Genes"
    Index_Genes.append(name)
    Sequencias[name] = []
    for u in list_index:
        li = []
        for v in u:
            index_Seq = np.where(np.array(tupla_SSRs[i]) == chaves[v])
            index_Seq = list(index_Seq)
            li.append(index_Seq)
        Sequencias[name].append(li)
Genes ={}
li =[]
for t in df_Seq:
    Genes_seq = "Genes_set_" + str(t)
    Genes_set = set(list(Sequencias[t]["Product"]))
    li.append(Genes_seq)
    Genes[Genes_seq] = Genes_set
Gen =[]

for k in li:
    Gen.append(list(Genes[k]))


freq = list(map(Counter, Gen))
gen_freq = {gen: [count[gen] for count in freq] for gen in {gen for count in freq for gen in count}}

##


abun = []
chav = list(gen_freq.keys())
li_index = []

for i in chav:
    abun.append(gen_freq[i])
for c in range(len(abun)):
    percent = ((sum(abun[c]))/(len(abun[c])))*100
    gen_freq[chav[c]] = percent

df_SSR_Seq = []
for i in df_Seq:
    l = i  + "_SSRs"
    df_SSR_Seq.append(l)

lista_arquivos =[]
for i in list_SEQ:
    filename= i[:-3]
    filename_saida = open(folder_output + "/" + i + "/" + str(filename) + "_SSRs_in_gene_abundance.txt", 'w')
    lista_arquivos.append(folder_output + "/" + i + "/" + str(filename) + "_SSRs_in_gene_abundance.txt")

cabecalho = ["Start", "End", "SSR", "Product", "SSR_count%", "Product%", "SSR_relat%", "SSR_Count_per_genome"]

for i in range(len(lista_arquivos)):
    filename_saida = open(lista_arquivos[i], 'w')
    linha = '\t'.join(cabecalho)
    filename_saida.write(linha + '\n')
    for x in range(len(list_index)):
        vec_lido = []

        for t in np.array(Sequencias[Index_Genes[i]][x][0][0]):
            soma_rel = 0
            posi_conc = str(Sequencias[df_Seq[i]]["Start"][t]) + str(Sequencias[df_Seq[i]]["End"][t])
            for m in res[Sequencias[df_SSR_Seq[i]][t]]:
                if  m != 0:
                    soma_rel += 1
            if(posi_conc not in vec_lido):
                vec_lido.append(posi_conc)
                corpo = [str(Sequencias[df_Seq[i]]["Start"][t]), str(Sequencias[df_Seq[i]]["End"][t]), str(Sequencias[df_Seq[i]]["SSR"][t]),
                        Sequencias[df_Seq[i]]["Product"][t], str(round((res[Sequencias[df_SSR_Seq[i]][t]][i]/sum(res[Sequencias[df_SSR_Seq[i]][t]]))*100, 2)) , str(gen_freq[Sequencias[df_Seq[i]]["Product"][t]]), str((soma_rel/len(res[Sequencias[df_SSR_Seq[i]][t]]))*100), str(res[Sequencias[df_SSR_Seq[i]][t]])]
                linha = '\t'.join(corpo)
                filename_saida.write(linha + '\n')
