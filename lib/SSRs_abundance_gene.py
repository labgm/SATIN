#!/usr/bin/ python

# -*- coding: utf-8 -*-
""" 

author: carlos willian

data:03/05/2023
proposito: Contar a abundancia de SSRs identificados por genomas em relacao as regioes genicas 
input: SATIN_OUTPUT
output: SSR_counting.txt (tsv)
"""


def lerarquivos_SATINgbk(folder_output):
    Sequencias = {}
    
    if folder_output[:-1] != "/":
        folder_output = str(folder_output) + "/"
    list_SEQ = os.listdir(folder_output) # lista com os nomes dos arquivos dentro da pasta "folder_output"
    #print(list_SEQ)
    for i in list_SEQ:
        df_name = i # retira  a terminacao do nome ".fasta.misa"
        df_data = pd.read_csv(folder_output + str(i) + "/" + str(i) + ".coding" , delimiter = ';',skiprows=1 , names= ("N", "Start", "End", "SSR","Gene", "locus_tag", "Product"), engine ='python') #abrir arquivos e definir os nomes das colunas
        Sequencias[df_name]=df_data
    return(Sequencias, list_SEQ)



def unique_genes_SATIN(Sequencias):
    genes =[]
    for i in Sequencias.keys():
        for t in Sequencias[i]['Gene']:
            genes.append(t)


    new_genes = set(genes)
    new_genes = list(new_genes)
    
    if '-' in new_genes:
        new_genes.remove('-')
    

    return new_genes


def selecionar_SSRs_SATIN(Sequencias, u):
    SSR = []
    genomas = []
    genes = u
    
    for i in Sequencias.keys():
        selected_gene_df= pd.DataFrame(Sequencias[i])
        
        selected_gene_df = selected_gene_df.loc[selected_gene_df['Gene'] == str(genes)]
        
        SSR.append(selected_gene_df['SSR'])
        genomas.append(i)
        
    return (SSR, genomas)


def contar_SSRs_entre_genomas(Sequencias, unique_genes):
    lista_SSRs = []
    df_Seq =[]
    Sequencias1 = {}
    SSR_count ={}
    
    for u in unique_genes:
        (SSR, genomas)= selecionar_SSRs_SATIN(Sequencias, u)
        #print(u)
        
        for k in range(len(list(Sequencias.keys()))): # Pegar index de cada genoma
            df_name = list(Sequencias.keys())[k]
#             print(df_name)
            df_Seq.append(df_name) # lista com os nomes dos genomas (df_Seq)
            df_name1 = df_name + "_SSRs" # lista com os nomes dos genomas + "_SSRs"
            Sequencias1[df_name1] = []  # Adicionar as SSRs obtidas por genoma
#             print(lista_SSRs)
            for j in SSR[k]: # lista de cada genoma (j é o SSR analisado)
                Sequencias1[df_name1].append(j)
                
            lista_SSRs.append(df_name1)
            
        for a in lista_SSRs:
            seq= "Sequencias1['" + str(a) + "']" # chamar cada Sequencias[df_name1]
            index = lista_SSRs.index(a)
            seq =  eval(seq) # chamar o valor a partir da string montada
            lista_SSRs[index] = seq
        lista_SSRs = tuple(lista_SSRs)
        freq = list(map(Counter, lista_SSRs))
        SSR_count[str(u)] = {motif: [count[motif] for count in freq] for motif in {motif for count in freq for motif in count}} 
        lista_SSRs = []
        df_Seq =[]
        Sequencias1 = {}
        
    return(SSR_count)



def salvar_contagem(SSR_count):
    cabecalho = ['Gene', 'SSR']
    for k in Sequencias.keys():
        cabecalho.append(k)
    with open('SSR_counting.txt', 'w') as f:
        f.write('\t'.join(cabecalho))
        f.write('\n')
        for w in SSR_count.keys():
            gene = w
            motivos = list(SSR_count[w].keys())
            contagem =  list(SSR_count[w].values())
            for i in range(len(motivos)):
                f.write(str(gene) + '\t')
                f.write(str(motivos[i]))
                for t in range(len(contagem[i])):
                    f.write('\t' + str(contagem[i][t]))
                f.write('\n')
                
        
        f.close()


if __name__=='__main__':
    from collections import Counter
    import pandas as pd
    import os
    import sys
    import numpy as np
    
    folder_output = sys.argv[1]
    (Sequencias, list_SEQ) = lerarquivos_SATINgbk(folder_output)# chamar a funcao e definir os valores sequencias(Dicionario) e list_Seq (lista dos arquivos dentro da pasta)

    print(list_SEQ)
    unique_genes = unique_genes_SATIN(Sequencias)
    SSR_count = contar_SSRs_entre_genomas(Sequencias, unique_genes)
    salvar_contagem(SSR_count)