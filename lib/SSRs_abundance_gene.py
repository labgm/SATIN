#!/usr/bin/ python

# -*- coding: utf-8 -*-
""" 

author: carlos willian

data:01/11/2023
proposito: Contar a abundancia de SSRs identificados por genomas em relacao as regioes genicas 
input: SATIN_OUTPUT
output: SSR_counting.txt (tsv)

To run:

"python lib/SSRs_abundance_gene.py OUTPUT_FOLDER/"

"""


from os import name


def lerarquivos_SATINgbk(folder_output):
    Sequencias = {}
    
    if folder_output[:-1] != "/":
        folder_output = str(folder_output) + "/"
    list_SEQ = os.listdir(folder_output) # lista com os nomes dos arquivos dentro da pasta "folder_output"
    #print(list_SEQ)
    for i in list_SEQ:
        df_name = i # retira  a terminacao do nome ".gbff"
        df_data = pd.read_csv(folder_output + str(i) + "/" + str(i) + ".coding" , delimiter = ';',skiprows=1 , names= ("N", "Start", "End", "SSR","Gene", "strand", "locus_tag", "Product", "ID"), engine ='python') #abrir arquivos e definir os nomes das colunas
        df_data['Gene'] = df_data['Gene'].astype(str).str.split('_').str[0] # Remove as terminações '_2' em genes que tenham o mesmo nome exceto estas terminações por indicar mutações
                
        Sequencias[df_name]=df_data
    return(Sequencias, list_SEQ)



def unique_genes_SATIN(Sequencias):
    genes =[]
    for i in Sequencias.keys():
        for t in Sequencias[i]['Gene']:
            genes.append(t)


    new_genes = set(genes)
    new_genes = list(new_genes)
    
    if str('nan') in new_genes:
        new_genes.remove('nan')


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

    ppp = 0
    
    for u in unique_genes:
        (SSR, genomas)= selecionar_SSRs_SATIN(Sequencias, u)
        # print(u)
        
        for k in range(len(list(Sequencias.keys()))): # Pegar index de cada genoma
            df_name = list(Sequencias.keys())[k]
#             print(df_name)
            df_Seq.append(df_name) # lista com os nomes dos genomas (df_Seq)
            df_name1 = df_name + "_SSRs" # lista com os nomes dos genomas + "_SSRs"
            Sequencias1[df_name1] = []  # Adicionar as SSRs obtidas por genoma
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




def contar_SSRs_entre_genomas_filtrado(SSR_count):
    SSR_count2 = {}

    for w in SSR_count.keys():
        gene = w
        
        motivos = list(SSR_count[w].keys())
        contagem =  list(SSR_count[w].values())
        SSR_count2[gene] = {}

        if len(motivos) > 1:

            motivos_em_genes = []
            
            only_alpha = ""
            numeros_motivos =[]
            for j in motivos: # lista de cada genoma (j é o SSR analisado)
                
                for m in j:
                    if ord(m) >=65 and ord(m)<=90:
                        only_alpha += m
                    elif ord(m) >= 97 and ord(m) <=122:
                        only_alpha += m
                
                number = re.findall(r'\d+', j)[0]
                numeros_motivos.append(number)
            
                motivos_em_genes.append(str(only_alpha))
                only_alpha = ""
            
            remover = []
            scanear_motivos_em_genes = motivos_em_genes.copy()
            for a in range(len(motivos_em_genes)-1):
                if len(scanear_motivos_em_genes) > a:
                    del scanear_motivos_em_genes[0] # remover a primeira sequencia que sera buscada para evitar que de hit em sequencias palindromas
                else:
                    continue
                motif = Seq(motivos_em_genes[a])
                if str(motif.complement()) in scanear_motivos_em_genes:
                    #print("COMPLEMENTO!!")
                    #print(scanear_motivos_em_genes, "scanear_motivos_em_genes")
                    indexes = int(len(motivos_em_genes)-1) - int(motivos_em_genes[::-1].index(str(motif.complement()))) # index da posicao do complemento da sequencia
                    #print(indexes)
                    #print(a, "valor do A")
                    
                    if numeros_motivos[a] == numeros_motivos[indexes]:
                        motivo_to_ad = motivos[a]
                        print(motivo_to_ad)
                        lista1= contagem[a]
                        lista2= contagem[indexes]
                        lista_sum = [x + y for x, y in zip(lista1, lista2)] # Somar duas listas
                        
                        contagem[a]= lista_sum
                        motivos[a]= motivo_to_ad

                        motivos_em_genes[a] = motivos_em_genes[indexes]
                        remover.append(indexes)
                        continue

                    else:
                        continue

                elif str(motif.complement())[::-1] in scanear_motivos_em_genes:
                    indexes = int(len(motivos_em_genes)-1) - int(motivos_em_genes[::-1].index(str(motif.complement())[::-1])) # index da posicao do inverso da sequencia
                    if numeros_motivos[a] == numeros_motivos[indexes]:

                        motivo_to_ad = motivos[a]
                        print(motivo_to_ad)
                        lista1= contagem[a]
                        lista2= contagem[indexes]
                        lista_sum = [x + y for x, y in zip(lista1, lista2)] # Somar duas listas
                        
                        contagem[a]= lista_sum
                        motivos[a]= motivo_to_ad
                        motivos_em_genes[a] = motivos_em_genes[indexes]
                        remover.append(indexes)
                        continue
                    else:
                        continue

                elif str(motif)[::-1] in scanear_motivos_em_genes:
                    indexes = int(len(motivos_em_genes)-1) - int(motivos_em_genes[::-1].index(str(motif)[::-1])) # index da posicao do inverso da sequencia
                    if numeros_motivos[a] == numeros_motivos[indexes]:

                        motivo_to_ad = motivos[a]
                        lista1= contagem[a]
                        lista2= contagem[indexes]
                        lista_sum = [x + y for x, y in zip(lista1, lista2)] # Somar duas listas
                        
                        contagem[a]= lista_sum
                        motivos[a]= motivo_to_ad
                        motivos_em_genes[a] = motivos_em_genes[indexes]
                        remover.append(indexes)
                        
                        continue
                    else:
                        continue
                elif str(motif[-1] + motif[0:-1]) in scanear_motivos_em_genes:
                    indexes = int(len(motivos_em_genes)-1) - int(motivos_em_genes[::-1].index(str(motif[-1] + motif[0:-1]))) # index da posicao do complemento da sequencia

                    if numeros_motivos[a] == numeros_motivos[indexes]:
                        motivo_to_ad = motivos[a]
                        lista1= contagem[a]
                        lista2= contagem[indexes]
                        lista_sum = [x + y for x, y in zip(lista1, lista2)] # Somar duas listas
                        
                        contagem[a]= lista_sum
                        motivos[a]= motivo_to_ad

                        
                        motivos_em_genes[a] = motivos_em_genes[indexes]
                        remover.append(indexes)
                        continue

                    else:
                        continue
                elif str(motif[-2:] + motif[:-2]) in scanear_motivos_em_genes:
                    indexes = int(len(motivos_em_genes)-1) - int(motivos_em_genes[::-1].index(str(motif[-2:] + motif[:-2]))) # index da posicao do complemento da sequencia
                    
                    if numeros_motivos[a] == numeros_motivos[indexes]:
                        motivo_to_ad = motivos[a]
                        lista1= contagem[a]
                        lista2= contagem[indexes]
                        lista_sum = [x + y for x, y in zip(lista1, lista2)] # Somar duas listas
                        
                        contagem[a]= lista_sum
                        motivos[a]= motivo_to_ad

                        
                        motivos_em_genes[a] = motivos_em_genes[indexes]
                        remover.append(indexes)
                        continue

                    else:
                        continue

                elif str(motif[-3:] + motif[:-3]) in scanear_motivos_em_genes:
                    indexes = int(len(motivos_em_genes)-1) - int(motivos_em_genes[::-1].index(str(motif[-3:] + motif[:-3]))) # index da posicao do complemento da sequencia
                    if numeros_motivos[a] == numeros_motivos[indexes]:
                        motivo_to_ad = motivos[a]
                        lista1= contagem[a]
                        lista2= contagem[indexes]
                        lista_sum = [x + y for x, y in zip(lista1, lista2)] # Somar duas listas
                        
                        contagem[a]= lista_sum
                        motivos[a]= motivo_to_ad

                        
                        motivos_em_genes[a] = motivos_em_genes[indexes]
                        remover.append(indexes)
                        continue

                    else:
                        continue
                else:

                    continue
            
            remover = sorted(remover, reverse= True)
            if len(remover) > 1:
                del remover[0] # selecionar o primeiro motivo para ser mantido em SSR_count2

            for motivoaremover in remover:
                del motivos[motivoaremover]
                del contagem[motivoaremover]
            
            for z in range(len(motivos)):

                SSR_count2[gene].update({motivos[z]:contagem[z]})
                
            
        else:
            SSR_count2[gene]= {motivos[0]:contagem[0]}
            continue

    return(SSR_count2)



def salvar_contagem2(SSR_count2):
    cabecalho = ['Gene', 'SSR']
    for k in Sequencias.keys():
        cabecalho.append(k)
    with open('SSR_counting.txt', 'w') as f:
        f.write('\t'.join(cabecalho))
        f.write('\n')
        for w in SSR_count2.keys():
            gene = w
            motivos = list(SSR_count2[w].keys())
            contagem =  list(SSR_count2[w].values())
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
    import re
    import numpy as np
    from Bio.Seq import Seq
    
    
    folder_output = sys.argv[1]
    (Sequencias, list_SEQ) = lerarquivos_SATINgbk(folder_output)# chamar a funcao e definir os valores sequencias(Dicionario) e list_Seq (lista dos arquivos dentro da pasta)

    unique_genes = unique_genes_SATIN(Sequencias)
    SSR_count = contar_SSRs_entre_genomas(Sequencias, unique_genes)
    # salvar_contagem(SSR_count)

    SSR_count2 = contar_SSRs_entre_genomas_filtrado(SSR_count)
    salvar_contagem2(SSR_count2)
    
