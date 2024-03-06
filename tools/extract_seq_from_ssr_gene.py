

""" 

author: SATIN

data:21/11/2023
proposito: EXTRACT THE FLANKING REGIONS OF THE SELECTED SSR 
input: Genome_types.txt, Folder_output, Selected_SSR and gene_name
output: SSR_gene_contigs_for_SSR.fasta

To run:

BEFORE YOU RUN FIRST PLEASE EDIT THE NAME OF THE GENOMES_GROUP FILE ON LINE 32 OF THIS SCRIPT 
THEN RUN THE COMAND ACCORDINGLY

''

python tools/extract_seq_from_ssr_gene.py output_folder/ "SSR" gene

''

"""


import pandas as pd
import os
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

file_output = "ecoli_types.txt"
folder_output = str(sys.argv[1])
SSR_sig = str(sys.argv[2])
GENE_sig = str(sys.argv[3])



def lerarquivos_SAT(folder_output):
    SAT = {}

    if folder_output[:-1] != "/":
        folder_output = str(folder_output) + "/"
    list_SEQ = os.listdir(folder_output) # lista com os nomes dos arquivos dentro da pasta "folder_output"
    
    for i in list_SEQ:
        df_name = i # retira  a terminacao do nome ".coding"  
        df_data = pd.read_csv(folder_output + str(i) + "/" + str(i) + ".coding" , delimiter = ';', skiprows=1 , names= ("INDEX", "START", "END", "SSR", "GENE", "STRAND", "LOCUS", "PRODUCT", "ID"), engine='python', quoting=3) #abrir arquivos e definir os nomes das colunas
        df_data['GENE'] = df_data['GENE'].astype(str).str.split('_').str[0] # Remove as terminações '_2' em genes que tenham o mesmo nome exceto estas terminações por indicar mutações
        SAT[df_name]=df_data
    return(SAT, list_SEQ)


(SAT, list_SEQ) = lerarquivos_SAT(folder_output)# chamar a funcao e definir os valores sequencias(Dicionario) e list_Seq (lista dos arquivos dentro da pasta)

def lerarquivos_Seq(folder_output):
    Seq = {}
    if folder_output[:-1] != "/":
        folder_output = str(folder_output) + "/"
    list_SEQ = os.listdir(folder_output) # lista com os nomes dos arquivos dentro da pasta "folder_output2"

    for i in list_SEQ:
        seq_name = i # retira  a terminacao do nome ".fasta"

        seq_data = list(SeqIO.parse(folder_output + str(i) + "/" + str(i) +  ".fasta"  ,"fasta"))
        Seq[seq_name]=seq_data
    return Seq

Seq = lerarquivos_Seq(folder_output)

def classificargrupos(file_output, list_SEQ):
    classes = pd.read_csv(file_output , delimiter = '\t', names= ("genomas", "grupos"), engine='python', quoting=3)

    return classes

classes = classificargrupos(file_output, list_SEQ)

def select_and_filter_SSRs(SAT, SSR_sig, GENE_sig):
    
    SAT_sig={}
    for i in SAT.keys():
        df_name = i
        df = SAT[i]
        df2 = df.loc[df["GENE"] == GENE_sig]
        df3 = df2.loc[df2["SSR"] == SSR_sig]
        SAT_sig[df_name]=df3

    
    return (SAT_sig)

(SAT_sig)= select_and_filter_SSRs(SAT, SSR_sig, GENE_sig)

intervalo = int(input("SIZE OF FLANKING SEQUENCE\n"))

def select_and_save_seq(SAT_sig, Seq, classes):
    Seq1= dict.copy(Seq)
    Seq_for_primer3 =[]
    for i in Seq.keys():
        grupo= ''
        genomas_SAT_sig = list(SAT_sig.keys())
        for z in range(len(classes["genomas"])):
            if classes["genomas"][z] == str(i):
                grupo = classes["grupos"][z]
            else:
                continue


        if i not in genomas_SAT_sig:
            continue
        else:

            
            if len(list(SAT_sig[i])) == 0:  # Caso o SAT_sig retorne nenhum resultado
                continue     # "continue" para buscar todas as sequencias. Se quiser testar o script use "break" para agilizar
                
            else:
                ID = list(SAT_sig[i]['ID'])
                for k in range(len(ID)): # Pega todos IDs do arquivo e separa em unIDades (para 1 genoma)
                    string2 = ID[k]
                    string = string2.split(" ")[0] # Divide as descrições das sequencias de acordo com o espaço " " e pega o primeiro valor (assumo que o primeiro é o seqid 
                    # print(string)
                    
                    for g in range(len(Seq1[i])): # Pega os indexes para cada SeqRecord() dentro do arquivo .fasta
   
                        if str(Seq1[i][g].id) == str(string):
                            # print(Seq1[i][g].id)
                                
                            tmp = Seq1[i][g]
                            START = pd.to_numeric(SAT_sig[i]['START']).to_list()  # Pega todos START do arquivo SAT_sig e cria uma lista (para 1 genoma)
                            END = pd.to_numeric(SAT_sig[i]['END']).to_list()     # Pega todos END do arquivo SAT_sig e cria uma lista (para 1 genoma)
                            END = END

                            START1 = int(START[k])
                            if (START1- intervalo) < 0:
                                START1 = 0
                            else:
                                START1 = int(START[k]) - intervalo

                            END1 = int(END[k])
                            if (END1 + intervalo) < len(Seq[i][g].seq):
                                END1 =  END1 + intervalo
                            else:
                                END1 = int(len(Seq[i][g].seq))

                            
                            tmp1 = SeqRecord(id = str(tmp.id) + "_" + str(k +1), seq = tmp.seq[START1:END1], description = str(i) + "  " + str(grupo) + "  " , name= tmp.name, )
                            Seq_for_primer3.append(tmp1)
                            del(tmp1)
                            del(tmp)
                        continue
                    else:
                        print("SSR found, but genome file " + str(i) +" not found!Skipped!")
                        continue
                        
                        
    return Seq_for_primer3

nome_editado = str(SSR_sig.replace(")", "_").replace("(", "")) + "_" + str(GENE_sig) + "_contigs_for_SSR.fasta"  
Seq_for_primer3 = select_and_save_seq(SAT_sig, Seq, classes)
SeqIO.write(Seq_for_primer3, nome_editado, "fasta")


print("File edited and saved as ", str(nome_editado))
