from controller.fasta import process_fasta
from controller.gbk import process_gbk
import os
from Bio import SeqIO

class MainController:

    def __init__(self, input_file_param: str):
        self.param, self.max_int = self.init_params(input_file_param)


    def init_params(self, input_file_param):

        PARAM, INT = '', ''
        try:
            with open(input_file_param, 'r') as arquivo:
                for linha in arquivo:
                    linha_ = linha.strip()
                    if linha_.split(' ')[0] == 'PARAM':
                        PARAM = str(linha.replace('PARAM ', '')).strip()
                    if linha_.split(' ')[0] == 'INT':
                        INT = str(linha.replace('INT ', '')).strip()
        except FileNotFoundError:
            print(f"O arquivo '{input_file_param}' não foi encontrado.")
        except Exception as e:
            print(f"Ocorreu um erro: {e}")

        if not PARAM or not INT:
            raise Exception(f"[ERROR] Parameters.")

        return PARAM, INT
    

    def read_fasta(self, input_fasta, output_path):
        process_fasta(input_fasta, output_path, self.param, self.max_int)


    def read_fastas(self, input_fasta_s, output_path):
        fasta_extensions = [".fasta", ".fa", ".fna"]
        for file in os.listdir(input_fasta_s):
            # Verifica se a extensão do arquivo está na lista de extensões FASTA
            if any(file.lower().endswith(ext) for ext in fasta_extensions):
                # Lê as sequências do arquivo FASTA
                local_file = os.path.join(input_fasta_s, file)
                name_fasta = local_file.rstrip('/').split('/')[-1]

                count_reads = len(list(SeqIO.parse(local_file, "fasta")))
                if count_reads < 2:
                    print(f"Fasta File: {name_fasta} ..")
                else:
                    print(f"MultiFasta File: {name_fasta} | ({count_reads}) seqs..")

                process_fasta(local_file, output_path, self.param, self.max_int)


    def read_gbk(self, input_gbk, output_path):
        process_gbk(input_gbk, output_path, self.param, self.max_int)


    def read_gbks(self, input_gbk_s, output_path):
        gbk_extensions = [".gbk", ".gb", ".gbff"]
        for file in os.listdir(input_gbk_s):
            if any(file.lower().endswith(ext) for ext in gbk_extensions):
                print(f"\nReading GBK file {file}")
                local_file = os.path.join(input_gbk_s, file)
                process_gbk(local_file, output_path, self.param, self.max_int)

