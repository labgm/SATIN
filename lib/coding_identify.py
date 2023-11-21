# SCRIPT QUE COMPARA OUTPUT.SATIN COM .PTT PARA DEFINIR MICROSATS QUE ESTÃO EM REGIOES CODIFICANTES (Gerando arquivo 'output.coding')
# OBS¹: AS PRIMEIRAS 4 LINHAS 'OUTPUT.SATIN' SERÁ IGNORADO (POIS SUBTENDE-SE QUE É HEADER)
# OBS²: A PRIMEIRA LINHA DO ARQUIVO '.ptt' também será ignorando pois será tratado como HEADER

import sys
import csv
import re
from Bio.Seq import Seq


def analyse_satin_ptt(input_ppt_file: str, input_satin_file: str, output_coding_file: str):

	with open(input_satin_file, "r", encoding="utf8") as file:

		line_satin_file = csv.reader(file, delimiter="\t")

		# Ignorando as 4 primeiras linhas por ser header
		next(line_satin_file); next(line_satin_file); next(line_satin_file); next(line_satin_file)

		with open(output_coding_file, "w") as handle:
		
			handle.write("#;Start;End;SSR;Gene;Strand;locus_tag;Product;ID\n")

			index_output = 0
			for row in line_satin_file:
				
				if not row: continue
				
				ssr, size, start, end, index, id = row

				start = int(start); end = int(end)

				# Lê arquivo ptt
				cont_header = 0
				file_i = open(input_ppt_file, 'r')

				for line_ptt_file in file_i:

					# IGNORANDO HEADER
					if cont_header == 0:
						cont_header += 1
						continue

					line_ptt_file = line_ptt_file.strip().split(";")
					
					start_ptt = int(line_ptt_file[0])
					end_ptt = int(line_ptt_file[1])
					

					# SE ESTÁ DENTRO DA FAIXA
					if start >= start_ptt and end <= end_ptt:
						index_output += 1
						gene = line_ptt_file[4].strip()
						strand = line_ptt_file[2].strip()
						synon = line_ptt_file[5].strip()
						prod = line_ptt_file[6].strip()
						ID = id.strip()
						if str(strand) == "-": #Caso a sequencia esteja em um strand negativo o inverso complementar do mesmo sera chamado
							only_alpha = ""
							for j in ssr: # lista de cada genoma (j é o SSR analisado)
								
								for m in j:
									if ord(m) >=65 and ord(m)<=90:
										only_alpha += m
									elif ord(m) >= 97 and ord(m) <=122:
										only_alpha += m
							number = re.findall(r'\d+', ssr)[0]
							motif = Seq(only_alpha[::-1])
							only_alpha= ""
							motif = "(" + str(motif.complement()) + ")" + str(number)
							ssr=motif
							
							handle.write(f"{index_output}.;{start};{end};{ssr};{gene};{strand};{synon};{prod};{ID}\n")
						else:
							handle.write(f"{index_output}.;{start};{end};{ssr};{gene};{strand};{synon};{prod};{ID}\n")



if __name__  == '__main__':

	
	input_ppt_file = sys.argv[1]
	input_satin_file = sys.argv[2]
	output_coding_file = sys.argv[3]
	
	analyse_satin_ptt(input_ppt_file, input_satin_file, output_coding_file)
