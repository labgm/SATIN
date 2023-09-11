import os
import subprocess
from lib.extract_fasta_from_gbk import extract_fasta
from lib.GBKtoPTT import convert_gbktoptt
from lib.coding_identify import analyse_satin_ptt

def process_gbk(input_gbk: str, output_path: str, config_params: str, max_distance_compound: str):
    """
        Processa arquivo GBK

        Args:
            - input_gbk: Local GBK file 
            - output_path: Local output
            - config_params: Configuration parameters in format -> Ex: 3:2 4:2 5:2
            - max_distance_compound: Maximum distance between SSRs to consider compound

        Returns:
            - Output Fasta file
            - Output PPT file
            - Output with result search SSRs
            - Output coding file
    """

    name = input_gbk.rstrip('/').split('/')[-1]
    local_output_new = f"{output_path}/{name}"
    try: os.mkdir(local_output_new)
    except: pass
    
    # Extract Fasta
    output_fasta = f"{local_output_new}/{name}.fasta"
    extract_fasta(input_gbk, output_fasta)

    # Gen PTT
    output_ptt_file = f"{local_output_new}/{name}.ptt"
    convert_gbktoptt(input_gbk, output_ptt_file)

    # Search SSRs
    output_satin_file = f"{local_output_new}/{name}.satin"
    execution_command = [f"./bin/search_ssr_mf"] + [output_fasta, output_satin_file, config_params, max_distance_compound]
    subprocess.run(execution_command)

    # Gen Coding
    output_coding_file = f"{local_output_new}/{name}.coding"
    analyse_satin_ptt(output_ptt_file, output_satin_file, output_coding_file)
