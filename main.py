import logging
import argparse
import sys
import os
import time
from controller.main import MainController

def parse_args(args):
    """Parse command line parameters

    Args:
      args ([str]): command line parameters as list of strings

    Returns:
      :obj:`argparse.Namespace`: command line parameters namespace
    """
    parser = argparse.ArgumentParser(
        description="Kafka2spark is a pyspark module which extract a offset messages from a Kafka topic and save them in a output path.")
    parser.add_argument(
        "--input-fasta", # 1 file fasta
        "-f",
        dest="input_fasta",
        required=False,
        help="Input fasta file",
        metavar="STRING")
    parser.add_argument(
        "--input-multi-fasta", # 1 file multi fasta
        "-mf",
        dest="input_multi_fasta",
        required=False,
        help="Input multi-fasta file",
        metavar="STRING")
    parser.add_argument(
        "--input-fastas", # Pasta com arquivos fasta ou multifasta
        "-fs",
        dest="input_fasta_s",
        required=False,
        help="Input folder fasta/multifasta files",
        metavar="STRING")
    parser.add_argument(
        "--input-gbk",
        "-g",
        dest="input_gbk",
        required=False,
        help="Input GBK file",
        metavar="STRING")
    parser.add_argument(
        "--input-gbks",
        "-gs",
        dest="input_gbk_s",
        required=False,
        help="Input folder GBK files",
        metavar="STRING")
    parser.add_argument(
        "--output",
        "-o",
        dest="output_path",
        required=True,
        help="Input folder to storage output",
        metavar="STRING")
    parser.add_argument(
        "-v",
        "--verbose",
        dest="log_level",
        help="Set loglevel to INFO",
        action="store_const",
        const=logging.INFO)
    parser.add_argument(
        "-vv",
        "--very-verbose",
        dest="log_level",
        help="Set loglevel to DEBUG",
        action="store_const",
        const=logging.DEBUG)

    return parser.parse_known_args(args)

def main(args):

    controller = MainController('Parameters.ini')

    init = time.time()

    args, unknown = parse_args(args)

    # Configuração básica do logging
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

    logging.info("###################################")
    logging.info("###################################")
    logging.info("Unknow Parameters")
    logging.info("{}".format(unknown))
    logging.info("###################################")
    logging.info("###################################")

    input_fasta = args.input_fasta
    input_multi_fasta = args.input_multi_fasta
    input_fasta_s = args.input_fasta_s
    input_gbk = args.input_gbk
    input_gbk_s = args.input_gbk_s
    output_path = args.output_path.rstrip('/')

    try: os.mkdir(f"{output_path}")
    except: pass


    if input_fasta: # -f          <file>
        controller.read_fasta(input_fasta, output_path)

    elif input_multi_fasta: # -mf <file>
        controller.read_fasta(input_multi_fasta, output_path)

    elif input_fasta_s: # -fs     <folder>
        controller.read_fastas(input_fasta_s, output_path)

    elif input_gbk: # -g         <file>
        controller.read_gbk(input_gbk, output_path)

    elif input_gbk_s: # -gs     <folder>
        controller.read_gbks(input_gbk_s, output_path)


    end = time.time()
    logging.info(f"Total execution: {(end-init)/60} minutes.")

if __name__ == '__main__':

    args = sys.argv[1:]
    main(args)