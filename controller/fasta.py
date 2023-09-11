import subprocess

def process_fasta(input_fasta: str, output_path: str, config_params: str, max_distance_compound: str):
    """Search SSRs from Fasta file
    
        Args:
            - input_fasta: Input Fasta file
            - output_path: Output path
            - config_params: Configuration parameters in format -> Ex: 3:2 4:2 5:2
            - max_distance_compound: Maximum distance between SSRs to consider compound

        Returns:
            - Output with result search SSRs
    """

    name_fasta = input_fasta.rstrip('/').split('/')[-1]
    output_satin_file = f"{output_path}/{name_fasta}.satin"

    execution_command = [f"./bin/search_ssr_mf"] + [input_fasta, output_satin_file, config_params, max_distance_compound]
    subprocess.run(execution_command)