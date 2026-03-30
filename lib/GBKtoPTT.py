#!/usr/bin/env python3
"""
Convert GenBank file to PTT (Protein Table) format
Extracts CDS features with Gene and Product information
"""

import sys
import warnings
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

def convert_gbktoptt(gbk_file, output_file):
    """
    Extract CDS features from GenBank file and write to PTT format
    
    PTT Format: 
    #Start;End;Strand;Length;Gene;locus_tag;Product
    """
    
    try:
        # Parse the GenBank file
        records = list(SeqIO.parse(gbk_file, "genbank"))
        
        if not records:
            print(f"Warning: No records found in {gbk_file}")
            return 0
        
        record = records[0]  # Use first record
        cds_features = []
        
        # Extract all CDS features
        for feature in record.features:
            if feature.type == "CDS":
                # Get location information
                start = int(feature.location.start) + 1  # 1-based coordinates
                end = int(feature.location.end)
                strand = "+" if feature.location.strand >= 0 else "-"
                length = end - start + 1
                
                # Extract qualifiers
                gene = feature.qualifiers.get("gene", [""])[0]
                locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
                product = feature.qualifiers.get("product", [""])[0]
                
                # Handle missing values
                gene = gene if gene else ""
                product = product if product else ""
                locus_tag = locus_tag if locus_tag else ""
                
                cds_features.append({
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "length": length,
                    "gene": gene,
                    "locus_tag": locus_tag,
                    "product": product
                })
        
        # Sort by start position
        cds_features.sort(key=lambda x: x["start"])
        
        # Write to PTT file
        with open(output_file, 'w') as f:
            # Write header
            f.write("#Start;End;Strand;Length;Gene;locus_tag;Product\n")
            
            # Write each CDS
            for cds in cds_features:
                line = f"{cds['start']};{cds['end']};{cds['strand']};{cds['length']};{cds['gene']};{cds['locus_tag']};{cds['product']}\n"
                f.write(line)
        
        print(f"Successfully extracted {len(cds_features)} CDS features to {output_file}")
        return len(cds_features)
        
    except ValueError as e:
        if "Premature end of line" in str(e):
            print(f"ERROR: The GenBank file '{gbk_file}' is corrupted or incomplete.")
            print("Please redownload the file using the instructions above.")
        else:
            print(f"ERROR: {e}")
        return 0
    except Exception as e:
        print(f"ERROR: Unexpected error: {e}")
        return 0

def main():
    if len(sys.argv) != 3:
        print("Usage: python GBKtoPTT.py <input.gbk> <output.ptt>")
        print("Example: python GBKtoPTT.py GEN/GCA_001951615.1.gbff output/GCA_001951615.1.ptt")
        sys.exit(1)
    
    gbk_file = sys.argv[1]
    output_file = sys.argv[2]
    
    # Create output directory if it doesn't exist
    import os
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Run extraction
    count = convert_gbktoptt(gbk_file, output_file)
    
    if count > 0:
        print(f"✅ Successfully created PTT file with {count} features")
        sys.exit(0)
    else:
        print(f"❌ Failed to create PTT file")
        sys.exit(1)

if __name__ == "__main__":
    main()