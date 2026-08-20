#!/usr/bin/env python3
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def extract_from_gbk(gbk_path, target, out_path):
    records = list(SeqIO.parse(gbk_path, "genbank"))
    out_records = []
    
    for idx, record in enumerate(records):
        if target.lower() == "complete":
            new_id = f"{record.id}_complete"
            out_records.append(SeqRecord(record.seq, id=new_id, description="Complete Genome"))
            continue
            
        for i, feature in enumerate(record.features):
            match = False
            search_terms = [target.lower()]
            if target.lower() == "ltr":
                search_terms.append("long_terminal_repeat")
                
            # Check feature type
            if any(term in feature.type.lower() for term in search_terms):
                match = True
            
            # Check qualifiers (gene, product, note, label, rpt_type)
            if not match:
                for qual_key in ['gene', 'product', 'note', 'label', 'rpt_type']:
                    if qual_key in feature.qualifiers:
                        qual_values = feature.qualifiers[qual_key]
                        if any(term in val.lower() for val in qual_values for term in search_terms):
                            match = True
                            break
                            
            if match:
                seq = feature.extract(record.seq)
                new_id = f"{record.id}_{target}_{i}"
                desc = f"Extracted {target} from {record.id}"
                out_records.append(SeqRecord(seq, id=new_id, description=desc))
                
    if not out_records:
        # Create an empty file to avoid snakemake complaining if no target found
        open(out_path, 'w').close()
        print(f"Warning: No features found for target '{target}' in {gbk_path}. Created empty file.")
    else:
        SeqIO.write(out_records, out_path, "fasta")
        print(f"Extracted {len(out_records)} sequences for target '{target}' to {out_path}.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract specific sequences from a GBK file.")
    parser.add_argument("--gbk", required=True, help="Input GenBank file")
    parser.add_argument("--target", required=True, help="Target feature (e.g., LTR, gag, complete)")
    parser.add_argument("--out", required=True, help="Output FASTA file")
    
    args = parser.parse_args()
    extract_from_gbk(args.gbk, args.target, args.out)
