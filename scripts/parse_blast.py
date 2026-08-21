#!/usr/bin/env python3
import argparse
import pandas as pd
import os

def parse_and_summarize(blast_files, out_path, ltr_range=500):
    all_hits = []
    
    # Columns we expect from our blastn command:
    # 6 qseqid sseqid sstart send sstrand pident qcovs
    colnames = ["qseqid", "sseqid", "sstart", "send", "sstrand", "pident", "qcovs"]
    
    for bf in blast_files:
        if not os.path.exists(bf) or os.path.getsize(bf) == 0:
            continue
            
        # Filename format expected: {genome}_{gbk}_{target}.tsv
        basename = os.path.basename(bf).replace('.tsv', '')
        parts = basename.split('_')
        
        # If the genome or gbk names have underscores, this might be tricky,
        # but we can try to guess or just use the filename as source.
        source_file = basename
        
        try:
            df = pd.read_csv(bf, sep='\t', header=None, names=colnames)
            df['source_file'] = source_file
            
            # Apply LTR range expansion (context)
            # Adjusting start and end to include the context margin
            df['adjusted_start'] = df.apply(
                lambda row: max(1, row['sstart'] - ltr_range) if row['sstrand'] == 'plus' else max(1, row['send'] - ltr_range), 
                axis=1
            )
            df['adjusted_end'] = df.apply(
                lambda row: row['send'] + ltr_range if row['sstrand'] == 'plus' else row['sstart'] + ltr_range, 
                axis=1
            )
            
            all_hits.append(df)
        except pd.errors.EmptyDataError:
            pass # File is empty, skip

    if all_hits:
        final_df = pd.concat(all_hits, ignore_index=True)
        # Sort by subject sequence and start position
        final_df = final_df.sort_values(by=["sseqid", "adjusted_start"])
        final_df.to_csv(out_path, sep='\t', index=False)
        print(f"Summarized {len(final_df)} hits into {out_path}.")
    else:
        # Create empty summary if no hits at all
        with open(out_path, 'w') as f:
            f.write("\t".join(colnames + ["source_file", "adjusted_start", "adjusted_end"]) + "\n")
        print(f"No hits found across provided BLAST files. Created empty summary at {out_path}.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Summarize BLAST results.")
    parser.add_argument("--blast_files", nargs='+', required=True, help="List of BLAST output TSV files")
    parser.add_argument("--out", required=True, help="Output summary TSV file")
    parser.add_argument("--ltr_range", type=int, default=500, help="Margin to add to hits (default: 500)")
    
    args = parser.parse_args()
    parse_and_summarize(args.blast_files, args.out, args.ltr_range)
