# scripts/merge_counts.py
import pandas as pd
from snakemake import input, output

# Merge all count files
counts = []
for f in input:
    sample = f.split('/')[-1].replace('_counts.txt', '')
    df = pd.read_csv(f, sep='\t', comment='#', index_col=0)  # featureCounts output
    counts.append(df.iloc[:, -1].rename(sample))

pd.concat(counts, axis=1).to_csv(output[0])
