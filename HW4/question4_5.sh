salloc -A ecoevo283 --ntasks=2 srun --pty /bin/bash -i

mamba activate deeptools

bamCoverage -b A4_filtered.bam -o A4_frag.bedgraph --binSize 1 --normalizeUsing RPKM --extendReads --outFileFormat bedgraph
bamCoverage -b A5_filtered.bam -o A5_frag.bedgraph --binSize 1 --normalizeUsing RPKM --extendReads --outFileFormat bedgraph

python

import pandas as pd
import matplotlib.pyplot as plt

df_A4 = pd.read_csv("A4_frag.bedgraph", sep="\t", header=None, names=["chr", "start", "end", "signal"])
df_A5 = pd.read_csv("A5_frag.bedgraph", sep="\t", header=None, names=["chr", "start", "end", "signal"])

region_start, region_end = 1904000, 1904500

df_A4_roi = df_A4[(df_A4["start"] >= region_start) & (df_A4["end"] <= region_end)]
df_A5_roi = df_A5[(df_A5["start"] >= region_start) & (df_A5["end"] <= region_end)]

plt.figure(figsize=(10, 5))
plt.step(df_A4_roi["start"], df_A4_roi["signal"], label="A4 Extended Reads", color="blue", where="post")
plt.step(df_A5_roi["start"], df_A5_roi["signal"], label="A5 Extended Reads", color="red", where="post")

plt.xlabel("Genomic Position on chrX")
plt.ylabel("Coverage (RPKM)")
plt.title("Coverage near chrX:1,904,042 (Extended Reads)")
plt.legend()

plt.xlim(region_start, region_end)

plt.savefig("coverage_profile_extended.png", dpi=300, bbox_inches='tight')
plt.show()

