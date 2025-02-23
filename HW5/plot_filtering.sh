srun -A CLASS-ECOEVO283 -c 1 --pty /bin/bash -i

python

import subprocess
import matplotlib.pyplot as plt

# Define VCF file paths
vcf_files = {
    "Raw SNPs": "all_variants.vcf.gz",
    "Step 1": "output1.vcf.gz",
    "Step 2": "output2.vcf.gz",
    "Final SNPs": "output3.vcf.gz"
}

# Function to count SNPs in a VCF file
def count_snps(vcf):
    cmd = f"bcftools view -v snps {vcf} | grep -vc '^#'"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return int(result.stdout.strip()) if result.stdout.strip().isdigit() else 0

# Count SNPs at each stage
snp_counts = {step: count_snps(vcf) for step, vcf in vcf_files.items()}

# Plot SNP loss
plt.figure(figsize=(7,5))
plt.bar(snp_counts.keys(), snp_counts.values(), color=['blue', 'orange', 'red', 'green'])
plt.ylabel("Number of SNPs")
plt.xlabel("Filtering Step")
plt.title("SNP Retention at Each Filtering Step")
plt.xticks(rotation=20)
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.show()

plt.savefig("snp_filtering_plot.png", dpi=300)
