import pandas as pd
import io
import gzip

def read_vcf(path):
    with gzip.open(path, 'rt') as f:
        lines = [l for l in f if not l.startswith('##')]
    df = pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})
    
    def extract_ad(sample):
        sample_fields = sample.split(':')
        return sample_fields[1]

    def extract_dp(sample):
        sample_fields = sample.split(':')
        return sample_fields[2]

    for sample in df.columns[9:]:
        df[f'AD_{sample}'] = df.apply(lambda row: extract_ad(row[sample]), axis=1)
        df[f'DP_{sample}'] = df.apply(lambda row: extract_dp(row[sample]), axis=1)
    
    return df

def compute_alt_allele_frequencies(df, samples):
    for sample in samples:
        df[[f"REF_{sample}", f"ALT_{sample}"]] = df[f"AD_{sample}"].str.split(',', expand=True).astype(float)
        df[f"AF_{sample}"] = df[f"ALT_{sample}"] / (df[f"REF_{sample}"] + df[f"ALT_{sample}"])
    return df

def output_per_sample(df, samples, output_csv=False):
    sample_data = {}
    for sample in samples:
        sample_df = df[['CHROM', 'POS', f'AF_{sample}', f'DP_{sample}']].copy()
        sample_data[sample] = sample_df
        if output_csv:
            sample_df.to_csv(f"{sample}_summary.csv", index=False)
    return sample_data

def main():
    path = '/dfs6/pub/nshiroon/EE283/DATA/DNAseq/output3.vcf.gz'
    samples = ["ADL06", "ADL10", "ADL14", "ALD09"]
    
    df = read_vcf(path)
    df = compute_alt_allele_frequencies(df, samples)
    sample_data = output_per_sample(df, samples, output_csv=True)
    
    print("Output complete. Summary files generated for each sample.")
    return sample_data

if __name__ == "__main__":
    main()

