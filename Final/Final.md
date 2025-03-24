# Writing a python program for important VCF calculations

### Finding a way to read VCF files in python and creating a mamba environment 

Through investigation I found two potential ways to read the VCF files and the necessary columns to do the calculations. The first way is through a package in Python called Pysam and the second is through a function somebody had created to read a VCF file as a pandas data frame. Given these options it is far easier to work with a pandas data frame. The following code is shown below for reading a vcf file.

```{python}
def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})
```

The function seems to read the vcf file line by line and filter out any metadata and keeping the column headers. The lines are then read in as a pandas csv which can the be reformatted and return a pandas data frame with the necessary information to perform the calculations in the problem.

Now that we have a way to read the vcf file the next step is to create a mamba environment with the correct packages to complete the task.

```{bash}
#Claim computing node (Interactive)
srun -A CLASS-ECOEVO283 -c 1 --pty /bin/bash -i

#Create mamba enviornment with python
mamba create -n vcf_analysis python

#Activate enviornment
mamba activate vcf_analysis

#Install necessary packages with pip 
pip install io os pandas matplotlib seaborn jupyter 
```

### Launching Jupyter and loading VCF

To allow for easy visualization and debugging within the HPC3, I employ Jupyter Notebook. Jupyter needs to be port-forwarded to the local machine as Jupyter is running on the cloud. To perform this we can use the following steps:

```{bash}
#Within the HPC3 create this easy shortcut to start an instance of jupyter with a random port
nano ~/.bashrc

#Create the function to launch jupyter notebook with a random acceptable port
jupyter_random() {
    local port=$((RANDOM % 3000 + 7000))
    jupyter notebook --no-browser --ip="$(hostname -s)" --port="$port"
}
#Create alias
alias nb=jupyter_random

#Launch alias
nb

#Should now see something like this
#http://127.0.0.1:7309/tree?token=57cbd012417142ac0b159707b051b3ffaea67437032074fe
```

Now that we have launched an instance of jupyter with port 7309 we need to connect the local computer to this port, we can do this with the following code locally

```{bash}
#Run from local command line, ensure computing node and generated port are accurate, 9090 represents the local port we will connect to
ssh -L 9090:hpc3-21-15:7309 nshiroon@hpc3.rcic.uci.edu
```

From here we can simply now launch [http://127.0.0.1:](http://127.0.0.1:7309){.uri}[9090](http://127.0.0.1:9090) to connect to an instance of jupyter allowing for interactive coding and visualization which will be important in terms of visualizing data frames, debugging, and plots.

The first step is to load all packages needed

```{python}
import pandas as pd
import io
import os
import matplotlib.pyplot as plt
import seaborn
```

Next I load in the function and define my vcf path and then check the data frame as a sanity check. What I noted was that the function did not work on .gz files, so to combat this and not have to have a uncompressed file in my storage, we can work with the gzip module in python to fix this.

```{python}
import gzip

def read_vcf(path):
#Added gzip to read the file
    with gzip.open(path, 'rt') as f: 
        lines = [l for l in f if not l.startswith('##')]
    
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

df = read_vcf('/dfs6/pub/nshiroon/EE283/DATA/DNAseq/output3.vcf.gz')

```

Checking the first row of the data frame we get the following

```{python}
df[:1]
CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ADL06	ADL10	ADL14	ALD09
0	2L	16483	.	T	A	1145.32	PASS	AC=4;AF=0.5;AN=8;BaseQRankSum=0.767;DP=255;Exc...	GT:AD:DP:GQ:PGT:PID:PL:PS	0|1:47,19:66:99:0|1:16479_A_AT:656,0,1901:16479	0|1:71,11:82:99:0|1:16479_A_AT:224,0,2916:16479	0|1:50,7:57:99:0|1:16479_A_AT:143,0,2070:16479	0/1:40,6:46:99:.:.:132,0,1636:.
```

Importing the VCF into python was successful and now further work can be done to write a program that can perform calculations.

### Calculating frequency of "ALT" allele 

The first important step here is that the AD (Allelic Depth) has not been properly parsed out when reading in the file. The format column within the data represents where we can find the AD depth per sample. Given this we need to further modify our read in function to account for AD depth for each sample in a few new columns. The easy fix here is to extract the second field from each sample since the format is GT:AD meaning the second item here should give us the Allelic Depth. By printing out the sample field and seeing the following: The second item here or sample_fields[1] will correspond to the AD for each sample.

```{python}
print(f"SAMPLE fields: {sample_fields}")
SAMPLE fields: ['0|1', '47,19', '66', '99', '0|1', '16479_A_AT', '656,0,1901', '16479']

```

Modifying the original function now gives us:

```{python}
def read_vcf(path):
    with gzip.open(path, 'rt') as f:
        lines = [l for l in f if not l.startswith('##')] 
    #Now read in df earlier to fix formatting in successive function
    df = pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

    #Define a new function to extract the AD per sample
    def extract_ad(sample):
        #The sample field that matches "AD" is the second based on sample_fields print out and FORMAT column in VCF
        sample_fields = sample.split(':')
        return sample_fields[1] 
    #Apply extract_ad to each sample 
    for sample in df.columns[9:]:
            df[f'AD_{sample}'] = df.apply(lambda row: extract_ad(row[sample]), axis=1)
    return df

df = read_vcf('/dfs6/pub/nshiroon/EE283/DATA/DNAseq/output3.vcf.gz')
```

Now that we have the AD information per sample in the data frame, we can now calculate the frequency of the ALT allele per sample per position with the following function:

```{python}
def compute_alt_allele_frequencies(df, samples):
    for sample in samples:
        #The two values in AD represent ref/alt, split them and then use them for frequency calculation 
        df[[f"REF_{sample}", f"ALT_{sample}"]] = df[f"AD_{sample}"].str.split(',', expand=True).astype(float)
        
        #Frequency of the ALT allele should be defined as ALT/(REF+ALT)
        df[f"AF_{sample}"] = df[f"ALT_{sample}"] / (df[f"REF_{sample}"] + df[f"ALT_{sample}"])

    return df

samples = ["ADL06", "ADL10", "ADL14", "ALD09"]
df = compute_alt_allele_frequencies(df, samples)
```

With this information we can now make a plot to view some results, in this case it might be best to visualize averages, the following plot shows the average ALT frequency per sample

```{python}
def plot_average_alt_allele_frequency(df, samples, save_path=None):
    avg_af = {f"AF_{sample}": df[f"AF_{sample}"].mean() for sample in samples}
    avg_af_df = pd.DataFrame(list(avg_af.items()), columns=['Sample', 'Average_AF'])
    
    plt.figure(figsize=(10, 6))
    sns.barplot(x='Sample', y='Average_AF', data=avg_af_df, palette='Set2')
    plt.title('Average ALT Allele Frequency per Sample')
    plt.xlabel('Sample')
    plt.ylabel('Average ALT Allele Frequency')
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path)  
    
    plt.show()

samples = ["ADL06", "ADL10", "ADL14", "ALD09"]
plot_average_alt_allele_frequency(df, samples, save_path="average_alt_allele_frequency.png")
```

The following plot is saved as "average_alt_allele_frequency.png" in the github Final Repo.

### Outputting Frequency and Coverage for each sample

The last step is to output the frequency and coverage for each sample. We can modify the current existing code to account for these and then as a final wrap-up package the code together into one nice python program for use. To get depth of coverage this is essentially the same process as allelic depth and the read-in function can be changed to index DP and save the value to a column for each row. The following code adjusts for this:

```{python}
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
    #Newly added code to the read function to read in depth of coverage and output in to a column
    def extract_dp(sample):
        sample_fields = sample.split(':')
        return sample_fields[2] 

    for sample in df.columns[9:]:
        df[f'AD_{sample}'] = df.apply(lambda row: extract_ad(row[sample]), axis=1)
        df[f'DP_{sample}'] = df.apply(lambda row: extract_dp(row[sample]), axis=1)
    
    return df

df = read_vcf('/dfs6/pub/nshiroon/EE283/DATA/DNAseq/output3.vcf.gz')

```

Now that we have coverage and also calculated frequency from before we can organize everything together per sample and output the results. Here I create a function to create a separate data frame for each sample with information of chromosome, position, frequency, and coverage as a .csv file. The following code computes this:

```{python}
def output_per_sample(df, samples, output_csv=False):
    sample_data = {}
    for sample in samples:
        sample_df = df[['CHROM', 'POS', f'AF_{sample}', f'DP_{sample}']].copy()
        sample_data[sample] = sample_df
        
        if output_csv:
            sample_df.to_csv(f"{sample}_summary.csv", index=False)
    return sample_data

sample_dfs = output_per_sample(df, samples, output_csv=True)
```

### Conclusion

Through tests and reiteration we now have a python package that can read a VCF file, create columns for the necessary information, complete a frequency calculation, then output a .csv summary for each sample. There are now two files within the Final repository, "vcf_summarize.py" and "vcf_summarize_job.sh" which include the finalized pipeline as well as the cluster job to submit for processing.
