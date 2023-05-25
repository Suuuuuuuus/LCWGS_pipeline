import pandas as pd
import sys

kmer1 = pd.read_csv('results/kmer/kmer_accuracy_read1.txt', sep = '\t', 
                   names = ['Code', 'Kmer_Error_Rate_Read1'], dtype = {'Code': 'string', 'Kmer_Error_Rate_Read1': 'float'})
kmer2 = pd.read_csv('results/kmer/kmer_accuracy_read2.txt', sep = '\t', 
                   names = ['Code', 'Kmer_Error_Rate_Read2'], dtype = {'Code': 'string', 'Kmer_Error_Rate_Read2': 'float'})
fastqc = pd.read_csv('results/fastqc/duplication_rate_fastqc.txt', sep = '\t', 
                   names = ['Code', 'Fastqc_Dup_Rate'], dtype = {'Code': 'string', 'Fastqc_Dup_Rate': 'float'})
samtools = pd.read_csv('results/dup_rate/duplication_rate_samtools.txt', sep = '\t', 
                   names = ['Code', 'Samtools_Dup_Rate'], dtype = {'Code': 'string', 'Samtools_Dup_Rate': 'float'})
coverage = pd.read_csv('results/coverage/per_sample_coverage.txt', sep = ' ', 
                   names = ['Code', 'Coverage'], dtype = {'Code': 'string', 'Coverage': 'float'})
uncoverage = pd.read_csv('results/coverage/per_chromosome_coverage/uncoverage_rate.txt', sep = ' ', 
                   names = ['Code', 'Uncoverage'], dtype = {'Code': 'string', 'Uncoverage': 'float'})
fragment_size = pd.read_csv('results/fragment_size/fragment_size.txt', sep = '\t', 
                   names = ['Code', 'Fragment_Size'], dtype = {'Code': 'string', 'Fragment_Size': 'float'})
whole = pd.merge(kmer1, kmer2, on='Code', how="left")
whole = pd.merge(whole, fastqc, on='Code', how="left")
whole = pd.merge(whole, samtools, on='Code', how="left")
whole = pd.merge(whole, coverage, on='Code', how="left")
whole = pd.merge(whole, uncoverage, on='Code', how="left")
whole = pd.merge(whole, fragment_size, on='Code', how="left")
whole = whole.sort_values(by=['Kmer_Error_Rate_Read1'])
whole.to_csv('results/lcwgs_results.csv', index=False)

