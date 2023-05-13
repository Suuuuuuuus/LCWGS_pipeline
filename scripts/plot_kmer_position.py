import io
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

code = sys.argv[1]

df1 = pd.read_csv("results/kmer/" + code + "/read1/" + code + "_position1.tsv.gz", sep='\t', compression='gzip')
df1['number_of_errors'] = pd.to_numeric(df1['number_of_errors'])
df1['number_of_reads'] = pd.to_numeric(df1['number_of_reads'])
df1['A'] = pd.to_numeric(df1['A'])
df1['G'] = pd.to_numeric(df1['G'])
df1['C'] = pd.to_numeric(df1['C'])
df1['T'] = pd.to_numeric(df1['T'])
df1['error_rate'] = df1['number_of_errors']/df1['number_of_reads']
df1['percentage A'] = df1['A']/df1['number_of_reads']
df1['percentage C'] = df1['C']/df1['number_of_reads']
df1['percentage G'] = df1['G']/df1['number_of_reads']
df1['percentage T'] = df1['T']/df1['number_of_reads']
df1 = df1.drop(columns = ['start_or_end', 'sum_of_base_qualities', 'number_of_reads', 'number_of_errors', 'A', 'G', 'C', 'T'])

df2 = pd.read_csv("results/kmer/" + code + "/read2/" + code + "_position2.tsv.gz", sep='\t', compression='gzip')
df2['number_of_errors'] = pd.to_numeric(df2['number_of_errors'])
df2['number_of_reads'] = pd.to_numeric(df2['number_of_reads'])
df2['A'] = pd.to_numeric(df2['A'])
df2['G'] = pd.to_numeric(df2['G'])
df2['C'] = pd.to_numeric(df2['C'])
df2['T'] = pd.to_numeric(df2['T'])
df2['error_rate'] = df2['number_of_errors']/df2['number_of_reads']
df2['percentage A'] = df2['A']/df2['number_of_reads']
df2['percentage C'] = df2['C']/df2['number_of_reads']
df2['percentage G'] = df2['G']/df2['number_of_reads']
df2['percentage T'] = df2['T']/df2['number_of_reads']
df2 = df2.drop(columns = ['start_or_end', 'sum_of_base_qualities', 'number_of_reads', 'number_of_errors', 'A', 'G', 'C', 'T'])

fig = plt.figure(figsize=(6,15))
ax1 = plt.subplot(3, 1, 1)
ax1.plot(df1['position'][:120], df1['error_rate'][:120]*100, label = 'read1')
ax1.plot(df2['position'][:120], df2['error_rate'][:120]*100, label = 'read2')
ax1.set_xlabel('Index on reads')
ax1.set_ylabel('Error rate (%)')
ax1.set_title('Error rate at each index')
ax1.set_xlim((-5, 155))
ax1.set_ylim((0, 30))
ax1.legend()

ax2 = plt.subplot(3, 1, 2)
ax2.plot(df1['position'], df1['percentage A']*100, label = 'A')
ax2.plot(df1['position'], df1['percentage C']*100, label = 'C')
ax2.plot(df1['position'], df1['percentage G']*100, label = 'G')
ax2.plot(df1['position'], df1['percentage T']*100, label = 'T')
ax2.set_xlabel('Index on reads')
ax2.set_ylabel('Base composition along read (%)')
ax2.set_title('A, T, C, G contents, read 1')
ax2.legend()
ax2.set_xlim((-5, 155))
ax2.set_ylim((10, 45))

ax3 = plt.subplot(3, 1, 3)
ax3.plot(df2['position'], df2['percentage A']*100, label = 'A')
ax3.plot(df2['position'], df2['percentage C']*100, label = 'C')
ax3.plot(df2['position'], df2['percentage G']*100, label = 'G')
ax3.plot(df2['position'], df2['percentage T']*100, label = 'T')
ax3.set_xlabel('Index on reads')
ax3.set_ylabel('Base composition along read (%)')
ax3.set_title('A, T, C, G contents, read 2')
ax3.legend()
ax3.set_xlim((-5, 155))
ax3.set_ylim((10, 45))

plt.savefig('graphs/kmer_position/' + code + '_kmer_position.png', bbox_inches = "tight", dpi=300)
