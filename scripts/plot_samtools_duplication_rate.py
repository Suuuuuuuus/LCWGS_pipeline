import io
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import random
from collections import Counter
import csv
import glob
import matplotlib.colors as mcolors

samtools_dup_rate = pd.read_csv("results/dup_rate/duplication_rate_samtools.txt", header = None, names = ['code', 'dup_rate'], sep = '\t',
                             dtype = {
        'code': 'string',
        'dup_rate': 'float'
    })

if samtools_dup_rate.shape[0] < 30:
    plt.figure(figsize=(5,4))
else:
    plt.figure(figsize=(10,4))
plt.bar(samtools_dup_rate['code'], samtools_dup_rate['dup_rate']*100)
plt.xticks(rotation=45)
plt.xlabel('Samples')
plt.ylabel('Samtools duplication rate (%)')
plt.title('Samtools duplication rate')

plt.savefig('graphs/samtools_duplication_rate.png', bbox_inches = "tight", dpi=300)
plt.show()