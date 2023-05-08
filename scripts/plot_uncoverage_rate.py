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

uncoverage_rate = pd.read_csv("results/coverage/per_chromosome_coverage/uncoverage_rate.txt", header = None, names = ['code', 'uncoverage_rate'], sep = ' ',
                             dtype = {
        'code': 'string',
        'uncoverage_rate': 'float'
    })

if uncoverage_rate.shape[0] < 30:
    plt.figure(figsize=(5,4))
else:
    plt.figure(figsize=(10,4))
plt.bar(uncoverage_rate['code'], uncoverage_rate['uncoverage_rate']*100)
plt.xticks(rotation=45)
plt.xlabel('Samples')
plt.ylabel('Uncoverage rate (%)')
plt.title('Uncoverage rate')

plt.savefig('graphs/uncoverage_rate.png', bbox_inches = "tight", dpi=300)
plt.show()