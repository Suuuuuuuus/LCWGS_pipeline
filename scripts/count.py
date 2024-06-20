import sys
sys.path.append('/well/band/users/rbx225/software/lcwgsus/')
import lcwgsus
from lcwgsus.variables import *
import pandas as pd

name = sys.argv[1]

dfs = [name + ".chr" + str(i) + ".vcf.gz" for i in CHROMOSOMES_ALL]
lcwgsus.get_n_variants_vcf(dfs)
