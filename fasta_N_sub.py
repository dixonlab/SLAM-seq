#!/usr/bin/env python

chrNameLength_path = '~/b38_decoy/chrNameLength.txt'
fasta_path = '~/b38_decoy/hs38d5.fa'
chr_snp_dir = '~/RPE1-SNPs/SNP_by_chr' #files in this directory need to be labeled as chr1.vcf, chr2.vcf, etc. It will ignore decoy chromosomes though.

out_dir = '~/RPE1-SNPs'

import pandas as pd
import math
import subprocess
from pathlib import Path
import os

length = pd.read_csv(chrNameLength_path,sep='\t',header=None)

os.makedirs(out_dir + '/fasta_N_sub_by_chr',exit_ok=True)

skiprows = 0
expected_rows = 0
chr_list = ''
for i in range(length.shape[0]):
	chr = length.loc[i,0]
	chr_length = length.loc[i,1]
	chr_check = chr.replace('chr','')
	print(chr)
	try: int(chr_check)
	except ValueError:
		print('here')
		if (chr_check != 'X') and (chr_check != 'Y') and (chr_check != 'M'):
			print('not X or Y or M')
			skiprows += math.ceil(chr_length/60)+1
			continue
	fasta = pd.read_csv(fasta_path,header=0,sep='\t',nrows=math.ceil(chr_length/60),skiprows=skiprows)
	if fasta.columns[0] != ('>' + chr): print('panic1')
	skiprows += math.ceil(chr_length/60)+1
	expected_rows += math.ceil(chr_length/60) + 1
	vcf_file = chr_snp_dir + '/' + chr + '.vcf'
	new_fasta = fasta.copy()
	if Path(vcf_file).exists():
		vcf = pd.read_csv(vcf_file,sep='\t',header=0)
		for row in range(vcf.shape[0]):
			zero_pos = vcf.loc[row,'POS'] - 1
			vcf_ref = vcf.loc[row,'REF']
			fasta_row = new_fasta.iloc[math.floor(zero_pos/60),0]
			fasta_row_pos = zero_pos % 60
			fasta_ref = fasta_row[fasta_row_pos]
			if vcf_ref != fasta_ref: print('panic2')
			new_fasta_row = fasta_row[:fasta_row_pos] + 'N' + fasta_row[fasta_row_pos+1:]
			new_fasta_ref = new_fasta_row[fasta_row_pos]
			if new_fasta_ref != 'N': print('panic3')
			new_fasta.iloc[math.floor(zero_pos/60),0] = new_fasta_row
	new_fasta.to_csv(out_dir + '/fasta_N_sub_by_chr/' + chr + '.fa',sep='\t',header=True,index=False)
	chr_list = chr_list + ' ' + out_dir + '/fasta_N_sub_by_chr/' + chr + '.fa'
print(chr_list)
print(expected_rows)
out_name = Path(fasta_path).stem
func = 'cat' + chr_list + ' > ' + out_dir + '/' + out_name + '_Nsub.fa'
subprocess.call(func, shell=True)
# will need to cat resulting chr files and then add in the decoy chr without subs