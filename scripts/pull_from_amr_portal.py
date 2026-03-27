#!/usr/bin/env python

import os, sys
import argparse
import pandas as pd
import urllib.request
#from contextlib import closing
#import yaml
import subprocess



def create_parser():
	parser = argparse.ArgumentParser(prog='pull_from_amr_portal.py', description="Script to pull completed assemblies from AMR Portal")
	parser.add_argument('-i', '--input', help="Input CSV from AMR Portal (manual download)")
	parser.add_argument('-o', '--output-dir', default=f"{os.getcwd()}", help="Output directory for GCA downloads")
	parser.add_argument('-f', '--filter', help="Use to filter columns")
	parser.add_argument('-c', '--include-count', action='store_true', help="Output counts as TSV")
	parser.add_argument('-t', '--threads', default=1, help="Number of threads")
	parser.add_argument('--pull-gca', action='store_true', help='Pull completed assembly FASTA output via the GCA identifiers')
	parser.add_argument('--pull-sra', action='store_true', help='Pull FASTQ sequencing read(s) from the SRA')

	return parser

def determine_latest():
	"""
	Pull from releases.yml to determine latest version to pull GCA
	"""
	url = "https://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/releases.yml"
	filename = releases.yml
	with closing(urllib.request.urlopen(url)) as r:
		with open(filename, 'wb') as f:
			version = yaml.load(filename)['latest'] 

	return version

def pull_GCA(accession, output, release=None):
	"""
	Given a GCA accession, pull using enaDataGet
	"""
	#url_path = f"https://ftp.ebi.ac.uk/pub/databases/amr_portal/"
	script_path = os.path.join(script_dir, 'enaBrowserTools','python3', 'enDataGet')
	if not os.path.exists(script_path):
		print("MISSING ENDATAGET!")
		sys.exit(1)
	else:
		params = ['python', script_path, '-f', 'fasta', f"{accession}"]
		subprocess.run(params)

def pull_SRA(accession, output, threads):
	"""
	Using custom download_sra.sh, download read files per accession
	"""
	script_path = os.path.join(script_dir, 'download_sra.sh')
	if not os.path.exists(script_path):
		print("MISSING DOWNLOAD_SRA.SH")
		sys.exit(1)
	else:
		params = ['bash', script_path, '-s', f"{accession}", '-o', start_dir, '-t', f"{threads}"]
		subprocess.run(params)

def parse_table(input, col=None):
	main_table = pd.read_csv(input, sep=',', header=0)

	if col is None:
		initial_filter = main_table
	else:
		initial_filter = main_table[col]

	return main_table, initial_filter

def counts(df, by, print_output=True):
	"""
	Given a column name, return grouped numbers of unique terms. If print, only returns to stdout, else returns as variable
	"""
	if print_output:
		print(df[by].value_counts(dropna=True))
	else:
		return(df[by].value_counts(dropna=True))

def main(args):
	# relevant columns
	relevant_col = ['pheno_geno_merged-assembly_ID',
				'pheno_geno_merged-SRA_accession',
				'pheno_geno_merged-antibiotic_name',
				'pheno_geno_merged-resistance_phenotype',
				'pheno_geno_merged-species',
				'pheno_geno_merged-host',
				'pheno_geno_merged-ast_standard',
				'pheno_geno_merged-isolation_source_category',
				'pheno_geno_merged-country']
	# initial tables (filtered by relevant columns)
	whole, subset = parse_table(args.input)

	print(whole)
	print(subset)

	if args.include_count:
		for item in relevant_col:
			with open(os.path.join(os.getcwd(), f'{item}.tsv'), 'w+') as out:
				out.write(pd.Series.to_csv(counts(subset, item, print_output=False), sep="\t"))

	# pull from GCA
	if args.pull_gca:
		for id in subset['pheno_geno_merged-assembly_ID'].tolist():
			pull_GCA(id, os.path.join(start_dir, 'gca_assemblies'))

	# pull from SRA
	if args.pull_sra:
		for id in subset['pheno_geno_merged-SRA_accession'].tolist():
			pull_SRA(id, os.path.join(start_dir, 'sra_reads'), args.threads)





pull_SRA()
if __name__ == "__main__":
	parser = create_parser()
	args = parser.parse_args()
	start_dir = os.getcwd()
	script_dir = os.path.abspath(sys.argv[0])

	main(args)


