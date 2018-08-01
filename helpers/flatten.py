#!/usr/bin/python3
import os
import argparse

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Flatten Array[Array[File]] into Array[File]')
	parser.add_argument('--tsv', type=str, help='result of write_tsv in cromwell')
	args = parser.parse_args()
	files = []
	with open(args.tsv) as tsv_file:
		for line in tsv_file:
			for file in line.split('\t'):
				if file:
					files.append(file)
	for file in files:
		print(file)