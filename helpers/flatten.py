#!/usr/bin/python3
import os
import argparse

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Flatten Array[Array[File]] into Array[File]')
	parser.add('--tsv', type=str, help='result of write_tsv in cromwell')
	args = parser.parse_args()
	with open(args.tsv) as file:
		for line in file:
			print(line)