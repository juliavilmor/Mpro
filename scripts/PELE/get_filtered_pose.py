import numpy as np
import pandas as pd
import argparse
import copy

def main(data,fnames,thresh_up,thresh_low):
	df = pd.read_csv(data)
	#print(df)
	dfaux = copy.copy(df)
	for i,fnames in enumerate(fnames):
		print('num filter low up')
		print(i,fnames,thresh_low[i],thresh_up[i])
		dfaux = dfaux[(dfaux[fnames] < thresh_up[i])]
		dfaux = dfaux[(dfaux[fnames] > thresh_low[i])]
		print(dfaux.shape)
	print(dfaux)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Filtering data.csv from results')
	parser.add_argument('-i', dest='data', help='Data csv file', required=True)
	parser.add_argument('--fnames', dest='fnames', nargs='+' , help='data.csv columns to use as filrers', required = True)
	parser.add_argument('--up_t', dest='thresh_up', nargs='+', type=float, help='list containing upper thresholds (it must be as long as fnames list)')
	parser.add_argument('--low_t', dest='thresh_low', nargs='+', type=float,help='list containing lower thresholds (it must be as long as fnames list)')
	
	args = parser.parse_args()
	data = args.data
	fnames = args.fnames
	thresh_up = args.thresh_up
	thresh_low = args.thresh_low
	if len(fnames) != len(thresh_up) or len(fnames) != len(thresh_low) or len(thresh_up) != len(thresh_low):
		raise ValueError('Filter names and thresholds must have the same length')
	main(data,fnames,thresh_up,thresh_low)
