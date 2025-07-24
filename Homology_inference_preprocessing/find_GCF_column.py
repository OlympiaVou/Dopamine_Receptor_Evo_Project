

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("table", type=str)	        # Save the input .csv file
args = parser.parse_args()

input = args.table                              # Save the input file path as input

import pandas as pd


#names=['Date','AgentName','Group','Direction'], skiprows=1, sep='\s+'
df = pd.read_csv(input, sep=',')

print(*df['Genome_Assembly_ID'].values, sep='\n')
