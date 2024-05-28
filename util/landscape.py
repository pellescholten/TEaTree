import os,sys,gzip,datetime,collections,argparse,errno
import matplotlib.pyplot as plt
import csv

parser=argparse.ArgumentParser()
parser.add_argument('-i', metavar='str', type=str, help='TSV file with the following columns: sequence order k80adjust trans transv defBases numCpGs insert delet alignLength fam', required=True)
parser.add_argument('-g', type = int, help="Genome size of species", required=True)
parser.add_argument('-o', metavar='str', type=str, help='Name of output', required=True)
args=parser.parse_args()

gsize=args.g
output_file=args.o
if not output_file.endswith('.tsv'):
    output_file += '.tsv'

dLandscape = {}
with open(args.i, "r") as infile:
   
    for line in infile:
        col = line.rstrip().split()
        length = int(col[5]) + int(col[8])
        order = col[1]
        k = int(float(col[2]))
        
        if dLandscape.get(k):
            if dLandscape.get(k).get(order):
                dLandscape[k][order] += length
            else:
                dLandscape[k][order] = length
        else:
            dLandscape[k] = {}
            dLandscape[k][order] = length


sorted_data = dict(sorted(dLandscape.items()))

# Write the nested dictionary to a CSV file
with open(output_file, 'w', newline='') as csvfile:
    fieldnames = ['Kimura', 'SuperFamily', 'BasePairs', 'Percentage of Genome']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')
    
    writer.writeheader()
    for key, subdict in sorted_data.items():
        for subkey, value in subdict.items():
            divided_value = value / gsize * 100
            writer.writerow({'Kimura': key, 'SuperFamily': subkey, 'BasePairs': value, 'Percentage of Genome': divided_value})


