import argparse
parser = argparse.ArgumentParser(description='select a group name, for example:"UNC-44-Proximal"')
parser.add_argument("ID", help="merge a group of loom files according to the group name")
args = parser.parse_args()
print(args.ID)
import os
import loompy

path="/athena/elementolab/scratch/yah2014/Projects/scRNAseq-Lung/data/velocyto"
os.chdir(path) # change current path
print(os.getcwd())
# List all filer folder's names.
file_folders=os.listdir(os.getcwd())  # list files
files=[s for s in file_folders if args.ID in s]
print(files)


# on the command line do: cp file1.loom merged.loom
output_filename=args.ID+"_merged.loom"
loompy.combine(files, output_filename, key="Accession")

