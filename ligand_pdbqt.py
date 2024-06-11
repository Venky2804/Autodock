#!/usr/bin/python3
import os, shutil
import argparse, sys
import subprocess
from joblib import Parallel, delayed
from progressbar import ProgressBar

parser = argparse.ArgumentParser(description="A program for converting ligand molecules to pdbqt")
parser.add_argument("-l", "--ligpath", help="Path of directory with ligand files.", action="store")
parser.add_argument("-o", "--output", help="Output directory to save pdbqt files [Default: ./ligands_pdbqt", action="store", default="ligands_pdbqt")
parser.add_argument("-p", "--psh", help="Path of mgltools pythonsh [Default: ~/Downloads/mgltools_x86_64Linux2_1.5.6/bin/pythonsh]", action="store", default="/home/haitoshailesh/Downloads/mgltools_x86_64Linux2_1.5.6/bin/pythonsh")
parser.add_argument("-t", "--threads", type=int, help="Number of threads to use. [Default: 1]", action="store", default="1")
parser.add_argument("-u", "--utildir", help="Path of mgltools utilities [Default: ~/Downloads/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/]", action="store", default="/home/haitoshailesh/Downloads/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24")
args = parser.parse_args()


if args.ligpath == None:
	print("\nArgument error. Please check your command.\nUse \"-h\" or \"--help\" for help\n")
	sys.exit(1)

if not os.path.isdir(args.utildir):
	print("\n!!!!!ERROR!!!!!\nPlease check path for mgltools utilities\n")
	sys.exit(1)

if not os.path.isfile(args.psh):
	print("\n!!!!!ERROR!!!!!\nPlease check path for mgltools pythonsh\n")
	sys.exit(1)

s1 = args.utildir+"/prepare_ligand4.py"

from datetime import datetime

time = datetime.now().strftime("%A %-d %B, %Y %I:%M:%S %p")
print("Conversion started at:", time+"\n")

done = open("finished.txt", "a+")
D = []
d = open("finished.txt").read().splitlines()
for x in d:
	D.append(x)

if args.threads > 1:
	print("Using", args.threads, "threads for Screening")

if not os.path.exists(args.output):
	os.mkdir(args.output)

pbar = ProgressBar().start()
count = 0

print("Converting ligand files to pdbqt")

def ligconv(ligand):
	global count
	done = open("finished.txt", "a+")
	if ligand.split('.')[0] not in D:
		c1 = args.psh+" "+s1+" -l "+args.ligpath+"/"+ligand+" -o "+args.output+"/"+ligand.split('.')[0]+".pdbqt"
		c1 = c1.split()
		st1 = subprocess.run(c1, stdout=subprocess.PIPE)
		if not st1.returncode == 0:
			print("Error converting ligand file.\nThere must be a problem with input files.\nCheck manually")
			sys.exit(1)
		print(ligand.split('.')[0], file  = done)
		count += 1
		pbar.update(count)
	done.close()

if __name__ == "__main__":
    	results = Parallel(n_jobs=args.threads, backend="multiprocessing")(delayed(ligconv)(i) for i in os.listdir(args.ligpath))

pbar.finish()

end = datetime.now().strftime("%A %-d %B, %Y %I:%M:%S %p")
print("\nConversion finished at:", end)
