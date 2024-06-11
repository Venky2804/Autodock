#!/usr/bin/python3
import os, shutil
import argparse, sys
import subprocess
import time
from joblib import Parallel, delayed

# patch joblib progress callback
class BatchCompletionCallBack(object):
    # Added code - start
    global total_n_jobs
    # Added code - end
    def __init__(self, dispatch_timestamp, batch_size, parallel):
        self.dispatch_timestamp = dispatch_timestamp
        self.batch_size = batch_size
        self.parallel = parallel

    def __call__(self, out):
        self.parallel.n_completed_tasks += self.batch_size
        this_batch_duration = time.time() - self.dispatch_timestamp

        self.parallel._backend.batch_completed(self.batch_size,
                                           this_batch_duration)
        self.parallel.print_progress()
        # Added code - start
        progress = self.parallel.n_completed_tasks / total_n_jobs
#        time_remaining = (this_batch_duration / self.batch_size) * (total_n_jobs - self.parallel.n_completed_tasks)
#        print( "\rProgress: [{0:50s}] {1:.1f}% est {2:1f}seconds left".format('#' * int(progress * 50), progress*100, time_remaining/60) , end="", flush=True)
        print("\rProgress: [{0:50s}] {1:.1f}%".format('#' * int(progress * 50), progress*100) , end="", flush=True)
        if self.parallel.n_completed_tasks == total_n_jobs:
            print('\n')
        # Added code - end
        if self.parallel._original_iterator is not None:
            self.parallel.dispatch_next()
      
import joblib.parallel
joblib.parallel.BatchCompletionCallBack = BatchCompletionCallBack

parser = argparse.ArgumentParser(description="A program for docking a ligand to protein")
parser.add_argument("-l", "--ligpath", help="Path of directory with ligand files (pdbqt).", action="store")
parser.add_argument("-R", "--receptor", help="Receptor/Protein (.pdb) file", action="store")
parser.add_argument("-p", "--psh", help="Path of mgltools pythonsh [Default: ~/Downloads/mgltools_x86_64Linux2_1.5.6/bin/pythonsh]", action="store", default="/home/haitoshailesh/Downloads/mgltools_x86_64Linux2_1.5.6/bin/pythonsh")
parser.add_argument("-t", "--threads", type=int, help="Number of threads to use. [Default: 1]", action="store", default="1")
parser.add_argument("-u", "--utildir", help="Path of mgltools utilities [Default: ~/Downloads/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/]", action="store_true", default="/home/haitoshailesh/Downloads/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/")
args = parser.parse_args()

if args.ligpath == None or args.receptor == None:
	print("\nArgument error. Please check your command.\nUse \"-h\" or \"--help\" for help\n")
	sys.exit(1)

if not os.path.isdir(args.utildir):
	print("\n!!!!!ERROR!!!!!\nPlease check path for mgltools utilities\n")
	sys.exit(1)
    
if not os.path.isfile(args.psh):
	print("\n!!!!!ERROR!!!!!\nPlease check path for mgltools pythonsh\n")
	sys.exit(1)

s2 = args.utildir+"prepare_receptor4.py"
s3 = args.utildir+"prepare_gpf4.py"
s5 = args.utildir+"prepare_dpf42.py"
s10 = args.utildir+"summarize_results4.py"
total_n_jobs = len(os.listdir(args.ligpath))

from datetime import datetime

strtime = datetime.now().strftime("%A %-d %B, %Y %I:%M:%S %p")
print("Docking started at:", strtime+"\n")

if args.threads > 1:
	print("Using", args.threads, "threads for Screening")

done = open("finished.txt", "a+")
D = []
d = open("finished.txt").read().splitlines()
for x in d:
	D.append(x)

workdir = os.getcwd()
if not os.path.exists("Results"):
	os.mkdir("Results")

c2 = args.psh+" "+s2+" -r "+args.receptor
print("Preparing receptor file", end="\t-----\t")# -->\t"+c2)
c2 = c2.split()
st2 = subprocess.run(c2, stdout = subprocess.DEVNULL)
if not st2.returncode == 0:
	print("Error preparing receptor file.\nThere must be a problem with input files.\nCheck manually")
	sys.exit(1)
print("Done\nDocking ligands to receptor")

def dock(ligand):
	if not ligand.split('.')[0] in D:
		done = open("finished.txt", "a+") #finished ligands log file
#		print("Docking", args.receptor.split('.')[0], "with", ligand.split('.')[0])
		if not os.path.exists("Results/"+ligand.split('.')[0]):
			os.mkdir("Results/"+ligand.split('.')[0])
		if ligand.endswith(".pdbqt"):
			shutil.copy(args.ligpath+"/"+ligand, "Results/"+ligand.split('.')[0])
			shutil.copy(args.receptor+"qt", "Results/"+ligand.split('.')[0])
			shutil.copy(args.receptor, "Results/"+ligand.split('.')[0])
		os.chdir("Results/"+ligand.split('.')[0])
		c3 = args.psh+" "+s3+" -l "+ligand.split(".")[0]+".pdbqt -r "+args.receptor+"qt -o "+args.receptor.split('.')[0]+"_"+ligand.split(".")[0]+".gpf"
#		print("Grid Parameter File generation")# -->\t"+c3)
		c3 = c3.split()
		st3 = subprocess.run(c3, stdout = subprocess.DEVNULL)
		if not st3.returncode == 0:
			print("Error generating Grid Parameter File.\nThere must be a problem with input files.\nCheck manually")
			sys.exit(1)
		c4 = "autogrid4 -p "+args.receptor.split('.')[0]+"_"+ligand.split(".")[0]+".gpf -l "+args.receptor.split('.')[0]+"_"+ligand.split(".")[0]+".glg"
#		print("GLG file generation")# -->\t"+c4)
		c4 = c4.split()
		st4 = subprocess.run(c4)
		if not st4.returncode == 0:
			print("Error generating GLG File.\nThere must be a problem with input files.\nCheck manually")
			sys.exit(1)
		c5 = args.psh+" "+s5+" -l "+ligand.split(".")[0]+".pdbqt -r "+args.receptor+"qt -p ga_num_evals=1750000 -p ga_pop_size=150 -p ga_run=1 -p rmstol=2.0 -o complex.dpf"
#		print("Docking Parameter File generation")# -->\t"+c5)
		c5 = c5.split()
		st5 = subprocess.run(c5, stdout = subprocess.DEVNULL)
		if not st5.returncode == 0:
			print("Error generating Docking Parameter File.\nThere must be a problem with input files.\nCheck manually")
			sys.exit(1)
		c6 = "autodock4 -p complex.dpf -l complex.dlg"
#		print("Docking")# -->\t"+c6)
		c6 = c6.split()
		st6 = subprocess.run(c6)
		if not st6.returncode == 0:
			print("Error Docking File.\nThere must be an error in previous steps or problem with input files.\nCheck manually")
			sys.exit(1)
#		print("Generating complex.pdb")
		dlg = open("complex.dlg").read().splitlines()
		f7 = open("docked-only-ligand.pdbqt", "w+")
		for line in dlg:
			if line.startswith("DOCKED: "):
				line = line.split("DOCKED: ")[1]
				print(line, file = f7)
		f7.close()
		f8 = open("docked-only-ligand.pdb", "w+")
		c8 = "cut -c-66 docked-only-ligand.pdbqt"
		c8 = c8.split()
		st8 = subprocess.run(c8, stdout=f8)
		f8.close()
		COMP = []
		f9 = open("complex.pdb", "w+")
		rec = open(args.receptor).read().splitlines()
		lig = open("docked-only-ligand.pdb").read().splitlines()
		COMP += rec
		COMP += lig
		for x in COMP:
			if not x == "END" or not x.startswith("END   "):
				print(x, file = f9)
		f9.close()
		os.chdir(workdir+"/Results")
		c10 = args.psh+" "+s10+" -d "+ligand.split('.')[0]+" -t 2 -L -a -o summary_2.0.txt"
		c10 = c10.split()
		st10 = subprocess.run(c10, stdout = subprocess.DEVNULL)
		if not st10.returncode == 0:
			print("Error Summarizing results.\nThere must be an error in previous steps or problem with input files.\nCheck manually")
			sys.exit(1)
		print(ligand.split('.')[0], file = done, flush=True)
		os.chdir(workdir)
		done.close()
#		print("Finished docking", args.receptor.split('.')[0], "with", ligand.split('.')[0])
#	else:
#		print(ligand.split('.')[0]+" docking already finished\nResults can be found at Results/"+ligand.split('.')[0])

if __name__ == "__main__":
	results = Parallel(n_jobs=args.threads, backend="multiprocessing")(delayed(dock)(i) for i in os.listdir(args.ligpath))

os.chdir(workdir+"/Results")
c11 = "sort summary_2.0.txt -k5n -t, -o summary_2.0.sort"
c11 = c11.split()
st11 = subprocess.run(c11)
if not st11.returncode == 0:
	print("Error Sorting results.\nThere must be an error in previous steps or problem with input files.\nCheck manually using\nsort summary_2.0.txt -k5n -t, -o summary_2.0.sort")

done.close()

end = datetime.now().strftime("%A %-d %B, %Y %I:%M:%S %p")
print("\nDocking finished at:", end)
