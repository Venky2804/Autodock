#!/usr/bin/python3
#script for running docking from s4
import os
import argparse, sys

parser = argparse.ArgumentParser(description="A program for docking a ligand to protein")
parser.add_argument("-l", "--ligand", help="Ligand file (.pdb).", action="store")
parser.add_argument("-R", "--receptor", help="Receptor/Protein (.pdb) file", action="store")
parser.add_argument("-u", "--utils", help="Path of mgltools utilities", action="store")
args = parser.parse_args()

psh = args.utils+"bin/pythonsh"
utildir = args.utils+"MGLToolsPckgs/AutoDockTools/Utilities24/"

if args.ligand == None or args.receptor == None:
    print("\nArgument error. Please check your command.\nUse \"-h\" or \"--help\" for help\n")
    sys.exit(1)

if not os.path.isdir(args.utils):
    print("\n!!!!!ERROR!!!!!\nPlease check path for mgltools utilities\n")
    sys.exit(1)

if not os.path.isdir(utildir):
    print("\n!!!!!ERROR!!!!!\nMgltools utilities not found at given location. Please check the location\n")
    sys.exit(1)
    
if not os.path.isfile(psh):
    print("\n!!!!!ERROR!!!!!\nPythonsh not found at given location. Please check the location\n")
    sys.exit(1)

s1 = utildir+"prepare_ligand4.py"
s2 = utildir+"prepare_receptor4.py"
s3 = utildir+"prepare_gpf4.py"
s5 = utildir+"prepare_dpf42.py"

from datetime import datetime

time = datetime.now().strftime("%A %-d %B, %Y %I:%M:%S %p")
print("Docking started at:", time+"\n")

print("Running AutoDock with", args.ligand.split('.')[0], "and", args.receptor.split('.')[0])
c1 = psh+" "+s1+" -l "+args.ligand
print("\nStep 1 -->\t"+c1)
os.system(c1)
c2 = psh+" "+s2+" -r "+args.receptor
print("Step 2 -->\t"+c2)
os.system(c2)
c3 = psh+" "+s3+" -l "+args.ligand+"qt -r "+args.receptor+"qt -o "+args.receptor.split('.')[0]+".gpf"
print("Step 3 -->\t"+c3)
os.system(c3)
c4 = "autogrid4 -p "+args.receptor.split('.')[0]+".gpf -l "+args.receptor.split('.')[0]+".glg"
print("Step 4 -->\t"+c4)
os.system(c4)
c5 = psh+" "+s5+" -l "+args.ligand+"qt -r "+args.receptor+"qt"+" \ "+" -p ga_num_evals=1750000 \ -p ga_pop_size=150 \ -p ga_run=1 \ -p rmstol=2.0 -o complex.dpf"
print("Step 5 -->\t"+c5)
os.system(c5)
c6 = "autodock4 -p complex.dpf -l complex.dlg"
print("Step 6 -->\t"+c6)
os.system(c6)
c7 = "grep '^DOCKED' complex.dlg | cut -c9- > docked-only-ligand.pdbqt"
print("Step 7 -->\t"+c7)
os.system(c7)
c8 = "cut -c-66 docked-only-ligand.pdbqt > docked-only-ligand.pdb"
print("Step 8 -->\t"+c8)
os.system(c8)
c9 = "cat "+args.receptor+" docked-only-ligand.pdb | grep -v '^END   ' | grep -v '^END$' > complex.pdb"
print("Step 9 -->\t"+c9)
os.system(c9)

end = datetime.now().strftime("%A %-d %B, %Y %I:%M:%S %p")
print("\nDocking finished at:", end)

print("\nFirst open complex.pdb in Pymol and export as PDB")
print("Then Open complex_pymol.pdb in Chimera\nReceptor color - Light blue\nLigand color - Light green")
print("Find Hydrogen bonds using Tools>Structural analysis>FindHBond")
print("\tCheck include intra-residue H-bond and If endpoint atom hidden, show > Apply>OK")
print("\tHydrogen bond color - Dark blue")
print("Check position and aminoacid at H-bond >> Tools>Sequence>Sequence >> Highlight the Aminoacid")
print("\tAction>Label>Residue>Name+specifier")
print("\tResidue label color - Black")
