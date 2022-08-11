import os 
import sys
import glob
import subprocess
import time
import argparse
import shutil
import errno
import re

shrodinger_path = '/data/general_software/schrodinger2019-1/'

def parserfunc():
    parser = argparse.ArgumentParser(
        description ='Given a directory with .pdb.gz/pdb files, extract all the ligands from them')
    
    parser.add_argument('-d', dest="pdbdir", help = "PDBs directory", required=True)
    parser.add_argument('-o', dest="outdir", help = "Output directory", required=True)
    
    args = parser.parse_args()
    return args

def main():
    arg = parserfunc()
    pdbdir = arg.pdbdir

    os.chdir(pdbdir)

    PDBs = glob.glob("*.pdb*")
    for PDB in PDBs:
        if '.pdb.gz' in PDB:
            cmd1 = "gunzip %s"%PDB
            os.system(cmd1)
            PDB = PDB.replace('.pdb.gz','.pdb')
            ID = PDB.replace(".pdb","")
        elif '.pdb' in PDB:
            ID = PDB.replace(".pdb","")
        print('Extracting %s ligands...' %(ID))
        cmd2 = '%srun split_structure.py -m pdb %s %s.pdb -many_files'%(shrodinger_path,PDB,ID)
        #print(cmd2)
        os.system(cmd2)

def organize_output_files():
    arg = parserfunc()
    outdir = arg.outdir
    pdbdir = arg.pdbdir

    try:
        os.mkdir(outdir)
    except OSError:
        shutil.rmtree(outdir)
        os.mkdir(outdir)
    
    os.chdir(outdir)

    try:
        os.mkdir("ligand")
        os.mkdir("receptor")
        os.mkdir("water")
        os.mkdir("cof_ion")
    except OSError:
        shutil.rmtree("ligand")
        shutil.rmtree("receptor")
        shutil.rmtree("water")
        shutil.rmtree("cof_ion")
        os.mkdir("ligand")
        os.mkdir("receptor")
        os.mkdir("water")
        os.mkdir("cof_ion")
    
    os.chdir('../')
    for filename in os.listdir(pdbdir):
        if re.search("ligand", filename):
            shutil.move("%s/%s"%(pdbdir, filename), "%s/ligand/%s"%(outdir, filename))
        if re.search("receptor", filename):
            shutil.move("%s/%s"%(pdbdir, filename), "%s/receptor/%s"%(outdir, filename))
        if re.search("water", filename):
            shutil.move("%s/%s"%(pdbdir, filename), "%s/water/%s"%(outdir, filename))
        if re.search("cof_ion", filename):
            shutil.move("%s/%s"%(pdbdir, filename), "%s/cof_ion/%s"%(outdir, filename))




if __name__ == '__main__':
    main()
    os.chdir('../')
    organize_output_files()
    
