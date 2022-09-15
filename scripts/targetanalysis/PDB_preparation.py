import os 
import sys
import glob
import subprocess
import time
import argparse
import shutil
import errno
import re
from prody import *

schrodinger_path = '/data/Sergi/general_software/schrodinger2019-1/'

def parserfunc():
    parser = argparse.ArgumentParser(
        description ='Given a directory with pdb files, redactar blablabla')

    parser.add_argument('-d', dest="pdbdir", help = "PDBs directory", required=True)
    parser.add_argument('-o', dest="outdir", help = "Output directory", required=True)
    
    args = parser.parse_args()
    return args

def pdb_extract():
    """This function extracts only the target structure (receptor) of the list
    of pdb files. It removes ligands, waters and ions/cofactors."""

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
        cmd2 = '%srun split_structure.py -m pdb %s %s.pdb -many_files'%(schrodinger_path,PDB,ID)
        #print(cmd2)
        os.system(cmd2)

def _organize_output_files():
    """This function organize all the output files obtained in the pdb extraction
    into their respective directories."""
    
    arg = parserfunc()
    outdir = arg.outdir
    pdbdir = arg.pdbdir
    
    os.chdir("../%s"%(outdir))

    os.mkdir("extraction")
    os.chdir("extraction")
    os.mkdir("ligand")
    os.mkdir("receptor")
    os.mkdir("water")
    os.mkdir("cof_ion")

    os.chdir('../../')
    for filename in os.listdir(pdbdir):
        if re.search("ligand", filename):
            shutil.move("%s/%s"%(pdbdir, filename), "%s/extraction/ligand/%s"%(outdir, filename))
        if re.search("receptor", filename):
            shutil.move("%s/%s"%(pdbdir, filename), "%s/extraction/receptor/%s"%(outdir, filename))
        if re.search("water", filename):
            shutil.move("%s/%s"%(pdbdir, filename), "%s/extraction/water/%s"%(outdir, filename))
        if re.search("cof_ion", filename):
            shutil.move("%s/%s"%(pdbdir, filename), "%s/extraction/cof_ion/%s"%(outdir, filename))

def get_protChains(pdbs,outname,delimiter=None,upresfilter=None,lowresfilter=None,verbose=True):
    """
    From a pdbs list retrieve all its chains (filtered if asked) into individual pdbs with only protein
    pdbs: 'list'. PDBs list
    outname: 'str'. Prefix to add to the new PDB chains
    delimiter: 'str'. Delimiter to obtain an identifier from each PDB name
    upresfilter: 'int'. Residue upper threhsold. Maximum number of residues allowed per chain
    lowresfilter: 'int'. Residue lower threshold. Minimum number of residues per chain
    """
    os.mkdir("split_chains")

    PDBnames = [os.path.basename(pdb) for pdb in pdbs]

    if delimiter != None:
        IDs = [ pdbname.replace(".pdb","").split(delimiter)[0] for pdbname in PDBnames]
    else:
        IDs = [ pdbname.replace(".pdb","") for pdbname in PDBnames]
    
    for i,pdb in enumerate(pdbs):
        pdbname = PDBnames[i]
        if verbose: print(pdbname)
        structurePDB = parsePDB('%s'%(pdb)) #From PDB to prody atom class
        hvPDB = structurePDB.getHierView() #Get prody hierarchical view
        chains = [chain.getChid() for chain in hvPDB] #Get chain identifiers
        for chain in chains:
            if verbose: print(chain)
            nresidues =  hvPDB[chain].numResidues()
            if upresfilter != None:
                if nresidues > upresfilter:
                    print("Too many residues for chain %s in %s"%(chain,pdbname))
                    continue
            if lowresfilter != None:
                if nresidues < lowresfilter:
                    print("Too few residues for chain %s in %s"%(chain,pdbname))
                    continue
            # save only protein
            writePDB('%s_%s%s.pdb'%(outname,IDs[i],chain),hvPDB[chain].select('protein'))

def pdbs_superimposition(pdbs,fix_pdb,outdir,verbose=True):
    """
    Given a list of mobile pdbs and a fix superimpose all mobile elements to the fix pdb
    pdbs: 'list'. PDBs list of mobile elements
    fix_pdb: 'str'. Single PDB which will be fixed during the multiple superimpositions (fix element)
    outdir: 'str'. Directory to store all the superimposed PDBs
    """
    os.mkdir("superimposition")

    PDBnames = [os.path.basename(pdb) for pdb in pdbs]
    IDs = [ pdbname.replace(".pdb","") for pdbname in PDBnames]

    for i, pdb in enumerate(pdbs):
        #Superimpose using TMalign
        pdbname = PDBnames[i]
        if verbose: print(pdbname)
        outname = '%s/%s'%(outdir,IDs[i])
        outname = outname.replace('//','/')
        
        TMalign_cmd = './../../bin/TMalign %s %s -o %s'%(pdb, fix_pdb, outname)
        print(TMalign_cmd)
        os.system(TMalign_cmd)
        os.system('rm %s %s_all %s_all_atm_lig %s_atm'%(outname,outname,outname,outname))
        os.system('mv %s_all_atm %s_all_atm.pdb'%(outname,outname))

        #Extract mobile pdb from the TMalign output
        structurePDB = parsePDB('%s_all_atm.pdb'%(outname))
        #Load Prody Hierarchical Views
        hvPDB = structurePDB.getHierView()
        #Write superimpose mobile element as a pdb
        aux = '%s/%s_super.pdb'%(outdir,IDs[i])
        aux = aux.replace("//","/")
        writePDB(aux,hvPDB['A'])
        os.system('rm %s_all_atm.pdb'%(outname))

def prepWizard(pdbs,outfmt='mae',max_processes=4):
    """
    Given a list of pdbs, use PrepWizard from the Schroodinger premium package to preapre the PDB(protonate, add missing side chains etc...)
    pdbs: 'list'. PDB list of elements to prepare
    outfmt: 'str'. Outfile format. Either .mae or .pdb
    max_processes: Number of processors used to paralalize the different executions
    """
    if outfmt != 'mae' and outfmt != 'pdb':
        raise ValueError('outfmt must be either mae or pdb')

    PDBnames = [os.path.basename(pdb) for pdb in pdbs]
    IDs = [ pdbname.replace(".pdb","") for pdbname in PDBnames]
    
    cmd_prepW = ['%s/utilities/prepwizard -fillsidechains -WAIT %s %s_prep.%s'%(schrodinger_path, pdb,IDs[i],outfmt) for i,pdb in enumerate(pdbs)]
    cmd_prepW = [cmd.replace('//','/') for cmd in cmd_prepW]
    processes = set()

    for cmd in cmd_prepW:
        print(cmd)
        processes.add(subprocess.Popen(cmd,shell=True))
    for p in processes:
        p.wait()
        if p.wait() != 0:
            print("There was an error")

def _clean_prepWizard(outdir,outfmt='mae'):
    """
    Move the PrepWizard output to an specified directory
    outdir: 'str'. Output directory
    outfmt: 'str'. Outfile format of the PrepWizard. Either .mae or .pdb
    """
    #os.mkdir("prepWizard")

    if outfmt != 'mae' and outfmt != 'pdb':
        raise ValueError('outfmt must be either mae or pdb')
    
    os.system('mv *_prep.%s %s'%(outfmt,outdir))
    logdir = '%s/logs'%outdir
    logdir = logdir.replace('//','/')
    if not os.path.isdir(logdir):
        os.system('mkdir %s'%logdir)
    os.system('mv *.mae %s'%logdir)
    os.system('mv *.log %s'%logdir)


def main():
    arg = parserfunc()
    outdir = arg.outdir

    # create output directory
    # try:
    #     os.mkdir(outdir)
    # except OSError:
    #     shutil.rmtree(outdir)
    #     os.mkdir(outdir)

    # PDB extraciton
    print("\n---------------PDB EXTRACTION--------------\n")
    pdb_extract()
    _organize_output_files()

    # Split chains
    print("\n---------------SPLIT CHAINS----------------\n")
    os.chdir(outdir)
    TARGs = glob.glob('extraction/receptor/*')
    get_protChains(pdbs=TARGs,outname='split_chains/Mpro',delimiter='_',upresfilter=310,lowresfilter=295)

    # # Superimposition
    # print("\n--------SUPERIMPOSE STRUCTURES-------------\n")
    # TARGs = glob.glob('split_chains/Mpro*')
    # pdbs_superimposition(pdbs=TARGs,fix_pdb='split_chains/Mpro_6lu7A.pdb',outdir='superimposition/')

    # # PrepWizard
    # print("\n------PROTEIN PREPARATION WIZARD------------\n")
    # TARGs = glob.glob('superimposition/*_super.pdb')
    # prepWizard(pdbs=TARGs,outfmt='mae',max_processes=10)
    # _clean_prepWizard(outdir='prepWizard',outfmt='mae')

if __name__ == '__main__':
    main()
