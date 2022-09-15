import glob
import numpy as np
import os
from prody import *
import shutil
import re

schrodinger_path = '/data/Sergi/general_software/schrodinger2019-1/'

# Define paths for ligands and receptors
datadir = '../data/SARS2'
receptors = glob.glob('../results/SARS2_analysis/extraction/receptor/*')
dirreceptor = os.path.dirname(receptors[0])
receptors = [os.path.basename(receptor) for receptor in receptors]
ligands = glob.glob('../results/SARS2_analysis/extraction/ligand/*')
dirligands = os.path.dirname(ligands[0]) 
ligands = [os.path.basename(ligand) for ligand in ligands]

def split_complex():
    """
    It splits the complex protein-ligand intro different chains,
    maintaining each ligand into each protein chain.
    It also deletes protein chains without ligand.
    """

    results_dir = '../results/SARS2_analysis/split_complex'
    if not os.path.isdir(results_dir):
        os.system('mkdir %s'%results_dir)

    # Dictionay with PDB IDs as keys and ligands as value
    dicCOMPLEX = {}
    for ligand in ligands:
        pdb = ligand.split('_')[0]
        if pdb not in dicCOMPLEX.keys():
            dicCOMPLEX[pdb]=[ligand]
        else:
            dicCOMPLEX[pdb].append(ligand)

    # Look for proteins with or without ligand
    print('Number of initial PDBs: ',str(len(receptors)))
    print('Number of PDBs with at least one LIG: ',str(len(dicCOMPLEX.keys())))

    # Match each ligand with its chain
    final_dict = {}
    for pdbComplex in dicCOMPLEX.keys():
        #print(pdbComplex)
        structureComplex = parsePDB('%s/%s.pdb'%(datadir,pdbComplex))
        hvComplex = structureComplex.select('protein').getHierView()
        chains = [chain.getChid() for chain in hvComplex]

        # coordinate of the binding site center of each chain
        chain_coords = {}
        for chain in chains:
            numres = hvComplex[chain].numResidues()
            if numres < 290 or numres > 330: continue
            resd = hvComplex[chain,145]
            if resd == None: continue
            atom = (resd['SG'])
            if atom == None: continue
            coord = atom.getCoords()
            chain_coords[chain] = coord
        
        # coordinates of each ligand centroid
        lig_coords = {}
        for lig in dicCOMPLEX[pdbComplex]:
            #print(lig)
            lig_struct = parsePDB('%s/%s'%(dirligands, lig))
            coords = calcCenter(lig_struct)
            lig_coords[lig] = coords

        # distance between coordinates + match under threshold of 10 A
        match_dict = {}
        auxLig = False
        for lig, coords in lig_coords.items():
            count = 1
            for chain, bs_coord in chain_coords.items():
                distance = np.linalg.norm(bs_coord-coords)
                if distance < 10:
                    auxLig = True
                    match_dict[chain] = lig
                    chain_struct = structureComplex.select('protein and chain %s'%(chain)).copy()
                    complex_structure = chain_struct + lig_struct
                    writePDB('%s/Mpro_complex_%s%s%s.pdb'%(results_dir, pdbComplex, chain, count), complex_structure)
                    count += 1

        if auxLig:
            final_dict[pdbComplex] = match_dict
    
    # print the dictionary for check the results                  
    print(final_dict)


# superimpose 
def pdbs_superimposition(pdbs,fix_pdb,outdir,verbose=True):
    """
    Given a list of mobile pdbs and a fix superimpose all mobile elements to the fix pdb
    pdbs: 'list'. PDBs list of mobile elements
    fix_pdb: 'str'. Single PDB which will be fixed during the multiple superimpositions (fix element)
    outdir: 'str'. Directory to store all the superimposed PDBs
    """

    PDBnames = [os.path.basename(pdb) for pdb in pdbs]
    IDs = [ pdbname.replace(".pdb","") for pdbname in PDBnames]

    for i, pdb in enumerate(pdbs):
        #Superimpose using TMalign
        pdbname = PDBnames[i]
        if verbose: print(pdbname)
        outname = '%s/%s'%(outdir,IDs[i])
        outname = outname.replace('//','/')

        TMalign_cmd = './../bin/TMalign %s %s -o %s'%(pdb, fix_pdb, outname)
        print(TMalign_cmd)
        os.system(TMalign_cmd)
        os.system('rm %s %s_all %s_atm'%(outname,outname,outname))
        os.system('mv %s_all_atm %s_all_atm.pdb'%(outname,outname))
        os.system('mv %s_all_atm_lig %s_all_atm_lig.pdb'%(outname,outname))

        #Extract mobile pdb from the TMalign output
        structurePDB = parsePDB('%s_all_atm_lig.pdb'%(outname))
        #Load Prody Hierarchical Views
        hvPDB = structurePDB.getHierView()
        #Write superimpose mobile element as a pdb
        aux = '%s/%s_super.pdb'%(outdir,IDs[i])
        aux = aux.replace("//","/")
        writePDB(aux,hvPDB['A'])
        os.system('rm %s_all_atm_lig.pdb'%(outname))

# check if the ligand affects the superimposition
def check_superimposition_coords(prot_pdbs, complex_pdbs):
    """"
    Check if the ligand does not affect the superimposition result
    by looking for difference between coordinates in the superimposed
    proteins and the superimposed ligand-protein complex.
    """
    for complex_pdb in complex_pdbs:
        complex_id = os.path.basename(complex_pdb).split('_')[2]
        complex_struct = parsePDB(complex_pdb)
        complex_struct.getHierView()
        complex_coords = complex_struct.select('resnum 145 and calpha').getCoords()

        for prot_pdb in prot_pdbs:
            prot_id = os.path.basename(prot_pdb).split('_')[1]
            prot_struct = parsePDB(prot_pdb)
            prot_struct.getHierView()
            prot_coords = prot_struct.select('resnum 145 and calpha').getCoords()

            dist = np.linalg.norm(complex_coords-prot_coords)
            if dist < 2:
                print("ok")

# extract ligands
def pdb_extract(pdb_dir):
    """This function extracts the target structure (receptor), the ligands,
    the waters and the ions/cofactors of the list of pdb files."""
    
    os.chdir(pdb_dir)
    PDBs = glob.glob("*_super.pdb*")
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

def _organize_extraction_files(pdb_dir, out_dir):
    """This function organize all the output files obtained in the pdb extraction
    into their respective directories."""
    
    cwd = os.getcwd()
    os.chdir(out_dir)
    os.mkdir("ligand")
    os.mkdir("receptor")
    os.mkdir("water")
    os.mkdir("cof_ion")
    os.chdir(cwd)

    for filename in os.listdir(pdb_dir):
        if re.search("ligand", filename):
            shutil.move("%s/%s"%(pdb_dir, filename), "%s/ligand/%s"%(out_dir, filename))
        if re.search("receptor", filename):
            shutil.move("%s/%s"%(pdb_dir, filename), "%s/receptor/%s"%(out_dir, filename))
        if re.search("water", filename):
            shutil.move("%s/%s"%(pdb_dir, filename), "%s/water/%s"%(out_dir, filename))
        if re.search("cof_ion", filename):
            shutil.move("%s/%s"%(pdb_dir, filename), "%s/cof_ion/%s"%(out_dir, filename))


if __name__ == '__main__':

    results_dir = '../results/SARS2_analysis'
    
    #split_complex()
    
    if not os.path.isdir('%s/superimpose_complex'%results_dir):
        os.system('mkdir %s/superimpose_complex'%results_dir)
    TARGs = glob.glob('%s/split_complex/Mpro_complex*'%results_dir)
    #pdbs_superimposition(pdbs=TARGs,fix_pdb='%s/split_chains/Mpro_6lu7A.pdb'%results_dir,outdir='%s/superimpose_complex/'%results_dir)

    TARGs_prot = glob.glob('%s/superimposition/*_super.pdb'%results_dir)
    TARGs_complex = glob.glob('%s/superimpose_complex/*_super.pdb'%results_dir)
    #check_superimposition_coords(prot_pdbs=TARGs_prot, complex_pdbs=TARGs_complex)

    if not os.path.isdir('%s/extract_complex'%results_dir):
        os.system('mkdir %s/extract_complex'%results_dir)
    #pdb_extract('%s/superimpose_complex'%results_dir)
    _organize_extraction_files(pdb_dir='%s/superimpose_complex'%results_dir, out_dir='%s/extract_complex/'%results_dir)

# git push
