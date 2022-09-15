import glob
import numpy as np
import os
from prody import *

receptors = glob.glob('../results/SARS2_analysis/extraction/receptor/*.pdb')
receptordir = os.path.dirname(receptors[0])
receptors = [os.path.basename(receptor) for receptor in receptors]
LIGs = glob.glob('../results/SARS2_analysis/extraction/ligand/*_ligand*')
liganddir = os.path.dirname(LIGs[0])
LIGs = [os.path.basename(LIG) for LIG in LIGs]
distthreshold = 12
resthresholds = [290,320]

dicLIGs = {}
for LIG in LIGs:
    LIGid = LIG.split("_")[0]
    if LIGid not in dicLIGs.keys():
        dicLIGs[LIGid] = [LIG]
    else:
        dicLIGs[LIGid].append(LIG)

print(len(receptors),len(dicLIGs.keys()))

#Check which LIGS are at the ATP binding site
for i,receptor in enumerate(dicLIGs.keys()):
    structureReceptor = parsePDB('%s/%s_receptor1.pdb'%(receptordir,receptor))
    hvReceptor = structureReceptor.getHierView()
    chains = [chain.getChid() for chain in hvReceptor]
    for chain in chains:
        numRes = hvReceptor[chain].numResidues()
        if numRes < resthresholds[0] or numRes > resthresholds[1]:
            continue
        resd = hvReceptor[chain,145]
        if resd == None: continue
        resd_ca = (resd['SG'])
        if resd_ca == None: continue
        CAcoord = resd_ca.getCoords()
        for lig in dicLIGs[receptor]:
            count=1
            structureLig = parsePDB('%s/%s'%(liganddir,lig))
            LIGcoord = calcCenter(structureLig)
            dist = np.linalg.norm(CAcoord-LIGcoord)
            if dist < 12:
                selection = 'chain %s'%chain
                structureChain = structureReceptor.select(selection).copy()
                structureComplex = structureChain + structureLig 
                print('%s%s_%d.pdb'%(receptor,chain,count))
                writePDB('%s%s_%d.pdb'%(receptor,chain,count),structureComplex)
                count+=1 
