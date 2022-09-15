import glob
from prody import *
import numpy as np
import os
import time
import subprocess
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, cophenet, fcluster 
from scipy.spatial.distance import pdist
import pickle

schrodinger_path = '/data/general_software/schrodinger2019-1'

def parserfunc():
    parser = argparse.ArgumentParser(
        description ='Given a directory with mae files, redactar blablabla')

    parser.add_argument('-d', dest="maedir", help = "MAEs directory", required=True)
    parser.add_argument('-o', dest="outdir", help = "Output directory", required=True)
    
    args = parser.parse_args()
    return args

def siteMap(maes,asl,delimiter=None,outfmt='mae',max_processes=4):
    """
    Run a SiteMap calculation for a list of MAEs (can't be pdbs).
    maes: 'list'. MAE list of elements
    asl: 'str'. ASL (atom specification Language)
    delimiter: 'str'. Delimiter to obtain an identifier from each MAE name
    outfmt: 'str'. Outfile format. Either .mae or .pdb
    max_processes: Number of processors used to paralalize the different executions
    """
    MAEnames = [os.path.basename(mae) for mae in maes]
    if delimiter != None:
        IDs = [ maename.replace(".mae","").split(delimiter)[0] for maename in MAEnames]
    else:
        IDs = [ maename.replace(".mae","") for maename in MAEnames]

    cmd_SiteMap = ['%s/sitemap -j %s -prot %s -sitebox 12 -resolution standard -reportsize 20 -writestructs no -writevis yes -maxsites 1 -siteasl "%s" -WAIT'%(schrodinger_path,IDs[i],mae,asl) for i,mae in enumerate(maes)]
    cmd_SiteMap = [cmd.replace("//","/") for cmd in cmd_SiteMap]

    processes = set()

    for cmd in cmd_SiteMap:
        print(cmd)
        processes.add(subprocess.Popen(cmd,shell=True))
    for p in processes:
        p.wait()
        if p.wait() != 0:
            print("There was an error")

def _clean_siteMap(outdir,outfmt='maegz'):
    """
    Move the SiteMap output to an specified directory
    outdir: 'str'. Output directory
    outfmt: 'str'. Outfile format of the PrepWizard. Either .mae or .pdb
    """

    if outfmt != 'maegz':
        raise ValueError('outfmt must be maegz')
    os.system('mv *.%s %s'%(outfmt,outdir))
    logdir = '%s/logs'%outdir
    logdir = logdir.replace('//','/')
    if not os.path.isdir(logdir):
        os.system('mkdir %s'%logdir)
    os.system('mv *.vis %s'%logdir)
    os.system('mv *.smap %s'%logdir)
    os.system('mv *.log %s'%logdir)

def _group_siteMap(sites,out,outdir):
    """
    Group all volume sites from SiteMap into a single mae file
    sites:
    out:
    outdir:
    """
    conc_sites = ''
    for site in sites:
        conc_sites = conc_sites + ' ' + site

    cmd = '%s/utilities/structcat -imae %s -omae %s'%(schrodinger_path,conc_sites,out)
    cmd = cmd.replace('//','/')
    print(cmd)
    os.system(cmd)
    try:
        aux = 'mv %s %s'%(out,outdir)
        os.system(aux)
    except:
        print('wrong outdir')

def _uncompress_maegz(inp):
    """
    inp:
    """
    out = inp.replace('.maegz','.mae')
    cmd = '%s/utilities/structcat -imae %s -omae %s'%(schrodinger_path,inp,out)
    cmd = cmd.replace('//','/')
    print(cmd)
    os.system(cmd)

def get_volumeOverlapMatrix(sites,out,max_processes=4):
    """
    Generate pairwise volume overlap matrix
    sites: 'str'. single file containing multiple SitMap files
    """
    cmd = '%s/run volume_cluster.py -j %s -HOST localhost:%d -sc -r 2 %s'%(schrodinger_path,out,max_processes,sites)
    cmd = cmd.replace('//','/')
    print(cmd)
    p = subprocess.Popen(cmd, shell=True)
    p.wait()
    if p.wait() != 0:
            print("There was an error")
    

def _clean_volumeMatrix(out,coutdir):
    """
    Move the VolumeMatrix output to an specified directory
    out: 'str'
    outdir: 'str'. Output directory
    """
    if not os.path.isdir(coutdir):
        os.system('mkdir %s'%coutdir)

    os.system('mv *.csv %s'%(coutdir))
    logdir = '%s/logs/'%coutdir
    logdir = logdir.replace('//','/')
    if not os.path.isdir(logdir):
        os.system('mkdir %s'%logdir)
    os.system('mv *.mae %s'%logdir)
    os.system('mv *.log %s'%logdir)



if __name__ == '__main__':
    arg = parserfunc()
    inpdir = arg.maedir
    outdir = arg.outdir

    # create output directory
    if not os.path.isdir('%s/siteMap'%outdir):
        os.system('mkdir %s/siteMap'%outdir)

    # Compute the volume of each target specific binding site
    print("\n-------------RUNNING SITEMAP----------------\n")
    TARGs = glob.glob('%s/*_prep.mae'%(inpdir))
    siteMap(maes=TARGs,asl = "(res.num 145) AND ((atom.ptype \' HB2 \'))",delimiter='_prep',outfmt='mae',max_processes=30)
    _clean_siteMap(outdir='%s/siteMap'%(outdir))
    sites = glob.glob('%s/siteMap/*_out.maegz'%outdir)
    _group_siteMap(sites=sites,out='Mpro_sites.maegz',outdir='%s/siteMap/'%outdir)
    _uncompress_maegz(inp='%s/siteMap/Mpro_sites.maegz'%outdir)

    # Find targets without binding site arround the specified atom
    print("\n--------------CHECK STEP--------------------\n")
    print("These targets do not have the binding site around the specified atom. Please, remove them for further analysis:")
    TARGs_IDs = [os.path.basename(TARG).split('_')[1] for TARG in TARGs]
    sites_IDs = [os.path.basename(site).split('_')[1] for site in sites]
    print(set(TARGs_IDs)-set(sites_IDs))

    # Get the volume overlapping matrix of the target sites
    print("\n-----------VOLUME OVERLAPPING MATRIX----------\n")
    get_volumeOverlapMatrix(sites='%s/siteMap/Mpro_sites.maegz'%(outdir),out='Mpro_volumeMatrix',max_processes=4)
    _clean_volumeMatrix(out='Mpro_volumeMatrix',coutdir='%s/volumeMatrix/'%(outdir))


