import glob
import argparse
from prody import *
import numpy as np
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, cophenet, fcluster
from scipy.spatial.distance import pdist
import pickle
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import SpectralClustering
from sklearn import manifold

from PDB_preparation import _organize_output_files

def parserfunc():
    parser = argparse.ArgumentParser(
        description ='Given volume matrix i tal, redactar blablabla')

    parser.add_argument('-i', dest="input", help = "Input volume overlapping matrix file", required=True)
    parser.add_argument('-o', dest="outdir", help = "Output directory", required=True)
    
    args = parser.parse_args()
    return args


class VolumeOverlappingMatrix(object):
    """""
    """""
    def __init__(self,csv,IDs=None):
        """
        """
        self.matrix = pd.read_csv(csv,delimiter=',',index_col=0)
        self.IDs = IDs
        if self.IDs!=None:
            self.matrix.set_axis(self.IDs, axis=1, inplace=True)
            self.matrix.set_axis(self.IDs, axis=0, inplace=True)

    def plot_hierarchical(self,out,fontsize=1):
        """
        """
        sns.set(font_scale=fontsize)
        cg = sns.clustermap(self.matrix,cmap="RdBu_r",yticklabels=True,xticklabels=True,vmin=0,vmax=1)
        plt.savefig(out,dpi=300)


    def get_dendrogram(self,verbose=True):
        """
        """
        dendrogram = linkage(self.matrix, 'average',metric='euclidean')
        self.dendrogram = dendrogram
        c, coph_dists = cophenet(self.dendrogram, pdist(self.matrix))
        self.c = c
        self.coph_dists = coph_dists
        if verbose:
            print('Cophenetic Correlation Coefficient: ', str(c))

    def fancy_dendrogram(self,*args, **kwargs):
        max_d = kwargs.pop('max_d', None)
        if max_d and 'color_threshold' not in kwargs:
            kwargs['color_threshold'] = max_d
        annotate_above = kwargs.pop('annotate_above', 0)
        ddata = dendrogram(*args, **kwargs)

        if not kwargs.get('no_plot', False):
            plt.title('Hierarchical Clustering Dendrogram (truncated)')
            plt.xlabel('sample index or (cluster size)')
            plt.ylabel('distance')
            for i, d, c in zip(ddata['icoord'], ddata['dcoord'], ddata['color_list']):
                x = 0.5 * sum(i[1:3])
                y = d[1]
                if y > annotate_above:
                    plt.plot(x, y, 'o', c=c)
                    plt.annotate("%.3g" % y, (x, y), xytext=(0, -5),
                                 textcoords='offset points',
                                 va='top', ha='center')
            if max_d:
                plt.axhline(y=max_d, c='k')
        return ddata

    def plot_dendrogram(self,out,max_d,annotate_above,p=None):
        """
        """
        plt.figure(figsize=(25, 10))
        plt.title('Hierarchical Clustering Dendrogram')
        plt.xlabel('sample index',fontsize=16)
        plt.ylabel('distance',fontsize=16)
        if p==None:
            self.fancy_dendrogram(self.dendrogram,leaf_rotation=90., leaf_font_size=8.,labels=self.IDs,max_d=max_d,annotate_above=annotate_above)
        else:
            self.fancy_dendrogram(self.dendrogram,leaf_rotation=90., leaf_font_size=8.,labels=self.IDs,max_d=max_d,annotate_above=annotate_above,p=p,truncate_mode='lastp',show_contracted=True)
        ax = plt.gca()
        ax.tick_params(axis='x', which='major', labelsize=5)
        ax.tick_params(axis='y', which='major', labelsize=14)
        plt.savefig(out,dpi=300)

    def get_dendrogram_clusters(self,max_d):
        """
        """
        clusters = fcluster(self.dendrogram, max_d, criterion='distance')
        self.clusters = clusters
        clusters_dic = {}
        for i,id in enumerate(self.IDs):
            cluster = self.clusters[i]
            if cluster not in clusters_dic:
                clusters_dic[cluster]=[id]
            else:
                clusters_dic[cluster].append(id)
        self.clusters_dic = clusters_dic

    def save_dendrogram_clusters(self,out,verbose=True):
        """
        """
        if verbose:
            print("Clusters obtained and its respective targets:")
            print(self.clusters_dic)
            print("\nLength of each cluster:")
            for cluster in self.clusters_dic.keys():
                print(cluster,len(self.clusters_dic[cluster]))

        with open(out+'.p', 'wb') as handle:
            pickle.dump(self.clusters_dic, handle)

    def get_dendrogram_clustersCenters(self):
        """
        """
        self.clustersCenters = {}
        for cluster in self.clusters_dic.keys():
            cluster_elements = self.clusters_dic[cluster]
            volume_matrix_cluster = self.matrix.loc[cluster_elements,cluster_elements]
            volume_matrix_cluster['mean'] = volume_matrix_cluster.mean(axis=1)
            print('center cluster%d: '%cluster, volume_matrix_cluster['mean'].idxmax())

            self.clustersCenters[cluster] = volume_matrix_cluster['mean'].idxmax()
        


    def remove_cluster_outliers(self, list_clusters):
        """
        """
        outlayer_clusters = list_clusters

        for cluster in self.clusters_dic.keys():
            if cluster in outlayer_clusters:
                #print(cluster)
                outlayers = self.clusters_dic[cluster]
                #print(outlayers)
                self.matrix.drop(labels=outlayers,axis=1,inplace=True)
                self.matrix.drop(labels=outlayers,axis=0,inplace=True)
                for outlayer in outlayers:
                    self.IDs.remove(outlayer)
        
    def organize_cluster_structures(self, results_dir, struct_dir, site_dir):
        """
        """
        if not os.path.isdir('%s/clusters'%results_dir):
            os.system('mkdir %s/clusters'%results_dir)

        for cluster in self.clusters_dic.keys():
            cluster_elements = self.clusters_dic[cluster]
            for struct in cluster_elements:
                element_site = '%s/Mpro_%s_super_out.maegz'%(site_dir, struct)
                element_prep = "%s/Mpro_%s_super_prep.mae"%(struct_dir, struct)
                if not os.path.isdir('%s/clusters/%d'%(results_dir, cluster)):
                    os.system('mkdir %s/clusters/%d'%(results_dir, cluster))
                cmd1 = 'cp %s %s/clusters/%d'%(element_site, results_dir, cluster)
                cmd2 = 'cp %s %s/clusters/%d'%(element_prep, results_dir, cluster)
                os.system(cmd1)
                os.system(cmd2)

    def get_closest_elements_cluster_center(self, num_elements):
        for cluster in self.clusters_dic.keys():
            print('Cluster %s:' %cluster)
            center = self.clustersCenters[cluster]
            elements = self.clusters_dic[cluster]
            submatrix = self.matrix.loc[elements,elements]
            print(submatrix.nlargest(num_elements, [center])[center])



    def get_MDS(self):
        """
        """
        data = 1 - self.matrix

        mds = manifold.MDS(n_components=2, random_state=1, dissimilarity="precomputed")
        mds.fit(data)
        self.points = mds.embedding_

        self.targets_dic = {}
        
        for i in range(len(self.points)):
            for target in data:
                self.targets_dic[target] = [self.points[[i],:]]
                
        #plt.scatter(self.points[:,0], self.points[:,1], color='silver', s=150)
        #plt.show()

    def get_agglomerative_clustering(self, num_clusters, out):
        """
        """
        mat = self.matrix
        X = self.points

        model = AgglomerativeClustering(affinity='euclidean', n_clusters=num_clusters,linkage='ward')
        labels = model.fit_predict(X)

        clusters_dic = {}
        clusters_transl = {}

        for i in range(num_clusters):
            T=X[labels==i]  
            for ind in T:
                if i not in clusters_dic:
                    clusters_dic[i] = [ind]
                else:
                    clusters_dic[i].append(ind)
                        
        #for cluster, point in clusters_dic.items():
        #    clusters_dic[cluster] = [target for target in self.targets_dic]
        
        print("Clusters obtained:")
        print(clusters_dic)
        
        plt.scatter(x=X[:,0], y=X[:,1], c=labels, cmap='rainbow', s=5)
        plt.title('Agglomerative clustering')
        plt.xlabel('First MDS dimension')
        plt.ylabel('Second MDS dimension')
        targets = []
        for target, point in self.targets_dic.items():
            targets.append(target)
        for i, txt in enumerate(targets):
            plt.annotate(txt, (X[:,0][i], X[:,1][i]), fontsize = 4)
        plt.savefig(out,dpi=5000)

    def get_spectral_clustering(self, num_clusters, out):
        """
        """
        mat = self.matrix
        X = self.points
        model = SpectralClustering(n_clusters=num_clusters, assign_labels='discretize', random_state=0)
        labels = model.fit_predict(X)

        clusters_dic = {}

        for i in range(num_clusters):
            T=X[labels==i]
            for ind in T:
                if i not in clusters_dic:
                    clusters_dic[i] = [ind]
                else:
                    clusters_dic[i].append(ind)
        
        print("Clusters obtained:")
        print(clusters_dic)
        
        plt.scatter(x=X[:,0], y=X[:,1], c=labels, cmap='rainbow', s=5)
        plt.title('Spectral clustering')
        plt.xlabel('First MDS dimension')
        plt.ylabel('Second MDS dimension')
        targets = []
        for target, point in self.targets_dic.items():
            targets.append(target)
        for i, txt in enumerate(targets):
            plt.annotate(txt, (X[:,0][i], X[:,1][i]), fontsize = 4)
        plt.savefig(out,dpi=5000)


if __name__ == '__main__':
    arg = parserfunc()
    input = arg.input
    outdir = arg.outdir
    
    # Analyse the volume overlapping matrix
    print("\n-----------VOLUME OVERLAPPING MATRIX------------\n")
    f = open('../results/SARS2_analysis/siteMap/Mpro_sites.mae','r')
    IDs = []
    for line in f:
        if 'Mpro' in line:
            ID = line.split('_')[1]
            if ID not in IDs:
                IDs.append(ID)
    volume_matrix = VolumeOverlappingMatrix(input,IDs)
    print(volume_matrix.matrix)

    print("\n-------------HIERARCHICAL CLUSTERING ANALYSIS------------\n")
    # Clustermap
    print("Saving clustermap at %s/\n"%(outdir))
    volume_matrix.plot_hierarchical(out='%s/Mpro_volumematrix_clustermap.pdf'%(outdir))
    print("Correlation metrics:")
    volume_matrix.get_dendrogram()

    # Dendograms
    print("\nSaving dendograms at %s/\n"%(outdir))
    volume_matrix.plot_dendrogram(out='%s/Mpro_volumematrix_dendogram1.pdf'%(outdir),max_d=3,annotate_above=20)
    volume_matrix.plot_dendrogram(out='%s/Mpro_volumematrix_dendogram2.pdf'%(outdir),max_d=3,annotate_above=20)
    volume_matrix.plot_dendrogram(out='%s/Mpro_volumematrix_dendogram_truncated.pdf'%(outdir),max_d=3,annotate_above=0.6,p=4)
    
    # Clusters
    print("Clustering.........\n")
    volume_matrix.get_dendrogram_clusters(max_d=3)
    volume_matrix.save_dendrogram_clusters(out='%s/Mpro_dendrogram_clusters'%(outdir))

    # Center of each cluster --> SELECTED TARGETS
    print("\nGetting the centre of each cluster.....\n")
    print("The selected targets are these ones:")
    with open('%s/Mpro_dendrogram_clusters.p'%(outdir), 'rb') as handle:
        clusters_dic = pickle.load(handle)
    volume_matrix.get_dendrogram_clustersCenters()

    ####### SARS2 ########
    # Remove outliers: insert the cluster number you want to remove in the list
    volume_matrix.remove_cluster_outliers([1, 2])
    # cluster again
    volume_matrix.plot_hierarchical(out='%s/Mpro_volumematrix_clustermap_out.pdf'%(outdir))
    volume_matrix.get_dendrogram()
    volume_matrix.plot_dendrogram(out='%s/Mpro_volumematrix_dendogram_out.pdf'%(outdir),max_d=2.4,annotate_above=20)
    volume_matrix.get_dendrogram_clusters(max_d=2.4)
    volume_matrix.save_dendrogram_clusters(out='%s/Mpro_dendrogram_clusters_out'%(outdir))
    # again
    volume_matrix.remove_cluster_outliers([3,4])
    volume_matrix.plot_hierarchical(out='%s/Mpro_volumematrix_clustermap_out2.pdf'%(outdir))
    volume_matrix.get_dendrogram()
    volume_matrix.plot_dendrogram(out='%s/Mpro_volumematrix_dendogram_out2.pdf'%(outdir),max_d=2.4,annotate_above=20)
    volume_matrix.get_dendrogram_clusters(max_d=2.4)
    volume_matrix.save_dendrogram_clusters(out='%s/Mpro_dendrogram_clusters_out2'%(outdir))
    with open('%s/Mpro_dendrogram_clusters_out2.p'%(outdir), 'rb') as handle:
        clusters_dic = pickle.load(handle)
    volume_matrix.get_dendrogram_clustersCenters()
    # organize cluster structures
    volume_matrix.organize_cluster_structures('../results/SARS2_analysis', struct_dir='../results/SARS2_analysis/prepWizard', site_dir='../results/SARS2_analysis/siteMap')
    # get 10 closest elements to the center of each cluster
    print("\nClosest elements to the center of the cluster:\n")
    volume_matrix.get_closest_elements_cluster_center(10)


    # ####### SARS ########
    # # Remove outliers: insert the cluster number you want to remove in the list
    # volume_matrix.remove_cluster_outliers([1, 4, 5])
    # # cluster again
    # volume_matrix.plot_hierarchical(out='%s/Mpro_volumematrix_clustermap_out.pdf'%(outdir))
    # volume_matrix.get_dendrogram()
    # volume_matrix.plot_dendrogram(out='%s/Mpro_volumematrix_dendogram_out.pdf'%(outdir),max_d=0.9,annotate_above=20)
    # volume_matrix.get_dendrogram_clusters(max_d=0.9)
    # volume_matrix.save_dendrogram_clusters(out='%s/Mpro_dendrogram_clusters_out'%(outdir))
    # with open('%s/Mpro_dendrogram_clusters_out.p'%(outdir), 'rb') as handle:
    #     clusters_dic = pickle.load(handle)
    # volume_matrix.get_dendrogram_clustersCenters()
    # # organize structures
    # volume_matrix.organize_cluster_structures('../results/SARS_analysis', struct_dir='../results/SARS_analysis/prepWizard', site_dir='../results/SARS_analysis/siteMap')
    # # get 10 closest elements to the center of each cluster
    # print("\nClosest elements to the center of the cluster:\n")
    # volume_matrix.get_closest_elements_cluster_center(10)


    # ####### MERS ########
    # print("\nRefinement....\n")
    # # Remove outliers: insert the cluster number you want to remove in the list
    # volume_matrix.remove_cluster_outliers([1])
    # # cluster again
    # volume_matrix.plot_hierarchical(out='%s/Mpro_volumematrix_clustermap_out.pdf'%(outdir))
    # volume_matrix.get_dendrogram()
    # volume_matrix.plot_dendrogram(out='%s/Mpro_volumematrix_dendogram_out.pdf'%(outdir),max_d=0.8,annotate_above=20)
    # volume_matrix.get_dendrogram_clusters(max_d=0.8)
    # volume_matrix.save_dendrogram_clusters(out='%s/Mpro_dendrogram_clusters_out'%(outdir))
    # with open('%s/Mpro_dendrogram_clusters_out.p'%(outdir), 'rb') as handle:
    #     clusters_dic = pickle.load(handle)
    # volume_matrix.get_dendrogram_clustersCenters()
    # # organize structures
    # volume_matrix.organize_cluster_structures('../results/MERS_analysis', struct_dir='../results/MERS_analysis/prepWizard', site_dir='../results/MERS_analysis/siteMap')
    # # get 10 closest elements to the center of each cluster
    # print("\nClosest elements to the center of the cluster:\n")
    # volume_matrix.get_closest_elements_cluster_center(10)







    # print("\n-----------OTHER CLUSTERING METHODS------------------\n")
    # # Multidimensional Scaling
    # volume_matrix.get_MDS()

    # # Agglomerative clustering
    # print("AGLOMERATIVE CLUSTERING\n")
    # volume_matrix.get_agglomerative_clustering(4, out='%s/Mpro_AgglomerativeClustering.pdf'%(outdir))

    # # Spectral clustering
    # print("\nSPECTRAL CLUSTERING\n")
    # volume_matrix.get_spectral_clustering(4, out='%s/Mpro_SpectralClustering.pdf'%(outdir))
