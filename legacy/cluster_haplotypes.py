#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# This script was originally authored by Taylor A Anderson. Please send questions or comments to tayloranthonyanderson@gmail.com
# verified for Python3 on 10-26-2020 (what a year eh?)

#INPUT: A VCF file for a single chromosome (imputed with no missing data) and 7 parameters
#OUTPUTS: A file containing PCs, hierarchically defined clusters for each sample*windowPosition into folder where run

#EXAMPLE COMMMAND
# The command "python3 HaplotypeAnalysis_Visualization/Haplotype_analysis_scripts/cluster_haplotypes.py 'HaplotypeAnalysis_Visualization/Example_files/SL4.0ch09_subset.vcf' ch09 250000 100000 10 2 80 10" would set the following parameters:
#   a chromosome 9 file [will need to change to your input vcf file]
#   an output filestem "ch09" [change to chromosome name you want as filestem]
#   a window size of 250 Kb [even integer]
#   a window step size of 100 Kb [even integer]
#   a minimum of 10 SNPs per window [integer]
#   a minimum d of 2 [even integer]
#   a maximum d of 80 [even integer]
#   a step size for d of 10 [integer]

#Load dependencies
import allel; print('imported scikit-allel', allel.__version__)
import numpy as np; print('imported numpy')
import pandas as pd; print('imported pandas')
from sklearn.cluster import AgglomerativeClustering
from sklearn.preprocessing import scale
from sklearn.metrics import silhouette_score
from sklearn.decomposition import PCA; print('imported scikit-learn')
import sys; print('imported system commands')

#Check input parameters
if(len(sys.argv[1:])!=8):
    print("\n")
    print("!! USER ERROR !!")
    print("Usage: python3 cluster_haplotypes.py [vcf_file] [chromosome_basename] [window_size] [window_step_size] [min_snps_cutoff] [min d] [max d] [step d]\n")
    print("arg1: Specify location of uncompressed chromosome level vcf file\n")
    print("arg2: Specify output file basename (typically the chromosome number)\n")
    print("arg3: Set window size to for iterating through the genome - MUST BE AN EVEN NUMBER\n")
    print("arg4: Set step size for window iterations\n")
    print("arg5: Set the minimum number of SNPs in each window in order to perform calculations, otherwise outputs NaN\n")
    print("arg6: Set the minimum distance threshold for merging clusters\n")
    print("arg7: Set the maximum distance threshold for merging clusters\n")
    print("arg8: Set the distance (d) step size\n")
    print("\n")
    sys.exit()

#Assign the user input parameters
gt,chromosome,window,step,cutoff,min_dist,max_dist,dist_step=sys.argv[1:]
gt=str(gt)
chromosome=str(chromosome)
window=int(window)
step=int(step)
cutoff=int(cutoff)
min_dist=int(min_dist)
max_dist=int(max_dist)
dist_step=int(dist_step)
n_components=2

#Report the parameters to the user
print("Beginning to classify haplotypes with user inputs:\n")
print("genotype (vcf) file: " + gt + "\n")
print("chromosome basename: " + chromosome + "\n")
print("sliding window size: " + str(window) + "\n")
print("sliding window step size: " + str(step) + "\n")
print("minimum number of SNPs in window for classification: " + str(cutoff) + "\n")
print("mergeby distance minimum (inclusive): " + str(min_dist) + "\n")
print("mergeby distance maximum (exclusive): " + str(max_dist) + "\n")
print("distance step size: " + str(dist_step) + "\n")
print("number of principal components (set internally): " + str(n_components) + "\n")

#Define a function that will look through windows of specified size and calculate the following:
def clusterf(gt, window, step, n_components, cutoff, chrom, min_dist, max_dist, dist_step):
    #Load genotype dataset
    callset = allel.read_vcf(gt)
    genos = allel.GenotypeArray(callset['calldata/GT'])
    samples = callset['samples']
    position = callset['variants/POS']

    #Initialize the data array
    window_positions = int( ((position[-1] - position[0])-window) / step)
    n_samples = len(samples)
    n_data_vectors = n_components + 2
    data = np.empty(shape=(window_positions, n_samples, n_data_vectors), dtype="float")

    #Create a vector of window midpoints for the dataframe below
    window_stop = window_start = position[0]
    midpoint=([])
    for i in range(window_positions):
        window_stop = window_start + window
        midpoint.append((window_start + window_stop) / 2)
        window_start = window_start + step

    #Perform the sliding window calculations, calculating PCs, and hierarchical clustering groups
    window_stop = window_start = position[0]
    pca = PCA(n_components=n_components, copy=True, whiten=False, svd_solver='auto', tol=0.0, iterated_power='auto', random_state=None)
    dlength = range(min_dist,max_dist,dist_step)
    silhouette_scores = np.zeros(len(dlength), dtype="float")
    for i in range(window_positions):
        window_stop = window_start + window
        index = (np.greater_equal(position,window_start) & np.less_equal(position,window_stop))
        count = np.sum(index)
        if count > cutoff:
            X = scale(np.transpose(genos[index,:,0] + genos[index,:,1]))
            #Find the distance threshold value for the window that best describes the data using the silhouette score
            j=0
            try:
                for d in range(min_dist,max_dist,dist_step):
                    silhouette_scores[j] = silhouette_score(X, AgglomerativeClustering(distance_threshold=d, n_clusters=None).fit(X).labels_)
                    j+=1
                    d_thresh = min_dist + (np.argmax(silhouette_scores)*dist_step)
                    dlist = np.repeat(np.array(d_thresh), n_samples)
            except ValueError:
                print("merge distance of " + str(d_thresh) + " results in 1 or nsample clusters, setting d to 10 and trying again")
                d_thresh = 10
                dlist = np.repeat(np.array(d_thresh), n_samples)

            hclusters = AgglomerativeClustering(distance_threshold=d_thresh, n_clusters=None).fit_predict(X)
            data[i,:,:] = np.column_stack((pca.fit_transform(X), dlist, hclusters))
        else:
            data[i,:,:] = np.nan
        print("Finished window: ", window_start, window_stop, " for chromosome ", chrom)
        window_start = window_start + step

    #Convert the data array to a pandas data frame, adding column names and passport information
    name_vector = pd.DataFrame(np.tile(samples, window_positions), columns=['Sample'])
    position_vector = pd.DataFrame(np.repeat(np.array(midpoint), n_samples), columns=['Positions'])
    chrom_ID = pd.DataFrame([[chrom] for i in range(len(position_vector))], columns=["Chromosome"])
    df = pd.DataFrame(data.reshape(window_positions*n_samples,n_data_vectors), columns=['PC1','PC2','dist','hclust']) #Need to change this if more than 2 PCs are chosen
    df = pd.concat([df, position_vector, name_vector, chrom_ID], axis=1)
    print("Finished Windowing Chromosome: " + chrom)
    filename = str(chromosome) + "_window" + str(window) + "_step" + str(step) + "_cutoff" + str(cutoff) + "_dmin" + str(min_dist) + "_dmax" + str(max_dist)  + "_dstep" + str(dist_step) + "_pcacomp" + str(n_components) + ".csv"
    print("Writing output to file: " + filename)
    return(df)

#Write out the data frame to file
out = clusterf(gt=gt, window=window, step=step, n_components=n_components, cutoff=cutoff, chrom=chromosome, min_dist=min_dist, max_dist=max_dist, dist_step=dist_step)
out.to_csv(str(chromosome) + "_window" + str(window) + "_step" + str(step) + "_cutoff" + str(cutoff) + "_dmin" + str(min_dist) + "_dmax" + str(max_dist)  + "_dstep" + str(dist_step) + "_pcacomp" + str(n_components) + ".csv")
