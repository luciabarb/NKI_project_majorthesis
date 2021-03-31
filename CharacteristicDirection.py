##Python test 

import h5py
import os
import numpy as np
import pandas as pd
import csv
import itertools
import sys 
import random
import math

from datetime import datetime

import warnings
from sklearn import decomposition
from sklearn.decomposition import PCA
from scipy.stats import chi2
from scipy.stats import spearmanr
from scipy.stats.mstats import zscore
from sklearn.metrics import pairwise_distances
from PyAstronomy import pyasl
from sklearn.metrics.pairwise import cosine_distances

def obtain_metadata(path_h5_directory, type='notnull'):
  
  metadata_sample_path = os.path.join('data','GEO','GSE70138','GSE70138_Broad_LINCS_inst_info_2017-03-06.txt')
  metadata_sample = pd.read_csv(metadata_sample_path,sep='\t')
  
  if type == 'null':  
    metadata_sample_path = os.path.join('data','GEO','GSE70138','metadata_6drugsDMSO_GSE70138_Broad_LINCS_inst_info_2017-03-06.txt')
    metadata_sample = pd.read_csv(metadata_sample_path, sep=' ')
  
  list_h5_file = os.listdir(path_h5_directory)
  
  metadata_sample['inst_id'] = metadata_sample['inst_id'].str.replace(':','_')
  
  all_inst_id = []
  for h5_file in list_h5_file: #Select the metadata found on the h5_files
    
    data_h5 = h5py.File(os.path.join(path_h5_directory, h5_file), 'r')
    
    all_inst_id.extend(np.char.decode(np.array(data_h5['rowid'])).tolist())
    
  
  
  metadata_sample = metadata_sample[metadata_sample['inst_id'].isin(all_inst_id)]
  
  
  
  metadata_sample['rep'] = metadata_sample.apply(lambda row: '_'.join([str(row.pert_iname),
                          str(row.pert_dose), str(row.pert_time), str(row.cell_id)]), axis=1)
  
  return(metadata_sample)

def select_control(path_h5_file, index_stat_rep, metadata_sample):
  
  data_h5 = h5py.File(path_h5_file, 'r')
  
  ##Substract data from data_h5
  data_inst_id = np.char.decode(np.array(data_h5['rowid'])).tolist()
  
  
  ##Filter metadata for the sample we have
  metadata_sample_filtered = metadata_sample[metadata_sample['inst_id'].isin(data_inst_id)]
  
  metadata_sample_DMSO = metadata_sample_filtered.loc[metadata_sample_filtered['pert_iname'] == 'DMSO']

  
  inst_id_sample_DMSO = metadata_sample_DMSO['inst_id'].to_numpy()
  
  ##Select h5 file that are DMSO -- control 
  #Check position of our matching ID
  
  index_DMSO = np.in1d(data_inst_id,inst_id_sample_DMSO).nonzero()
  
  level3_data_DMSO =  np.array(data_h5['data'])[index_DMSO,index_stat_rep,:][0].tolist() 
  
  return(level3_data_DMSO)

def select_perturbed_id(path_h5_file,metadata_sample):
  data_h5 = h5py.File(path_h5_file, 'r')
  
  ##Substract data from data_h5
  data_inst_id = np.char.decode(np.array(data_h5['rowid'])).tolist()
  
  ##Filter metadata for the sample we have
  metadata_sample_filtered = metadata_sample[metadata_sample['inst_id'].isin(data_inst_id)]
  
  
  
  metadata_sample_perturbed = metadata_sample_filtered.loc[metadata_sample_filtered['pert_iname'] != 'DMSO']
  
  unique_rep_inst_id = {} #Save in a dictionary the name of the ID + characteristics
  
  for rep in metadata_sample_perturbed.rep.unique(): #Check unique characteristics
    metadata_sample_rep = metadata_sample_perturbed.loc[metadata_sample_perturbed['rep'] == rep]
    
    unique_rep_inst_id[rep] = metadata_sample_rep['inst_id'].to_numpy()[0] #None of the plates have more than one experiment with equal conditions in the same plate
  
  #Returns dicitonary, being keys the characteristic and values the inst_id (sample id)
  return(unique_rep_inst_id)

def select_perturbed(path_h5_file,id_perturbed,index_stat_rep):
  data_h5 = h5py.File(path_h5_file, 'r')
  
  ##Substract data from data_h5
  data_inst_id = np.char.decode(np.array(data_h5['rowid'])).tolist()
  
  
  #Check position of our matching ID
  inst_id = [i for i,x in enumerate(data_inst_id) if x == id_perturbed]
  
  level3_data_perturbed =  np.array(data_h5['data'])[inst_id,index_stat_rep,:][0].tolist() ##TODO: Make this more efficient
  
  
  #Returns dicitonary, being keys the characteristic and values the inst_id (sample id)
  return(level3_data_perturbed)

def dictionaries_ID(path_h5_directory, metadata):
  
  e=0
  unique_plate_ID = {} #First level of the dictionary is all the plates (non replicated)
  
  
  list_h5_file = os.listdir(path_h5_directory) ## All plates ID
  
  list_h5_file = [s for s in list_h5_file if "null_mito" in s]
  
  h5_files_processed = [] #List to append processed h5 files to not repeat them later
  for h5_file in list_h5_file: # Loop through all h5 files PLATES
    
    if h5_file not in h5_files_processed: #Make sure we haven't seen this plate before
        
      replicate_plate_ID = {} #Second level of the dictionary is all the replicates id
      replicate_plate_ID[h5_file] = {}
      
      for matching_file in list_h5_file: #Check which one matches
        if matching_file != h5_file:
            
          if h5_file.startswith(matching_file[:matching_file.find('X')]): #Remove the `_array_QN`+ 6 last characters as this is the cell id
              
            replicate_plate_ID[matching_file] = {}
              
            h5_files_processed.append(matching_file)
            h5_files_processed.append(h5_file)
              
        
      for replicate in replicate_plate_ID.keys():  ## Loop through all replicate plate files
        unique_rep_inst_id = select_perturbed_id(os.path.join(path_h5_directory,replicate), metadata) #Obtain 3rd level of dictionary
        replicate_plate_ID[replicate] = unique_rep_inst_id
        e+=len(unique_rep_inst_id)
        
        
      unique_plate_ID[h5_file[:h5_file.find('X')]] = replicate_plate_ID
      # Returns dictionary of diciontaries: 1st level plate replicates, 2nd level plate IDs (file name), 3rd level the type of experiment and in values the pert_iname
  
  
  
  
  print('Total experiments to be run:', e)
  return(unique_plate_ID) 

def nipals(X,a,it=10,tol=1e-4):
    # Nipals algorithm for Principal Component Analysis
    # This function is written largely based on nipals function from R chemometrics package.

    Xh = np.array(X)
    
    (obsCount,varCount) = np.shape(Xh)
    #Xh = X - np.tile(np.mean(X,axis=0),(obsCount,1)) #I think they are centralizing but this has already been made
    

    T = np.zeros((obsCount,a))
    
    P = np.zeros((varCount,a))
    
    
    pcvar = np.zeros(varCount)
    varTotal = np.sum(np.var(Xh,axis=0))
    
    currVar = varTotal
    nr = 0
    
    for h in range(a):
        th = np.reshape(Xh[:,0],(obsCount,-1)) #Select a column vector xi of matrix X (Scores)
        
        ende = False
        
        while not ende:
            nr = nr + 1
           
            if (np.dot(Xh.T,th) == 0).any():
              print('np.dot(Xh.T,th)',np.dot(Xh.T,th))
              print('(Xh.T)',Xh.T)
              print('th',th)
              sys.exit()
              
            if np.isnan(np.dot(th.T,th)).any():
              print(nr)
              sys.exit()
              
              
            ph = np.dot(Xh.T,th)/np.dot(th.T,th) #Project matrix X onto Th (scores) in order to find the loadings (ph)
            
            ph = ph/np.sqrt(np.dot(ph.T,ph)) # Normalize the loading vector v to length 1
            
            thnew = np.dot(Xh,ph)/np.dot(ph.T,ph) #Project X onto loadings (ph) in order to find new score vector (Thnew)
            
            prec = np.dot((thnew-th).T,(thnew-th)) #Check convergence: calculate difference between the previous scores and the current scores.
            
            
            th = thnew
            
            if prec <= np.power(tol,2): #If the difference prec is larger than a pre-defined threshold, repeat
                ende = True
                
            if it <= nr:
                ende = True
                print('Iteration stops without convergence')
                print(prec)
                sys.exit()
                
                #https://cran.r-project.org/web/packages/chemometrics/vignettes/chemometrics-vignette.pdf
        
        Xh = Xh - np.dot(th,ph.T) #Error
        T[:,h] = th[:,0]
        P[:,h] = ph[:,0]
        oldVar = currVar
        currVar = np.sum(np.var(Xh,axis=0))
        pcvar[h] = ( oldVar - currVar )/varTotal
        nr = 0
        
    return T, P, pcvar
  
def PCAoutliers(X, m1):
  '''Based on https://github.com/sorgerlab/L1000chDir/blob/master/chDirCode/pcaOutliers.m
  X: Control samples matrix where rows are genes and columns are samples.
  m1: boolean vector, where control samples are found in original matrix.
  '''
  
  X = X.transpose()
  
  dist = pairwise_distances(X, metric='euclidean')
  
  means = dist.mean(axis=0)#mean value of each column (samples)
  
  
  max_outlier = round(means.shape[0]*0.5) #TODO
  
  r = pyasl.generalizedESD(means,maxOLs=max_outlier,alpha=0.01) #TODO: Double check the maxOLs
    
  indices_outliers = r[1] 
    
  indices_notoutliers = list(set(np.arange(dist.shape[0])) - set(indices_outliers))
  print('------------ Number outliers control:',len(indices_outliers))
    
  m1_position = np.where(m1 == True)[0]
  m1_outliers = m1_position[indices_notoutliers]
    
  m1_final = np.isin(np.arange(0,len(m1)),  m1_outliers)
  
  return(m1_final)

def removeConstantGenes(X, genes_id):
  '''Based on https://github.com/sorgerlab/L1000chDir/blob/master/chDirCode/pcaOutliers.m
  remove genes with significant low variance using the outlier function
  X: Matrix where rows are genes and columns are samples.
  '''
  
  logStd = np.log(X.std(1)) #log, Std for row
  Idx_X = np.arange(len(logStd))
  Ιdx = np.array(genes_id) #Original index 978 genes
  
  if len(Idx_X) != len(Ιdx):
    raise ValueError("List genes and matrix dimensions do not conicide")
  
  ids = dict(zip(Ιdx,Idx_X))
  
  InfΒo, nonInfΒo  = np.isinf(logStd), ~np.isinf(logStd) #Look for the inf
  InfIdx, nonInfIdx = Ιdx[InfΒo], Ιdx[nonInfΒo]
  
  max_outliers = len(nonInfIdx) - 2
  
  r  = pyasl.generalizedESD(logStd[nonInfΒo],maxOLs=max_outliers,alpha=0.01) #TODO: Double check the maxOLs
  
  pre_outlierIdx = np.array(r[1]) #Take into account that this is the index for logStd[nonInfΒo], not the original genes
  outlierIdx = nonInfIdx[pre_outlierIdx] #Original indices

  outlierIdx_pos = [ids[Ιdx_genes] for Ιdx_genes in outlierIdx]
  
  constantGenesBo = logStd[outlierIdx_pos] < np.mean(logStd) #if logSD is lower than the mean SD
  constantGenes = outlierIdx[constantGenesBo]
  
  constantGenes = np.concatenate((InfIdx, constantGenes))
  
  print('------------ Number constant genes: ',len(constantGenes))
  
  return(constantGenes)

def chdir(data, sampleclass, genes, gamma=1., sort=True, calculate_sig=False, nnull=10, sig_only=False, norm_vector=True):
  """
  Calculate the characteristic direction for a gene expression dataset
	
	Input:
		data: numpy.array, is the data matrix of gene expression where rows correspond to genes and columns correspond to samples
		sampleclass: list or numpy.array, labels of the samples, it has to be consist of 0, 1 and 2, with 0 being columns to be excluded, 1 being control and 2 being perturbation
				example: sampleclass = [1,1,1,2,2,2]
		genes: list or numpy.array, row labels for genes 
		gamma: float, regulaized term. A parameter that smooths the covariance matrix and reduces potential noise in the dataset
		sort: bool, whether to sort the output by the absolute value of chdir
		calculate_sig: bool, whether to calculate the significance of characteristic directions
		nnull: int, number of null characteristic directions to calculate for significance
		sig_only: bool, whether to return only significant genes; active only when calculate_sig is True
		norm_vector: bool, whether to return a characteristic direction vector normalized to unit vector
	Output:
		A list of tuples sorted by the absolute value in descending order characteristic directions of genes.
			If calculate_sig is set to True, each tuple contains a third element which is the ratio of characteristic directions to null ChDir
	"""
	
  ## check input
  data.astype(float)
  sampleclass = np.array(list(map(int, sampleclass)))
  # masks
  m_non0 = sampleclass != 0
  m1 = sampleclass[m_non0] == 1
  m2 = sampleclass[m_non0] == 2
  

  if type(gamma) not in [float, int]:
    raise ValueError("gamma has to be a numeric number")
  if set(sampleclass) != set([1,2]) and set(sampleclass) != set([0,1,2]):
    raise ValueError("sampleclass has to be a list whose elements are in only 0, 1 or 2")
  # if m1.sum()<2 or m2.sum()<2:
  # 	raise ValueError("Too few samples to calculate characteristic directions")
  if len(genes) != data.shape[0]:
    raise ValueError("Number of genes does not match the demension of the expression matrix")

  data = data[:, m_non0]
  
  
  ############################ Detect outliers control samples
  
  m1 = PCAoutliers(data[:,m1],m1) #Outliers control
  
  data_m1, data_m2 = data[:,m1], data[:,m2]
  
  data = np.concatenate((data_m1,data_m2), axis=1)
  
  m1 = np.repeat([True,False],[data_m1.shape[1],data_m2.shape[1]])
  m2 = ~m1
  
  ############################Remove constant genes
  ids = dict(zip(genes,np.arange(data.shape[0]))) #Correlate the original genes id with the number of rows
  
  constantGenesIdx = removeConstantGenes(data,genes)
  notconstantGenesIdx = np.setdiff1d(np.array(genes),constantGenesIdx)
  genes =  np.concatenate((notconstantGenesIdx,constantGenesIdx))
  notconstantGenesIdx_pos = [ids[Ιdx_genes] for Ιdx_genes in notconstantGenesIdx]
  
  data = data[notconstantGenesIdx_pos, :] #Remove constant genes
  ############################
  
  #Normalize data
  #data = zscore(data) # standardize for each genes across samples
  data = (data.T - data.mean(1)).T # Centralize the matrix by the mean of each row, which is a requirement of PCA.
  
  
  ## start to compute
  n1 = m1.sum() # number of controls
  n2 = m2.sum() # number of experiments
  
  
  ## the difference between experiment mean vector and control mean vector.
  meanvec = data[:,m2].mean(axis=1) - data[:,m1].mean(axis=1) 
  
  
   #PCA can be performed before LDA to regularize the problem and avoid over-fitting.
  ############################FOR REGULAR PCA
  '''# initialize the pca object
  pca = PCA()
  pca.fit(data.T)
  
  cumsum = pca.explained_variance_ratio_ # explained variance of each PC
  print('cumsum',cumsum)
  keepPC = len(cumsum[cumsum > 0.001])
  print('Number keepPC: ',keepPC)
  
  v = pca.components_[0:keepPC].T # rotated data 
  r = pca.transform(data.T)[:,0:keepPC] # transformed data'''
  ########################################
  
  
  ############################FOR NIPALS IMPLEMENTED FUNCTION ---
  #  the number of output components desired from PCA. We only want to calculate
  #  the chdir in a subspace that capture most variance in order to save computation 
  #  workload. The number is set 20 because considering the number of genes usually 
  #  present in an expression matrix 20 components would  capture most of the variance.
  (rowCount,colCount) = np.shape(data)
  componentsCount = np.min([colCount-1,30]);
  
  # use the nipals PCA algorithm to calculate scores, loadings, and explained_var. 
  # explained_var are the variances captured by each component 
  
  scores, loadings, explained_var = nipals(data.T,componentsCount,1e5,1e-4)
  scores = scores.T
  loadings = loadings.T
  
  
  ## compute the number of PCs to keep
  
  # We only want components that cpature 95% of the total variance or a little above.
  captured_variance = 0
  for i in range(len(explained_var)):
      captured_variance += explained_var[i]
      if captured_variance > 0.999:
          break
  keepPC = i+1
  scores = scores[0:keepPC] # R in Neil's algorithm
  loadings = loadings[0:keepPC] # V in Neil's algorithm

  r = scores.T
  v = loadings.T
  ############################
  
  
  ############################FOR NIPALS FUNC.
  '''pca = PCA(data.T, method='nipals')
  cumsum =  pca.eigenvals/sum(pca.eigenvals)
  
  #keepPC = len(cumsum[cumsum > 0.001])
  keepPC = len(cumsum)
  print('cumsum',cumsum)
  pca_components = np.asarray(pca.loadings).T
  v = pca_components[0:keepPC].T
  r = pca.scores.T[:,0:keepPC] ##pca.scores is ncomp x nobs, select ncomp of interest
  '''
  ############################

  #print('Shape of loadings: ',v.shape) #978xkeepPC
  #print('Shape of scores: ',r.shape) #obsxkeepPC
  
  
  dd = ( np.dot(r[m1].T,r[m1]) + np.dot(r[m2].T,r[m2]) ) / float(n1+n2-2) # covariance
  sigma = np.mean(np.diag(dd)) # the scalar covariance

  shrunkMats = np.linalg.inv(gamma*dd + sigma*(1-gamma)*np.eye(keepPC))
  
  b = np.dot(v, np.dot(np.dot(v.T, meanvec), shrunkMats))

  if norm_vector:
    b /= np.linalg.norm(b) # normalize b to unit vector
  
  b = np.concatenate((b, np.repeat(0,len(constantGenesIdx) ))) #Add constant genes to 0 
  
  grouped = zip([abs(item) for item in b],b,genes)
  
  if sort:
    grouped = sorted(grouped,key=lambda x: x[0], reverse=True)
  
  if not calculate_sig: # return sorted b and genes.
    #res = [(item[1],item[2]) for item in grouped]
    
    res = pd.DataFrame({'CD':[item[1] for item in grouped]},index=[item[2] for item in grouped]).sort_index(axis=0)
    
    
    return res
  else: # generate a null distribution of chdirs
    nu = n1 + n2 - 2
    y1 = np.random.multivariate_normal(np.zeros(keepPC), dd, nnull).T * np.sqrt(nu / chi2.rvs(nu,size=nnull))
    y2 = np.random.multivariate_normal(np.zeros(keepPC), dd, nnull).T * np.sqrt(nu / chi2.rvs(nu,size=nnull))
    y = y2 - y1 ## y is the null of v
    

    nullchdirs = []
    for col in y.T:
      bn = np.dot(np.dot(np.dot(v,shrunkMats), v.T), np.dot(col,v.T))
      bn /= np.linalg.norm(bn)
      bn = bn ** 2
      bn.sort()
      bn = bn[::-1] ## sort in decending order
      bn = np.concatenate((bn, np.repeat(0,len(constantGenesIdx) ))) #Add constant genes to 0 
      nullchdirs.append(bn)

    nullchdirs = np.array(nullchdirs).T
    
    
    nullchdirs = nullchdirs.mean(axis=1)
    b_s = b ** 2 
    b_s.sort()
    b_s = b_s[::-1] # sorted b in decending order
    
    print('b_s',b_s)
    print('nullchdirs', nullchdirs)
    relerr = b_s / nullchdirs ## relative error
    
    # ratio_to_null
    ratios = np.cumsum(relerr)/np.sum(relerr)- np.linspace(1./len(genes),1,len(genes))
    #Lucía note: first part of ratios is cumulative sum and normalized to 1, then a vector
    ## that goes from 1/978, to 1, in step of 978 i subtracted (expected vector if each gene = contribution)
    ## if the results were not significant then it should be 0
    
    ##Lucía note: The # significant genes: b is ordered (components which contributes the most) 
    ## and then the ratios are caculated, the # genes are cut when max. ratio is found
    
    
    res = [(item[1],item[2], ratio) for item, ratio in zip(grouped, ratios)] 
    print('Number of significant genes: %s'%(np.argmax(ratios)+1))
    lol
    if sig_only:
      return res[0:np.argmax(ratios)+1]
    else:
      return res, res[0:np.argmax(ratios)+1]

def final_CD_pval(CD_replicates):
  #CD_replicates: pandas data frame, rows n_genes and columns replicates 
  null_cosdis_path = os.path.join('data','L1000_bayesian','Bayesian_GSE70138_Level5','null_distribution_cosdis_distribution.txt')  
  null_cosdis = pd.read_csv(null_cosdis_path, sep=',',index_col=0)
  
  genes = CD_replicates.index
  
  if CD_replicates.shape[1] == 1:
    #If there's only one replicate
    Di = np.array(CD_replicates['CD'])
    
  else:
    #More than one replicate
    sum_CDi, CDi_div = 0, 0
    print('------------------- Correlation plates replicates')
    
    
    for (CDi_id, CDi) in CD_replicates.iteritems():
      CDi = CDi.to_numpy()
      wi = 0
      cos_distance_CDi = 0
      for (CDj_id, CDj) in CD_replicates.drop(CDi_id, axis=1).iteritems():
        #Compute cos distance between replicates
        CDj = np.array(CDj)
        cos_dist = cos_distance(CDi, CDj)
        cos_distance_CDi += cos_dist #sum the cosine distances between CDi and all its replicates
        print('-------------------', CDi_id, ' and ', CDj_id,'Cos dist: ',cos_dist)
      
      cos_distance_CDi = cos_distance_CDi/(len(CD_replicates.columns)-1) #Normalize cos distance by number of replicates
      p_value = distribution2pval(np.array(null_cosdis['0']), cos_distance_CDi) #Compute p.value, compared to null distribution
      print('-------------------', CDi_id, 'p-value: ',p_value)
      
      
      if p_value <= 0.05: #If p.value significative, add CDi
        CDi_div += 1
        sum_CDi += CDi
      
    if CDi_div == 1: #If only one CDi was significative, add that one
      Di= sum_CDi
      
    elif CDi_div == 0: #If none of the replicates CDi were significative, set everything to 0
      Di = np.array([0]*978)
      
    else: #Else, normalize over the number of replicates which p.value <0.05
      Di = sum_CDi/CDi_div
      Di = Di/np.linalg.norm(Di)
  
  grouped = list(zip(Di,genes))
  
  return(grouped)
               
def save_documents_plate(data,file_name,output_path):
  
  output_path = os.path.join(output_path,file_name)
  
  f = open(output_path, 'w')
  
  if len(data[0]) == 3:
    line = ' '.join(('CD','gene','p_value'))
  else:
    line = ' '.join(('CD','gene'))
    
  f.write(line + '\n')
  for t in data:
    line = ' '.join(str(x) for x in t)
    f.write(line + '\n')
  f.close()

def chdir_plate(path_h5_directory,metadata,output_path):
  
 
  genes_ID = list(range(1,979))
  
  
  unique_plate_ID = dictionaries_ID(path_h5_directory,metadata)
  
  
  for main_plate in  list(unique_plate_ID.keys()): #Loop through 1st dictionary 
    
    
    processed_experiments = []
    name_file_plates_dict =  unique_plate_ID[main_plate] #These plates have the same experiments
    
    for rep_plate in list(name_file_plates_dict.keys())[::-1]: #Loop through 2nd dictionary
      experiments_dict = name_file_plates_dict[rep_plate]
      
      
      for experiment in list(experiments_dict.keys())[::-1]:  #Loop through 3rd dictionary 
        
        
        if experiment not in processed_experiments:
          processed_experiments.append(experiment)
          
          print('--Experiment ID:', experiment)
          print('##########################################################################')
          for index_stat_rep in range(0,100): #Loop through the 100 statitical replicates
            
            print('New replicate ID: ', index_stat_rep  ,'(experiment:',experiment,')' )
            print('------------------------------------------------------------')
            
            
            diff_analysis = pd.DataFrame(index=genes_ID)
            file_name = ''.join(('chdir_landmark','_',main_plate,'exp_',experiment,'_rep_',str(index_stat_rep),'.txt'))
            
            output_path_files = os.listdir(output_path)#Plates already processed
            if file_name not in output_path_files: #If it hasn't been processed yet (in case I have to re-run it)
              
              for n,rep_name_file_plates in enumerate(name_file_plates_dict.keys()):#Loop through 2nd dictionary
              
                data_level3 = []
              
                print('----Plate ID:',rep_name_file_plates)
                
                if experiment in name_file_plates_dict[rep_name_file_plates].keys():
                  
                  path_file = os.path.join(path_h5_directory,rep_name_file_plates)
                  index_experiment = name_file_plates_dict[rep_name_file_plates][experiment]
                    
                    
                  data_level3.append(select_perturbed(path_file,index_experiment,index_stat_rep))
                  
                  
                  data_level3_control = select_control(path_file,index_stat_rep,metadata)
                  data_level3.extend(data_level3_control)
                    
                    
                  list_sample = [2]+[1]*len(data_level3_control)  #List samples first 2 is the perturbed samples, 1 is for control samples
                  
                    
                  print('------ Num perturbed:', len(data_level3)-len(data_level3_control))
                  print('------ Num control: ',len(data_level3_control))
                    
                  data_level3 = np.transpose(np.array(data_level3))
                
                  diff_analysis_rep = chdir(data_level3, list_sample, genes_ID, calculate_sig=False, nnull=500)
                  
                  diff_analysis = diff_analysis.join(diff_analysis_rep, rsuffix=n )
                  
                  
                else: 
                  
                  print('------------------------------------------------------------')
                  print('For ', rep_name_file_plates, 'experiment: ', experiment, 'does not exist.')
                  print('------------------------------------------------------------')
                  
              diff_analysis = final_CD_pval(diff_analysis)
              save_documents_plate(diff_analysis,file_name, output_path)
              print('------------------------------------------------------------')
          
          print('##########################################################################')
          print('New experiment')
          
def cos_distance(CD1,CD2):
  
  CD1 = CD1.reshape(1, -1)
  CD2 = CD2.reshape(1, -1)
  cos_distance = cosine_distances(CD1,CD2)[0][0]
  
  return(cos_distance)

def PermuPmetric(path_h5_directory,final_output_path,metadata_sample,c,k, genes):
  ''''generate null distribution of projected distances by only sampling control
  replicates. H0= variation between technical replicates = variation between the CDs
                  for an equal number of perturnabgens (randomly selected from =/= cell line + treatment)
  
  Input:
    data is a matrix(978*n) of normlized 978-gene control replicates 
    c is the prevalent control replicate count in a plate
    k is the prevalent experiment replicate count (if most experiments have 4 replicates, k will be 4).
  Output:
   pMetrics is a vector of projected distances that represent the null distribution.
   unitVTotal store all the computed null chdirs.
   
  p-value of an experiment could be easily computed by compare its projected 
  distance to this pMetrics, the null distriubtion of projected distances.
  '''
  print('Start null distribution')
  def random_combination(iterable, r):
    "Random selection from itertools.combinations(iterable, r)"
    pool = tuple(iterable)
    n = len(pool)
    indices = sorted(random.sample(range(n), r))
    return tuple(pool[i] for i in indices)
  
  def find_data(numeric_id):
    numeric_id = str("{:.2f}".format(round(numeric_id/100, 2)))
    index_rep = int(numeric_id[-2:])
    id = int(numeric_id[:(len(numeric_id)-3)])
    
    metadata = metadata_sample.iloc[[id]]
    
    file = [s for s in list_h5_file if metadata['det_plate'].values[0] in s][0]
    path_h5_file = os.path.join(path_h5_directory, file)
    
    data_h5 = h5py.File(path_h5_file, 'r')
    data_inst_id = np.char.decode(np.array(data_h5['rowid'])).tolist()
    
    index_exp = [i for i, s in enumerate(data_inst_id) if metadata['inst_id'].values[0] in s][0]
    
    level3_data =  np.array(data_h5['data'])[index_exp,index_rep,:].tolist()
    
    return(level3_data)
  
  list_h5_file = os.listdir(path_h5_directory) ## All plates ID
  
  #Perform 10000 permutations at most.
  maxi = 10000
  
  sampleCount = metadata_sample.shape[0]*100
  sampleIdx = range(sampleCount)
  print('Combinations....')
  combinations = np.array([list(random_combination(sampleIdx,c)) for _ in range(maxi)]) #Only maxi random combinations of the vector split in c
  combinationsCount = maxi  #All possible combinations
  print('Combinations done')
  
  pMetrics = np.zeros((combinationsCount,1)) #Initialize distance vectors
  unitVTotal = np.zeros((len(genes),combinationsCount)) #Initialize chdir matrix 978xcombinationcounts 
  
  for i in range(combinationsCount):
    print('---New combination:', i)
    combination = combinations[i]  #Chose only one
    
    pseudoExpIdx = np.setxor1d(sampleIdx,combination) #Make sure ctrl and exp are differnt samples
    
    pseudoCtrl = []
    for comb in combination:
      pseudoCtrl.append(find_data(comb))
    
    pseudoExp = [find_data(random.choice(pseudoExpIdx))] 
    list_sample = [2]*len(pseudoExp)+[1]*len(pseudoCtrl)  #List samples first 2 is the perturbed samples, 1 is for control samples
    
    data = np.concatenate((np.transpose(np.array(pseudoExp)),np.transpose(np.array(pseudoCtrl))), axis=1)
    
    unitV = chdir(data,list_sample, genes)
    
    
    unitVTotal[:,i] = unitV.loc[:,'CD'].values
  
  pd.DataFrame(unitVTotal).to_csv(os.path.join(final_output_path,"null_Chdir_distribution.txt"))
  for i in range(combinationsCount):
    paring_sample = list(range(combinationsCount))
    paring_sample.remove(i)
    
    randomi = random.sample(paring_sample, k)
    sum_cos_distance = 0
    
    for e in randomi:
      sum_cos_distance += cos_distance(unitVTotal[:,i],unitVTotal[:,e])
      
    pMetrics[i] = sum_cos_distance/len(randomi)
    
    pd.DataFrame(pMetrics).to_csv(os.path.join(final_output_path,"null_cosdis_distribution.txt"))
  return(pMetrics, unitVTotal)

def distribution2pval(distribution,obsDist):
  '''% calculate the pval of an experiment
  Input:
  distribution: output of permuall.m function.
  obsDist: the average cosine distance among chdir replicates of an 
          experiment, calculated using: mean(pdist(expmChdirReps','cosine'))
  Output: p value.'''
  
  pval = len(np.where((distribution >= obsDist)==False)[0])/distribution.size
    
  return(pval)

def CD_diff_labs(final_output_path):
  
  output_path_files = sorted(os.listdir(final_output_path)) #Plates processed
  
  output_path_files = [s for s in output_path_files if "null_mito" in s]
  
  processed_same_exp = []
  
  genes_ID = list(range(1,979))
  
  
  for level5_file in output_path_files:
    same_exp = level5_file[level5_file.find('LJP')+6:]
    print(same_exp)
    
    file_name = ''.join(('final_chdir_landmark_GSE70138_DPEAK',same_exp))
    
    if same_exp not in processed_same_exp:
          print(same_exp)
          print('########################################')
          print('-ID: ', same_exp)
          processed_same_exp.append(same_exp)
          diff_analysis = pd.DataFrame(index=genes_ID)
          
          name_files_same_exp = [x for x in output_path_files if same_exp in x]
          
          for n, name_file in enumerate(name_files_same_exp):
            print('------File name:', name_file)
            
            CD = pd.read_csv(os.path.join(final_output_path,name_file),sep=" ") 
            CD = CD.set_index(CD.gene).drop('gene', axis=1)
            diff_analysis = diff_analysis.join(CD, rsuffix=n)
            
            
          diff_analysis = final_CD_pval(diff_analysis)
          
          save_documents_plate(diff_analysis,file_name, os.path.join(final_output_path,'final'))
          
  
  return(cos_distance)

def main():
  ##h5 file
  path_h5_file = os.path.join('data','L1000_bayesian','Bayesian_GSE70138_Level3_samples_QN')
  
  final_output_path = os.path.join('data','L1000_bayesian','Bayesian_GSE70138_Level5')
  
  metadata = obtain_metadata(path_h5_file)
  
  PermuPmetric( path_h5_file,final_output_path,metadata,30,4, list(range(1,979)))

  chdir_plate(path_h5_file,metadata,final_output_path)
  
  
  CD_diff_labs(final_output_path)

if __name__ == '__main__':
    main()


