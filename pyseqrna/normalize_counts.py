#!/usr/bin/env python
'''

:Title: This script contains quality trimming tools for pySeqRNA

:Author: Naveen Duhan

:Version: 0.1

'''
import re
import logging
import pandas as pd
import numpy as np
import itertools
import matplotlib.pyplot as plt
from scipy import stats
from collections import defaultdict
from pyseqrna import pyseqrna_utils as pu
from pyseqrna.pyseqrna_utils import PyseqrnaLogger
pLog=PyseqrnaLogger(mode="a", log="ncount")



def boxplot(data =None, countType=None, colors=None, figsize=(20,10), **kwargs):
    """
    This function make a boxplot with boxes colored according to the countType they belong to

    Args:
        data ([float], optional): [list of array-like of float]. Defaults to None.
        countType ([str], optional): [list of string, same length as `data`]. Defaults to None.
        colors ([type], optional): [description]. Defaults to None.
    Other Args:
        kwargs ([dict], optional): [keyword arguments for `plt.boxplot`]
    """
    allTypes = sorted(set(countType))

    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    countType_color = dict(zip(allTypes, itertools.cycle(colors)))

    countType_data = defaultdict(list)


    for d, t in zip(data, countType):

        for c in allTypes:

            countType_data[c].append([])

        countType_data[t][-1] = d

    fig, ax = plt.subplots(figsize=figsize)

    lines = []

    for at in allTypes:
        
        for key in ['boxprops', 'whiskerprops', 'flierprops']:

            kwargs.setdefault(key, {}).update(color=countType_color[at])
       
        box = ax.boxplot(countType_data[at], **kwargs)

        lines.append(box['whiskers'][0])

    ax.legend(lines, allTypes)

    return fig, ax





class Normalization():
    """
    This class is for calculation of normalized counts from raw counts
    """

    def __init__(self, countFile=None, featureFile=None, typeFile='GFF', keyType='ncbi', attribute='ID', feature='gene', geneColumn='Gene'):
        """[summary]

        Args:
            countFile ([type], optional): [description]. Defaults to None.
            featureFile ([type], optional): [description]. Defaults to None.
            type (str, optional): [description]. Defaults to 'GFF'.
            attribute (str, optional): [description]. Defaults to 'ID'.
        """

        self.countFile = countFile
        self.featureFile = featureFile
        self.typeFile = typeFile
        self.attribute = attribute
        self.geneColumn = geneColumn
        self.feature = feature
        self.keyType = keyType

        return
    def getGeneLength(self, file, feature='gene', typeFile='GFF', attribute='ID'):
        """[summary]

        Args:
            file ([type], optional): [description]. Defaults to None.
            feature (str, optional): [description]. Defaults to 'gene'.
            type (str, optional): [description]. Defaults to 'GFF'.
            attribute (str, optional): [description]. Defaults to 'ID'.
        """

        gtf = pd.read_csv(file, sep="\t", header=None, comment="#")

        gtf.columns = ['seqname', 'source', 'feature', 'start', 'end', 's1', 'strand', 's2', 'identifier']

        gtf = gtf[gtf['feature'] == feature]

        if typeFile.upper() == 'GFF':
            
            try:

                if self.keyType.lower() == 'ncbi':

                        gtf['Gene'] = list(map(lambda x: re.search(attribute+'[=](.*?)[;]',x,re.MULTILINE).group(1), gtf['identifier'].values.tolist()))
                    
                if self.keyType.lower() =='ensembl':

                        gtf['Gene'] = list(map(lambda x: re.search(attribute+'[=](.*?)[;]',x,re.MULTILINE).group(1), gtf['identifier'].values.tolist()))
            
            except Exception:

                pLog.exception("Attribute not present in GFF file")

        if typeFile.upper() == 'GTF':

            try:
                gtf['Gene'] = list(map(lambda x:  re.search(attribute+'\s+["](.*?)["]',x,re.MULTILINE).group(1), gtf['identifier'].values.tolist()))

            except Exception:

                pLog.exception("Attribute not present in GTF file")

        gene_list = gtf.values.tolist()

        gtf['Length'] = list(map(lambda x: x[4]-x[3]+1,gene_list))

        final = gtf [['Gene', 'Length']].copy()
        final['Gene']= final['Gene'].str.replace("gene:", "")
        final['Gene']= final['Gene'].str.replace("gene-", "")
        final = final.drop_duplicates(subset='Gene')
        finalDF = final.set_index('Gene')

        return finalDF

    def CPM(self, plot=True, figsize=(20,10)):
        """
        This function convert counts to counts per million
        """

        df = pd.read_excel(self.countFile)
        gene_names = df[self.geneColumn]

        countDF = df.set_index(self.geneColumn)

        counts = np.asarray(countDF,dtype=int)

        cpm = (counts * 1e6) / counts.sum(axis=0) 

        cpmDF = pd.DataFrame(data=cpm,index=countDF.index,columns=countDF.columns)

        if plot:
            
            logCounts = list(np.log(counts.T+1))

            logNorm_counts = list(np.log(cpm.T+1))

            fig, ax = boxplot(data=logCounts + logNorm_counts,countType=['Raw counts']*len(countDF.columns)+['CPM counts']*len(countDF.columns),
                         labels= list(countDF.columns)+list(cpmDF.columns), figsize=figsize )

            ax.set_xlabel('Sample Name')

            ax.set_ylabel('log counts')
            
            ax.set_xticklabels(list(countDF.columns)+list(cpmDF.columns), rotation = 90)
        
        cpmDF.insert(0, 'Gene', gene_names)
        return cpmDF, fig, ax


    def RPKM(self, plot=True, figsize=(20,10)):
        """
        This function convert counts to reads per killobase per million
        Returns:
        [type]: [description]
        """

        df = pd.read_excel(self.countFile)

        countDF = df.set_index(self.geneColumn)

        geneDF = self.getGeneLength(self.featureFile,self.feature,self.typeFile, self.attribute)

        match_index = pd.Index.intersection(countDF.index, geneDF.index)

        counts = np.asarray(countDF.loc[match_index], dtype=int)

        gene_lengths = np.asarray(geneDF.loc[match_index]['Length'],dtype=int)

        gene_names = np.array(match_index)

        total_read_per_sample = counts.sum(axis=0)  

        rpkm =  1e9 * counts / (total_read_per_sample[np.newaxis, :] * gene_lengths[:, np.newaxis])

        rpkmDF = pd.DataFrame(data=rpkm,columns=countDF.columns)
        

        if plot:
            
            logCounts = list(np.log(counts.T+1))

            logNorm_counts = list(np.log(rpkm.T+1))

            fig, ax = boxplot(data=logCounts + logNorm_counts,countType=['Raw counts']*len(countDF.columns)+['RPKM counts']*len(countDF.columns),
                         labels= list(countDF.columns)+list(rpkmDF.columns), figsize=figsize )

            ax.set_xlabel('Sample Name')

            plt.xticks(rotation=90)

            ax.set_ylabel('log counts')
        rpkmDF.insert(0, 'Gene', gene_names)
       
        return rpkmDF, fig


    def TPM(self, plot=True, figsize=(20,10)):

        """
        This function convert counts to reads per killobase per million
        Returns:
        [type]: [description]
        """

        df = pd.read_excel(self.countFile)

        countDF = df.set_index(self.geneColumn)

        geneDF = self.getGeneLength(self.featureFile,self.feature,self.typeFile, self.attribute)

        match_index = pd.Index.intersection(countDF.index, geneDF.index)

        counts = np.asarray(countDF.loc[match_index], dtype=int)

        gene_lengths = np.asarray(geneDF.loc[match_index]['Length'],dtype=int)

        gene_names = np.array(match_index)

        total_read_per_sample = counts.sum(axis=0)  

        data =  1e3 * counts / (total_read_per_sample[np.newaxis, :] * gene_lengths[:, np.newaxis])

        tpm = (data * 1e6) / data.sum()

        tpmDF = pd.DataFrame(data=tpm,index=gene_names,columns=countDF.columns)

        if plot:
            
            logCounts = list(np.log(counts.T+1))

            logNorm_counts = list(np.log(tpm.T+1))

            fig, ax = boxplot(data=logCounts + logNorm_counts,countType=['Raw counts']*len(countDF.columns)+['TPM counts']*len(countDF.columns),
                         labels= list(countDF.columns)+list(tpmDF.columns), figsize=figsize )

            ax.set_xlabel('Sample Name')

            ax.set_ylabel('log counts')

        tpmDF.insert(0, 'Gene', gene_names)

        return tpmDF, fig, ax


    def meanRatioCount(self):

        df = pd.read_excel(self.countFile)
        gene_names = df['Gene']
        df2 = df.set_index('Gene')
        col = df2.columns.tolist()
        df3 = df2[df2[col] !=0]
        sqrt = stats.gmean(df3[col],axis=1)
        
        df4 = df3.div(sqrt,axis=0)
        
        m = df4.median()
        
        res = df2.div(m,axis=1)

        res.insert(0, 'Gene', gene_names)
        
        return res


















  

    

