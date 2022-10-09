#!/usr/bin/env python
'''

:Title: This module converts raw read counts to normalize counts

:Created: August 3, 2021

:Author: Naveen Duhan

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



def boxplot(data =None, countType=None,  figsize=(20,10), **kwargs):
    """
    This function make a boxplot with boxes colored according to the countType they belong to

    :param data: List containg log data of raw and normalize read counts.

    :param countType: Columns data of raw and normalize read counts.

    :param figsize: Figure size.
   
    :param kwargs: Other optional arguments for boxplot.
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
        """
        
        :param countFile: Raw read counts file.

        :param featureFile: Gene feature file.

        :param type: Type of gene feature file GFF/GTF. Default is GFF.

        :param keyType: If feature file is from NCBI or ENSEMBL.

        :param attribute: Attribute to consider from gene feature file. Defaults to 'ID'.

        :param feature: type of feature.

        :param geneColumn: Gene column name in raw read count file.

        """

        self.countFile = countFile
        self.featureFile = featureFile
        self.typeFile = typeFile
        self.attribute = attribute
        self.geneColumn = geneColumn
        self.feature = feature
        self.keyType = keyType

        return

    def _getGeneLength(self):
        """
        This function calculate the gene length from gene feature file.
        
        """

        gtf = pd.read_csv(self.featureFile, sep="\t", header=None, comment="#")

        gtf.columns = ['seqname', 'source', 'feature', 'start', 'end', 's1', 'strand', 's2', 'identifier']

        gtf = gtf[gtf['feature'] == self.feature]

        if self.typeFile.upper() == 'GFF':
            
            try:

                if self.keyType.lower() == 'ncbi':

                        gtf['Gene'] = list(map(lambda x: re.search(self.attribute+'[=](.*?)[;]',x,re.MULTILINE).group(1), gtf['identifier'].values.tolist()))
                    
                if self.keyType.lower() =='ensembl':

                        gtf['Gene'] = list(map(lambda x: re.search(self.attribute+'[=](.*?)[;]',x,re.MULTILINE).group(1), gtf['identifier'].values.tolist()))
            
            except Exception:

                pLog.exception("Attribute not present in GFF file")

        if self.typeFile.upper() == 'GTF':

            try:
                gtf['Gene'] = list(map(lambda x:  re.search(self.attribute+'\s+["](.*?)["]',x,re.MULTILINE).group(1), gtf['identifier'].values.tolist()))

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
        This function convert counts to counts per million (CPM)

        :param plot: True if to plot log raw counts and log CPM counts on boxplot.

        :param figsize: Figure size.
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
        This function convert counts to reads per killobase per million (RPKM)
       
        :param plot: True if to plot log raw counts and log RPKM counts on boxplot.

        :param figsize: Figure size.
        """

        df = pd.read_excel(self.countFile)

        countDF = df.set_index(self.geneColumn)

        geneDF = self._getGeneLength()

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
        
        :param plot: True if to plot log raw counts and log TPM counts on boxplot.

        :param figsize: Figure size.
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


    def meanRatioCount(self, plot=True, figsize=(20,10)):

        """
        This function convert counts to medianRatio count
        
        :param plot: True if to plot log raw counts and log TPM counts on boxplot.

        :param figsize: Figure size.
        """

        df = pd.read_excel(self.countFile)
        gene_names = df['Gene']
        df2 = df.set_index('Gene')
        col = df2.columns.tolist()
        df3 = df2[df2[col] !=0]
        sqrt = stats.gmean(df3[col],axis=1)
        
        df4 = df3.div(sqrt,axis=0)
        
        m = df4.median()
        
        res = df2.div(m,axis=1)

        if plot:
            
            logCounts = list(np.log(df2.T+1))

            logNorm_counts = list(np.log(res.T+1))

            fig, ax = boxplot(data=logCounts + logNorm_counts,countType=['Raw counts']*len(df2.columns)+['medianRatio counts']*len(df2.columns),
                         labels= list(df2.columns)+list(res.columns), figsize=figsize )

            ax.set_xlabel('Sample Name')

            ax.set_ylabel('log counts')

        res.insert(0, 'Gene', gene_names)

        
        return res, fig, ax


















  

    

