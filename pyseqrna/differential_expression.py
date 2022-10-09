#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
:Title: This module finds differentially expressed genes from raw read counts

:Created: October 22, 2021

:Author: Naveen Duhan
'''

from itertools import combinations
import pandas as pd
import numpy as np
import logging
import requests
from io import StringIO, TextIOWrapper
from xml.etree import ElementTree
from future.utils import native_str
from pyseqrna.pyseqrna_utils import PyseqrnaLogger
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri, numpy2ri, Formula
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from rpy2.rinterface_lib.callbacks import logger as rpy2_logger

rpy2_logger.setLevel(logging.ERROR)

log = PyseqrnaLogger(mode='a', log="diff")

to_dataframe = robjects.r('function(x) data.frame(x)')

numpy2ri.activate()



def runDESeq2(countDF=None, targetFile=None, design=None,combination=None,  gene_column='Gene',  mmg =False, subset=False, lib=None):
    """
    This function is a wrapper to DESeq2 package in R for differeantial expression analysis from raw read counts.

    :param countFile: Raw read count file. 
    
    :param targetFile: Tab-delimited target file with replication and sample name

    :param design: [description]. Defaults to None.

    :param combination: Comparison list contaning samples to compare.

    :param gene_column: First column in raw read count file. Defaults to 'Gene'.

    :param mmg: True if raw read counts are from multimapped gene groups.

    :param subset: If runDESeq2 subset raw read count according to comparison.

    :param lib: library path of DESeq2 to use.

    :returns: DataFrame

    :rtype: A datafram containing all gene differantial expression in all combinations.
        

    """
    try:
        if lib is None:
            deseq = importr('DESeq2')
        else:
            deseq = importr('DESeq2', lib_loc=lib)
    except Exception:
        log.error("DESeq2 installation was not found ")

    if mmg:
        gene_id = countDF[['MMG', gene_column]].values
        countDF.set_index(['MMG', gene_column], inplace=True)
    else:

        gene_id = countDF[[gene_column]].values
        countDF.set_index(gene_column, inplace=True)

    

    if mmg:

        deseq_results=pd.DataFrame(gene_id,columns=['MMG',gene_column]) # Intialize a dataframe with gene names as first column for deseq2 results

    else: 
        deseq_results=pd.DataFrame(gene_id,columns=[gene_column])

    if subset:

        loopSamples = np.array(range(0,len(combination))) #  if subset is True: Create an array of combination to subset count data

    else:

        loopSamples = np.array([1]) # If subset is False: then consider all count data as single set
    
    # Now loopthrough combination in loopSamples 

    for j in loopSamples:

        if subset:

            c1,c2 = combination[j].split("-")

            subDF = countDF.filter(regex='|'.join([c1,c2] )) # Subset the count data for one combination

            subTF = targetFile[targetFile['sample'].str.contains('|'.join([c1,c2] ))] # Subset the target data for one combination

            with localconverter(robjects.default_converter + pandas2ri.converter):

                count_matrix = robjects.conversion.py2rpy(subDF)

                design_matrix = robjects.conversion.py2rpy(subTF)

            designFormula="~ "+design

            design_formula = Formula(designFormula)

            dds=deseq.DESeqDataSetFromMatrix(countData=count_matrix, colData=design_matrix, design= design_formula)

            comb = [combination[j]] # put the combination in comb for constrast

        else:

            with localconverter(robjects.default_converter + pandas2ri.converter):

                count_matrix = robjects.conversion.py2rpy(countDF)

                design_matrix = robjects.conversion.py2rpy(targetFile)

            designFormula="~ "+design

            design_formula = Formula(designFormula)

            dds=deseq.DESeqDataSetFromMatrix(countData=count_matrix, colData=design_matrix,  design= design_formula)

            comb=combination
    
        dds1 = deseq.DESeq(dds, quiet=True) # run deseq2


        # Itrate through all the given combination for DEG comparision 

        for co in comb:

            c1,c2=co.split("-")

            R_contrast = robjects.vectors.StrVector(np.array([design,c1,c2]))

            result = deseq.results(dds1, contrast=R_contrast)

            result=to_dataframe(result)

            with localconverter(robjects.default_converter + pandas2ri.converter):

                result = robjects.conversion.rpy2py(result)

            result=pd.DataFrame(result)

            result['padj']=result['padj'].replace(np.nan,1)

            result['log2FoldChange'] = result['log2FoldChange'].replace(np.nan, 0)
            result.columns = ['baseMean','logFC','lfcSE','stat','pvalue','FDR']
            result.reset_index(drop=True, inplace=True)

            result.columns=[s+"("+co+")" for s in result.columns]  # Add combination names to the column names 

            deseq_results.reset_index(drop=True, inplace=True)

            deseq_results = pd.concat([deseq_results,result],axis=1)

    return deseq_results

def run_edgeR(countDF=None, targetFile=None, combination=None,  gene_column='Gene',  mmg=False, subset=False, replicate=True, bcv= 0.4, lib=None):
    """
    This function is a wrapper to edgeR package in R for differeantial expression analysis from raw read counts.

    :param countFile: Raw read count file. 
    
    :param targetFile: Tab-delimited target file with replication and sample name

    :param design: [description]. Defaults to None.

    :param combination: Comparison list contaning samples to compare.

    :param gene_column: First column in raw read count file. Defaults to 'Gene'.

    :param mmg: True if raw read counts are from multimapped gene groups.

    :param subset: If runDESeq2 subset raw read count according to comparison.

    :param replicate: False if there are no replicates.

    :param bcv:  Biological coefficient of variation if there are no replicate.

    :param lib: library path of DESeq2 to use.

    :returns: DataFrame

    :rtype: A datafram containing all gene differantial expression in all combinations.
    """
    try:
        if lib is None:
            edgeR = importr('edgeR')
            limma = importr('limma')
        else:
            edgeR = importr('edgeR', lib_loc=lib)
            limma = importr('limma', lib_loc=lib)
            
    except Exception:
        log.error("edgeR installation not found")

    if mmg:
        gene_id = countDF[['MMG', gene_column]].values
        countDF.set_index(['MMG', gene_column], inplace=True)
    else:

        gene_id = countDF[[gene_column]].values
        countDF.set_index(gene_column, inplace=True)

    if mmg:

        edgeR_results=pd.DataFrame(gene_id,columns=['MMG',gene_column]) # Intialize a dataframe with gene names as first column for edgeR results

    else: 
        edgeR_results=pd.DataFrame(gene_id,columns=[gene_column])

    if subset:

        loopSamples = np.array(range(0,len(combination))) #  if subset is True: Create an array of combination to subset count data

    else:

        loopSamples = np.array([1]) # If subset is False: then consider all count data as single set
    
    # Now loopthrough combination in loopSamples 

    for j in loopSamples:
        

        if subset:

            c1,c2 = combination[j].split("-")

            subDF = countDF.filter(regex='|'.join([c1,c2] )) # Subset the count data for one combination

            subTF = targetFile[targetFile['sample'].str.contains('|'.join([c1,c2] ))] # Subset the target data for one combination

            groups = subTF['sample']

            with localconverter(robjects.default_converter + pandas2ri.converter):

                count_matrix = robjects.conversion.py2rpy(subDF)

                group = robjects.conversion.py2rpy(groups)

            dds = edgeR.DGEList(counts=count_matrix,group=group)
    
            dds = edgeR.calcNormFactors(dds)
            
            robjects.r.assign('group', group)

            robjects.r.assign('dds', dds)

            design = robjects.r('design<-model.matrix(~0+dds$samples$group,data=dds$samples)')

            design = pd.DataFrame(design)

            col= robjects.r('levels(dds$samples$group)')

            design.columns= col

            design =design.to_records(index=False)

            cont= robjects.vectors.StrVector([combination[j]])
           
            contrasts = limma.makeContrasts(contrasts=cont,levels=design)

        else:

            counts = np.asarray(countDF,dtype=int)

            cpm = (counts * 1e6) / counts.sum(axis=0) 

            cpm =pd.DataFrame(data=cpm,index=countDF.index,columns=countDF.columns)

            cpm = countDF[cpm.sum(axis=1)>1]

            countDF = countDF[countDF.sum(axis=1)>2]

            groups = targetFile['sample']

            with localconverter(robjects.default_converter + pandas2ri.converter):

                count_matrix = robjects.conversion.py2rpy(countDF)

                group = robjects.conversion.py2rpy(groups)

            dds = edgeR.DGEList(counts=count_matrix,group=group)
    
            dds = edgeR.calcNormFactors(dds)
            
            robjects.r.assign('group', group)

            robjects.r.assign('dds', dds)

            design = robjects.r('design<-model.matrix(~0+dds$samples$group,data=dds$samples)')

            design = pd.DataFrame(design)

            col= robjects.r('levels(dds$samples$group)')

            design.columns= col

            design =design.to_records(index=False)
            
            cont= robjects.vectors.StrVector(combination)
           
            contrasts = limma.makeContrasts(contrasts=cont,levels=design)
        
        # Itrate through all the given combination for DEG comparision 
    

    # counts = np.asarray(countDF,dtype=int)
        if replicate:
            dds = edgeR.estimateGLMCommonDisp(dds, design)

            dds = edgeR.estimateGLMTrendedDisp(dds, design)

            dds = edgeR.estimateGLMTagwiseDisp(dds, design)

            fit = edgeR.glmFit(dds, design)
        else:

            fit = edgeR.glmFit(dds, design, dispersion=float(bcv**2))

        for i in range(len(contrasts.T)):

            lrt = edgeR.glmLRT(fit, contrast=contrasts.T[i])
            
            deg = edgeR.topTags(lrt,countDF.shape[0] )
            
            result=to_dataframe(deg)

            with localconverter(robjects.default_converter + pandas2ri.converter):

                result = robjects.conversion.rpy2py(result)
            
            result=pd.DataFrame(result)

            result.columns = ['logFC','logCPM','LR','pvalue','FDR']

            result.reset_index(drop=True, inplace=True)

            result.columns=[s+"("+combination[i]+")" for s in result.columns]  # Add combination names to the column names 

            edgeR_results.reset_index(drop=True, inplace=True)

            edgeR_results = pd.concat([edgeR_results,result],axis=1)
    
    return edgeR_results

def degFilter(degDF=None, CompareList=None, FDR=0.05, FOLD=2, plot=True, figsize=(10,6), replicate=True, extraColumns=False):
    """
    This function filter all gene expression file based on given FOLD and FDR

    :param degDF: A datafram containing all gene differantial expression in all combinations.

    :param CompareList: A list of all the sample comparison.

    :param FDR: False Discovery Rate for filtering DEGs. Defaults to 0.05.

    :param FOLD: Fold change value. The log2 of the value will be calculated. Defaults to 2.

    :param plot: True if want to plot DEGs per sample on barplot. Defaults to True.
    """
    Up = []
    Down = []
    Total = []
    DEGs = {}
    Ups = {}
    Downs = {}
    # summary = pd.DataFrame()

    if extraColumns:

        degDF = degDF.set_index(['Gene', 'Name', 'Description'])

    else:

        degDF = degDF.set_index('Gene')
        
    for c in CompareList:

        dk = degDF.filter(regex=c, axis=1)

        FDRR = "FDR("+c+")"

        LFC = "logFC("+c+")"

        if replicate:

            fdr = dk[dk[FDRR]<=FDR].dropna()

            upDF = fdr[fdr[LFC]>=np.log2(FOLD)]

            downDF = fdr[fdr[LFC]<=-np.log2(FOLD)]
        else:
            upDF = dk[dk[LFC]>=np.log2(FOLD)]

            downDF = dk[dk[LFC]<=-np.log2(FOLD)]

        Up.append(upDF.shape[0])

        Down.append(downDF.shape[0])

        Total.append(upDF.shape[0]+downDF.shape[0])
        
        upDF =upDF.reset_index()
        downDF =downDF.reset_index()
        final = pd.concat([upDF, downDF], axis=0)
      
        DEGs[c] = final
        Ups[c] = upDF
        Downs[c] = downDF

    summary =pd.DataFrame({"Comparisons": CompareList, "Total_DEGs": Total, "Up_DEGs": Up, "Down_DEGs": Down})

    if plot == True:

        category_names= ['UP', 'Down']
        labels= summary['Comparisons'].values.tolist()

        if len(labels)>10:
            hg = 15
        else:
            hg = 10
            

        updata= summary['Up_DEGs'].values.tolist()
        downdata = summary['Down_DEGs'].values.tolist()
        my_range=list(range(1,len(summary.index)+1))
        fig, ax = plt.subplots(figsize=(hg,hg),dpi=300)
        ax.barh(labels,updata, color='mediumseagreen')
        ax.barh(labels,downdata, left=updata, color='salmon')
        plt.xticks(fontsize=10)
        plt.yticks(fontsize=10)
        plt.xlabel("Number of Genes", fontsize=10)
        plt.ylabel("Comparisons", fontsize=10)
        plt.legend(['Up-regulated', 'Down-regulated'],  ncol =2, loc='center', bbox_to_anchor=(0.5, 1.1))


        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_bounds((0, len(my_range)))
        # add some space between the axis and the plot
        ax.spines['left'].set_position(('outward', 8))
        ax.spines['bottom'].set_position(('outward', 5))
        # plt.savefig('deg.png', dpi=300, bbox_inches='tight')
        if replicate:
            plt.title(f'Filter DEGs (Fold:{FOLD} and FDR:{FDR})', loc='center', pad=40)
        else:
            plt.title(f'Filter DEGs (Fold:{FOLD} )', loc='center', pad=40)
        fig.tight_layout()
        
    return {'summary': summary, "filtered": DEGs,"filteredup":Ups, "filtereddown":Downs, "plot": fig}

class Gene_Description:

    """
    This class fetch gene name and description for genes.
    """

    def __init__(self, species, type, combinations=None, degFile= None, filtered=True) :
        self.species = species
        self.type = type
        self.combinations = combinations
        self.filtered = filtered
        self.degFile = degFile
        self.names = self._query(species=self.species, type=self.type)
        
        return


    def _get_request(self, url,  **params):

        if params:
            r = requests.get(url, params=params, stream=True)
        else:
            r = requests.get(url)
        r.raise_for_status()

        return r

    def _add_attr_node(self, root, attr):
        attr_el = ElementTree.SubElement(root, 'Attribute')
        attr_el.set('name', attr)

        return


    def _query(self, species, type):

        # first need to check if the species is animal or plant

        if type == 'animals':
            uri= "https://ensembl.org/biomart/martservice"
            scheme = 'default'
            fspecies = species+"_gene_ensembl"
        if type == 'plants':
            uri = "https://plants.ensembl.org/biomart/martservice"
            scheme = 'plants_mart'
            fspecies = species+"_eg_gene"

        # build query

        root = ElementTree.Element('Query')
        root.set('virtualSchemaName', scheme)
        root.set('formatter', 'TSV')
        root.set('header', '1')
        root.set('uniqueRows', native_str(int(True)))
        root.set('datasetConfigVersion', '0.6')

        dataset = ElementTree.SubElement(root, 'Dataset')
        dataset.set('name', fspecies)
        dataset.set('interface', 'default')
        attributes = ["ensembl_gene_id", "external_gene_name","description"]
        for attr in attributes:
            self._add_attr_node(dataset, attr)

        response = self._get_request(
        uri , query = ElementTree.tostring(root))
        result = pd.read_csv(StringIO(response.text), sep='\t')
        result.columns = ['Gene', 'Name','Description']
        
        return result


    def add_names(self):
        """
        This function add gene name and description in DEGs.

        :returns: DataFrame
        :rtype: DEGs file with gene name and description.
        """
        

        if self.filtered:

            file = self.degFile.split("_")[0]

            wd = pd.ExcelWriter(f"{file}_DEGs.xlsx")

            for c in self.combinations:

                df = pd.read_excel(self.degFile, sheet_name=c)

                col =df.columns

                col = col.insert(1, 'Name')

                col = col.insert(2, 'Description')

                final = df.merge(self.names, on='Gene')

                final = final[col]

                final.to_excel(wd, sheet_name=c, index=False)

            wd.save()
        else:

            if isinstance(self.degFile, pd.DataFrame):

                deg = self.degFile

            else:

                deg = pd.read_excel(self.degFile)
            
            col =deg.columns

            col = col.insert(1, 'Name')

            col = col.insert(2, 'Description')

            final = deg.merge(self.names, on='Gene')

            final = final[col]

            print("please provide the filtered DEGs file")

        return final

    def add_names_annotation(self, file):

        """
        This function add gene name and description in functional annotation.

        :returns: DataFrame
        :rtype: Functional annotation file with gene name and description.
        """

        pl = self.names.values.tolist()
        pp = {}
        for p in pl:
            pp[p[0]]=p[1]
        try:
            gene_names =[]
            for i,row in file.iterrows():
                Genes = str(file.at[i,'Genes']).split(",")
                result = []
                for gene in Genes:
                    result.append(pp[gene].upper())
                res = ",".join(result)

                gene_names.append([res])
                
            file['Gene_Name'] = gene_names


            for i, row in file.iterrows():
                
                res = ','.join(file.at[i,'Gene_Name'])

                file.at[i, 'Gene_Name'] = res

        except Exception:
            file = file

        return file
       



        