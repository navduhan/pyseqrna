

from itertools import combinations
import pandas as pd
import numpy as np
import logging
from pyseqrna.pyseqrna_utils import PyseqrnaLogger
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri, numpy2ri, Formula
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr
import matplotlib.pyplot as plt



log = PyseqrnaLogger(mode='a', log="diff")



to_dataframe = robjects.r('function(x) data.frame(x)')

numpy2ri.activate()

def runDESeq2(countDF=None, targetFile=None, design=None,combination=None,  gene_column='Gene', subset=True, lib=None):
    """[summary]

    Args:
        countFile ([type], optional): [description]. Defaults to None.
        targetFile ([type], optional): [description]. Defaults to None.
        design ([type], optional): [description]. Defaults to None.
        combination ([type], optional): [description]. Defaults to None.
        gene_column (str, optional): [description]. Defaults to 'Gene'.
        subset (bool, optional): [description]. Defaults to False.

    Returns:
        [type]: [description]
    """
    try:
        if lib is None:
            deseq = importr('DESeq2')
        else:
            deseq = importr('DESeq2', lib_loc="/home/nav/R/x86_64-pc-linux-gnu-library/4.0")
    except Exception:
        log.error("DESeq2 installation was not found ")
    gene_id = countDF[[gene_column]].values

    countDF.set_index(gene_column, inplace=True)

    deseq_results=pd.DataFrame(gene_id,columns=[gene_column]) # Intialize a dataframe with gene names as first column for deseq2 results

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

def run_edgeR(countDF=None, targetFile=None, combination=None,  gene_column='Gene', subset=False, replicate=True, bcv= 0.4):
    """[summary]

    Args:
        countFile ([type], optional): [description]. Defaults to None.
        targetFile ([type], optional): [description]. Defaults to None.
        design ([type], optional): [description]. Defaults to None.
        combination ([type], optional): [description]. Defaults to None.
        gene_column (str, optional): [description]. Defaults to 'Gene'.
        subset (bool, optional): [description]. Defaults to False.

    Returns:
        [type]: [description]
    """
    try:
        edgeR = importr('edgeR')
        limma = importr('limma')
    except Exception:
        log.error("edgeR installation not found")
    
    gene_id = countDF[[gene_column]].values

    countDF.set_index(gene_column, inplace=True)

    edgeR_results=pd.DataFrame(gene_id,columns=[gene_column]) # Intialize a dataframe with gene names as first column for deseq2 results

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

def degFilter(degDF=None, CompareList=None, FDR=0.05, FOLD=2, plot=True, figsize=(10,6), replicate=True):
    """[summary]

    Args:
        degDF ([type], optional): [description]. Defaults to None.
        CompareList ([type], optional): [description]. Defaults to None.
        FDR (float, optional): [description]. Defaults to 0.05.
        FOLD (int, optional): [description]. Defaults to 2.
        plot (bool, optional): [description]. Defaults to True.
    """
    Up = []
    Down = []
    Total = []
    DEGs = {}
    Ups = {}
    Downs = {}
    summary = pd.DataFrame()

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

        final = pd.concat([upDF, downDF], axis=0)

        DEGs[c] = final
        Ups[c] = upDF
        Downs[c] = downDF

    summary =summary.append(pd.DataFrame({"Comparisons": CompareList, "Total_DEGs": Total, "UP_DEGs": Up, "Down_DEGs": Down}))

    if plot == True:

        category_names= ['UP', 'Down']
        labels= summary['Comparisons'].values.tolist()
        data = np.array(summary[['UP_DEGs', 'Down_DEGs']].values.tolist())
        data_cum = data.cumsum(axis=1)
        category_colors = plt.get_cmap('RdBu_r')(
                np.linspace(0.15, 0.85, data.shape[1]))

        fig, ax = plt.subplots(figsize=figsize)
        fig = plt.gcf()
        fig.set_size_inches(18.5, 10.5)
        ax.invert_yaxis()
        # ax.xaxis.set_visible(False)
        ax.set_xlim(0, np.sum(data, axis=1).max())

        for i, (colname, color) in enumerate(zip(category_names, category_colors)):
            widths = data[:, i]
            
            starts = data_cum[:, i] - widths
      
            ax.barh(labels, widths, left=starts, height=0.8,
                    label=colname, color=color)
            xcenters = starts + widths / 2
            if 0 in xcenters:
                pass
            else:
                r, g, b, _ = color
                text_color = 'white' if r * g * b < 0.5 else 'darkgrey'
                for y, (x, c) in enumerate(zip(xcenters, widths)):
                    ax.text(x, y, str(int(c)), ha='center', va='center',
                            color=text_color)
        ax.legend( bbox_to_anchor=([1.2,0.6,0,0 ]),
                    loc='center right', fontsize='small', frameon=False)
        if replicate:
            ax.set_title(f'Filter DEGs (Fold:{FOLD} and FDR:{FDR})', loc='center')
        else:
            ax.set_title(f'Filter DEGs (Fold:{FOLD} )', loc='center')
        for spine in ax.spines:
            ax.spines[spine].set_visible(False)
        ax.margins(0.01)
        fig.tight_layout()
        
    return {'summary': summary, "filtered": DEGs,"filteredup":Ups, "filtereddown":Downs, "plot": fig}