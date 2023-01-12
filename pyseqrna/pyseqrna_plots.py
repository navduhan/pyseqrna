#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
:Title: This module function generate visualization 

:Created : November 2, 2021

:Author : Naveen Duhan
"""


import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import seaborn as sns
import matplotlib.patches as patches
from itertools import chain

def plotVolcano(degDF=None, comp=None,FOLD=2,pValue=0.05,color=('red','grey','green'), dim=(8,5), dotsize=8, markerType='o', alpha=0.5):
    """
    This function plots a Volcano plot. 

    :param degDF: All gene expression file.

    :param comp: Sample comparison.

    :param FOLD: FOLD change. Defaults to 2.

    :param pValue: Pvalues. Defaults to 0.05.

    :param color: Colors to be used in plot. Defaults to ('red','grey','green').

    :param dim: Dimensions of the plot. Defaults to (8,5).

    :param dotsize: Dotsize on plot. Defaults to 8.

    :param markeType: Shape to use. Defaults to 'o'.

    :param alpha: Transparency of plot. Defaults to 0.5.
    """

    PVAL = "pvalue("+comp+")"

    LFC = "logFC("+comp+")"
    _y = r'$ log_{2}(Fold Change)$'
    _x = r'$ -log_{10}(P-value)$'

    dk = degDF.filter(regex=comp, axis=1)
    
    try:
        
        final = dk[dk[PVAL]<=pValue].copy()

        final.loc[final[LFC]>=np.log2(FOLD),'colorADD'] = color[2]

        final.loc[final[LFC]<=-np.log2(FOLD), 'colorADD'] = color[0]

        final['colorADD'].fillna(color[1], inplace=True)  

        final['log(10)_pvalue'] = -(np.log10(final[PVAL]))
        
        color_values = {col: i for i, col in enumerate(color)}

        color_num = [color_values[i] for i in final['colorADD']]
        
        legendlabels=['Significant down', 'Not significant', 'Significant up']  

        fig, ax = plt.subplots(figsize=dim)

        a = ax.scatter(final[LFC], final['log(10)_pvalue'], c=color_num, cmap=ListedColormap(color), alpha=alpha,s=dotsize, marker=markerType)
        ax.legend(handles=a.legend_elements()[0], labels=legendlabels,   bbox_to_anchor=(1.46,1,0,0))
        ax.set_xlabel(_x)
        ax.set_ylabel(_y)
        fig.tight_layout()

    except Exception:
        return 'No Volcano'

    return fig 


def plotMA(degDF=None, countDF=None, comp=None,FOLD=2, FDR=0.05, color=('red','grey','green'), dim=(8,5), dotsize=8, markerType='o', alpha=0.5):
    """
    This function plots a MA plot. 

    :param degDF: All gene expression file.

    :param comp: Sample comparison.

    :param FOLD: FOLD change. Defaults to 2.

    :param FDR: FDR value. Defaults to 0.05.

    :param color: Colors to be used in plot. Defaults to ('red','grey','green').

    :param dim: Dimensions of the plot. Defaults to (8,5).

    :param dotsize: Dotsize on plot. Defaults to 8.

    :param markeType: Shape to use. Defaults to 'o'.

    :param alpha: Transparency of plot. Defaults to 0.5.
    """

    LFC = "logFC("+comp+")"
    PVAL = "FDR("+comp+")"
    # baseMean ="baseMean("+comp+")"
    _y = r'$ log_{2}(Fold Change)$'
    _x = r'$ log_{2}(Mean Count)$'

    dk = degDF.filter(regex=comp, axis=1)

    if dk.empty:
        pass
    else:
    

        cdf1 = countDF.filter(regex=comp.split("-")[0], axis=1)
        cdf2 = countDF.filter(regex=comp.split("-")[1], axis=1)

        cdf_mean1=cdf1.mean(axis=1)
        cdf_mean2=cdf2.mean(axis=1)

        counts = pd.concat([cdf_mean1,cdf_mean2], axis=1)

        counts.columns = ['first', 'second']

        final = pd.concat([dk,counts], axis=1)

        final = final.fillna(0)

        final.loc[(final[LFC]>=np.log2(FOLD)) & (final[PVAL]<=FDR),'colorADD'] = color[2]

        final.loc[(final[LFC]<=-np.log2(FOLD)) & (final[PVAL]<=FDR), 'colorADD'] = color[0]

        final['colorADD'].fillna(color[1], inplace=True)  

        final['mean'] = final[['first', 'second']].mean(axis=1)

        finalDF = final.loc[final['mean']>0].copy()

        np.seterr(divide = 'ignore') 

        finalDF['log(2)_Mean'] = np.log2(finalDF['first']) + np.log2(finalDF['second']) / 2
        
        color_values = {col: i for i, col in enumerate(color)}

        color_num = [color_values[i] for i in finalDF['colorADD']]
        
        legendlabels=['Significant down', 'Not significant', 'Significant up']  

        fig,ax = plt.subplots(figsize=dim)

        a = ax.scatter(finalDF['log(2)_Mean'], finalDF[LFC], c=color_num, cmap=ListedColormap(color),
                            alpha=alpha, s=dotsize, marker=markerType)
        try:
            ax.legend(handles=a.legend_elements()[0], labels=legendlabels, bbox_to_anchor=(1.46,1,0,0))
        except ValueError:  #raised if `y` is empty.
            pass
        plt.axhline(y=0, color='#7d7d7d', linestyle='--')
        
        ax.set_xlabel(_x)
        ax.set_ylabel(_y)

        fig.tight_layout()
 
   
    return fig 




def plotHeatmap(degDF= None, combinations=None, num=50, figdim=(12,10), type='counts'):

    """
    This function plots a heatmap based on FOLD change or counts.

    :param degDF: All gene expression file or Counts file.

    :param combinations: All sample combinations. 

    :param num: Total number of genes to plot. Default is 50 (25 up and 25 down)

    :param figdim: Figure dimensions.

    :param type: Heatmap to create on counts/degs. 
    """

    # Need to add argument for all deg hetamp.
    degDF = degDF.set_index('Gene')
    if type.lower() =='degs':
        com2=[]
        labelsd = []
        for i in combinations:
            fc="logFC("+i+")"
            com2.append(fc)
            labelsd.append(i)
        topGene=degDF.nlargest(int(num/2),com2)
        botGene=degDF.nsmallest(int(num/2),com2)
        final=pd.concat([topGene,botGene])
        fin = final[com2]
        fin.columns= labelsd
    elif type.lower() =='counts':
        fin=degDF.nlargest(int(num), degDF.columns)
        
    fig, ax = plt.subplots(figsize=figdim)

    sns.heatmap(fin, cmap="seismic",ax=ax)

    fig.tight_layout()

    return fig,ax


defaultColors = [
    # r, g, b, a
    [92, 192, 98, 0.5],
    [90, 155, 212, 0.5],
    [241, 90, 96, 0.4],
    [255, 255,102,0.3],
    [255, 117, 0, 0.3],
]
defaultColors = [
    [i[0] / 255.0, i[1] / 255.0, i[2] / 255.0, i[3]]
    for i in defaultColors
]
def _insertEllipse(fig, ax, x, y, w, h, a,  fillcolor):
    e = patches.Ellipse(xy=(x, y),width=w, height=h,angle=a,
        fill="blue",linewidth=2, color=fillcolor)
    ax.add_patch(e)
def _insertText(fig, ax, x, y, text, fontsize=None, col="black", ha="center", va="center", fontweight= 600):
    ax.text(x, y, text, horizontalalignment=ha,
        verticalalignment=va,fontsize=fontsize, fontweight=fontweight,
        color=col)

def _GenerateCollection(data=None):
    N = len(data)
    GenesData = [set(data[i]) for i in range(N)]  # sets for separate groups
    AllGenes = set(chain(*data))                     # union of all sets
    GeneCollections = {}
    for n in range(1, 2**N):
        key = bin(n).split('0b')[-1].zfill(N)
        value = AllGenes
        sets_for_intersection = [GenesData[i] for i in range(N) if  key[i] == '1']
        sets_for_difference = [GenesData[i] for i in range(N) if  key[i] == '0']
        for s in sets_for_intersection:
            value = value & s
        for s in sets_for_difference:
            value = value - s
        GeneCollections[key] = value
    return GeneCollections

def plotVenn(DEGFile=None, FOLD=2, comparisons=None, degLabel="",  fontsize=14, figsize=(12,12),dpi=300):
    """
    This function plots a Venn diagram for filtered degs in samples.

    :param DEGFile: Filtered deg excel file containg samples sheet-wise.

    :param FOLD: FOLD change. Defaults to 2.

    :param comparisons: Comparison list. Defaults to None.

    :param degLabel: How to put labes either total/ up-down. Defaults to "" i.e. up-down.

    :param fontsize: Font size. Defaults to 14.
    
    :param figsize: Figure size. Defaults to (12,12).

    :param dpi: Figure DPI resolution. Defaults to 300.
    """
    global labelsUp, labelsDown, labels, fig
    data = []
    Up = []
    Down = []

    for com in comparisons:
        df = pd.read_excel(DEGFile, sheet_name=com)
        id = "logFC("+com+")"
        down = df[df[id] <= -np.log(FOLD)]
        up = df[df[id] >= np.log(FOLD)]
        data.append(df['Gene'])
        Up.append(up['Gene'])
        Down.append(down['Gene'])

    if degLabel == "total":
        GeneCollections = _GenerateCollection(data)

        labels = {k: "" for k in GeneCollections}

        for k in GeneCollections:
            labels[k] += str(len(GeneCollections[k]))
        
    else:
        upGene = _GenerateCollection(Up)

        downGene = _GenerateCollection(Down)

        labelsUp = {k: "" for k in upGene}

        for k in upGene:
            labelsUp[k] += str(len(upGene[k]))
        labelsDown = {k: "" for k in downGene}
        for k in upGene:
            labelsDown[k] += str(len(downGene[k]))

    fig = plt.figure(0, figsize=figsize, dpi=dpi)
    ax = fig.add_subplot(111, aspect='equal')
    ax.set_axis_off()
    ax.set_ylim(bottom=0.0, top=1.0)
    ax.set_xlim(left=0.0, right=1.0)

    if len(comparisons)== 4:
        colors = [defaultColors[i] for i in range(4)]

        ellipsePoint =[(0.350, 0.400, 0.72, 0.45, 140.0),(0.450, 0.500, 0.72, 0.45, 140.0),(0.544, 0.500, 0.72, 0.45, 40.0),(0.644, 0.400, 0.72, 0.45, 40.0)]
        for e, c in zip(ellipsePoint, colors):

            _insertEllipse(fig, ax, e[0],e[1],e[2],e[3],e[4], c)

        # dataPoints=[(0.85,0.42),(0.68, 0.72),(0.77, 0.59),(0.32, 0.72),(0.71, 0.30),(0.50, 0.66),(0.65, 0.50),(0.14, 0.42),
        #             (0.50, 0.17),(0.29, 0.30),(0.37, 0.26),(0.23, 0.59),(0.63, 0.26),(0.35, 0.50),(0.50, 0.38)]
        dataPoints=[(0.85,0.42),(0.66, 0.72),(0.77, 0.59),(0.32, 0.72),(0.69, 0.30),(0.50, 0.66),(0.65, 0.50),(0.14, 0.42),
                    (0.50, 0.17),(0.30, 0.30),(0.38, 0.25),(0.23, 0.59),(0.63, 0.26),(0.35, 0.50),(0.50, 0.38)]
        labelPoints=[('0001', ''),( '0010', ''),( '0011', ''),( '0100', ''),( '0101', ''),( '0110', ''),( '0111', ''),
                     ( '1000', ''),( '1001', ''),( '1010', ''),( '1011', ''),( '1100', ''),( '1101', ''),( '1110', ''),( '1111', '')]

        for d, l in zip(dataPoints,labelPoints):
            if degLabel != "total":
                _insertText(fig, ax, d[0], d[1], labelsUp.get(l[0],''),col="blue", fontsize=fontsize, fontweight=600)
                _insertText(fig, ax, d[0], (d[1]-0.03), labelsDown.get(l[0],''), col="red", fontsize=fontsize, fontweight=600)
            else:
                _insertText(fig, ax, d[0], d[1], labels.get(l[0], ''), col="blue", fontsize=fontsize, fontweight=600)
        # legend
        _insertText(fig, ax, 0.13, 0.18, comparisons[0],  fontsize=fontsize,fontweight=600,  ha="right")
        _insertText(fig, ax, 0.18, 0.83, comparisons[1],  fontsize=fontsize,fontweight=600,  ha="right", va="bottom")
        _insertText(fig, ax, 0.82, 0.83, comparisons[2],  fontsize=fontsize,fontweight=600,  ha="left", va="bottom")
        _insertText(fig, ax, 0.87, 0.18, comparisons[3],  fontsize=fontsize, fontweight=600, ha="left", va="top")

    elif len(comparisons)== 3:
        colors = [defaultColors[i] for i in range(3)]
        ellipsePoint = [( 0.333, 0.633, 0.5, 0.5, 0.0),(0.666, 0.633, 0.5, 0.5, 0.0),(0.500, 0.310, 0.5, 0.5, 0.0)]
        for e, c in zip(ellipsePoint, colors):
            _insertEllipse(fig, ax, e[0],e[1],e[2],e[3],e[4], c)
        dataPoints = [(0.50, 0.27), (0.73, 0.65), (0.61, 0.46), (0.27, 0.65), (0.39, 0.46), ( 0.50, 0.65), (0.50, 0.51)]
        labelPoints = [('001', ''), ('010', ''), ('011', ''), ('100', ''), ('101', ''), ('110', ''), ('111', '')]

        for d, l in zip(dataPoints, labelPoints):
            if degLabel != "total":
                _insertText(fig, ax, d[0], d[1], labelsUp.get(l[0], ''), col="blue", fontsize=fontsize, fontweight='normal')
                _insertText(fig, ax, d[0], (d[1] - 0.03), labelsDown.get(l[0], ''), col="red", fontsize=fontsize, fontweight='normal')
            else:
                _insertText(fig, ax, d[0], d[1], labels.get(l[0], ''), col="blue", fontsize=fontsize, fontweight='normal')
        # legend
        _insertText(fig, ax, 0.15, 0.87, comparisons[0], fontsize=fontsize, fontweight=600, ha="right", va="bottom")
        _insertText(fig, ax, 0.85, 0.87, comparisons[1], fontsize=fontsize, fontweight=600, ha="left", va="bottom")
        _insertText(fig, ax, 0.50, 0.02, comparisons[2], fontsize=fontsize, fontweight=600, ha="left", va="top")

    elif len(comparisons)== 2:
        colors = [defaultColors[i] for i in range(2)]
        ellipsePoint = [( 0.375, 0.3, 0.5, 0.5, 0.0),(0.625, 0.3, 0.5, 0.5, 0.0)]
        for e, c in zip(ellipsePoint, colors):
            _insertEllipse(fig, ax, e[0],e[1],e[2],e[3],e[4], c)
        dataPoints = [(0.74, 0.30), (0.26, 0.30), (0.50, 0.30)]
        labelPoints = [('01', ''), ('10', ''), ('11', '')]

        for d, l in zip(dataPoints, labelPoints):
            if degLabel != "total":
                _insertText(fig, ax, d[0], d[1], labelsUp.get(l[0], ''), col="blue", fontsize=fontsize)
                _insertText(fig, ax, d[0], (d[1] - 0.03), labelsDown.get(l[0], ''), col="red", fontsize=fontsize)
            else:
                _insertText(fig, ax, d[0], d[1], labels.get(l[0], ''), col="blue", fontsize=fontsize)
        # legend
        _insertText(fig, ax, 0.20, 0.56, comparisons[0], fontsize=fontsize, ha="right", va="bottom")
        _insertText(fig, ax, 0.80, 0.56, comparisons[1], fontsize=fontsize, ha="left", va="bottom")

    else:
        print("Please provide combination between 2-4")

    return fig


