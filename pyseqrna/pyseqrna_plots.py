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
import matplotlib.lines as mlines
from matplotlib.colors import ListedColormap
import seaborn as sns
import matplotlib.patches as patches
from itertools import chain
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

def plotVolcano(degDF=None, comp=None, FOLD=2, pValue=0.05, color=('red', 'grey', 'green'), 
                dim=(8, 5), font=14, dotsize=10, markerType='o', alpha=0.5):
    """
    This function plots a Volcano plot.

    :param degDF: DataFrame containing differential expression results.
    :param comp: Sample comparison (string).
    :param FOLD: Fold change threshold. Defaults to 2.
    :param pValue: P-values threshold. Defaults to 0.05.
    :param color: Colors to be used in the plot. Defaults to ('red', 'grey', 'green').
    :param dim: Dimensions of the plot. Defaults to (8, 5).
    :param dotsize: Dot size on the plot. Defaults to 10.
    :param markerType: Shape to use for markers. Defaults to 'o'.
    :param alpha: Transparency of the plot. Defaults to 0.5.
    :return: The figure object.
    """
    
    PVAL = "pvalue(" + comp + ")"
    LFC = "logFC(" + comp + ")"
    _x = r'$ log_{2}(Fold Change)$'
    _y = r'$ -log_{10}(P-value)$'
    
    # Filter differential expression results based on comparison
    dk = degDF.filter(like=comp, axis=1)

    if dk.empty:
        print("No data found for the specified comparison.")
        return None

    try:
        final = dk[dk[PVAL] <= pValue].copy()

        final.loc[final[LFC] >= np.log2(FOLD), 'colorADD'] = color[2]  # Significant up
        final.loc[final[LFC] <= -np.log2(FOLD), 'colorADD'] = color[0]  # Significant down
        final.fillna({'colorADD': color[1]}, inplace=True)  # Not significant

        final['log(10)_pvalue'] = -(np.log10(final[PVAL]))

        # Map colors to numerical values
        color_values = {col: i for i, col in enumerate(color)}
        color_num = [color_values[i] for i in final['colorADD']]

        # Create the plot
        fig, ax = plt.subplots(figsize=dim)
        scatter = ax.scatter(final[LFC], final['log(10)_pvalue'], c=color_num, cmap=ListedColormap(color), 
                             alpha=alpha, s=dotsize, marker=markerType)

        # Create legend handles
        handles, _ = scatter.legend_elements()
        unique_colors = list(set(color_num))
        legend_labels = [label for i, label in enumerate(['Significant down', 'Not significant', 'Significant up']) if i in unique_colors]

        if len(handles) == len(legend_labels):
            ax.legend(handles, legend_labels, bbox_to_anchor=(1.46, 1, 0, 0))
        else:
            print("Mismatch in handles and labels, adjusting legend.")

        ax.set_xlabel(_x, fontsize=font, fontweight='bold')
        ax.set_ylabel(_y, fontsize=font, fontweight='bold')
        fig.tight_layout()

        return fig

    except Exception as e:
        print(f"An error occurred: {e}")
        return 'No Volcano'


def plotMA(degDF=None, countDF=None, comp=None, FOLD=2, FDR=0.05, 
           color=('red', 'grey', 'green'), font=14, dim=(8, 5), 
           dotsize=8, markerType='o', alpha=0.5):
    """
    This function plots a MA plot. 

    :param degDF: DataFrame containing differential expression results.
    :param countDF: DataFrame containing raw counts.
    :param comp: Sample comparison (string).
    :param FOLD: Fold change threshold. Defaults to 2.
    :param FDR: False discovery rate threshold. Defaults to 0.05.
    :param color: Colors to be used in plot. Defaults to ('red', 'grey', 'green').
    :param dim: Dimensions of the plot. Defaults to (8, 5).
    :param dotsize: Dot size on plot. Defaults to 8.
    :param markerType: Shape to use for markers. Defaults to 'o'.
    :param alpha: Transparency of plot. Defaults to 0.5.
    :return: The figure object.
    """
    
    LFC = "logFC(" + comp + ")"
    PVAL = "FDR(" + comp + ")"
    _y = r'$ log_{2}(Fold Change)$'
    _x = r'$ log_{2}(Mean Count)$'

    try:
        # Filter the differential expression DataFrame
        dk = degDF.filter(regex=comp, axis=1)

        if dk.empty:
            print("No data found for the specified comparison.")
            return None

        cdf1 = countDF.filter(regex=comp.split("-")[0], axis=1)
        cdf2 = countDF.filter(regex=comp.split("-")[1], axis=1)

        cdf_mean1 = cdf1.mean(axis=1)
        cdf_mean2 = cdf2.mean(axis=1)

        counts = pd.concat([cdf_mean1, cdf_mean2], axis=1)
        counts.columns = ['first', 'second']

        final = pd.concat([dk, counts], axis=1)
        final = final.fillna(0)

        # Color assignment based on thresholds
        final.loc[(final[LFC] >= np.log2(FOLD)) & (final[PVAL] <= FDR), 'colorADD'] = color[2]
        final.loc[(final[LFC] <= -np.log2(FOLD)) & (final[PVAL] <= FDR), 'colorADD'] = color[0]
        final.fillna({'colorADD': color[1]}, inplace=True)

        final['mean'] = final[['first', 'second']].mean(axis=1)
        finalDF = final.loc[final['mean'] > 0].copy()

        np.seterr(divide='ignore')

        # Calculate log2 mean
        finalDF['log(2)_Mean'] = (np.log2(finalDF['first']) + np.log2(finalDF['second'])) / 2

        # Map colors to numerical values
        color_values = {col: i for i, col in enumerate(color)}
        color_num = [color_values[i] for i in finalDF['colorADD']]

        # Create the plot
        fig, ax = plt.subplots(figsize=dim)
        a = ax.scatter(finalDF['log(2)_Mean'], finalDF[LFC], c=color_num, cmap=ListedColormap(color),
                       alpha=alpha, s=dotsize, marker=markerType)

        # Create legend handles
        legend_elements = []
        unique_colors = list(set(color_num))

        if 0 in unique_colors:
             legend_elements.append(mlines.Line2D([0], [0], marker=markerType, color='w', label='Significant down',
                          markerfacecolor=color[0], markersize=dotsize))

        if 1 in unique_colors:
             legend_elements.append(mlines.Line2D([0], [0], marker=markerType, color='w', label='Not significant',
                          markerfacecolor=color[1], markersize=dotsize))

        if 2 in unique_colors:
             legend_elements.append(mlines.Line2D([0], [0], marker=markerType, color='w', label='Significant up',
                          markerfacecolor=color[2], markersize=dotsize))

        if legend_elements:
            ax.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left')

        plt.axhline(y=0, color='#7d7d7d', linestyle='--')
        ax.set_xlabel(_x, fontsize=font, fontweight='bold')
        ax.set_ylabel(_y, fontsize=font, fontweight='bold')
        fig.tight_layout()

        return fig

    except Exception as e:
        print(f"An error occurred: {e}")
        return None

def plotHeatmap(degDF=None, combinations=None, num=50, figdim=(15, 10), extraColumns=False, type='counts',
                color_map='seismic', annot=False, scale=True, rowclus=True, colclus=True, zscore=None,
                xlabel=True, ylabel=True, tickfont=(10, 10), theme=None):

    """
    This function plots a clustered heatmap of the top expressed and downregulated genes.
    """

    if extraColumns:
         degDF = degDF.drop(columns=['Name', 'Description'])
         degDF = degDF.set_index('Gene')
    else:
        if 'Gene' in degDF.columns:
            degDF = degDF.set_index('Gene')

    heatmap_data = pd.DataFrame()

    if type.lower() == 'degs':
        degDF = degDF.fillna(0)
        
        if theme == 'dark':
            plt.style.use('dark_background')
        else:
            plt.style.use('default')

        top_genes = set()
        
        for combination in combinations:
            col_name = f'logFC({combination})'
            if col_name in degDF.columns:
                top_expressed = degDF.nlargest(int(num), col_name)
                top_downregulated = degDF.nsmallest(int(num), col_name)
                top_genes.update(top_expressed.index.tolist())
                top_genes.update(top_downregulated.index.tolist())

        if not top_genes:
            return None, None

        # Filter data
        valid_combinations = []
        for c in combinations:
            col_name = f'logFC({c})'
            if col_name in degDF.columns:
                heatmap_data[c] = degDF.loc[list(top_genes), col_name]
                valid_combinations.append(c) # keep original name for plot columns
        
        # Rename columns to just combination name (removing logFC wrapper) if desired, 
        # but user code kept 'c' as key which presumably works if we assign straight to it.
        # Actually user code: heatmap_data[c] = ... so column name is 'c'.
        
    elif type.lower() == 'counts':
        
        degDF = degDF.replace(0, np.nan).dropna(axis=1, how="all")
        
        # Logic from old code for counts
        # It used degDF.columns because counts df only has sample columns (after index set)
        
        topGene = degDF.nlargest(int(num), degDF.columns) # Old code used num/2 for up/down but this is counts? 
        # Old code: topGene=degDF.nlargest(int(num/2),degDF.columns)
        #           botGene=degDF.nsmallest(int(num/2),degDF.columns)
        #           final=pd.concat([topGene,botGene])
        # Preserving that logic:
        
        topGene = degDF.nlargest(int(num/2), degDF.columns)
        botGene = degDF.nsmallest(int(num/2), degDF.columns)
        heatmap_data = pd.concat([topGene, botGene])
        
        # Should we remove duplicates if any?
        heatmap_data = heatmap_data[~heatmap_data.index.duplicated(keep='first')]

    if heatmap_data.empty:
        return None, None

    # Dynamic figure size (optional override or usage)
    # The user snippet calculated figdim but we have a default/arg. 
    # Let's use the arg figdim if provided, or dynamic?
    # User's logic:
    # fig_width = max(10, len(heatmap_data.columns) * 1.2)
    # fig_height = max(10, len(heatmap_data) * 0.2)
    # figdim = (fig_width, fig_height)
    # But function signature has figdim arg. I will honor the arg if it matches old behavior, 
    # but the user code snippet *ignored* the input figdim and recalculated it.
    # To be safe and "according to this", I'll use the dynamic calc as default if it looks reasonable?
    # Since existing calls expect a specific size, maybe I stick to figdim?
    # The user passed figdim=(15,10) in their snippet but then recalculated it.
    # I'll stick to dynamic calc as it's often better for heatmaps.

    fig_width = max(10, len(heatmap_data.columns) * 1.2)
    fig_height = max(10, len(heatmap_data) * 0.2)
    figdim = (fig_width, fig_height)

    # Plotting
    if rowclus or colclus:
        # zscore: 0 for rows, 1 for columns, None for no scaling. 
        # sns.clustermap z_score param takes 0 or 1.
        # The user snippet passed zscore=None.
        
        g = sns.clustermap(
            heatmap_data,
            row_cluster=rowclus,
            col_cluster=colclus,
            cmap=color_map,
            annot=annot,
            cbar=scale,
            z_score=zscore,
            xticklabels=xlabel,
            yticklabels=ylabel,
            figsize=figdim,
            linewidths=0.5,
            linecolor='black' # removed theme support as requested
        )
        return g, None

    else:
        fig, ax = plt.subplots(figsize=figdim, dpi=300)
        sns.heatmap(
            heatmap_data,
            ax=ax,
            cmap=color_map,
            annot=annot,
            cbar=scale,
            linewidths=0.5,
            linecolor='black'
        )
        return fig, ax


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


def pcaPlot(ncountdf=None ,legends=False, fontsize=14, figsize=(12,12),dpi=300):

    """This function plots PCA plot based on normalized counts

    :param ncountdf: Normalized counts File.

    :param legenes: show legends.

    :param fontsize: Font size. Defaults to 14.
    
    :param figsize: Figure size. Defaults to (12,12).

    :param dpi: Figure DPI resolution. Defaults to 300.

    """

    # Transpose the normalize counts to make sample as rows and genes as columns

    ncountdf = ncountdf.set_index(['Gene'])

    ndf = ncountdf.T

    pca = make_pipeline(StandardScaler(), PCA())
    # Perform PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(ndf)

    # Create a DataFrame with the principal components
    pca_df = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'], index=ndf.index)

    # Assign unique colors to each sample using seaborn color palette
    colors = sns.color_palette("husl", n_colors=len(pca_df.index))

    # Plot the PCA results with sample labels and color
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    legend_labels = []

    for i, sample in enumerate(pca_df.index):
        scatter = ax.scatter(pca_df['PC1'].iloc[i], pca_df['PC2'].iloc[i], color=colors[i])
        legend_labels.append(sample)

    # Add labels to data points with adjustText
    texts = [plt.text(pca_df['PC1'].iloc[i], pca_df['PC2'].iloc[i], sample, ha='right', va='bottom') for i, sample in enumerate(pca_df.index)]

    # Adjust text positions to avoid overlap
    adjust_text(texts, arrowprops=dict(arrowstyle="-", color='k', lw=0.5))

    if legends:
        # Organize legends in columns with a maximum of 10 legends per column
        num_legends_per_column = 15
        num_columns = len(legend_labels) // num_legends_per_column + 1

        legend = ax.legend(legend_labels, loc='center left', bbox_to_anchor=(1, 0.5),
                        title='Sample Legends', ncol=num_columns)

    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


    ax.xlabel('Principal Component 1 (PC1)', fontsize = fontsize)
    ax.ylabel('Principal Component 2 (PC2)', fontsize = fontsize)

    return fig
