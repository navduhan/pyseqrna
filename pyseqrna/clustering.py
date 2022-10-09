#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
:Title: This module function generate a similarity dendrogram between samples

:Created : August 2, 2021

:Author : Naveen Duhan
"""

from itertools import count
import seaborn as sns
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as shc

def leaf_label(temp):
    """Funtion to generate leaf label for dendrogram

    :param temp: Temp leaves arrangment for dendrogram

    :return: Return leaf labels for dendrogram.
    """
    label =[]
    for i in range(len(temp)):
        label.append(temp[i])
    return label
def clusterSample(countDF = None):
    """Function to cluster samples based on similarity

    :param countDF: A dataframe of read counts or normalized read counts.

    :return: A clustered dendrogram of samples
    """
    if isinstance(countDF, pd.DataFrame):
        counts = countDF
    if not isinstance(countDF, pd.DataFrame):
        if countDF.endswith(".xlsx"):
            counts = pd.read_excel(countDF)
          
        if countDF.endswith(".csv"):
            counts = pd.read_csv(countDF)
        if countDF.endswith(".txt"):
            counts = pd.read_csv(countDF, sep="\t")
   
    counts = counts.set_index('Gene') 
    corrCount = counts.corr()
    linked = shc.linkage(corrCount, 'ward')
    R = shc.dendrogram(
                    linked,
                    truncate_mode='lastp',  # show only the last p merged clusters
                    p=len(counts.columns),  # show only the last p merged clusters
                    no_plot=True,
                    )

    temp = {R["leaves"][ii]: counts.columns[ii] for ii in range(len(R["leaves"]))}
    
    height = len(counts.columns)/2
    fig, ax = plt.subplots(figsize=(height, 10))

    shc.dendrogram(
            linked,
            truncate_mode='lastp',  # show only the last p merged clusters
            p=len(counts.columns),  # show only the last p merged clusters
            # leaf_label_func=leaf_label(),
            labels=leaf_label(temp),
            orientation='left',
            leaf_font_size=8.,
            show_contracted=True,
            ax= ax  # to get a distribution impression in truncated branches
            )
    ax = plt.gca()

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    for xlabel_i in ax.get_xticklabels():
        xlabel_i.set_visible(False)
        xlabel_i.set_fontsize(0.0)
    for tick in ax.get_xticklines():
        tick.set_visible(False)

    return fig
