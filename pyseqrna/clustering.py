import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as shc

def llf(temp,xx):
   
    return "{}".format(temp[xx])
def clusterSample(countDF = None):

    corrCount = countDF.corr()
    linked = shc.linkage(corrCount, 'ward')
    R = shc.dendrogram(
                    linked,
                    truncate_mode='lastp',  # show only the last p merged clusters
                    p=len(r.columns),  # show only the last p merged clusters
                    no_plot=True,
                    )

    temp = {R["leaves"][ii]: countDF.columns[ii] for ii in range(len(R["leaves"]))}
    
   

    shc.dendrogram(
            linked,
            truncate_mode='lastp',  # show only the last p merged clusters
            p=len(r.columns),  # show only the last p merged clusters
            leaf_label_func=llf(temp),
            orientation='left',
           
            leaf_font_size=8.,
            show_contracted=True,  # to get a distribution impression in truncated branches
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

    return plt
