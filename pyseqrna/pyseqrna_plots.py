import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import seaborn as sns



def plotVolcano(degDF=None, comp=None,FOLD=2,pValue=0.05,color=('red','grey','green'), dim=(8,5), dotsize=8, markerType='o', alpha=0.5):
    """[summary]

    Args:
        degDF ([type], optional): [description]. Defaults to None.
        comp ([type], optional): [description]. Defaults to None.
        FOLD (int, optional): [description]. Defaults to 2.
        pValue (float, optional): [description]. Defaults to 0.05.
        color (tuple, optional): [description]. Defaults to ('red','grey','green').
        dim (tuple, optional): [description]. Defaults to (8,5).
        dotsize (int, optional): [description]. Defaults to 8.
        markeType (str, optional): [description]. Defaults to 'o'.
        alpha (float, optional): [description]. Defaults to 0.5.
    """

    PVAL = "pvalue("+comp+")"

    LFC = "logFC("+comp+")"
    _y = r'$ log_{2}(Fold Change)$'
    _x = r'$ -log_{10}(P-value)$'

    dk = degDF.filter(regex=comp, axis=1)
    
    if dk.empty:
        pass
    else:
        
        final = dk[dk[PVAL]<=pValue].copy()

        final.loc[final[LFC]>=np.log2(FOLD),'colorADD'] = color[2]

        final.loc[final[LFC]<=-np.log2(FOLD), 'colorADD'] = color[0]

        final['colorADD'].fillna(color[1], inplace=True)  

        final['log(10)_pvalue'] = -(np.log10(dk[PVAL]))
        
        color_values = {col: i for i, col in enumerate(color)}

        color_num = [color_values[i] for i in final['colorADD']]
        
        legendlabels=['Significant down', 'Not significant', 'Significant up']  

        fig, ax = plt.subplots(figsize=dim)

        a = ax.scatter(final[LFC], final['log(10)_pvalue'], c=color_num, cmap=ListedColormap(color), alpha=alpha,s=dotsize, marker=markerType)
        ax.legend(handles=a.legend_elements()[0], labels=legendlabels,   bbox_to_anchor=(1.46,1,0,0))
        ax.set_xlabel(_x)
        ax.set_ylabel(_y)
        fig.tight_layout()
    

    return fig ,ax

# degDF = pd.read_excel("./test_pyseq/Raw_DEGs_all.xlsx")
# comp = ['M1-V1', 'M1-M6', 'M1-V6', 'M1-A1', 'V6-A12']
# countDF = pd.read_csv("./test_pyseq/Counts_final.txt", sep="\t")

def plotMA(degDF=None, countDF=None, comp=None,FOLD=2,color=('green','grey','red'), dim=(8,5), dotsize=8, markerType='o', alpha=0.5):
    """
    

    Args:
        degDF ([type], optional): [description]. Defaults to None.
        comp ([type], optional): [description]. Defaults to None.
        FOLD (int, optional): [description]. Defaults to 2.
        color (tuple, optional): [description]. Defaults to ('red','grey','green').
        dim (tuple, optional): [description]. Defaults to (8,5).
        dotsize (int, optional): [description]. Defaults to 8.
        markeType (str, optional): [description]. Defaults to 'o'.
        alpha (float, optional): [description]. Defaults to 0.5.
    """

    LFC = "logFC("+comp+")"
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

        final.loc[final[LFC]>=np.log2(FOLD),'colorADD'] = color[2]

        final.loc[final[LFC]<=-np.log2(FOLD), 'colorADD'] = color[0]

        final['colorADD'].fillna(color[1], inplace=True)  
        final['mean'] = final[['first', 'second']].mean(axis=1)
        finalDF = final.loc[final['mean']>0].copy()

        np.seterr(divide = 'ignore') 

        finalDF['log(2)_Mean'] = np.log2(finalDF['first']) + np.log2(finalDF['second']) / 2
        
        color_values = {col: i for i, col in enumerate(color)}

        color_num = [color_values[i] for i in finalDF['colorADD']]
        
        legendlabels=['Significant up', 'Not significant', 'Significant down']  

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
   

    return fig ,ax






def plotHeatmap(degDF= None, combinations=None, num=50, figdim=(12,10), type='counts'):

    # Need to add argument for all deg hetamp.
    degDF = degDF.set_index('Gene')
    if type=='deg':
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
    elif type =='counts':
        fin=degDF.nlargest(int(num), degDF.columns)
        
    fig, ax = plt.subplots(figsize=figdim)

    sns.heatmap(fin, cmap="seismic",ax=ax)

    fig.tight_layout()

    return fig,ax
