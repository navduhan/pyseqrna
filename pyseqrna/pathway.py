#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
:Title: This module is for performing KEGG pathway enrichment analysis of differentially expressed genes

:Created: January 10, 2022

:Author: Naveen Duhan
'''
import io
import os
import re
import requests
import math
from urllib.request import urlopen, urlretrieve
from statsmodels.stats.multitest import multipletests
import scipy.stats as stats
import numpy as np
import pandas as pd
from io import StringIO
import matplotlib.pyplot as plt
from pyseqrna import pyseqrna_utils as pu
from matplotlib.cm import ScalarMappable
from pyseqrna.pyseqrna_utils import PyseqrnaLogger


log = PyseqrnaLogger(mode='a', log="pathway")


class Pathway:


    """
    This class is for KEGG Pathway enrichment 

    :param species: Species name. Ex. for Arabidopsis thaliana it is athaliana

    :param keyType: Genes are from NCBI or ENSEMBL. Default is ENSEMBL.

    :param gff: Gene feature file.
    
    """

    def __init__(self, species=None, keyType=None, gff=None):

        self.species = species 

        self.keyType = keyType  

        self.gff = gff

        if self.keyType.lower() == 'ncbi':

             log.info("Fetching Pathways from KEGG API")

             self.df, self.background_count = self._kegg_list()

             self.idmapping = pu.parse_gff(self.gff)

        if self.keyType.lower() == 'ensembl':

             log.info("Fetching Pathways from pyseqrna API")

             self.df, self.background_count = self._get_pathways()
             

    def _q(self,op, arg1, arg2=None, arg3=None):

        URL = "http://rest.kegg.jp/%s"

        if arg2 and arg3:

            args = "%s/%s/%s/%s" % (op, arg1, arg2, arg3)

        elif arg2:

            args = "%s/%s/%s" % (op, arg1, arg2)

        else:

            args = "%s/%s" % (op, arg1)

        resp = urlopen(URL % (args))

        if "image" == arg2:

            return resp

        handle = io.TextIOWrapper(resp, encoding="UTF-8")

        handle.url = resp.url

        return handle

    def _kegg_organism(self):

        resp =  urlopen("http://rest.kegg.jp/list/organism")

        handle = io.TextIOWrapper(resp, encoding="UTF-8")

        handle.url = resp.url

        df = pd.read_csv(handle, sep="\t", names=['Organism ID', 'Organism Code', 'Organism Name', 'Taxonomy'])

        organisms = df[['Organism ID', 'Organism Code', 'Organism Name']]

        return organisms

    def _kegg_list(self):

        resp= self._q("list","pathway",self.species)

        data =[]

        for r in resp:

            a, b=r.split("\t")

            data.append([a.split(":")[1],b.rstrip()])



        pathway={}
        parse=None
        nad=None
        bg=[]
        for d in data:
            resp=self._q("get",d[0])
            # resp=_q("get",'dosa04626')
            genes = []
            for line in resp:
                line = line.strip()
                # print(line)


                if not line.startswith("/"):
                    if not line.startswith(" "):
                        first_word = line.split(" ")[0]
                        if first_word.isupper() and first_word.isalpha():
                            parse = first_word    
                        if parse == "NAME":
                            nad = line.replace(parse,"").strip()
                            desc = nad.split(" - ")[0]       
                        if parse== "GENE":
                            gened = line.replace(parse,"").strip().split(" ")[0]
                            genes.append(gened)
            bg.extend(genes)
            pathway[d[0]]= [d[0],desc,genes,len(genes) ]

        df = pd.DataFrame(pathway).T
        df.columns= ['ID', 'Term', 'Gene', 'Gene_length']

        x=np.array(bg)
        bg_count = len(np.unique(x))

        return df, bg_count

    def _get_pathways(self):

        r = requests.get(f"http://bioinfo.usu.edu/pyseqrna-api/list/pathways/{self.species}")
        m = re.sub('<[^<]+?>', '', r.text)
        df = pd.read_csv(io.StringIO(m), sep="\t", names =['Species', 'Gene', 'ID', 'Term'])
        dd= df.values.tolist()
        kk = {}
        kd = {}
        genes = []

        for d in dd:

            if d[2] in kk:

                kk[d[2]].append(d[1].upper())
            
            else:

                kk[d[2]]= [d[1].upper()]

                kd[d[2]] = d[3]

            genes.append(d[1])

        pathways = []

        for k in kk:

            pathways.append([k, kd[k],kk[k], len(kk[k])])

        dp = pd.DataFrame(pathways, columns=['ID', 'Term', 'Gene', 'Gene_length'])

        x= np.array(genes)

        bg_count = len(np.unique(x))

        return dp, bg_count

    def _fdr_calc(self, x):
        """
        Assumes a list or numpy array x which contains p-values for multiple tests
        Copied from p.adjust function from R  
        """
        p_vals = pd.Series(x)
        from scipy.stats import rankdata
        ranked_p_values = rankdata(p_vals)
        fdr = p_vals * len(p_vals) / ranked_p_values
        fdr[fdr > 1] = 1
        return fdr


    def dotplotKEGG(self, df=None, nrows=20, colorBy='logPvalues'):
        """
        This function creates a dotplot for KEGG pathway enrichment.

        :param df: KEGG pathway enrichment file from enrichKEGG function.

        :param nrows: Number of rows to plot. Default to 20 rows.

        :param colorBy: Color dot on plots with logPvalues / FDR. Defaults to 'logPvalues'.
        
        :returns: a dotplot
        """

        if colorBy=='logPvalues':
            df['logPvalues'] = round(-np.log10(df['Pvalues']),2)
            title = '-log10(Pvalues)'
        
        if colorBy=='FDR':
            df = df
            title = 'FDR'

        df =df.sort_values('Counts', ascending=False)
        df = df.head(nrows)
        df =df.sort_values('Counts', ascending=True)
        fsize = math.ceil(nrows/2)
        fig, ax = plt.subplots(figsize=(10, 10), dpi=300)
        scatter = ax.scatter(x=df['Counts'], y= df['Description'], s=df['Counts'], c=df[colorBy])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_bounds((0, nrows))
        # add some space between the axis and the plot
        ax.spines['left'].set_position(('outward', 8))
        ax.spines['bottom'].set_position(('outward', 5))
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.xlabel("Counts", fontsize=12, fontweight='bold')
        plt.ylabel("Description", fontsize=12, fontweight='bold')

        cbar = plt.colorbar(scatter,shrink=.25, pad=.2, aspect=10)
        cbar.ax.set_title(title,pad=20, fontweight='bold')
        fig.tight_layout()

        return fig

    def barplotKEGG(self, df=None,nrows=20, colorBy='logPvalues' ):

        """
        This function creates a barplot for KEGG pathway enrichment.

        :param df: KEGG pathway enrichment file from enrichKEGG function.

        :param nrows: Number of rows to plot. Default to 20 rows.

        :param colorBy: Color bar on plots with logPvalues / FDR. Defaults to 'logPvalues'.
        
        :returns: a barplot
        """
        
        if colorBy=='logPvalues':
            df['logPvalues'] = round(-np.log10(df['Pvalues']),2)
            title = '-log10(Pvalues)'
        
        if colorBy=='FDR':
            df = df
            title = 'FDR'

        df =df.sort_values('Counts', ascending=False)
        df = df.head(nrows)
        df = df.sort_values('Counts', ascending=True)
        counts = df['Counts'].values.tolist()
        terms = df['Description'].values.tolist()

        data_color_normalized = [x / max(df[colorBy]) for x in df[colorBy]]
        
        fig, ax = plt.subplots(figsize=(10, 10), dpi=300)

        my_cmap = plt.cm.get_cmap('RdYlBu')
        colors = my_cmap(data_color_normalized)

        rects = ax.barh(terms, counts, color=colors)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_bounds((0, nrows))
        # add some space between the axis and the plot
        ax.spines['left'].set_position(('outward', 8))
        ax.spines['bottom'].set_position(('outward', 5))
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.xlabel("Counts", fontsize=12, fontweight='bold')
        plt.ylabel("Description", fontsize=12, fontweight='bold')
        sm = ScalarMappable(cmap=my_cmap, norm=plt.Normalize(0,max(df[colorBy])))

        sm.set_array([])

        cbar = plt.colorbar(sm, shrink=0.25,pad=.02, aspect=10)
        cbar.ax.set_title(title,pad=20,fontweight='bold')
        fig.tight_layout()

        return fig

    def enrichKEGG(self,file, pvalueCutoff=0.05, plot=True, plotType= 'dotplot', nrows=20,colorBy='logPvalues'):

                # Need to add support for different database IDs
        """
        This function performs KEGG pathway enrichment of DEGs.

        :param file: Differentially expressed genes in a sample.

        :param pvalueCutoff: P-value cutoff for enrichment. Default is 0.05.

        :param plot: True if a plot is needed. Default is True.

        :param plotType: KEGG pathway enrichment visualization on dotplot/barplot. Default is dotplot.

        :param nrows: Number of rows to plot. Default to 20 rows.

        :param colorBy: Color dot on plots with logPvalues / FDR. Defaults to 'logPvalues'.
        
        :returns:  a dictionary

        :rtype results: KEGG pathway enrichment results.

        :rtype plot: a dotplot/barplot
        """

        log.info(f"Performing KEGG enrichment analysis on {file}")  

        df_keggList = self.df[['ID', 'Gene']].values.tolist()

        kegg_dict = {}

        for value in df_keggList:
            kegg_dict[value[0]] = value[1]

        count = self.df[['ID', 'Gene_length']].values.tolist()
        kegg_count = {}

        for c in count:
            kegg_count[c[0]] = c[1]
        
        df_List = self.df[['ID', 'Term']].values.tolist()
        kegg_description = {}

        for line in df_List:

            kegg_description[line[0]] = line[1]
        

        user_gene_ids = dict()
        gene_kegg_count = dict()

        userID_count_kegg = dict()

        user_unique_gene_id = dict()
        user_genecount = []

        for item in kegg_dict:

            user_gene_ids[item] = []

            gene_kegg_count[item] = kegg_count[item]

            userID_count_kegg[item] = 0
            # GO terms
        bg_gene_count = int(self.background_count)

        if self.keyType.lower() == 'ncbi':

            ufile = pd.read_csv(file)

            id_intermediate = ufile.merge(self.idmapping, on='Gene').drop_duplicates()

            read_id_file = id_intermediate['entrez'].values.tolist()
           

            for gene_id in read_id_file:
               
                gene_id = gene_id.strip().upper()
        # remove the duplicate ids and keep unique
                user_unique_gene_id[gene_id] = 0
               

        if self.keyType.lower()== 'ensembl':

            read_id_file = open(file, 'r')
            for gene_id in read_id_file:
                
                gene_id = gene_id.strip().upper()
            # remove the duplicate ids and keep unique
                user_unique_gene_id[gene_id] = 0
            read_id_file.close()



        pathway_annotation_count = 0
        for k1 in kegg_dict:
            for k2 in user_unique_gene_id:
                if k2 in kegg_dict[k1]:
                    # if the user input id present in df_dict_glist increment count
                    user_gene_ids[k1].append(k2)
                    userID_count_kegg[k1] += 1
                    pathway_annotation_count += 1
                    if k2 not in user_genecount:
                        user_genecount.append(k2)

        
        pvalues = []
        enrichment_result = []
        # get total mapped genes from user list
        # mapped_user_ids = sum(userID_count_kegg.values())
        mapped_user_ids = len(user_genecount)

        for k in userID_count_kegg:
            gene_in_pathway = userID_count_kegg[k]

            # gene_not_in_pathway_but_in_query = mapped_user_ids - gene_in_pathway
            gene_not_in_catgory_but_in_genome = gene_kegg_count[k] - gene_in_pathway
            # bg_gene_kegg_ids = gene_kegg_count[k]
            bg_in_genome = bg_gene_count - mapped_user_ids - (gene_in_pathway + gene_not_in_catgory_but_in_genome) \
                + gene_in_pathway
            gene_ids = user_gene_ids[k]
            gID = ""

            for g in gene_ids:
                gID += g+","

            gID = gID.rsplit(",", 1)[0]
            pvalue = stats.hypergeom.sf(
                gene_in_pathway - 1, bg_gene_count, gene_kegg_count[k], mapped_user_ids)

            if gene_in_pathway > 0:
                pvalues.append(pvalue)

                enrichment_result.append([k, kegg_description[k], 
                                        f"{gene_in_pathway}/{mapped_user_ids}", f"{kegg_count[k]}/{bg_gene_count}", pvalue,len(gene_ids), gID])

        # fdr = list(multipletests(pvals=pvalues, method='fdr_bh')[1])
        fdr = list(self._fdr_calc(pvalues))

        end = pd.DataFrame(enrichment_result)
        if end.shape[0]>1:
            end.columns = ['Pathway_ID', 'Description',  'GeneRatio', 'BgRatio','Pvalues', 'Counts', 'Genes' ]
            end.insert(5, 'FDR', fdr)
            end = end[end['Pvalues']<=pvalueCutoff]

            if self.keyType.lower()=='ncbi':

                 results = pu.change_ids(id_intermediate, end)

            if self.keyType.lower()== 'ensembl':

                results = end

            if results.shape[0]>=nrows:
                nrows = nrows
            else:
                nrows = results.shape[0]

            if plot:
                if plotType == 'dotplot':
                    fig = self.dotplotKEGG(results,nrows,colorBy)
                if plotType == 'barplot':
                    fig = self.barplotKEGG(results,nrows,colorBy)

            return {'result': results, 'plot': fig}

        
        return 'No Pathway'


        

    

