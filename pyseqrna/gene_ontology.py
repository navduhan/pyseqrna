#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
:Title: gene_ontology module is for performing gene ontology enrichment analysis of differentially expressed genes

:Created: January 5, 2022

:Author: Naveen Duhan
'''

from __future__ import division
from statsmodels.stats.multitest import multipletests
import scipy.stats as stats
import numpy as np
import requests
import pandas as pd
import math
from io import StringIO, TextIOWrapper
from xml.etree import ElementTree
from future.utils import native_str
from pyseqrna.pyseqrna_utils import PyseqrnaLogger
from pyseqrna import pyseqrna_utils as pu
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from urllib.request import urlopen, urlretrieve
import textwrap
import json
log = PyseqrnaLogger(mode='a', log="go")


class GeneOntology:

    """
    This class is for Gene Ontology enrichment 

    :param species: Species name. Ex. for Arabidopsis thaliana it is athaliana

    :param type: Species is from plants or animals.

    :param keyType: Genes are from NCBI or ENSEMBL. Default is ENSEMBL.

    :param taxid: Taxonomy ID if keyType is NCBI.

    :param gff: Gene feature file.
    
    """
    def __init__(self, species=None, type=None,  keyType = 'ensembl', taxid=None, gff= None):

        self.species = species 
        self.type = type 
        self.keyType = keyType
        self.taxid = taxid
        self.gff = gff

        if self.keyType.lower() == 'ensembl':
            log.info("Fetching Gene Ontology from Biomart")
            self.gdata = self._query(self.species,self.type)
            log.info("Processing Gene Ontology data from Biomart")
            self.df, self.background_count = self._preprocessBioMart(self.gdata)
        if self.keyType.lower() == 'ncbi':
             log.info("Fetching Gene Ontology from pyseqrnautils API")
             self.df, self.background_count = self._parse_GO(self.taxid)
             self.idmapping = pu.parse_gff(self.gff)
        
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

    def _go_organism(self, type = "plants"):
            if type == "plants":
                resp =  urlopen("https://plants.ensembl.org/biomart/martservice?type=datasets&requestid=biomaRt&mart=plants_mart")

            if type == 'animals':
                resp =  urlopen("https://www.ensembl.org/biomart/martservice?type=datasets&requestid=biomaRt&mart=ENSEMBL_MART_ENSEMBL")

            handle = TextIOWrapper(resp, encoding="UTF-8")
            handle.url = resp.url

            df =pd.read_csv(handle, sep="\t", names=['Table', 'Code', 'Organism', 'A', 'Assembly', 'B', 'C', 'Default', 'Date'])

            organism = df[['Code', 'Organism']]

            return organism
            
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
        attributes = ["ensembl_gene_id", "ensembl_transcript_id",
                    "go_id", "name_1006", "namespace_1003", "definition_1006"]
        for attr in attributes:
            self._add_attr_node(dataset, attr)

        response = self._get_request(
        uri , query=ElementTree.tostring(root))
        result = pd.read_csv(StringIO(response.text), sep='\t')
        result.columns = ['Gene', 'Transcript', 'GO_ID',
                    'GO_term', 'GO_ontology', 'GO_def']
        
        return result


    def _preprocessBioMart(self,data):

        df = data
        
        df2 = df[df['GO_ID'].notna()]
        gg = list(df2['Gene'])
        x = np.array(gg)

        bg_count = len(np.unique(x))

        lines = df2.values.tolist()
        GeneID = {}

        for line in lines:

            if line[2] not in GeneID:

                GeneID[line[2]] = [line[0]]

            else:
                GeneID[line[2]].append(line[0])

        GO_rest = {}

        for line in lines:

            if line[2] not in GO_rest:

                GO_rest[line[2]] = [
                                    line[2], line[3], line[4]]

        ds = [GO_rest, GeneID]
        d = {}
        for k in GO_rest.keys():
            d[k] = list(d[k] for d in ds)

        dd = []
        
        for k, v in d.items():
            v[1] = [i for i in v[1] if str(i) != 'NaN']
            v[1] =[x.upper() for x in v[1]]
            if v[0][2] == 'cellular_component':
                v[0][2] = 'CC'
            if v[0][2] == 'molecular_function':
                v[0][2] = 'MF'
            if v[0][2] == 'biological_process':
                v[0][2] = 'BP'

            dd.append([v[0][0], v[0][1], v[0][2], 
                        v[1], len(v[1])])

        finalDF = pd.DataFrame(dd, columns=[
                            'ID', 'Term', 'Ontology',  'Gene', 'Gene_length'])

        return finalDF, bg_count

    def _parse_GO(self, taxid ):

        resp = requests.get(f"http://bioinfo.usu.edu/pyseqrnautils/api/go/?taxid={taxid}")

        df = pd.DataFrame(json.loads(resp.text))

        gg = list(df['entrez'])
        x = np.array(gg)

        bg_count = len(np.unique(x))

        lines = df.values.tolist()
        GeneID = {}

        for line in lines:

            if line[3] not in GeneID:

                GeneID[line[3]] = [line[2]]

            else:
                GeneID[line[3]].append(line[2])

        GO_rest = {}

        for line in lines:

            if line[3] not in GO_rest:

                GO_rest[line[3]] = [
                                    line[3], line[5], line[7]]

        ds = [GO_rest, GeneID]
        d = {}
        for k in GO_rest.keys():
            d[k] = list(d[k] for d in ds)

        dd = []

        for k, v in d.items():
            dd.append([v[0][0], v[0][1], v[0][2], 
                                v[1], len(v[1])])

        finalDF = pd.DataFrame(dd, columns=[
                                    'ID', 'Term', 'Ontology',  'Gene', 'Gene_length'])
        
        return finalDF, bg_count

    def _fdr_calc(self, x):
        """
        Assumes a list or numpy array x which contains p-values for multiple tests
        Copied from p.adjust function from R  
        """
        if len(x)!=0:
            p_vals = pd.Series(x)
        else:
            p_vals = pd.Series(x, dtype=object)
        from scipy.stats import rankdata
        ranked_p_values = rankdata(p_vals)
        fdr = p_vals * len(p_vals) / ranked_p_values
        fdr[fdr > 1] = 1
        return fdr

    def dotplotGO(self, df=None, nrows=20, colorBy='logPvalues'):
        """
        This function creates a dotplot for Gene Ontology enrichment.

        :param df: Gene Ontology enrichment file from enrichGO function.

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

        fig, ax = plt.subplots(figsize=(10,10), dpi=300)
        scatter = ax.scatter(x=df['Counts'], y= df['GO Term'], s=df['Counts'], c=df[colorBy])
        ax.xaxis.get_major_locator().set_params(integer=True)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_bounds((0, nrows))
        # add some space between the axis and the plot
        ax.spines['left'].set_position(('outward', 8))
        ax.spines['bottom'].set_position(('outward', 5))
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.xlabel("Counts", fontsize=12, fontweight='bold')
        plt.ylabel("GO Description", fontsize=12, fontweight='bold')

        cbar = plt.colorbar(scatter,shrink=.25, pad=.2, aspect=10)
        cbar.ax.set_title(title,pad=20, fontweight='bold')
        

        return fig

    def barplotGO(self, df=None,nrows=20, colorBy='logPvalues' ):

        """
        This function creates a barplot for Gene Ontology enrichment.

        :param df: Gene Ontology enrichment file from enrichGO function.

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
        df =df.sort_values('Counts', ascending=True)
        counts = df['Counts'].values.tolist()
        terms = df['GO Term'].values.tolist()

        data_color_normalized = [x / max(df[colorBy]) for x in df[colorBy]]
        fsize = math.ceil(nrows/2)
        fig, ax = plt.subplots(figsize=(10, 10), dpi=300)

        my_cmap = plt.cm.get_cmap('RdYlBu')
        colors = my_cmap(data_color_normalized)

        ax.barh(terms, counts, color=colors)
        ax.xaxis.get_major_locator().set_params(integer=True)
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_bounds((0, nrows))
        # add some space between the axis and the plot
        ax.spines['left'].set_position(('outward', 8))
        ax.spines['bottom'].set_position(('outward', 5))
        plt.xticks(fontsize=10)
        plt.yticks(fontsize=10, wrap=True)

        plt.xlabel("Counts", fontsize=12, fontweight='bold')
        plt.ylabel("GO Description", fontsize=12, fontweight='bold')
        sm = ScalarMappable(cmap=my_cmap, norm=plt.Normalize(0,max(df[colorBy])))

        sm.set_array([])

        cbar = plt.colorbar(sm, shrink=0.25,pad=.02, aspect=10)
        cbar.ax.set_title(title,pad=20,fontweight='bold')
        

        return fig

    def enrichGO( self, file= None, pvalueCutoff=0.05, plot=True, plotType= 'dotplot', nrows=20,colorBy='logPvalues'):

        # Need to add support for different database IDs
        """
        This function performs Gene Ontology enrichment of DEGs.

        :param file: Differentially expressed genes in a sample.

        :param pvalueCutoff: P-value cutoff for enrichment. Default is 0.05.

        :param plot: True if a plot is needed. Default is True.

        :param plotType: Gene Ontology enrichment visualization on dotplot/barplot. Default is dotplot.

        :param nrows: Number of rows to plot. Default to 20 rows.

        :param colorBy: Color dot on plots with logPvalues / FDR. Defaults to 'logPvalues'.
        
        :returns:  a dictionary

        :rtype results: Gene Ontology enrichment results.

        :rtype plot: a dotplot/barplot
        """

        log.info(f"Performing GO enrichment analysis on {file}")  
        

        df_goList = self.df[['ID', 'Gene']].values.tolist()
        count = self.df[['ID', 'Gene_length']].values.tolist()
        df_List = self.df[['ID', 'Term', 'Ontology']].values.tolist()

        go_dict = {}

        for value in df_goList:
            go_dict[value[0]] = str(value[1]).upper()

        
        go_count = {}

        for c in count:
            go_count[c[0]] = c[1]

       
        KOdescription = {}

        for line in df_List:

            KOdescription[line[0]] = [line[1], line[2]]

        get_gene_ids_from_user = dict()
        gene_GO_count = dict()

        get_user_id_count_for_GO = dict()

        user_provided_uniq_ids = dict()
        user_genecount = []
        for item in go_dict:

            get_gene_ids_from_user[item] = []

            gene_GO_count[item] = go_count[item]

            get_user_id_count_for_GO[item] = 0
            # GO terms
        bg_gene_count = self.background_count

        if self.keyType.lower() == 'ncbi':

            ufile = pd.read_csv(file)

            id_intermediate = ufile.merge(self.idmapping, on='Gene').drop_duplicates()

            read_id_file = id_intermediate['entrez'].values.tolist()

            for gene_id in read_id_file:

                gene_id = gene_id.strip().upper()
        # remove the duplicate ids and keep unique
                user_provided_uniq_ids[gene_id] = 0

        if self.keyType.lower()== 'ensembl':

            read_id_file = open(file, 'r')
            for gene_id in read_id_file:
                gene_id = gene_id.strip().upper()
            # remove the duplicate ids and keep unique
                user_provided_uniq_ids[gene_id] = 0
            read_id_file.close()


        anot_count = 0
        for k1 in go_dict:
            for k2 in user_provided_uniq_ids:
                if k2 in go_dict[k1]:
                    # if the user input id present in df_dict_glist increment count
                    get_gene_ids_from_user[k1].append(k2)
                    get_user_id_count_for_GO[k1] += 1
                    anot_count += 1
                    if k2 not in user_genecount:
                        user_genecount.append(k2)
        pvalues = []
        enrichment_result = []
        # get total mapped genes from user list
        # mapped_query_ids = sum(get_user_id_count_for_GO.values())
        mapped_query_ids = len(user_genecount)

        for k in get_user_id_count_for_GO:
            gene_in_category = get_user_id_count_for_GO[k]

            gene_not_in_category_but_in_sample = mapped_query_ids - gene_in_category
            gene_not_in_catgory_but_in_genome = gene_GO_count[k] - gene_in_category
            # bg_gene_GO_ids = gene_GO_count[k]
            # bg_in_genome = bg_gene_count - mapped_query_ids - (gene_in_category + gene_not_in_catgory_but_in_genome) \
            #     + gene_in_category
            gene_ids = get_gene_ids_from_user[k]
            gID = ""

            for g in gene_ids:
                gID += g+","

            gID = gID.rsplit(",", 1)[0]
            pvalue = stats.hypergeom.sf(
                gene_in_category -1, bg_gene_count, gene_GO_count[k], mapped_query_ids)
            

            if gene_in_category > 0:
                pvalues.append(pvalue)

                enrichment_result.append([k, KOdescription[k][0], KOdescription[k][1], 
                                        f"{gene_in_category}/{mapped_query_ids}", f"{go_count[k]}/{bg_gene_count}", pvalue, len(gene_ids), gID])
        
        fdr = list(self._fdr_calc(pvalues))

        # a = [i for i in fdr if i <= 0.05]

        end = pd.DataFrame(enrichment_result)

        if end.shape[0]>1:
            end.columns = ['GO ID', 'GO Term', 'Ontology',  'GeneRatio', 'BgRatio','Pvalues', 'Counts', 'Genes' ]
        
            end.insert(7, 'FDR', fdr)

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
                    fig = self.dotplotGO(results,nrows,colorBy)
                if plotType == 'barplot':
                    fig = self.barplotGO(results,nrows,colorBy)
                
            return {'result': results, 'plot': fig}

        return "No Gene Ontology"

