import io
from urllib.request import urlopen
from statsmodels.stats.multitest import multipletests
import scipy.stats as stats
import numpy as np
import pandas as pd
from io import StringIO
from pyseqrna.pyseqrna_utils import PyseqrnaLogger

log = PyseqrnaLogger(mode='a', log="pathway")


def kegg_list(species):
    URL = "http://rest.kegg.jp/list/pathway/%s"
    resp = urlopen(URL % (species))
    handle = io.TextIOWrapper(resp, encoding="UTF-8")
    handle.url = resp.url
    result = pd.read_csv(StringIO(handle.read()), sep='\t', names=['Pathway', 'Description'])
    path = [i.split(':', 1)[1] for i in list(result['Pathway'])]
    desc = [i.split(' - ', 1)[0] for i in list(result['Description'])]
    final = list(zip(path,desc))
    kegg_description = {}

    for f in final:
        kegg_description[f[0]]=f[1]
    
    return kegg_description

def kegg_gene(species):
    URL = "http://rest.kegg.jp/link/%s/pathway"
    resp = urlopen(URL % (species))
    handle = io.TextIOWrapper(resp, encoding="UTF-8")
    handle.url = resp.url
    result = pd.read_csv(StringIO(handle.read()), sep='\t', names=['Pathway', 'Accession'])
    path = [i.split(':', 1)[1] for i in list(result['Pathway'])]
    desc = [i.split(':', 1)[1] for i in list(result['Accession'])]
    pp = list(zip(path,desc))

    pathway = {}
    for p in pp:
        if p[0] not in pathway:
            pathway[p[0]]= [p[1]]
        else:
            pathway[p[0]].append(p[1])
    
    x = np.array(desc)

    bg_count = len(np.unique(x))

    finalp = {}

    for k, v in pathway.items():
        finalp[k]=[k,v,len(v)]
    
    final = pd.DataFrame.from_dict(finalp).T
    final.columns = ['ID', 'Gene', 'Gene_length']

    return final, bg_count

def enrichKEGG(file, species):

    log.info("Reading annotation from KEGG")
    
    df, background_count = kegg_gene(species)
    kegg_description = kegg_list(species)
    log.info(f"Performing KEGG enrichment analysis on {file}")  

    df_keggList = df[['ID', 'Gene']].values.tolist()

    kegg_dict = {}

    for value in df_keggList:
        kegg_dict[value[0]] = value[1]

    count = df[['ID', 'Gene_length']].values.tolist()
    kegg_count = {}

    for c in count:
        kegg_count[c[0]] = c[1]
    

    get_gene_ids_from_user = dict()
    gene_kegg_count = dict()

    get_user_id_count_for_kegg = dict()

    user_provided_uniq_ids = dict()

    for item in kegg_dict:

        get_gene_ids_from_user[item] = []

        gene_kegg_count[item] = kegg_count[item]

        get_user_id_count_for_kegg[item] = 0
        # GO terms
    bg_gene_count = background_count

    read_id_file = open(file, 'r')
    for gene_id in read_id_file:
        gene_id = gene_id.strip().upper()
    # remove the duplicate ids and keep unique
        user_provided_uniq_ids[gene_id] = 0
    read_id_file.close()


    anot_count = 0
    for k1 in kegg_dict:
        for k2 in user_provided_uniq_ids:
            if k2 in kegg_dict[k1]:
                # if the user input id present in df_dict_glist increment count
                get_gene_ids_from_user[k1].append(k2)
                get_user_id_count_for_kegg[k1] += 1
                anot_count += 1
    
    pvalues = []
    enrichment_result = []
    # get total mapped genes from user list
    mapped_query_ids = sum(get_user_id_count_for_kegg.values())

    for k in get_user_id_count_for_kegg:
        gene_in_category = get_user_id_count_for_kegg[k]

        gene_not_in_category_but_in_sample = mapped_query_ids - gene_in_category
        gene_not_in_catgory_but_in_genome = gene_kegg_count[k] - gene_in_category
        bg_gene_kegg_ids = gene_kegg_count[k]
        bg_in_genome = bg_gene_count - mapped_query_ids - (gene_in_category + gene_not_in_catgory_but_in_genome) \
            + gene_in_category
        gene_ids = get_gene_ids_from_user[k]
        gID = ""

        for g in gene_ids:
            gID += g+","

        gID = gID.rsplit(",", 1)[0]
        pvalue = stats.hypergeom.sf(
            gene_in_category - 1, bg_gene_count, gene_kegg_count[k], mapped_query_ids)

        if gene_in_category > 0:
            pvalues.append(pvalue)

            enrichment_result.append([k, kegg_description[k], 
                                    f"{gene_in_category}/{mapped_query_ids}", f"{kegg_count[k]}/{bg_gene_count}", pvalue,len(gene_ids), gID])

    fdr = list(multipletests(pvals=pvalues, method='fdr_bh')[1])

    a = [i for i in fdr if i <= 0.05]

    end = pd.DataFrame(enrichment_result)
    end.columns = ['Pathway_ID', 'Description',  'GeneRatio', 'BgRatio','Pvalues', 'Count', 'Genes' ]
    end.insert(5, 'FDR', fdr)


    return end

  

