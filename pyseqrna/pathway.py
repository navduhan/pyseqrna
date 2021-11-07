import pandas as pd

import io
from urllib.request import urlopen


from statsmodels.stats.multitest import multipletests
import scipy.stats as stats
import numpy as np
import requests
import pandas as pd
from io import StringIO
from xml.etree import ElementTree
from future.utils import native_str

def _q(op, arg1, arg2=None, arg3=None):
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


resp= _q("list","pathway","ath")

data =[]
for r in resp:
    a, b=r.split("\t")
    data.append([a.split(":")[1],b.rstrip()])
    
    
genes = []
pathway={}
parse=None
nad=None
for d in data:
    resp=_q("get",d[0])
    # resp=_q("get",'dosa04626')

    for line in resp:
        line = line.strip()
        print(line)


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
               
    pathway[d[0]]= [d[0],desc,genes,len(genes) ]





df = pd.DataFrame(pathway).T
df.columns= ['ID', 'Term', 'Gene', 'Gene_length']


def enrichKEGG(df, file):
    background_count = df['Gene_length'].sum()


    df_goList = df[['ID', 'Gene']].values.tolist()

    go_dict = {}

    for value in df_goList:
        go_dict[value[0]] = value[1]

    count = df[['ID', 'Gene_length']].values.tolist()
    go_count = {}

    for c in count:
        go_count[c[0]] = c[1]

    df_List = df[['ID', 'Term']].values.tolist()
    KOdescription = {}

    for line in df_List:

        KOdescription[line[0]] = [line[1]]

    get_gene_ids_from_user = dict()
    gene_GO_count = dict()

    get_user_id_count_for_GO = dict()

    user_provided_uniq_ids = dict()

    for item in go_dict:

        get_gene_ids_from_user[item] = []

        gene_GO_count[item] = go_count[item]

        get_user_id_count_for_GO[item] = 0
        # GO terms
    bg_gene_count = background_count

    read_id_file = open("pySeqRNA_results.1/diff_genes/M1-V6.txt", 'r')
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


    pvalues = []
    enrichment_result = []
    # get total mapped genes from user list
    mapped_query_ids = sum(get_user_id_count_for_GO.values())

    for k in get_user_id_count_for_GO:
        gene_in_category = get_user_id_count_for_GO[k]

        gene_not_in_category_but_in_sample = mapped_query_ids - gene_in_category
        gene_not_in_catgory_but_in_genome = gene_GO_count[k] - gene_in_category
        bg_gene_GO_ids = gene_GO_count[k]
        bg_in_genome = bg_gene_count - mapped_query_ids - (gene_in_category + gene_not_in_catgory_but_in_genome) \
            + gene_in_category
        gene_ids = get_gene_ids_from_user[k]
        gID = ""

        for g in gene_ids:
            gID += g+"/"

        gID = gID.rsplit("/", 1)[0]
        pvalue = stats.hypergeom.sf(
            gene_in_category - 1, bg_gene_count, gene_GO_count[k], mapped_query_ids)

        if gene_in_category > 0:
            pvalues.append(pvalue)

            enrichment_result.append([k, KOdescription[k][0], 
                                    f"{gene_in_category}/{mapped_query_ids}", f"{go_count[k]}/{bg_gene_count}", pvalue,len(gene_ids), gID])

    fdr = list(multipletests(pvals=pvalues, method='fdr_bh')[1])

    a = [i for i in fdr if i <= 0.05]

    end = pd.DataFrame(enrichment_result)
    end.columns = ['Pathway_ID', 'Description',  'GeneRatio', 'BgRatio','Pvalues', 'Count', 'Genes' ]
    end.insert(5, 'FDR', fdr)

    return end

