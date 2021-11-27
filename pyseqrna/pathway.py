import io
from urllib.request import urlopen
from statsmodels.stats.multitest import multipletests
import scipy.stats as stats
import numpy as np
import pandas as pd
from io import StringIO
from pyseqrna.pyseqrna_utils import PyseqrnaLogger

log = PyseqrnaLogger(mode='a', log="pathway")


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

def kegg_list(sp):
    resp= _q("list","pathway",sp)

    data =[]
    for r in resp:
        a, b=r.split("\t")
        data.append([a.split(":")[1],b.rstrip()])



    pathway={}
    parse=None
    nad=None
    bg=[]
    for d in data:
        resp=_q("get",d[0])
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


def enrichKEGG(file, df, background_count):

   
    
    
    log.info(f"Performing KEGG enrichment analysis on {file}")  

    df_keggList = df[['ID', 'Gene']].values.tolist()

    kegg_dict = {}

    for value in df_keggList:
        kegg_dict[value[0]] = value[1]

    count = df[['ID', 'Gene_length']].values.tolist()
    kegg_count = {}

    for c in count:
        kegg_count[c[0]] = c[1]
    
    df_List = df[['ID', 'Term']].values.tolist()
    kegg_description = {}

    for line in df_List:

        kegg_description[line[0]] = line[1]
    

    get_gene_ids_from_user = dict()
    gene_kegg_count = dict()

    get_user_id_count_for_kegg = dict()

    user_provided_uniq_ids = []

    for item in kegg_dict:

        get_gene_ids_from_user[item] = []

        gene_kegg_count[item] = kegg_count[item]

        get_user_id_count_for_kegg[item] = 0
        # GO terms
    bg_gene_count = int(background_count)

    read_id_file = open(file, 'r')
    for gene_id in read_id_file:
        gene_id = gene_id.strip().upper()
    # remove the duplicate ids and keep unique
        user_provided_uniq_ids.append(gene_id)
    read_id_file.close()


    anot_count = 0
    for k1 in kegg_dict:
        intersection = set(user_provided_uniq_ids). intersection(kegg_dict[k1]) 
                # if the user input id present in df_dict_glist increment count
        get_gene_ids_from_user[k1]=list(intersection)
        get_user_id_count_for_kegg[k1] = len(intersection)
        anot_count += len(intersection)
    
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

  

