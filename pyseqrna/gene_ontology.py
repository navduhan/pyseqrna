from __future__ import division
from statsmodels.stats.multitest import multipletests
import scipy.stats as stats
import numpy as np
import requests
import pandas as pd
from io import StringIO
from xml.etree import ElementTree
from future.utils import native_str
from pyseqrna.pyseqrna_utils import PyseqrnaLogger


log = PyseqrnaLogger(mode='a', log="go")

def get_request(url,  **params):

    if params:
        r = requests.get(url, params=params, stream=True)
    else:
        r = requests.get(url)
    r.raise_for_status()

    return r


def _add_attr_node(root, attr):
    attr_el = ElementTree.SubElement(root, 'Attribute')
    attr_el.set('name', attr)


def query(species):
    root = ElementTree.Element('Query')
    root.set('virtualSchemaName', 'plants_mart')
    root.set('formatter', 'TSV')
    root.set('header', '1')
    root.set('uniqueRows', native_str(int(True)))
    root.set('datasetConfigVersion', '0.6')

    dataset = ElementTree.SubElement(root, 'Dataset')
    dataset.set('name', species+"_eg_gene")
    dataset.set('interface', 'default')
    attributes = ["ensembl_gene_id", "ensembl_transcript_id",
                  "go_id", "name_1006", "namespace_1003", "definition_1006"]
    for attr in attributes:
        _add_attr_node(dataset, attr)

    response = get_request(
        "https://plants.ensembl.org/biomart/martservice", query=ElementTree.tostring(root))
    result = pd.read_csv(StringIO(response.text), sep='\t')
    result.columns = ['Gene', 'Transcript', 'GO_ID',
                  'GO_term', 'GO_ontology', 'GO_def']
    
    return result


def preprocessBioMart(data):
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
                                line[2], line[3], line[4], line[5]]

    ds = [GO_rest, GeneID]
    d = {}
    for k in GO_rest.keys():
        d[k] = list(d[k] for d in ds)

    dd = []
    
    for k, v in d.items():
        v[1] = [i for i in v[1] if str(i) != 'NaN']

        if v[0][2] == 'cellular_component':
            v[0][2] = 'CC'
        if v[0][2] == 'molecular_function':
            v[0][2] = 'MF'
        if v[0][2] == 'biological_process':
            v[0][2] = 'BP'

        dd.append([v[0][0], v[0][1], v[0][2], v[0][3],
                    v[1], len(v[1])])

    finalDF = pd.DataFrame(dd, columns=[
                           'ID', 'Term', 'Ontology', 'Function', 'Gene', 'Gene_length'])

    return finalDF, bg_count


def fdr_calc(x):
    """
    Assumes a list or numpy array x which contains p-values for multiple tests
    Copied from p.adjust function from R  
    """
    o = [i[0] for i in sorted(enumerate(x), key=lambda v:v[1],reverse=True)]
    ro = [i[0] for i in sorted(enumerate(o), key=lambda v:v[1])]
    q = sum([1.0/i for i in range(1,len(x)+1)])
    l = [q*len(x)/i*x[j] for i,j in zip(reversed(range(1,len(x)+1)),o)]
    l = [l[k] if l[k] < 1.0 else 1.0 for k in ro]
    return l

def enrichGO(gdata, file):

    log.info("Fetching Gene Ontology from Biomart")
    

    df, background_count = preprocessBioMart(gdata)
    log.info(f"Performing GO enrichment analysis on {file}")  

    df_goList = df[['ID', 'Gene']].values.tolist()

    go_dict = {}

    for value in df_goList:
        go_dict[value[0]] = str(value[1]).upper()

    count = df[['ID', 'Gene_length']].values.tolist()
    go_count = {}

    for c in count:
        go_count[c[0]] = c[1]

    df_List = df[['ID', 'Term', 'Ontology', 'Function']].values.tolist()
    KOdescription = {}

    for line in df_List:

        KOdescription[line[0]] = [line[1], line[2], line[3]]

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
    bg_gene_count = background_count

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
        bg_gene_GO_ids = gene_GO_count[k]
        bg_in_genome = bg_gene_count - mapped_query_ids - (gene_in_category + gene_not_in_catgory_but_in_genome) \
            + gene_in_category
        gene_ids = get_gene_ids_from_user[k]
        gID = ""

        for g in gene_ids:
            gID += g+"/"

        gID = gID.rsplit("/", 1)[0]
        pvalue = stats.hypergeom.sf(
            gene_in_category -1, bg_gene_count, gene_GO_count[k], mapped_query_ids)
        

        if gene_in_category > 0:
            pvalues.append(pvalue)

            enrichment_result.append([k, KOdescription[k][0], KOdescription[k][1], KOdescription[k][2],
                                    f"{gene_in_category}/{mapped_query_ids}", f"{go_count[k]}/{bg_gene_count}", pvalue, len(gene_ids), gID])
    
    fdr = list(fdr_calc(pvalues))

    a = [i for i in fdr if i <= 0.05]
    end = pd.DataFrame(enrichment_result)
    end.columns = ['GO ID', 'GO Term', 'Ontology', 'Definition', 'GeneRatio', 'BgRatio','Pvalues', 'Counts', 'Genes' ]
    
    end.insert(7, 'FDR', fdr)

    return end

