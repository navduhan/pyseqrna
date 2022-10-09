#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
:Title: This module count multimapped read groups in aligned files.

:Created: June 5, 2022

:Author: Naveen Duhan
'''
from numpy import append
import pysam
import re
import functools
import logging
import pandas as pd
import subprocess
from pyseqrna.pyseqrna_utils import PyseqrnaLogger
from pyseqrna import pyseqrna_utils as pu

log = PyseqrnaLogger(mode='a', log='mmgg')

def _Bam2Bed(bamDict):

    outBed = []
    alignFiles = []
    for k , v in bamDict.items():
        alignFiles.append(v[2])

    for sample in alignFiles:

        out = sample.split(".bam")[0]+".bed"

        outBed.append(out)
        input = sample.split('.bam')[0]+"_sorted.bam"
        # pysam.sort("-o",f"{sample.split('.bam')[0]}_sorted.bam", sample)
        reads = pysam.AlignmentFile(input,'rb')
       

        with open(out, 'w') as fp:

            for read in reads:

                try:
               
                    NH = read.get_tag('NH')
                    
                except:
                    pass
               
                if read.is_reverse:

                    strand= "-"

                else:

                    strand= "+"

                if NH >1:

                    fp.write(read.reference_name+"\t"+str(read.pos)+"\t"+str(read.aend)+"\t"+read.query_name+"\t"+str(read.flag)+"\t"+strand+"\n")
    return outBed

def _gffToBed(file=None, feature='gene'):

    ext = pu.get_file_extension(file)

    gtf = pd.read_csv(file,sep="\t", header=None, comment="#")

    gtf.columns = ['seqname', 'source', 'feature', 'start', 'end', 's1','strand', 's2', 'identifier']

    gtf = gtf[gtf['feature'] == feature]


    out=file+".bed"

    try:
  
        if ext =='.gff3' or ext =='.gff':

            gtf['gene'] = list(map(lambda x: re.search(r'ID=(.*?);',x,re.MULTILINE).group(1),gtf['identifier'].values.tolist()))

        elif ext =='.gtf':

            gtf['gene'] = list(map(lambda x: re.search(r'gene_id\s+["](.*?)["]',x,re.MULTILINE).group(1),gtf['identifier'].values.tolist()))
    except Exception:

        log.error(f"Feature file is invalid Attribute should start with gene_id or ID")
  
    gene_list = gtf.values.tolist()


    with open(out, 'w') as fp:

        for gene in gene_list:
            if gene[0] not in ("Mt", "Pt"):
                fp.write(f"{gene[0]}\t{gene[3]-1}\t{gene[4]}\t{gene[9]}\t.\t{gene[6]}\n")     
    return out

def _filterGenes(mmg=None):
    col = mmg.columns
    df = mmg['Gene'].values.tolist()
    substrings = {w1 for w1 in df for w2 in df if w1 in w2 and w1 != w2}
    rgenes = list(set(df) - substrings)
    fgene = []
    for r in rgenes:
        if "-" in r:
            fgene.append(r)
    fgenes = pd.DataFrame(fgene, columns=['Gene'])
    final = fgenes.merge(mmg, on='Gene')

    final = final[col]
    mmgid = []
   
    for m in range(0,final.shape[0]):
        mmgid.append("MMG_"+str(m))
    final.MMG = mmgid
    return final

def countMMG(sampleDict=None,bamDict=None, gff=None, feature="gene",minCount=100, percentSample=0.5 ):

    """
    This function calculates multimapped gene groups.

    :param sampleDict: a dictionary containing samples information.

    :param bamDict: a dictionary containing BAM files.

    :param gff: gene feature file.

    :param feature: feature type.

    :param minCount: minimum number of reads per sample.

    :param percentSample: minimum number of reads in percent sample.

    :returns: DataFrame

    :rtype: A DataFrame containing multimapped read groups counts.
    
    """

    bedFiles = _Bam2Bed(bamDict=bamDict)
    gffBED = _gffToBed(file=gff, feature=feature)
    intbed=[]

    Cmd =  'sort -k 1,1 -k2,2n '+ gffBED+' >tmp && mv tmp '+ gffBED
    
    job_id = subprocess.call(Cmd, shell=True)

    for sample in bedFiles:

        out= sample.split(".bed")[0]+".intersect.bed"

        intbed.append(out)
        scmd =  'sort -k 1,1 -k2,2n '+ sample+' >tmp && mv tmp '+ sample
        job_id = subprocess.call(scmd, shell=True)
        cmd = "intersectBed -sorted -abam {} -bed -wb -b {} >{}".format(sample,gffBED, out)
        job_id = subprocess.call(cmd, shell=True)

        
    listDF=[]

    for b in intbed:

        dd = {}
        with open(b, 'r') as fp:
            lines = [line.rstrip() for line in fp]
            for line in lines:
                data = re.split("\s+", line)
                if data[3] not in dd:
                    dd[data[3]] = [data[9]]
                else:
                    if data[9] != dd[data[3]]:
                        dd[data[3]].append(data[9])

        geneKey = []
        for key in dd:
            geneKey.append('-'.join(sorted(tuple((list(dict.fromkeys(dd[key])))))))  # remove duplicated gene name from list


        df=pd.DataFrame(geneKey)

    
        df.columns=[b]
        df= pd.DataFrame(df[b].value_counts(dropna=False))
        df.reset_index(level=0, inplace=True)
        df.columns=['Gene', b]
        listDF.append(df)

    df_final = functools.reduce(lambda left,right: pd.merge(left,right, on='Gene', how='outer'), listDF).fillna(0)
    
    final=df_final.values.tolist()

    colna=df_final.columns
    newNames = ['MMG', 'Gene']
    for key, value in sampleDict.items():
        
        for c in colna:
           
            if c.split("/")[-1].startswith(value[0]):
                newNames.append(key)


    kk=[]
    count=0
    
    for f in final:
       
        a=sum(i > minCount for i in f[1:])
        k=len(f[1:])*percentSample
        if a >= k:
            count += 1
            f.insert(0,"MMG_"+str(count))
            kk.append(f)

    try:
        df = pd.DataFrame(kk).astype(str).replace('\.0', '', regex=True)
        df.columns= newNames

        df['Gene']= df['Gene'].str.replace("gene:", "")
        df['Gene']= df['Gene'].str.replace("gene-", "")

        dk = _filterGenes(mmg=df)

    except Exception:
        log.warning("No multimapped groups found")
    

    return dk
    
 