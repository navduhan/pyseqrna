#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
:Title: This modules contains feature counts in aligned reads for pySeqRNA

:Created: July 29, 2021

:Author: Naveen Duhan
'''

import os
import re
import sys
import shutil
import subprocess
import pkg_resources
from pyseqrna.pyseqrna_utils import PyseqrnaLogger
from pyseqrna import pyseqrna_utils as pu
import pandas as pd

log = PyseqrnaLogger(mode="a", log="quant")


def featureCount(configFile=None, bamDict=None, gff=None, slurm=False,  mem=8, cpu=8, tasks=1, outDir=".", dep=''):
    """
    This function counts feature in the aligned BAM files using featureCounts tool.

    :param configFile: Paramters file for featureCounts.
    
    :param bamDict: A dictionary containing all the aligned BAM files.

    :param gff: Gene feature file GFF/GTF

    :param slurm: True if SLURM scheduling is available. Defaults to False

    :param mem: Memory in GB. Defaults to 10

    :param cpu: Total number of CPU to use per task. Defaults to 8

    :param task: Number of tasks per job. Defaults to 1

    :param outDir: Output directory. Defaults to present working directory.
     
    :param dep: Slurm job id if depends on other job. Defaults to ''

    :returns: DataFrame

    :rtype: A DataFrame containing read counts per feature per sample.

    """
    if configFile != None:

        try:

            config = pu.parse_config_file(configFile)
            
        except Exception:

            log.error("Please provide a valid config file")

    else:
        stream = pkg_resources.resource_stream('pyseqrna', "param/featureCount.ini")
        config = pu.parse_config_file(stream.name)

        log.info("Using default config file featureCount.ini")

    ext = pu.get_file_extension(gff)

    if ext == 'GTF':
        config = pu.change_attribute(config)
    else:
        config = config

   
    inputFiles=""
    col = ['Gene']

    for k, sample in bamDict.items():
        col.append(sample[0])

        inputFiles += sample[2] + " "

    for key, args in config.items():

        outFile = os.path.join(outDir,"Counts.txt")

        arg = ' '.join(args[0:])


    execPATH = shutil.which('featureCounts')  # get absolute path of featureCounts
    
    if execPATH is None:

        log.error("featureCount not found in PATH")
        sys.exit()
    else:

        featureCountCmd = f"{execPATH}  -a {gff} -o {outFile} {arg} {inputFiles}"
        # print(featureCountCmd)
        if slurm:
            try:
                    job_id = pu.clusterRun(job_name='featureCount', sout=os.path.join(outDir, "featureCount.out"), serror=os.path.join(outDir, "featureCount.err"), command= featureCountCmd, mem=mem, cpu=cpu, tasks=tasks, dep=dep)

                    log.info("Job successfully submited  for feature count with jobid {} ".format (job_id))

            except Exception:

                    log.error("Slurm job sumission failed")

        else:

            try:
                with open(os.path.join(outDir, "featureCount.out"), 'w+') as fout:
                    with open(os.path.join(outDir, "featureCount.err"), 'w+') as ferr:
                        job_id = subprocess.call(featureCountCmd, shell=True,stdout=fout,stderr=ferr)

                        log.info("Job successfully completed for feature count")

            except Exception:
                
                log.exception("Job sumission failed")
    
        df= pd.read_csv(outFile, sep="\t", comment="#", low_memory=False)
        newCountDF = df.drop(columns=["Chr", "Start", "End", "Strand","Length"])
        newCountDF.columns = col
        newCountDF['Gene']= newCountDF['Gene'].str.replace("gene:", "")
        newCountDF['Gene']= newCountDF['Gene'].str.replace("gene-", "")
        newCountDF.to_excel(os.path.join(outDir,"Raw_Counts.xlsx"), index=False)
        os.remove(outFile)
    return job_id


def htseqCount(configFile=None, bamDict=None, gff=None, slurm=False, mem=8, cpu=8, tasks=1, outDir=".", dep=''):
    """
    This function counts feature in the aligned BAM files using HTSeq.

    :param configFile: Paramters file for featureCounts.
    
    :param bamDict: A dictionary containing all the aligned BAM files.

    :param gff: Gene feature file GFF/GTF

    :param slurm: True if SLURM scheduling is available. Defaults to False

    :param mem: Memory in GB. Defaults to 10

    :param cpu: Total number of CPU to use per task. Defaults to 8

    :param task: Number of tasks per job. Defaults to 1

    :param outDir: Output directory. Defaults to present working directory.
     
    :param dep: Slurm job id if depends on other job. Defaults to ''

    :returns: DataFrame

    :rtype: A DataFrame containing read counts per feature per sample.

    """
    if configFile != None:

        try:

            config = pu.parse_config_file(configFile)
            
        except Exception:

            log.error("Please provide a valid config file")

    else:

        stream = pkg_resources.resource_stream('pyseqrna', "param/htseq.ini")
        config = pu.parse_config_file(stream.name)
        log.info("Using default config file htseq.ini")

    ext = pu.get_file_extension(gff)

    if ext == 'GTF':
        config = pu.change_attribute(config)
    else:
        config = config

    inputFiles = ""
    col = ['Gene']
    for k, sample in bamDict.items():

        inputFiles += sample[2] + " "
        col.append(k)

    # df = pd.DataFrame(columns=col)

    for key, args in config.items():

        outFile = os.path.join(outDir,"Counts.txt")

        arg = ' '.join(args[0:])

       
    try:

        execPATH = os.popen('which htseq-count').read().rstrip()  # get absolute path of htseq-count
    
    except Exception:

        log.error("htseq-count not found in PATH")

    
    htseqCountCmd = f"{execPATH} {arg}  {inputFiles} {gff} > {outFile}"


    if slurm:

        try:
            job_id = pu.clusterRun(job_name='htseqCount', sout=os.path.join(outDir, "htseq.out"), serror=os.path.join(outDir, "htseq.err"), command= htseqCountCmd, mem=mem, cpu=cpu, tasks=tasks, dep=dep)

            log.info("Job successfully submited  for feature count with jobid {} ".format (job_id))

        except Exception:

            log.error("Slurm job sumission failed")

    else:

        try:
            with open(os.path.join(outDir, "htseq.out"), 'w+') as fout:
                with open(os.path.join(outDir, "htseq.err"), 'w+') as ferr:
                    job_id = subprocess.call(htseqCountCmd, shell=True,stdout=fout,stderr=ferr)

                    log.info("Job successfully completed for feature count")

        except Exception:
            
            log.exception("Job sumission failed")
    df = pd.read_csv(outFile, sep="\t")
    df.columns = col    
    newCountDF = df.iloc[:-5]
    newCountDF['Gene']= newCountDF['Gene'].str.replace("gene:", "")
    newCountDF['Gene']= newCountDF['Gene'].str.replace("gene-", "")
    newCountDF.to_excel(os.path.join(outDir,"Raw_Counts.xlsx"), index=False)
    os.remove(outFile)
    return job_id
