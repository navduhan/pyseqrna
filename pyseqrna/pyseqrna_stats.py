#!/usr/bin/env python
'''
Title: This script contains quality trimming tools for pySeqRNA
Author: Naveen Duhan
Version: 0.1
'''

import pandas as pd 
import pysam
import pyfastx
import pandas as pd
from pyseqrna.pyseqrna_utils import PyseqrnaLogger
from pyseqrna import pyseqrna_utils as pu
import matplotlib as plt

log = PyseqrnaLogger(mode='a', log="stats")


def getNreads(file):
    """
        Get total number of reads in fastq file
    Args:
        file ([type]): fastq file

    Returns:
        [type]: total number of reads
    """
    result = len(pyfastx.Fastq(file))
    # result = subprocess.check_output(f"gzcat {file} | echo $((`wc -l`/4))", shell=True).decode('utf-8').rstrip()
    # logger.info(f"{result} input reads in {file}")
    
    return int(result)

def getAligned_reads(file):

    aligned=0

    for read in pysam.AlignmentFile(file,'rb'):
        if read.is_unmapped == False and read.is_secondary==False:
            aligned += 1

    return aligned


def getUniquely_mapped(file):

    uniq_mapped = 0
    for read in pysam.AlignmentFile(file,'rb'):
    
        try:
            if read.get_tag('NH')==1 and read.is_secondary==False:
                uniq_mapped += 1
        except:
            pass

    return uniq_mapped

def getMulti_mapped(file):

    multi_mapped = 0 

    for read in pysam.AlignmentFile(file,'rb'):
        try:
            if read.get_tag('NH') > 1 and read.is_secondary==False:
                multi_mapped += 1
        except:
            pass

    return multi_mapped


def align_stats(sampleDict=None,trimDict=None, bamDict=None,ribodict=None, pairedEND=False):

    Ireads = {}
    Nreads = {}
    Rreads = {}
    Areads = {}
    Ureads = {}
    Mreads = {}
    
    for sp in sampleDict:
        try:
            if pairedEND:
                Ireads[sp] = getNreads(sampleDict[sp][2])*2
            else:
                Ireads[sp] = getNreads(sampleDict[sp][2])
        except Exception:
            log.error(f"Not able to count Input read number in {sp}")
    # for rr in ribodict:
    #     try:
    #         if pairedEND:
    #             Rreads[sp] = getNreads(sampleDict[rr][2])*2
    #         else:
    #             Rreads[sp] = getNreads(sampleDict[rr][2])
    #     except Exception:
    #         log.error(f"Not able to count Input read number in {rr}")
    
    for tf in trimDict:
        try:
            if pairedEND:

                Nreads[tf] = getNreads(trimDict[tf][2])*2
            else:
                Nreads[tf] = getNreads(trimDict[tf][2])
        except Exception:
            log.error(f"Not able to count Trim read number in {tf}")

    for bf in bamDict:

        pysam.sort(bamDict[bf][2], "-o", bamDict[bf][2].split(".bam")[0] + "_sorted.bam")
        
        pysam.index(bamDict[bf][2].split(".bam")[0] + "_sorted.bam")
        
        try:
            Areads[bf]= getAligned_reads(bamDict[bf][2].split(".bam")[0] + "_sorted.bam")
        except Exception:
            log.error(f"Not able to count Aligned read number in {bf}")
        try:
            Ureads[bf] = getUniquely_mapped(bamDict[bf][2].split(".bam")[0] + "_sorted.bam")
        except Exception:
            log.error(f"Not able to count Uniquely mapped read number in {bf}")
        try:
            Mreads[bf] = getMulti_mapped(bamDict[bf][2].split(".bam")[0] + "_sorted.bam")
        except Exception:
            log.error(f"Not able to count Multi mapped read number in {bf}")
    try:
        total = []
        for k in Ireads:

            total.append([k,Ireads[k],Nreads[k],round(Nreads[k]/Ireads[k]*100,2),
                    Areads[k],round(Areads[k]/Nreads[k]*100,2),Ureads[k],round(Ureads[k]/Areads[k]*100,2),Mreads[k], round(Mreads[k]/Areads[k]*100,2)])
        if pairedEND:

            totalDF = pd.DataFrame(total,columns=['Sample','Input_reads2x','Cleaned2x', '%_Cleaned2x','Aligned','%_Aligned','Uniquely_mapped','%_Uniquely_mapped','Multi_mapped', '%_Multi_mapped'])
        else:
            totalDF = pd.DataFrame(total,columns=['Sample','Input_reads','Cleaned', '%_Cleaned','Aligned','%_Aligned','Uniquely_mapped','%_Uniquely_mapped','Multi_mapped', '%_Multi_mapped'])
    except Exception:
        log.error(f"Not able to generate align stats")


    return totalDF


