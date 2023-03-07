#!/usr/bin/env python
'''
:Title: This module generate read alignment statistics

:Created: Spetember 10, 2021

:Author: Naveen Duhan
'''

from tkinter import E
from matplotlib.colors import same_color
import pandas as pd 
import pysam
import pyfastx
import pandas as pd
from pyseqrna.pyseqrna_utils import PyseqrnaLogger
from pyseqrna import pyseqrna_utils as pu
import matplotlib as plt
from multiprocessing import Pool
import multiprocessing
import subprocess 
log = PyseqrnaLogger(mode='a', log="stats")


def _getNreads(file, rdict, sp, paired=False):

    """
    
    Get total number of reads in fastq file

    file ([type]): fastq file

    [type]: total number of reads

    """
    result = len(pyfastx.Fastq(file))
    
    if paired:
        rdict[sp] = int(result)*2
        
        log.info(f"{result*2} input reads in {sp}")
    else: 

        rdict[sp] = int(result)
        
        log.info(f"{result} input reads in {sp}")
    
    return rdict

def _getTreads(file, rdict, sp, paired=False):

    """
    
    Get total number of reads in fastq file

    file ([type]): fastq file

    [type]: total number of reads

    """
    result = len(pyfastx.Fastq(file))
    
    if paired:
        rdict[sp] = int(result)*2
        
        log.info(f"{result*2} cleaned reads in {sp}")
    else: 

        rdict[sp] = int(result)
        
        log.info(f"{result} cleaned reads in {sp}")
    
    return rdict

def _getAligned_reads(file, rdict, sp):

    aligned=0

    for read in pysam.AlignmentFile(file,'rb'):
        if read.is_unmapped == False and read.is_secondary==False:
            aligned += 1
    rdict[sp] = int(aligned)
    log.info(f"{aligned} input reads aligned {sp}")
    return rdict


def _getUniquely_mapped(file, rdict, sp):

    uniq_mapped = 0
    for read in pysam.AlignmentFile(file,'rb'):
    
        try:
            if read.get_tag('NH')==1 and read.is_secondary==False:
                uniq_mapped += 1
        except:
            pass
    rdict[sp] = int(uniq_mapped)
    log.info(f"{uniq_mapped} input reads uniquely mapped in {sp}")
    return rdict
    

def _getMulti_mapped(file,rdict, sp):

    multi_mapped = 0 

    for read in pysam.AlignmentFile(file,'rb'):
        try:
            if read.get_tag('NH') > 1 and read.is_secondary==False:
                multi_mapped += 1
        except:
            pass
    rdict[sp] = int(multi_mapped)
    log.info(f"{multi_mapped} input reads multi mapped in {sp}")
    return rdict
    
def _sort_bam(file):
    outfile = file.split(".bam")[0] + "_sorted.bam"
    samtools_cmd = f'samtools sort {file} > {outfile}'
    try:
        with open("bamsort.out", 'w+') as fout:
            with open("bamsort.err", 'w+') as ferr:
                job = subprocess.call(
                    samtools_cmd, shell=True, stdout=fout, stderr=ferr)
                
                log.info(
                    "Sorting bam completed for {} ".format(file))

    except Exception:

        log.error("Bam sorting failed")

def _index_bam(file):
   
    samtools_cmd = f'samtools index -c {file}'
    try:
        with open("bam_index.out", 'w+') as fout:
            with open("bam_index.err", 'w+') as ferr:
                job = subprocess.call(
                    samtools_cmd, shell=True, stdout=fout, stderr=ferr)
                
                log.info(
                    "Bam Indexing completed for {} ".format(file))

    except Exception:

        log.error("Bam indexing failed")


def align_stats(sampleDict=None,trimDict=None, bamDict=None,riboDict=None, pairedEND=False):

    """
    This function calculates the alignment statistics

    :param sampleDict: Raw Reads sample dictionary containing all samples

    :param trinDict: Dictionary containing trimmed samples

    :param bamDict: Dictionary containing all samples bam files.

    :param riboDict: Dictionary containing filtered reads.

    :returns: DataFrame
    :rtype: A DataFrame containg alignment statistics
    """
        
    manager = multiprocessing.Manager()
    Ireads = manager.dict()
    Nreads = manager.dict()
    Rreads = manager.dict()
    Areads = manager.dict()
    Ureads = manager.dict()
    Mreads = manager.dict()
    nprocesses = []
    tprocesses = []
    bprocesses = []
    sprocesses = []
    iprocesses = []
    uprocesses = []
    mprocesses = []
    
    try:
        for sp in sampleDict:
            
            if pairedEND:
                p=multiprocessing.Process(target= _getNreads, args=(sampleDict[sp][2],Ireads, sp,True,))
            else:
                p=multiprocessing.Process(target= _getNreads, args=(sampleDict[sp][2],Ireads, sp,))
        
            nprocesses.append(p)
           
        for process in nprocesses:
            process.start()
            
        for process in nprocesses:
            
            process.join()
            
    except Exception:
        log.error(f"Not able to count Input read number in {sp}")
   
    try:
        for tf in trimDict:
        
            if pairedEND:
                p=multiprocessing.Process(target= _getTreads, args=(trimDict[tf][2],Nreads, tf,True,))
            else:
                p=multiprocessing.Process(target= _getTreads, args=(trimDict[tf][2],Nreads, tf, ))
        
            tprocesses.append(p)

            

        for process in tprocesses:
            process.start()
            
        for process in tprocesses:
            
            process.join()
       
    except Exception:
        log.error(f"Not able to count Trim read number in {tf}")

    for bf in bamDict:
        p=multiprocessing.Process(target = _sort_bam, args=(bamDict[bf][2],))
        
        sprocesses.append(p)

    for process in sprocesses:
            process.start()
            
    for process in sprocesses:
        
        process.join()
        
    for bf in bamDict:   
        file = bamDict[bf][2].split(".bam")[0] + "_sorted.bam"
       
        p=multiprocessing.Process(target = _index_bam, args=(file,))
        iprocesses.append(p)

    for process in iprocesses:
            process.start()
            
    for process in iprocesses:
        
        process.join()

    try:
        for bf in bamDict:
        
            p=multiprocessing.Process(target= _getAligned_reads, args=(bamDict[bf][2].split(".bam")[0] + "_sorted.bam",Areads, bf,))
        
            bprocesses.append(p)

        for process in bprocesses:
            process.start()
            
        for process in bprocesses:
            
            process.join()

    except Exception:
        log.error(f"Not able to count Aligned read number in {bf}")
    try:
        for bf in bamDict:
            p=multiprocessing.Process(target= _getUniquely_mapped, args=(bamDict[bf][2].split(".bam")[0] + "_sorted.bam",Ureads, bf,))
        
            uprocesses.append(p)
            p.start()
            
        for process in uprocesses:
            process.start()
            
        for process in uprocesses:
            
            process.join()
                
    except Exception:
        log.error(f"Not able to count Uniquely mapped read number in {bf}")
    try:
        for bf in bamDict:
            p=multiprocessing.Process(target= _getMulti_mapped, args=(bamDict[bf][2].split(".bam")[0] + "_sorted.bam",Mreads, bf,))
        
            mprocesses.append(p)
            p.start()
            
        for process in mprocesses:
            process.start()
            
        for process in mprocesses:
            
            process.join()       
            
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


