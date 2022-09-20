#!/usr/bin/env python
'''
Title: This script contains quality trimming tools for pySeqRNA
Author: Naveen Duhan
Version: 0.1
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


def getNreads(file, rdict, sp, paired=False):

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

def getAligned_reads(file, rdict, sp):

    aligned=0

    for read in pysam.AlignmentFile(file,'rb'):
        if read.is_unmapped == False and read.is_secondary==False:
            aligned += 1
    rdict[sp] = int(aligned)
    log.info(f"{aligned} input reads aligned {sp}")
    return rdict


def getUniquely_mapped(file, rdict, sp):

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
    

def getMulti_mapped(file,rdict, sp):

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
    
def sort_bam(file):
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

def index_bam(file):
   
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


def align_stats(sampleDict=None,trimDict=None, bamDict=None,ribodict=None, pairedEND=False):
        
    manager = multiprocessing.Manager()
    Ireads = manager.dict()
    Nreads = manager.dict()
    Rreads = manager.dict()
    Areads = manager.dict()
    Ureads = manager.dict()
    Mreads = manager.dict()
    processes = []
    
    for sp in sampleDict:
        try:
            if pairedEND:
                p=multiprocessing.Process(target= getNreads, args=(sampleDict[sp][2],Ireads, sp,True,))
            else:
                p=multiprocessing.Process(target= getNreads, args=(sampleDict[sp][2],Ireads, sp,))
        
            processes.append(p)
            p.start()
            
            for process in processes:
                
                process.join()
            
        except Exception:
            log.error(f"Not able to count Input read number in {sp}")
   
    
    for tf in trimDict:
        try:
            if pairedEND:
                p=multiprocessing.Process(target= getNreads, args=(trimDict[tf][2],Nreads, tf,True))
            else:
                p=multiprocessing.Process(target= getNreads, args=(trimDict[tf][2],Nreads, tf,))
        
            processes.append(p)
            p.start()
            
            for process in processes:
                process.join()
       
        except Exception:
            log.error(f"Not able to count Trim read number in {tf}")

    for bf in bamDict:
            p=multiprocessing.Process(target=sort_bam, args=(bamDict[bf][2],))
        
            processes.append(p)
            p.start()
            
            for process in processes:
                process.join()
        
    for bf in bamDict:   
        file = bamDict[bf][2].split(".bam")[0] + "_sorted.bam"
       
        p=multiprocessing.Process(target=index_bam, args=(file,))
        processes.append(p)
        p.start()
        
        for process in processes:
            process.join()

    for bf in bamDict:
        try:
            p=multiprocessing.Process(target= getAligned_reads, args=(bamDict[bf][2].split(".bam")[0] + "_sorted.bam",Areads, bf,))
        
            processes.append(p)
            p.start()
            
            for process in processes:
                process.join()

        except Exception:
            log.error(f"Not able to count Aligned read number in {bf}")
        try:
            p=multiprocessing.Process(target= getUniquely_mapped, args=(bamDict[bf][2].split(".bam")[0] + "_sorted.bam",Ureads, bf,))
        
            processes.append(p)
            p.start()
            
            for process in processes:
                process.join()
                
        except Exception:
            log.error(f"Not able to count Uniquely mapped read number in {bf}")
        try:
            p=multiprocessing.Process(target= getMulti_mapped, args=(bamDict[bf][2].split(".bam")[0] + "_sorted.bam",Mreads, bf,))
        
            processes.append(p)
            p.start()
            
            for process in processes:
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


