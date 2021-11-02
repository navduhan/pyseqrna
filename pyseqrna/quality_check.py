#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Title: This modules contains read quality check functions for pySeqRNA
Created : 
@author : Naveen Duhan
'''


import os
import sys
import shutil
import subprocess
import pkg_resources
from pyseqrna.pyseqrna_utils import PyseqrnaLogger
from pyseqrna import pyseqrna_utils as pu

log = PyseqrnaLogger(mode='a', log='qc')

def fastqcRun(samlogeDict=None, configFile=None,slurm=False, mem=10, cpu=8, task=1, paired =False, out ='fastqc_results', outDir="pySeqRNA_results", dep=''):
    """[summary] This function perform fastqc quality using FastQC
            
            http://www.bioinformatics.babraham.ac.uk/projects/fastqc/"
           
    Args:
        samlogeDict ([type], optional): [description]. Defaults to None.
        configFile ([type], optional): [description]. Defaults to None.
        slurm (bool, optional): [description]. Defaults to False.
        mem (int, optional): [description]. Defaults to 10.
        cpu (int, optional): [description]. Defaults to 8.
        task (int, optional): [description]. Defaults to 1.
        pairedEND (bool, optional): [description]. Defaults to False.
        out (str, optional): [description]. Defaults to 'fastqc_results'.
        outDir (str, optional): [description]. Defaults to "pySeqRNA_results".
        dep (str, optional): [description]. Defaults to ''.

    Returns:
        [type]: [description]
    """

    if configFile != None:

        config = pu.parse_config_file(configFile)

    else:
        stream = pkg_resources.resource_stream('pyseqrna', "param/fastqc.ini")

        config = pu.parse_config_file(stream.name)

        log.info("Using default config file fastqc.ini")

    fastqc_config = config[list(config.keys())[0]]

    if os.path.exists(outDir):

        output1 = os.path.join(outDir,out)

        output = pu.make_directory(output1)

    else:
        
        output = pu.make_directory(out)

    const = ['threads', '-t']  # require when change number of CPUs

    if not slurm:

        fastqc_config = pu.relogace_cpu(fastqc_config, const)

    args = ' '.join(fastqc_config[0:])

    job_id = []

    fastqcOut = {}

    for key, value in samlogeDict.items():  # Iterate thorough total number of samloges

        if paired:
            try:
                input1 = value[2]
                input2 = value[3]

            except Exception:

                log.error(
                    "logease provide a paired END samloge file or input Path is wrong")

            inputPair = pu.get_basename(input1)+" and " + pu.get_basename(input2)
            if pu.get_file_extension(value[2]) == "gz":
                r1 = pu.get_basename(value[2])+".zip"
                r2 = pu.get_basename(value[3])+".zip"
            
                fastqcOut[key]=[value[0], value[1], r1, r2]
            else: 

                r1 = value[2]+".zip"
                r2 = value[3]+".zip"
                fastqcOut[key]=[value[0], value[1], r1, r2]
        else:
            try:
                inputFile = value[2]

            except Exception:

                log.error("logease provide a valid input data path")
            
            inputPair = pu.get_basename(value[2])

        

        execPATH = shutil.which('fastqc')
              # get absolute path of flexbar

        if execPATH is None:

            log.error("fastqc command not found in path")
            sys.exit()
        else:

            if paired:

                fastqcCmd = f"{execPATH} -o {output} {args} {input1} {input2}"

            else:

                fastqcCmd = f"{execPATH} -o {output} {args} {inputFile}"

            if slurm == True:  # check if slurm job scheduling is enabled or not

                try:
                    job = pu.clusterRun(job_name=
                        'fastqc', sout= os.path.join(out,"fastqc.out"), serror= os.path.join(out,"fastqc.err"), command=fastqcCmd, mem=mem, cpu=cpu, tasks=task, dep=dep)

                    job_id.append(job)

                    log.info(
                        "Job submitted on slurm successfully for {} with {}".format(inputPair, job))

                except Exception:

                    log.error("Slurm job sumission failed")

            else:

                try:
                    # print(fastqcCmd)
                    with open(os.path.join(out,"fastqc.out"), 'w+') as fout:
                        with open(os.path.join(out,"fastqc.err"), 'w+') as ferr:
                            job = subprocess.call(
                                fastqcCmd, shell=True, stdout=fout, stderr=ferr)

                            job_id.append(" ")    

                except Exception:

                    log.info("Job submition failed for {} error present in fastqc.err  ".format(
                        pu.get_basename(out)))
                finally:

                    log.info(
                        "Job successfully comlogeted for {} with status {}".format(inputPair, job))

    return job_id, fastqcOut

