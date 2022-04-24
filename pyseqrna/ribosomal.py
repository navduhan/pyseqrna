#!/usr/bin/env python
"""

Title: This script contains sortMeRNA function for removing ribosomal RNA from reads.
Author: Naveen Duhan
Version: 0.1

"""

import os
import shutil
import glob
import sys
import subprocess
from pyseqrna.pyseqrna_utils import PyseqrnaLogger
from pyseqrna import pyseqrna_utils as pu

# Initialize logger

log = PyseqrnaLogger(mode="a", log='ribosomal')


def sortmernaRun(sampleDict=None, outDir="pySeqRNA_results", rnaDatabases=None, pairedEND=False, mem=10, cpu=8, task=1, slurm=False, dep=''):

    """
    This Function execute sortMeRNA to remove ribosomal RNA from fastq reads.

    Args:
    sampleDict ([type], optional): [description]. Defaults to None.
    outDir (str, optional): [description]. Defaults to "pySeqRNA_results".
    pairedEND (bool, optional): [description]. Defaults to False.
    mem (int, optional): [description]. Defaults to 10.
    cpu (int, optional): [description]. Defaults to 8.
    task (int, optional): [description]. Defaults to 1.
    slurm (bool, optional): [description]. Defaults to False.
    dep (str, optional): [description]. Defaults to ''.

    """

    output = "sortMeRNA_results"

    if os.path.exists(outDir):

        output1 = os.path.join(outDir, output)

        output = pu.make_directory(output1)

    else:

        output = pu.make_directory(output)

    try:

        # refDB = glob.glob("pyseqrna/example/data/sortmerna_db/rRNA_databases/*.fasta")

        refDB = glob.glob(f"{rnaDatabases}/*.fasta")

        sortmernaREF = ""

        for ref in refDB:
            sortmernaREF += ' '.join(["--ref", ref, " "])

    except Exception:

        log.error("Please provide valid fasta files for reference")

    outsortmeRNA = {}

    job_id = []

    for key, sample in sampleDict.items():

        aligned_out = os.path.join(output, ''.join([sample[0], "_aligned"]))

        filterd_out = os.path.join(output, sample[0])

        workdir = os.path.join(output, "workdir_"+sample[0])

        pu.make_directory(workdir)

        if pairedEND:

            aligned_out = os.path.join(
                output, ''.join([sample[0], "_aligned"]))

            filterd_out = os.path.join(output, sample[0])

            outsortmeRNA[key] = [filterd_out +
                                 "_fwd.fastq", filterd_out+"_rev.fastq"]

            # workdir = os.path.join(output,"workdir_"+sample[0])

            # pu.make_directory(workdir)

        else:

            outsortmeRNA[key] = filterd_out+".fastq"

        execPATH = shutil.which('sortmerna')  # get absolute path of sortmerna

        if execPATH is None:

            log.error("sortmerna not found in path")
            sys.exit()

        else:

            if pairedEND:

                sortmernaCmd = f"{execPATH} {sortmernaREF} --reads {sample[2]} --reads {sample[3]} --aligned {aligned_out} --other {filterd_out}  --fastx true --out2 true -threads {cpu} -v"

            else:

                sortmernaCmd = f"{execPATH} {sortmernaREF} --reads {sample[2]} --aligned {aligned_out} --other {filterd_out}  --fastx true -threads {cpu} -v"

            # print(sortmernaCmd)
            if slurm:

                try:
                    job = pu.clusterRun(job_name='sortmeRNA', sout=os.path.join(output, "sortmeRNA.out"), serror=os.path.join(
                        output, "sortmeRNA.err"), command=sortmernaCmd, mem=mem, cpu=cpu, tasks=task, dep=dep)

                    job_id.append(job)
                    log.info("Job successfully submited for {} with {}".format(
                        sample[0], job_id))

                except Exception:

                    log.error("Slurm job sumission failed")

            else:

                try:
                    with open(os.path.join(output, "sortmeRNA.out"), 'w+') as fout:
                        with open(os.path.join(output, "sortmeRNA.err"), 'w+') as ferr:
                            job = subprocess.call(
                                sortmernaCmd, shell=True, stdout=fout, stderr=ferr)
                            job_id.append('')
                            log.info(
                                "Job successfully submited for {} ".format(sample[0]))

                except Exception:

                    log.error("Job sumission failed")

    return outsortmeRNA, job_id
