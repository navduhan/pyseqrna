
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
:Title: This modules contains read align class functions for pySeqRNA

:Created : July 21, 2021

:Author : Naveen Duhan
"""

from cmath import e
import os
from posixpath import join
import sys
import shutil
import subprocess
import pkg_resources
from pyseqrna.pyseqrna_utils import PyseqrnaLogger
from pyseqrna import pyseqrna_utils as pu


## Initialize logger

log = PyseqrnaLogger(mode="a",log='aligner')

class STAR_Aligner:
    """Class for STAR alignment program

    :param configFile: Path to STAR config file. This file will used to get the parameters for STAR alignment program

    :param slurm: To run commands with slurm task-scheduler.
    """

    def __init__(self, genome=None, configFile=None, outDir=None, slurm=False):
        """Default arguments for class STAR_Aligner

            :param genome: Reference genome file. 
            :param config_file: Parameter file for STAR. Default is from default parameters.
            :param outdir: Output directory to store alignement results. Default current working directory.
            :param slurm: To execute jobs with slurm scheduler.  
        """

        self.genome = genome

        self.slurm = slurm

        self.outDir = outDir

        if configFile != None:

            try:

                self.config = pu.parse_config_file(configFile)

                log.info(f"Using  config file {configFile}")
                
            except Exception:

                log.error("Please provide a valid config file")
        else:
            stream = pkg_resources.resource_stream('pyseqrna', "param/STAR.ini")

            self.config = pu.parse_config_file(stream.name)
           
            log.info("Using default config file STAR.ini")

        return
    
    def build_index(self,  mem= 20, tasks = 1, cpu= 8 , gff=None,  dep= ""):

        """This function build geneome index for read alignment

            :param mem: Provide memory in GB to use. Defaults to 20.

            :param tasks: Number of cpu-tasks to run. Defaults to 1.

            :param gff: Gene feature file to index with genome. Defaults to None.

            :param cpu: Total number of threads to use. Default 8.

            :param dep: slurm job id on which this job depends.  Defaults to ''.
        """
        const = ['runThreadN', '-n'] #require when change number of CPUs

        if not self.slurm:
            
            config = pu.replace_cpu(self.config['index'], const)

        else:

            config = self.config['index']

        directory = str(config[0]).split(" ")[1]

        if os.path.exists(self.outDir):

            output1 =  os.path.join(self.outDir,directory)

            output = pu.make_directory(output1)

        else:
            
            output = pu.make_directory(directory)
        # check genome file extension
        genome_extensions = tuple(['.fa', '.fasta', '.fna', '.fa.gz', 'fna.gz', '.fasta.gz'])
        
        if self.genome.endswith(genome_extensions):

            if self.genome.endswith(".gz"):

                os.system(' '.join(["gunzip", self.genome]))

            os.system(' '.join(["cp", self.genome, output]))

            log.info(f"{pu.get_basename(self.genome)} copied successfully in {directory}")

        else:

            log.error(f"Please provide a valid genome fasta file with these {genome_extensions} extensions")

            sys.exit(1)

        GenomeFasta = os.path.join(output, pu.get_basename(self.genome))

        arg = ' '.join(config[1:]) 

        execPATH = shutil.which('STAR') # get absolute path of STAR

        if execPATH is None:

            log.error("STAR aligner not found in path")

            sys.exit(1)
        
        else:

            if gff != None:

                star_command = f"{execPATH} --genomeDir {output} {arg}  --genomeFastaFiles {GenomeFasta} --sjdbGTFfile {gff}"

            else:

                star_command = f"{execPATH} --genomeDir {output} {arg}  --genomeFastaFiles {GenomeFasta} "

            if self.slurm:

                try:
                
                    job_id = pu.clusterRun(job_name='star_index', sout=os.path.join(output, "star_index"), serror=os.path.join(output, "star_index"), command= star_command, mem=mem, cpu=cpu, tasks=tasks, dep=dep)

                    log.info("Job successfully submited for {} with {} for indexing".format(GenomeFasta, job_id))

                except Exception:

                    log.error("Slurm job sumission failed")
            else:

                try:
                    with open(os.path.join(output, "star_index.out"), 'w+') as fout:
                        with open(os.path.join(output, "star_index.err"), 'w+') as ferr:
                            job_id = subprocess.call(star_command, shell=True,stdout=fout,stderr=ferr)

                    log.info("Job successfully completed for {} for indexing".format(GenomeFasta))

                except Exception:
                    
                    log.error("Job sumission failed")
    
        return job_id

           
    def check_index(self):

        """Function to check if star index is valid and exists.

           :return: Return true if genome index is valid.
           """
        for k, args in self.config.items():

            if k =='index':

                directory = str(args[0]).split(" ")[1]

                output =  os.path.join(self.outDir,directory)

        files = ['chrLength.txt', 'chrNameLength.txt', 'chrName.txt', 'chrStart.txt', 'genomeParameters.txt', 'Genome'] 

        if (os.path.exists(output) and os.path.isdir(output)):

            for f in files:

                if not(os.path.isfile(os.path.join(output, f))):

                    return False

            return True

        return  False

    def run_Alignment(self, target=None, pairedEND=False, mem= 20,cpu=8, tasks=1,  dep='' ):
        """This function align reads against indexed reference genome.

            :param target: target dictionary containing sample information.

            :param pairedEND: True if samples are paired.

            :param mem: Provide memory in GB to use. Default 20 GB.

            :param tasks: Number of cpu-tasks to run. Defaults to 1.

            :param cpu: Total number of threads to use. Default 8.

            :param dep: slurm job id on which this job depends.  Defaults to ''.
        """
        const = ['runThreadN', '-n'] #require when change number of CPUs

        if not self.slurm:
            
            config = pu.replace_cpu(self.config['alignment'], const)

        else:

            config = self.config['alignment']

        directory = str(config[0]).split(" ")[1]

        genomeIndex=os.path.join(self.outDir,directory)

        output1 = os.path.join(self.outDir, "star_results" )

        if os.path.exists(self.outDir):

            output = pu.make_directory(output1)

        else:
            parent, base = os.path.split(output1)
            output = pu.make_directory(base)

        arg = ' '.join(config[1:])

        outstarLog = {}

        job_id =[]
        
        for key, sample in target.items():  # Iterate thorough total number of samples

            try:

                if pairedEND:

                    outPrefix = os.path.join(output, sample[0])

                    outstarLog[key] = [sample[0],sample[1],outPrefix + "Aligned.out.bam"]

                    inputPair = f'{pu.get_basename(sample[2])} and {pu.get_basename(sample[3])}'

                else:

                    outPrefix = os.path.join(output, sample[0])

                    outstarLog[key] = [sample[0],sample[1],outPrefix + "Aligned.out.bam"]

                    inputPair = f'{pu.get_basename(sample[2])}'

            except Exception:

                    log.error(f'Please provide a valid dictionary with samples')

            execPATH = shutil.which('STAR') # get absolute path of STAR

            if execPATH is None:

                log.error("STAR aligner not found in path")

                sys.exit(1)
            
            else:
            
                if pairedEND:

                    star_command =f"{execPATH}  --genomeDir {genomeIndex} {arg} --outFileNamePrefix {outPrefix} --readFilesIn {sample[2]} {sample[3]}"

                else:

                    star_command =f"{execPATH}  --genomeDir {genomeIndex} {arg} --outFileNamePrefix {outPrefix}  --readFilesIn {sample[2]}"
                
                if self.slurm:
                    
                    try:
                            job = pu.clusterRun(job_name='star_align',sout=os.path.join(output, "star_index.out"), serror=os.path.join(output, "star_index.err"), command=star_command, mem=mem, cpu=cpu, tasks=tasks, dep=dep)

                            job_id.append(job)

                            log.info("Job successfully submited for {} with {} for alignment".format(inputPair, job))

                    except Exception:

                            log.error("Slurm job sumission failed")

                else:

                    try:

                        with open(os.path.join(output, "star_index.out"), 'w+') as fout:
                            with open(os.path.join(output, "star_index.err"), 'w+') as ferr:
                                job = subprocess.call(star_command, shell=True,stdout=fout,stderr=ferr)

                                job_id.append(job)

                                log.info("Job successfully completed for {} for alignment".format(inputPair))

                    except Exception:
                        
                        log.exception("Job sumission failed")

        return outstarLog, job_id



class hisat2_Aligner():

    """Class for HISAT2 alignment program

    :param configFile: Path to HISAT2 config file. This file will used to get the parameters for HISAT2 alignment program.

    :param slurm: To run commands with slurm task-scheduler.
    """
    def __init__(self, genome=None, configFile=None,  outDir='pySeqRNA_results', slurm=False):

        """Default arguments for class hisat2_Aligner

            :param genome: Reference genome file. 
            :param config_file: Parameter file for HISAT2. Default is from default parameters.
            :param outdir: Output directory to store alignement results. Default current working directory.
            :param slurm: To execute jobs with slurm scheduler.  
        """

        self.genome = genome

        self.slurm = slurm

        self.outDir = outDir

        if configFile != None:

            try:

                self.config = pu.parse_config_file(configFile)

                log.info(f"Using  config file {configFile}")
                
            except Exception:

                log.error("Please provide a valid config file")
        else:
            stream = pkg_resources.resource_stream('pyseqrna', "param/hisat2.ini")

            self.config = pu.parse_config_file(stream.name)
           
            log.info("Using default config file hisat2.ini")
        return
    
    def build_index(self, mem=8, tasks=1, cpu= 8, dep=''):

        """This function build geneome index for read alignment

            :param mem: Provide memory in GB to use. Defaults to 20.

            :param tasks: Number of cpu-tasks to run. Defaults to 1.

            :param gff: Gene feature file to index with genome. Defaults to None.

            :param cpu: Total number of threads to use. Default 8.

            :param dep: slurm job id on which this job depends.  Defaults to ''.
        """
        const = ['--threads','-p']

        if not self.slurm:
            
            config = pu.replace_cpu(self.config['index'], const)

        else:

            config = self.config['index']

        directory = str(config[0]).split(" ")[0]

        indexName = str(config[1]).split(" ")[0]

        if os.path.exists(self.outDir):

            output1 =  os.path.join(self.outDir,directory)

            output = pu.make_directory(output1)

        else:
            
            output = pu.make_directory(directory)

        genome_extensions = tuple(['.fa', '.fasta', '.fna', '.fa.gz', 'fna.gz', '.fasta.gz'])
        
        if self.genome.endswith(genome_extensions):

            if self.genome.endswith(".gz"):

                os.system(' '.join(["gunzip", self.genome]))

            os.system(' '.join(["cp", self.genome, output]))

            log.info(f"{pu.get_basename(self.genome)} copied successfully in {directory}")

        else:

            log.error(f"Please provide a valid genome fasta file with these {genome_extensions} extensions")

            sys.exit(1)


        if indexName != 'NA':

            basename = os.path.join(output, indexName)

        else:

            basename = os.path.join(output,pu.get_basename(self.genome))

        GenomeFasta = os.path.join(output, pu.get_basename(self.genome))

        arg = ' '.join(config[2:])

        execPATH = shutil.which('hisat2-build') # get absolute path of hisat2

        if execPATH is None:

            log.error("hisat2-build not found in path")
            sys.exit(1)
        else:
            hisat2_command = f"{execPATH} {arg} {GenomeFasta} {basename} "

            if self.slurm:
                try:

                    job_id = pu.clusterRun(job_name='hisat2-build',sout=os.path.join(output, "hisat2_build.out") , serror= os.path.join(output, "hisat2_build.err"), command= hisat2_command, mem=mem, cpu=cpu, tasks=tasks, dep=dep)

                    log.info("Job successfully submited for {} with {} for indexing".format(GenomeFasta, job_id))

                except Exception:

                    log.error("Slurm job sumission failed")
            else:

                try:

                    with open(os.path.join(output, "hisat2_build.out"), 'w+') as fout:
                            with open(os.path.join(output, "hisat2_build.err"), 'w+') as ferr:
                                job_id = subprocess.call(hisat2_command, shell=True,stdout=fout,stderr=ferr)

                                log.info("Job successfully submited for {} for indexing".format(GenomeFasta))

                except Exception:
                    
                    log.error("Job sumission failed")


        return job_id

    
    def check_index(self, largeIndex=False):

        """Function to check if star index is valid and exists.

            :param: True if genome indexed with large index.

            :return: Return true if genome index is valid.
        """
        for k, args in self.config.items():

            if k =='index':

                directory = str(args[0]).split(" ")[0]

                output =  os.path.join(self.outDir , directory)
        
        files = [r'.1.ht', r'.2.ht', r'.3.ht', r'.4.ht', r'.5.ht', r'.6.ht', r'.7.ht', r'.8.ht']

        
        if (os.path.exists(output) and os.path.isdir(output)):
            
            if largeIndex:

                for f in files:

                    if not(os.path.join(output, _) for _ in os.listdir(output) if _.endswith(f+"l")):

                        return False
            else:

                for f in files:

                    if not(os.path.join(output, _) for _ in os.listdir(output) if _.endswith(f)):

                        return False

            return True

        return  False

    
    def run_Alignment(self, target=None, pairedEND=False,  mem= 20, cpu=8, tasks=1, dep=''):

        """This function align reads against indexed reference genome.

            :param target: target dictionary containing sample information.

            :param pairedEND: True if samples are paired.

            :param mem: Provide memory in GB to use. Default 20 GB.

            :param tasks: Number of cpu-tasks to run. Defaults to 1.

            :param cpu: Total number of threads to use. Default 8.

            :param dep: slurm job id on which this job depends.  Defaults to ''.
        """

        consta = ['--threads', '-p']

        if self.slurm:

            config = self.config['alignment']
            
        else:
            config = pu.replace_cpu(self.config['alignment'], consta)
            
        
        directory = str(config[0])

        reference = str(config[1])

        genomeIndex = os.path.join(self.outDir,directory,reference)

        if os.path.exists(self.outDir):

            output1 =  os.path.join(self.outDir,"hisat2_results")

            output = pu.make_directory(output1)

        else:
            
            output = pu.make_directory("hisat2_results")

        arg = ' '.join(config[2:])

        outhisat2 = {}
        outsummary = {}
        job_id = []

        for key, sample in target.items():  # Iterate thorough total number of samples

            if pairedEND:

                outPrefix = os.path.join(output, sample[0])

                outBAM = outPrefix + "hisat2.bam"

                outhisat2[key] = [sample[0], sample[1], outBAM]

                summary = outPrefix + "hisat2.txt"

                outsummary [key] = summary

            else:
              
                outPrefix = os.path.join(output, sample[0])

                outBAM = outPrefix + "hisat2.bam"

                outhisat2[key] = [sample[0], sample[1], outBAM]

                summary = outPrefix + "hisat2.txt"

                outsummary [key] = summary

            
            execPATH = shutil.which('hisat2')  # get absolute path of hisat2
            
            if execPATH is None:

                log.error("hisat2 aligner not found in path")
                sys.exit(1)
            
            else:
            

                if pairedEND:

                    hisat2_command = f"{execPATH} -x {genomeIndex}  {arg} -1 {sample[2]} -2 {sample[3]} --summary-file {summary} | samtools view -Sbh >{outBAM}"

                else:

                    hisat2_command = f"{execPATH} -x {genomeIndex} {arg} -U {sample[2]} --summary-file {summary} | samtools view -Sbh >{outBAM}"

                
                print(hisat2_command)
            
                if self.slurm:
                    try:
                        job = pu.clusterRun(job_name='hisat2_align', sout=os.path.join(output, "hisat2.out") , serror=os.path.join(output, "hisat2.err") ,command= hisat2_command, mem=mem, cpu=cpu, tasks=tasks, dep=dep)
                        job_id.append(job)
                        log.info("Job successfully submited for {} with {} for alignment".format(outPrefix, job))

                    except Exception:

                        log.error("Slurm job sumission failed")

                else:

                    try:
                        with open(os.path.join(output, "hisat2.out"), 'w+') as fout:
                            with open(os.path.join(output, "hisat2.err"), 'w+') as ferr:
                                job = subprocess.call(hisat2_command, shell=True,stdout=fout,stderr=ferr)
                                job_id.append(job)
                                log.info("Job successfully completed for {} for alignment".format(outPrefix))

                    except Exception:
                        
                        log.exception("Job sumission failed")
            

        return  outhisat2, job_id
