from pyseqrna import pyseqrna_utils as pu
from pyseqrna import quality_check as qc
from pyseqrna import quality_trimming as qt
from pyseqrna import ribosomal as rb
from pyseqrna import aligners as al


data = pu.read_input_file("pyseqrna/example/input_Sample_PE.txt", "pyseqrna/example/data/" , paired=True)
samples= data['samples']


a = rb.sortmernaRun(sampleDict=samples,pairedEND=True)


