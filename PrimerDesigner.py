from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
import Bio.SeqUtils as sequtils
import RNA
import subprocess

def primer_design(sequence:str,start:int,end:int):
    """creates a list of forward and reverse primers for a given sequence"""
    seq_length=len(sequence)
    #forward primer
    forward_primer=[]
    if start>20:
        primer_start=start-20
    else:
        primer_start=0 
    for i in range(primer_start,start):
        forward_primer.append(Seq(sequence[i:i+20]))
    #reverse primer
    reverse_primer=[]
    if seq_length-end>20:
        primer_end=end+20

    else:
        primer_end=seq_length

    rev_comp=Seq(sequence[0:primer_end]).reverse_complement()
    for i in range(0,20):
        reverse_primer.append(rev_comp[i:i+20])
    return forward_primer,reverse_primer

def List_to_Primes(prime_list:list):
    """converts a list of sequences to a list of Primer objects"""
    return [Primer(prime) for prime in prime_list]

class Primer:
    def __init__(self,sequence:str):
        self.sequence=sequence
        self.length=len(sequence)
        self.gc_content=sequtils.gc_fraction(sequence)
        self.tm=mt.Tm_NN(sequence)
        self.secondary_structure,self.dg=RNA.fold(str(sequence))
        
    def __str__(self):
        return {"sequence":self.sequence,
                "length":self.length,
                "gc_content":self.gc_content,
                "tm":self.tm,
                "secondary_structure":self.secondary_structure,
                "dg":self.dg}
        
    def produce_secondary_structure(self):
        with open("sequence_and_structure.dbn", "w") as file:
            file.write(f">{self.sequence}\n")
            file.write(f"{primer_seq}\n")
            file.write(f"{structure}\n")
