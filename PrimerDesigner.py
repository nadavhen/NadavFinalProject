from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
import Bio.SeqUtils as sequtils
import RNA

def primer_design(sequence:str,start:int,end:int,p_length:int=20):
    """creates a list of forward and reverse primers for a given sequence"""
    seq_length=len(sequence)
    #forward primer
    forward_primer=[]
    if start>p_length:
        primer_start=start-p_length
    else:
        primer_start=0
        start=p_length
    for i in range(primer_start,start):
        forward_primer.append(Seq(sequence[i:i+p_length]))
    #reverse primer
    reverse_primer=[]
    if seq_length-end>p_length:
        primer_end=end+p_length

    else:
        primer_end=seq_length

    rev_comp=Seq(sequence[0:primer_end]).reverse_complement()
    for i in range(0,p_length):
        reverse_primer.append(rev_comp[i:i+p_length])
    return forward_primer,reverse_primer

def ListToPrimes(prime_list:list):
    """converts a list of sequences to a list of Primer objects"""
    return [Primer(prime) for prime in prime_list]


def checkBindingSites(primer_df, sequence):
    """Returns a DataFrame with the count of binding sites for the primers."""
    primer_df["binding_sites"] = 0
    rev_comp = str(Seq(sequence).reverse_complement())

    for i in primer_df.index:
        primer_seq = str(Seq(primer_df.loc[i, "sequence"]))
        binding_sites = sequence.count(primer_seq) + rev_comp.count(primer_seq)
        
        primer_df.at[i, "binding_sites"] = binding_sites

    return primer_df

def filterPrimers(primer_df, max_binding_sites=1, max_tm=60, min_tm=50, min_dg=-5, min_gc=35, max_gc=60):
    """Filters the primers based on the given parameters."""
    primer_df = primer_df[primer_df["binding_sites"] <= max_binding_sites]
    primer_df = primer_df[primer_df["tm"] <= max_tm]
    primer_df = primer_df[min_tm <= primer_df["tm"]]
    primer_df = primer_df[min_dg <= primer_df["dg"]]
    primer_df = primer_df[primer_df["gc_content"] <= max_gc]
    primer_df = primer_df[min_gc <= primer_df["gc_content"]]

    return primer_df
    

class Primer:
    def __init__(self,sequence:str):
        self.sequence=sequence
        self.length=len(sequence)
        self.gc_content=sequtils.gc_fraction(sequence)*100
        self.tm=mt.Tm_NN(sequence)
        self.secondary_structure,self.dg=RNA.fold(str(sequence))
        
    def __str__(self):
        string={"sequence":self.sequence,
                "length":self.length,
                "gc_content":self.gc_content,
                "tm":self.tm,
                "secondary_structure":self.secondary_structure,
                "dg":self.dg}
        return str(string)
    def to_dict(self):
        return {"sequence":str(self.sequence),
                "length":self.length,
                "gc_content":self.gc_content,
                "tm":self.tm,
                "secondary_structure":self.secondary_structure,
                "dg":self.dg}
    
