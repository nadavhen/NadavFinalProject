from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
import Bio.SeqUtils as sequtils
import RNA

def primer_design(sequence: str, start: int, end: int, p_length: int = 20):
    """
    Creates a list of forward and reverse primers for a given sequence.
    - sequence: The DNA sequence for which to design primers.
    - start: The start position in the sequence.
    - end: The end position in the sequence.
    - p_length: The desired length of the primers (default 20).
    Returns a tuple containing two lists: forward primers and reverse primers.
    """
    seq_length = len(sequence)
    forward_primer = []
    # Calculate the starting point for forward primer generation
    if start > p_length:
        primer_start = start - p_length
    else:
        primer_start = 0
        start = p_length
    # Generate forward primers
    for i in range(primer_start, start):
        forward_primer.append(Seq(sequence[i:i + p_length]))

    reverse_primer = []
    # Calculate the ending point for reverse primer generation
    if seq_length - end > p_length:
        primer_end = end + p_length
    else:
        primer_end = seq_length
    # Generate reverse primers
    rev_comp = Seq(sequence[0:primer_end]).reverse_complement()
    for i in range(0, p_length):
        reverse_primer.append(rev_comp[i:i + p_length])

    return forward_primer, reverse_primer

def ListToPrimes(prime_list: list):
    """
    Converts a list of sequences to a list of Primer objects.
    - prime_list: List of primer sequences (str).
    Returns a list of Primer objects.
    """
    return [Primer(prime) for prime in prime_list]

def checkBindingSites(primer_df, sequence):
    """
    Returns a DataFrame with the count of binding sites for the primers in the given sequence.
    - primer_df: DataFrame containing primer sequences.
    - sequence: The DNA sequence to check for binding sites.
    Returns the input DataFrame with an added column for binding site counts.
    """
    primer_df["binding_sites"] = 0
    rev_comp = str(Seq(sequence).reverse_complement())

    for i in primer_df.index:
        primer_seq = str(Seq(primer_df.loc[i, "sequence"]))
        # Count occurrences in both sequence and its reverse complement
        binding_sites = sequence.count(primer_seq) + rev_comp.count(primer_seq)
        primer_df.at[i, "binding_sites"] = binding_sites

    return primer_df

def filterPrimers(primer_df, max_binding_sites=1, max_tm=60, min_tm=50, min_dg=-5, min_gc=35, max_gc=60):
    """
    Filters primers based on binding sites, melting temperature (Tm), Gibbs free energy (ΔG), and GC content criteria.
    - primer_df: DataFrame containing primers and their properties.
    - max_binding_sites: Maximum allowed binding sites (default 1).
    - max_tm: Maximum allowed melting temperature (default 60°C).
    - min_tm: Minimum allowed melting temperature (default 50°C).
    - min_dg: Minimum allowed Gibbs free energy (default -5 kcal/mol).
    - min_gc: Minimum allowed GC content (default 35%).
    - max_gc: Maximum allowed GC content (default 60%).
    Returns a filtered DataFrame based on the specified criteria.
    """
    primer_df = primer_df[primer_df["binding_sites"] <= max_binding_sites]
    primer_df = primer_df[primer_df["tm"] <= max_tm]
    primer_df = primer_df[min_tm <= primer_df["tm"]]
    primer_df = primer_df[min_dg <= primer_df["dg"]]
    primer_df = primer_df[primer_df["gc_content"] <= max_gc]
    primer_df = primer_df[min_gc <= primer_df["gc_content"]]

    return primer_df

class Primer:
    def __init__(self, sequence: str):
        """
        Initializes a Primer object with properties derived from the given sequence.
        - sequence: The primer DNA sequence (str).
        Properties include sequence length, GC content, melting temperature (Tm), secondary structure, and Gibbs free energy (ΔG).
        """
        self.sequence = sequence
        self.length = len(sequence)
        self.gc_content = sequtils.gc_fraction(sequence) * 100
        self.tm = mt.Tm_NN(sequence)
        self.secondary_structure, self.dg = RNA.fold(str(sequence))
        
    def __str__(self):
        """
        Returns a string representation of the Primer object's properties.
        """
        string = {"sequence": self.sequence,
                  "length": self.length,
                  "gc_content": self.gc_content,
                  "tm": self.tm,
                  "secondary_structure": self.secondary_structure,
                  "dg": self.dg}
        return str(string)
    
    def to_dict(self):
        """
        Converts the Primer object's properties to a dictionary.
        Returns a dictionary representation of the Primer object.
        """
        return {"sequence": str(self.sequence),
                "length": self.length,
                "gc_content": self.gc_content,
                "tm": self.tm,
                "secondary_structure": self.secondary_structure,
                "dg": self.dg}