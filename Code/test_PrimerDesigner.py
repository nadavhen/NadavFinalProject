import PrimerDesigner
import pandas as pd
from Bio.Seq import Seq

sequence="ATATCCAGCTTACCTATTAGCCCCTGATCACTTCTTTCTTCACAGCCGATGATAGTAACTAACTCAGGGCACGGGGGTTACGACTTGATCTTTCACGCAG"
primer_df=pd.DataFrame(columns=["sequence","length","gc_content","tm","secondary_structure","dg"])
res=PrimerDesigner.primer_design(sequence,0,100,20)
def test_primer_design():
    # Test if the primer_design function returns the correct forward and reverse primers
    assert(res==([Seq('ATATCCAGCTTACCTATTAG'),
                Seq('TATCCAGCTTACCTATTAGC'),
                Seq('ATCCAGCTTACCTATTAGCC'),
                Seq('TCCAGCTTACCTATTAGCCC'),
                Seq('CCAGCTTACCTATTAGCCCC'),
                Seq('CAGCTTACCTATTAGCCCCT'),
                Seq('AGCTTACCTATTAGCCCCTG'),
                Seq('GCTTACCTATTAGCCCCTGA'),
                Seq('CTTACCTATTAGCCCCTGAT'),
                Seq('TTACCTATTAGCCCCTGATC'),
                Seq('TACCTATTAGCCCCTGATCA'),
                Seq('ACCTATTAGCCCCTGATCAC'),
                Seq('CCTATTAGCCCCTGATCACT'),
                Seq('CTATTAGCCCCTGATCACTT'),
                Seq('TATTAGCCCCTGATCACTTC'),
                Seq('ATTAGCCCCTGATCACTTCT'),
                Seq('TTAGCCCCTGATCACTTCTT'),
                Seq('TAGCCCCTGATCACTTCTTT'),
                Seq('AGCCCCTGATCACTTCTTTC'),
                Seq('GCCCCTGATCACTTCTTTCT')],
                [Seq('CTGCGTGAAAGATCAAGTCG'),
                Seq('TGCGTGAAAGATCAAGTCGT'),
                Seq('GCGTGAAAGATCAAGTCGTA'),
                Seq('CGTGAAAGATCAAGTCGTAA'),
                Seq('GTGAAAGATCAAGTCGTAAC'),
                Seq('TGAAAGATCAAGTCGTAACC'),
                Seq('GAAAGATCAAGTCGTAACCC'),
                Seq('AAAGATCAAGTCGTAACCCC'),
                Seq('AAGATCAAGTCGTAACCCCC'),
                Seq('AGATCAAGTCGTAACCCCCG'),
                Seq('GATCAAGTCGTAACCCCCGT'),
                Seq('ATCAAGTCGTAACCCCCGTG'),
                Seq('TCAAGTCGTAACCCCCGTGC'),
                Seq('CAAGTCGTAACCCCCGTGCC'),
                Seq('AAGTCGTAACCCCCGTGCCC'),
                Seq('AGTCGTAACCCCCGTGCCCT'),
                Seq('GTCGTAACCCCCGTGCCCTG'),
                Seq('TCGTAACCCCCGTGCCCTGA'),
                Seq('CGTAACCCCCGTGCCCTGAG'),
                Seq('GTAACCCCCGTGCCCTGAGT')]))
f=PrimerDesigner.ListToPrimes(res[0])
sequences=[str(i) for i in f]

def test_Primers():
    # Test if ListToPrimes correctly converts sequences to Primer objects and matches expected output
    assert(sequences==["{'sequence': Seq('ATATCCAGCTTACCTATTAG'), 'length': 20, 'gc_content': 35.0, 'tm': 42.57933707458011, 'secondary_structure': '....................', 'dg': 0.0}",
 "{'sequence': Seq('TATCCAGCTTACCTATTAGC'), 'length': 20, 'gc_content': 40.0, 'tm': 45.44796883890268, 'secondary_structure': '......(((........)))', 'dg': -1.100000023841858}",
 "{'sequence': Seq('ATCCAGCTTACCTATTAGCC'), 'length': 20, 'gc_content': 45.0, 'tm': 48.12932163415337, 'secondary_structure': '.....(((........))).', 'dg': -1.899999976158142}",
 "{'sequence': Seq('TCCAGCTTACCTATTAGCCC'), 'length': 20, 'gc_content': 50.0, 'tm': 50.19875292218006, 'secondary_structure': '....(((........)))..', 'dg': -1.899999976158142}",
 "{'sequence': Seq('CCAGCTTACCTATTAGCCCC'), 'length': 20, 'gc_content': 55.00000000000001, 'tm': 51.29216297807886, 'secondary_structure': '...(((........)))...', 'dg': -1.899999976158142}",
 "{'sequence': Seq('CAGCTTACCTATTAGCCCCT'), 'length': 20, 'gc_content': 50.0, 'tm': 50.17287522080596, 'secondary_structure': '..(((........)))....', 'dg': -1.899999976158142}",
 "{'sequence': Seq('AGCTTACCTATTAGCCCCTG'), 'length': 20, 'gc_content': 50.0, 'tm': 50.17287522080602, 'secondary_structure': '.(((........))).....', 'dg': -1.899999976158142}",
 "{'sequence': Seq('GCTTACCTATTAGCCCCTGA'), 'length': 20, 'gc_content': 50.0, 'tm': 50.198752922180006, 'secondary_structure': '(((........)))......', 'dg': -1.7000000476837158}",
 "{'sequence': Seq('CTTACCTATTAGCCCCTGAT'), 'length': 20, 'gc_content': 45.0, 'tm': 47.3595705245603, 'secondary_structure': '....................', 'dg': 0.0}",
 "{'sequence': Seq('TTACCTATTAGCCCCTGATC'), 'length': 20, 'gc_content': 45.0, 'tm': 47.39292197609336, 'secondary_structure': '......(((((...))))).', 'dg': -0.10000000149011612}",
 "{'sequence': Seq('TACCTATTAGCCCCTGATCA'), 'length': 20, 'gc_content': 45.0, 'tm': 48.38516388904134, 'secondary_structure': '.....(((((...)))))..', 'dg': -0.10000000149011612}",
 "{'sequence': Seq('ACCTATTAGCCCCTGATCAC'), 'length': 20, 'gc_content': 50.0, 'tm': 50.17287522080602, 'secondary_structure': '....(((((...)))))...', 'dg': -0.10000000149011612}",
 "{'sequence': Seq('CCTATTAGCCCCTGATCACT'), 'length': 20, 'gc_content': 50.0, 'tm': 49.85366706413515, 'secondary_structure': '...(((((...)))))....', 'dg': -0.10000000149011612}",
 "{'sequence': Seq('CTATTAGCCCCTGATCACTT'), 'length': 20, 'gc_content': 45.0, 'tm': 48.03668074088381, 'secondary_structure': '..(((((...))))).....', 'dg': -0.10000000149011612}",
 "{'sequence': Seq('TATTAGCCCCTGATCACTTC'), 'length': 20, 'gc_content': 45.0, 'tm': 48.06801988702034, 'secondary_structure': '.(((((...)))))......', 'dg': -0.10000000149011612}",
 "{'sequence': Seq('ATTAGCCCCTGATCACTTCT'), 'length': 20, 'gc_content': 45.0, 'tm': 49.62460970886957, 'secondary_structure': '....................', 'dg': 0.0}",
 "{'sequence': Seq('TTAGCCCCTGATCACTTCTT'), 'length': 20, 'gc_content': 45.0, 'tm': 49.88348396272676, 'secondary_structure': '....................', 'dg': 0.0}",
 "{'sequence': Seq('TAGCCCCTGATCACTTCTTT'), 'length': 20, 'gc_content': 45.0, 'tm': 49.88348396272676, 'secondary_structure': '....................', 'dg': 0.0}",
 "{'sequence': Seq('AGCCCCTGATCACTTCTTTC'), 'length': 20, 'gc_content': 50.0, 'tm': 51.33873965085621, 'secondary_structure': '....................', 'dg': 0.0}",
 "{'sequence': Seq('GCCCCTGATCACTTCTTTCT'), 'length': 20, 'gc_content': 50.0, 'tm': 51.338739650856155, 'secondary_structure': '....................', 'dg': 0.0}"])
    


def test_checkBindingSites():
    # Test if checkBindingSites correctly counts binding sites for each primer in the sequence
    primer_df=pd.DataFrame(columns=["sequence","length","gc_content","tm","secondary_structure","dg"])
    for primer in f:
        primer_df.loc[len(primer_df)]=primer.to_dict()
    primer_df=PrimerDesigner.checkBindingSites(primer_df,sequence)
    assert(primer_df["binding_sites"].tolist()==[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

def test_filterPrimers():
    # Test if filterPrimers correctly filters primers based on specified criteria
    primer_df=pd.DataFrame(columns=["sequence","length","gc_content","tm","secondary_structure","dg"])
    for primer in f:
        primer_df.loc[len(primer_df)]=primer.to_dict()
    primer_df=PrimerDesigner.checkBindingSites(primer_df,sequence)
    primer_df=PrimerDesigner.filterPrimers(primer_df, max_binding_sites=1, max_tm=60, min_tm=50, min_dg=-5, min_gc=35, max_gc=60)
    primer_df=primer_df.values.tolist()
    assert(primer_df==[['TCCAGCTTACCTATTAGCCC',
  20,
  50.0,
  50.19875292218006,
  '....(((........)))..',
  -1.899999976158142,
  1],
 ['CCAGCTTACCTATTAGCCCC',
  20,
  55.00000000000001,
  51.29216297807886,
  '...(((........)))...',
  -1.899999976158142,
  1],
 ['CAGCTTACCTATTAGCCCCT',
  20,
  50.0,
  50.17287522080596,
  '..(((........)))....',
  -1.899999976158142,
  1],
 ['AGCTTACCTATTAGCCCCTG',
  20,
  50.0,
  50.17287522080602,
  '.(((........))).....',
  -1.899999976158142,
  1],
 ['GCTTACCTATTAGCCCCTGA',
  20,
  50.0,
  50.198752922180006,
  '(((........)))......',
  -1.7000000476837158,
  1],
 ['ACCTATTAGCCCCTGATCAC',
  20,
  50.0,
  50.17287522080602,
  '....(((((...)))))...',
  -0.10000000149011612,
  1],
 ['AGCCCCTGATCACTTCTTTC',
  20,
  50.0,
  51.33873965085621,
  '....................',
  0.0,
  1],
 ['GCCCCTGATCACTTCTTTCT',
  20,
  50.0,
  51.338739650856155,
  '....................',
  0.0,
  1]]) 