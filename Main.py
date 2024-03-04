import GUI as gui
import PrimerDesigner
import pandas as pd


res=gui.run_gui()
if res[0]=="file":
    sequence=""
    with open(res[1]) as file:
        for line in file:
            if line[0]!=">":
                sequence+=line.strip()
elif res[0]=="text":
    sequence=res[1]

seq_length=len(sequence)
parameters=gui.para_gui(seq_length)

#breaking parameters dictionary into variables
length=parameters["length"]
start=parameters["start"]
end=parameters["end"]
gc=parameters["gc"]
tm=parameters["tm"]
dg=parameters["dg"]
pl=parameters["pl"]
amplicon_option=parameters["amplicon_option"]

#defining the amplicon
if amplicon_option==1:
    start=0
    end=length-1
    amplicon=sequence[start:end]
elif amplicon_option==2:
    amplicon=sequence[start-1:end-1]
elif amplicon_option==3:
    start=0
    end=seq_length-1
    amplicon=sequence[start:end]

#designing the primers
final_primers=pd.DataFrame(columns=["sequence","position","length","gc_content","tm","secondary_structure","dg"])
while True:
    forward_primer,reverse_primer=PrimerDesigner.primer_design(sequence,start,end,p_length=pl)
    forward_primer=PrimerDesigner.ListToPrimes(forward_primer)
    reverse_primer=PrimerDesigner.ListToPrimes(reverse_primer)
    forward_primer_df=pd.DataFrame(columns=["sequence","length","gc_content","tm","secondary_structure","dg"])
    reverse_primer_df=pd.DataFrame(columns=["sequence","length","gc_content","tm","secondary_structure","dg"])


    for primer in forward_primer:
        forward_primer_df.loc[len(forward_primer_df)]=primer.to_dict()
    for primer in reverse_primer:
        reverse_primer_df.loc[len(reverse_primer_df)]=primer.to_dict()

    forward_primer_df=PrimerDesigner.checkBindingSites(forward_primer_df,sequence)
    reverse_primer_df=PrimerDesigner.checkBindingSites(reverse_primer_df,sequence)
    forward_primer_df=PrimerDesigner.filterPrimers(forward_primer_df, max_binding_sites=1, max_tm=tm+5, min_tm=tm-5, min_dg=dg, min_gc=gc-10, max_gc=gc+10)
    reverse_primer_df=PrimerDesigner.filterPrimers(reverse_primer_df, max_binding_sites=1, max_tm=tm+5, min_tm=tm-5, min_dg=dg, min_gc=gc-10, max_gc=gc+10)
    forward_primer_df['position'] = 'forward'
    reverse_primer_df['position'] = 'reverse'
    forward_primer_df = forward_primer_df[final_primers.columns]
    reverse_primer_df = reverse_primer_df[final_primers.columns]
    final_primers = pd.concat([final_primers, forward_primer_df, reverse_primer_df], ignore_index=True)

    if len(forward_primer_df)>0 and len(reverse_primer_df)>0:
        break
    else:
        pl+=1
gui.displayAndSave_gui(final_primers)
