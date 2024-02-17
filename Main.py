import GUI as gui
from Bio.Seq import Seq


res=gui.run_gui()
if res[0]=="file":
    with open(res[1]) as file:
        sequence=file.read()
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
    




