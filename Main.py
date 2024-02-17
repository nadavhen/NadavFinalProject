import GUI as gui

res=gui.run_gui()
if res[0]=="file":
    with open(res[1]) as file:
        sequence=file.read()
elif res[0]=="text":
    sequence=res[1]

seq_length=len(sequence)
parameters=gui.para_gui(seq_length)

#