# Import necessary libraries and modules
import GUI as gui  # GUI module for graphical user interface interactions
import PrimerDesigner  # Module for designing primers
import pandas as pd  # Pandas library for data manipulation

# Run the GUI to get the input method and sequence or filepath
res = gui.run_gui()
# If the input method was a file, read the sequence from the file, excluding header lines
if res[0] == "file":
    sequence = ""
    with open(res[1]) as file:
        for line in file:
            if line[0] != ">":  # Skip header lines starting with ">"
                sequence += line.strip()  # Concatenate sequence lines, removing whitespace
# If the input was provided as text, directly use the provided sequence
elif res[0] == "text":
    sequence = res[1]

# Calculate the length of the sequence
seq_length = len(sequence)
# Run another part of the GUI to get parameters for primer design based on sequence length
parameters = gui.para_gui(seq_length)

# Extract primer design parameters from the user input
length = parameters["length"]
start = parameters["start"]
end = parameters["end"]
gc = parameters["gc"]
tm = parameters["tm"]
dg = parameters["dg"]
pl = parameters["pl"]
amplicon_option = parameters["amplicon_option"]

# Define the amplicon (target sequence segment for primer design) based on user selection
if amplicon_option == 1:
    start = 0
    end = length - 1
    amplicon = sequence[start:end]
elif amplicon_option == 2:
    amplicon = sequence[start - 1:end - 1]
elif amplicon_option == 3:
    start = 0
    end = seq_length - 1
    amplicon = sequence[start:end]

# Prepare a DataFrame to store final primer designs
final_primers = pd.DataFrame(columns=["sequence", "position", "length", "gc_content", "tm", "secondary_structure", "dg"])

# Loop to design primers, adjusting primer length if necessary
while True:
    # Design forward and reverse primers based on the sequence and specified parameters
    forward_primer, reverse_primer = PrimerDesigner.primer_design(sequence, start, end, p_length=pl)
    forward_primer = PrimerDesigner.ListToPrimes(forward_primer)
    reverse_primer = PrimerDesigner.ListToPrimes(reverse_primer)
    
    # Prepare DataFrames to store primer properties
    forward_primer_df = pd.DataFrame(columns=["sequence", "length", "gc_content", "tm", "secondary_structure", "dg"])
    reverse_primer_df = pd.DataFrame(columns=["sequence", "length", "gc_content", "tm", "secondary_structure", "dg"])

    # Populate the DataFrames with primer properties
    for primer in forward_primer:
        forward_primer_df.loc[len(forward_primer_df)] = primer.to_dict()
    for primer in reverse_primer:
        reverse_primer_df.loc[len(reverse_primer_df)] = primer.to_dict()

    # Check binding sites and filter primers based on specified criteria
    forward_primer_df = PrimerDesigner.checkBindingSites(forward_primer_df, sequence)
    reverse_primer_df = PrimerDesigner.checkBindingSites(reverse_primer_df, sequence)
    forward_primer_df = PrimerDesigner.filterPrimers(forward_primer_df, max_binding_sites=1, max_tm=tm+5, min_tm=tm-5, min_dg=dg, min_gc=gc-10, max_gc=gc+10)
    reverse_primer_df = PrimerDesigner.filterPrimers(reverse_primer_df, max_binding_sites=1, max_tm=tm+5, min_tm=tm-5, min_dg=dg, min_gc=gc-10, max_gc=gc+10)
    
    # Assign primer positions and align DataFrame columns with final_primers
    forward_primer_df['position'] = 'forward'
    reverse_primer_df['position'] = 'reverse'
    forward_primer_df = forward_primer_df[final_primers.columns]
    reverse_primer_df = reverse_primer_df[final_primers.columns]
    # Concatenate the new primer data into the final_primers DataFrame
    final_primers = pd.concat([final_primers, forward_primer_df, reverse_primer_df], ignore_index=True)

    # Break the loop if both forward and reverse primers are found or if primer length exceeds 30
    if len(forward_primer_df) > 0 and len(reverse_primer_df) > 0:
        break
    else:
        if pl >= 30:
            break
        pl += 1  # Increment primer length and retry

# Display the final primers in the GUI and offer to save them
gui.displayAndSave_gui(final_primers)