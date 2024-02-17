import tkinter as tk
from tkinter import ttk
from tkinter import filedialog

def run_gui():
    root = tk.Tk()
    root.title("GUI with Drop-down Menu")

    # Initialize return_value as an attribute of root to ensure it is accessible outside nested functions
    root.return_value = None

    def validate_input(event):
        input_text = text_box.get("1.0", "end-1c")
        if all(char in "gtacGTAC" for char in input_text):
            sub_button.grid(row=2, column=1, sticky='W')
            label_invalid.grid_forget()
        else:
            label_invalid.grid(row=2, column=1, sticky="W")
            sub_button.grid_forget()

    def on_select(event):
        selected_option = dropdown.get()
        if selected_option == "raw text":
            label_tb.grid(row=1, column=0, sticky="N")
            text_box.grid(row=1, column=1, sticky="E", columnspan=3)
            browse_button.grid_forget()
            label_browse.grid_forget()
            path_label.grid_forget()
            sub_button.grid(row=2, column=1, sticky='W')
            root.geometry("500x200")
        elif selected_option == "FASTA or Genebank file":
            label_browse.grid(row=1, column=0)
            browse_button.grid(row=1, column=1)
            text_box.grid_forget()
            label_tb.grid_forget()
            sub_button.grid_forget()
            root.geometry("350x100")

    def browse_file():
        file_path = filedialog.askopenfilename()
        if file_path:  # Ensure a file was selected
            print("Selected file:", file_path)
            path_label.config(text=file_path)
            path_label.grid(row=2, column=1, sticky="W")
            sub_button.grid(row=3, column=1, sticky='W')

    def submit():
        selected_option = dropdown.get()
        if selected_option == "raw text":
            root.return_value = ['text', text_box.get("1.0", "end-1c")]
        elif selected_option == "FASTA or Genebank file":
            root.return_value = ['file', path_label.cget("text")]
        root.destroy()

    # GUI layout and configuration
    dropdown = ttk.Combobox(root, values=["-select input type-", "raw text", "FASTA or Genebank file"])
    dropdown.current(0)  # Set the default selected option
    dropdown.grid(row=0, column=1, sticky="w")
    dropdown.bind("<<ComboboxSelected>>", on_select)

    label_select = tk.Label(root, text="Select input:")
    label_select.grid(row=0, column=0, sticky="e")

    text_box = tk.Text(root, height=10, width=50)
    text_box.bind("<KeyRelease>", validate_input)

    label_tb = tk.Label(root, text="Enter DNA sequence:")
    label_browse = tk.Label(root, text="Select:")
    path_label = tk.Label(root, width=30)  # Adjust the width of the label as needed
    label_invalid = tk.Label(root, text="Invalid input", fg="red")

    browse_button = ttk.Button(root, text="Browse", command=browse_file)
    sub_button = ttk.Button(root, text="Submit", command=submit)

    # Start the GUI event loop
    root.mainloop()

    # Return the value set on root object after the GUI is closed
    return root.return_value


def para_gui(length_default):
    def on_validate(P):
        if P.strip() == "":
            return True
        try:
            value = int(P)
            return 1 <= value <= length_default
        except ValueError:
            return False

    def on_validate_range(P, start=True):
        if P.strip() == "":
            return True
        try:
            value = int(P)
            if start:
                end_value = int(end_entry.get() or 0)
                return 1 <= value <= length_default and value < end_value
            else:
                start_value = int(start_entry.get() or 0)
                return 1 <= value <= length_default and value > start_value
        except ValueError:
            return False

    def update_input_fields():
        length_entry.config(state='disabled')
        start_entry.config(state='disabled')
        end_entry.config(state='disabled')
        if radio_value.get() == 1:
            length_entry.config(state='normal')
        elif radio_value.get() == 2:
            start_entry.config(state='normal')
            end_entry.config(state='normal')
    
    def toggle_gc():
        gc_spinbox.config(state='normal' if activate_gc_var.get() else 'disabled')

    def toggle_tm():
        tm_spinbox.config(state='normal' if activate_tm_var.get() else 'disabled')

    def toggle_dg():
        dg_spinbox.config(state='normal' if activate_dg_var.get() else 'disabled')

    def submit_and_close():
        """Collect input values, store in a dict, and close the window."""
        result['length'] = length_entry.get() if radio_value.get() == 1 else None
        result['range_start'] = start_entry.get() if radio_value.get() == 2 else None
        result['range_end'] = end_entry.get() if radio_value.get() == 2 else None
        result['amplicon_option'] = radio_value.get()
        result['gc'] = gc_spinbox.get() if activate_gc_var.get() else None
        result['tm'] = tm_spinbox.get() if activate_tm_var.get() else None
        result['dg'] = dg_spinbox.get() if activate_dg_var.get() else None
        root.destroy()

    root = tk.Tk()
    root.title("Nucleotide Options")

    result = {}  # Dictionary to store the results

    vcmd = (root.register(on_validate), '%P')
    vcmd_range_start = (root.register(lambda P: on_validate_range(P, start=True)), '%P')
    vcmd_range_end = (root.register(lambda P: on_validate_range(P, start=False)), '%P')

    font = ('Arial', 10)
    entry_width = 10
    spinbox_width = 5
    padding = {'padx': 5, 'pady': 5}

    radio_value = tk.IntVar(value=1)
    activate_gc_var = tk.BooleanVar(value=False)
    activate_tm_var = tk.BooleanVar(value=False)
    activate_dg_var = tk.BooleanVar(value=False)

    # Amplicon Placement
    amplicon_placement_label = tk.Label(root, text="Amplicon Placement", font=('Arial', 12, 'bold'))
    amplicon_placement_label.pack(fill='x', **padding)

    # Length option
    length_frame = tk.Frame(root)
    length_frame.pack(fill='x', **padding)
    tk.Radiobutton(length_frame, text="Length", variable=radio_value, value=1, command=update_input_fields, font=font).pack(side=tk.LEFT, **padding)
    length_entry = tk.Entry(length_frame, width=entry_width, font=font, validate='key', validatecommand=vcmd)
    length_entry.pack(side=tk.LEFT, **padding)
    length_entry.insert(0, str(length_default))

    # Range option
    range_frame = tk.Frame(root)
    range_frame.pack(fill='x', **padding)
    tk.Radiobutton(range_frame, text="Range", variable=radio_value, value=2, command=update_input_fields, font=font).pack(side=tk.LEFT, **padding)
    start_entry = tk.Entry(range_frame, width=entry_width, font=font, validate='key', validatecommand=vcmd_range_start)
    dash_label = tk.Label(range_frame, text="-", font=font)
    end_entry = tk.Entry(range_frame, width=entry_width, font=font, validate='key', validatecommand=vcmd_range_end)
    start_entry.pack(side=tk.LEFT, **padding)
    dash_label.pack(side=tk.LEFT)
    end_entry.pack(side=tk.LEFT, **padding)
    start_entry.insert(0, "1")
    end_entry.insert(0, str(length_default))
    start_entry.insert(0, "1")

    # All option
    all_frame = tk.Frame(root)
    all_frame.pack(fill='x', **padding)
    tk.Radiobutton(all_frame, text="All", variable=radio_value, value=3, command=update_input_fields, font=font).pack(side=tk.LEFT, **padding)

    # Other Parameters
    other_parameters_label = tk.Label(root, text="Other Parameters", font=('Arial', 12, 'bold'))
    other_parameters_label.pack(fill='x', **padding)

    # GC%
    gc_frame = tk.Frame(root)
    gc_frame.pack(fill='x', **padding)
    activate_gc_checkbtn = tk.Checkbutton(gc_frame, text="", variable=activate_gc_var, onvalue=True, offvalue=False, command=toggle_gc, font=font)
    gc_label = tk.Label(gc_frame, text="GC%:", font=font)
    gc_spinbox = tk.Spinbox(gc_frame, from_=0, to=100, wrap=True, textvariable=tk.IntVar(value=40), width=spinbox_width, font=font, state='disabled')
    activate_gc_checkbtn.pack(side=tk.LEFT, **padding)
    gc_label.pack(side=tk.LEFT, **padding)
    gc_spinbox.pack(side=tk.LEFT, **padding)

    # Tm
    tm_frame = tk.Frame(root)
    tm_frame.pack(fill='x', **padding)
    activate_tm_checkbtn = tk.Checkbutton(tm_frame, text="", variable=activate_tm_var, onvalue=True, offvalue=False, command=toggle_tm, font=font)
    tm_label = tk.Label(tm_frame, text="Tm:", font=font)
    tm_spinbox = tk.Spinbox(tm_frame, from_=0, to=100, increment=0.1, format="%.1f", wrap=True, textvariable=tk.DoubleVar(value=55), width=spinbox_width, font=font, state='disabled')
    activate_tm_checkbtn.pack(side=tk.LEFT, **padding)
    tm_label.pack(side=tk.LEFT, **padding)
    tm_spinbox.pack(side=tk.LEFT, **padding)

    # ΔG
    dg_frame = tk.Frame(root)
    dg_frame.pack(fill='x', **padding)
    activate_dg_checkbtn = tk.Checkbutton(dg_frame, text="", variable=activate_dg_var, onvalue=True, offvalue=False, command=toggle_dg, font=font)
    dg_label = tk.Label(dg_frame, text="MaximalΔG:", font=font)
    dg_spinbox = tk.Spinbox(dg_frame, from_=-100, to=0, increment=0.1, textvariable=tk.DoubleVar(value=-10), width=spinbox_width, font=font, state='disabled')
    activate_dg_checkbtn.pack(side=tk.LEFT, **padding)
    dg_label.pack(side=tk.LEFT, **padding)
    dg_spinbox.pack(side=tk.LEFT, **padding)

    # Submit button
    submit_button = tk.Button(root, text="Submit", command=submit_and_close, font=font)
    submit_button.pack(**padding)

    # Initially disable start and end entries
    start_entry.config(state='disabled')
    end_entry.config(state='disabled')

    root.mainloop()
    return result  # Return the collected input values


