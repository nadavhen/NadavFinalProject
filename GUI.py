import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
import pandas as pd

def run_gui():
    root = tk.Tk()
    root.title("Enter DNA sequence")

    # Initialize return_value as an attribute of root to ensure it is accessible outside nested functions
    root.return_value = None

    def validate_input(event=None):
        input_text = text_box.get("1.0", "end-1c")
        if all(char in "gtacGTAC" for char in input_text):
            sub_button.pack(side='left', padx=(10, 0))
            label_invalid.pack_forget()
        else:
            label_invalid.pack(side='left', padx=(10, 0))
            sub_button.pack_forget()

    def browse_file():
        file_path = filedialog.askopenfilename(filetypes=(("FASTA files", "*.fasta"), ("FASTA files", "*.fa"), ("All files", "*.*")))
        if file_path:  # Ensure a file was selected
            path_label.config(text=file_path)
            path_label.pack(side='left', padx=(10, 0))
            sub_button.pack(side='left', padx=(10, 0))
            
    def submit():
        selected_option = dropdown.get()
        if selected_option == "raw text":
            root.return_value = ['text', text_box.get("1.0", "end-1c")]
        elif selected_option == "FASTA file":
            root.return_value = ['file', path_label.cget("text")]
        print(root.return_value)  # For demonstration, print the return value
        root.destroy()

    # Dropdown for selecting input type
    dropdown = ttk.Combobox(root, values=["-select input type-", "raw text", "FASTA file"], state="readonly")
    dropdown.current(0)  # Set default selection
    dropdown.pack(fill='x', padx=10, pady=10)
    dropdown.bind("<<ComboboxSelected>>", lambda e: update_layout())

    # Text input frame setup
    text_frame = tk.Frame(root)
    label_tb = tk.Label(text_frame, text="Enter DNA sequence:")
    label_tb.pack(side='top', fill='x', padx=10, pady=(10, 0))
    text_box = tk.Text(text_frame, height=5, width=50)
    text_box.pack(side='top', fill='both', expand=True, padx=10)
    text_box.bind("<KeyRelease>", validate_input)
    label_invalid = tk.Label(text_frame, text="Invalid input", fg="red")

    # File input frame setup
    file_frame = tk.Frame(root)
    label_browse = tk.Label(file_frame, text="Select FASTA file:")
    label_browse.pack(side='left', padx=10, pady=10)
    browse_button = ttk.Button(file_frame, text="Browse", command=browse_file)
    browse_button.pack(side='left', pady=10)
    path_label = tk.Label(file_frame, width=30)

    # Submit button
    sub_button = ttk.Button(root, text="Submit", command=submit)

    # Function to update layout based on dropdown selection
    def update_layout():
        if dropdown.get() == "raw text":
            file_frame.pack_forget()
            text_frame.pack(fill='both', expand=True)
            sub_button.pack_forget()  # Reset submit button state
            validate_input()  # Validate input to decide if submit button should be shown
        elif dropdown.get() == "FASTA file":
            text_frame.pack_forget()
            file_frame.pack(fill='both', expand=True)
            label_invalid.pack_forget()
            sub_button.pack(side='left', padx=(10, 0))

    root.mainloop()
    return root.return_value  # Return the collected input values



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
    
    def toggle_pl():
        pl_spinbox.config(state='normal' if activate_pl_var.get() else 'disabled')

    def submit_and_close():
        """Collect input values, store in a dict, and close the window."""
        result['length'] = int(length_entry.get()) if radio_value.get() == 1 else 0
        result['start'] = int(start_entry.get()) if radio_value.get() == 2 else 0
        result['end'] = int(end_entry.get()) if radio_value.get() == 2 else 0
        result['amplicon_option'] = radio_value.get()
        result['pl'] = int(pl_spinbox.get()) if activate_pl_var.get() else 20
        result['gc'] = int(gc_spinbox.get()) if activate_gc_var.get() else 40
        result['tm'] = int(tm_spinbox.get()) if activate_tm_var.get() else 55
        result['dg'] = int(dg_spinbox.get()) if activate_dg_var.get() else -10
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
    activate_pl_var = tk.BooleanVar(value=False)

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
    
    # Primer Length
    pl_frame = tk.Frame(root)
    pl_frame.pack(fill='x', **padding)
    activate_pl_checkbtn = tk.Checkbutton(pl_frame, text="", variable=activate_pl_var, onvalue=True, offvalue=False, command=toggle_pl, font=font)
    pl_label = tk.Label(pl_frame, text="Primer\nLenght:", font=font)
    pl_spinbox = tk.Spinbox(pl_frame, from_=15, to=50, textvariable=tk.DoubleVar(value=20), width=spinbox_width, font=font, state='disabled')
    activate_pl_checkbtn.pack(side=tk.LEFT, **padding)
    pl_label.pack(side=tk.LEFT, **padding)
    pl_spinbox.pack(side=tk.LEFT, **padding)

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
    dg_label = tk.Label(dg_frame, text="Minimal\nΔG:", font=font)
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


def displayAndSave_gui(df: pd.DataFrame):
    def save_file():
        file_name = filedialog.asksaveasfilename(defaultextension=".txt",
                                                filetypes=[("Text documents", "*.txt"), ("All files", "*.*")])
        if file_name:
            # Assuming saving the DataFrame content
            df.to_csv(file_name, index=False)
            status_label.config(text="Saved!", fg="green")
            status_label.pack(pady=5)

    def display_dataframe(df):
        for i in tree.get_children():
            tree.delete(i)
        tree["columns"] = list(df.columns)
        tree.column("#0", width=0, stretch=tk.NO)
        for col in df.columns:
            col_width = max(df[col].astype(str).map(len).max(), len(col)) * 10
            tree.column(col, anchor=tk.CENTER, width=col_width)  # Center column values
            tree.heading(col, text=col, anchor=tk.CENTER)  # Center column headings
        for row in df.itertuples(index=False, name=None):
            tree.insert("", tk.END, values=row)  # Insert row values


    # Create the main window
    root = tk.Tk()
    root.title("Results")

    # Styling for the Treeview
    style = ttk.Style()
    style.configure("Treeview", rowheight=25, background="#D3D3D3",
                    fieldbackground="#D3D3D3", foreground="black")
    style.configure("Treeview.Heading", font=('Calibri', 13, 'bold'), background="lightblue")
    style.map('Treeview', background=[('selected', '#A9A9A9')])  # Change selection color to dark gray

    # Adjust alignment and background for cells
    style.configure("Treeview", rowheight=30)  # Adjust row height for better visibility
    style.layout("Treeview", [('Treeview.treearea', {'sticky': 'nswe'})])  # Layout adjustment

    # Create and pack the Treeview widget
    tree = ttk.Treeview(root)
    display_dataframe(df)
    tree.pack(pady=20, fill='x', expand=True)

    # Create a 'Save' button
    save_button = tk.Button(root, text="Save", command=save_file)
    save_button.pack(pady=5)

    # Create a status label (initially hidden or empty)
    status_label = tk.Label(root, text="", fg="green")

    # Start the GUI event loop
    root.mainloop()
