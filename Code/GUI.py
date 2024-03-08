import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
import pandas as pd

def run_gui():
    # Create the main window
    root = tk.Tk()
    root.title("Enter DNA sequence")  # Set window title

    # Initialize return_value as an attribute of root to store user input
    root.return_value = None

    # Function to validate DNA sequence input
    def validate_input(event=None):
        input_text = text_box.get("1.0", "end-1c")  # Get input from text box
        if all(char in "gtacGTAC" for char in input_text):  # Check if input is a valid DNA sequence
            sub_button.pack(side='left', padx=(10, 0))  # Show submit button if input is valid
            label_invalid.pack_forget()  # Hide invalid input warning
        else:
            label_invalid.pack(side='left', padx=(10, 0))  # Show invalid input warning
            sub_button.pack_forget()  # Hide submit button

    # Function to browse and select a file
    def browse_file():
        # Open file dialog to select a FASTA file
        file_path = filedialog.askopenfilename(filetypes=(("FASTA files", "*.fasta"), ("FASTA files", "*.fa"), ("All files", "*.*")))
        if file_path:  # If a file is selected
            path_label.config(text=file_path)  # Display selected file path
            path_label.pack(side='left', padx=(10, 0))
            sub_button.pack(side='left', padx=(10, 0))  # Show submit button
            
    # Function to submit user input
    def submit():
        selected_option = dropdown.get()  # Get selected input type from dropdown
        if selected_option == "raw text":
            root.return_value = ['text', text_box.get("1.0", "end-1c")]  # Store text input
        elif selected_option == "FASTA file":
            root.return_value = ['file', path_label.cget("text")]  # Store file path
        print(root.return_value)  # For demonstration, print the return value
        root.destroy()  # Close the GUI window

    # Setup dropdown for selecting input type
    dropdown = ttk.Combobox(root, values=["-select input type-", "raw text", "FASTA file"], state="readonly")
    dropdown.current(0)  # Set default selection
    dropdown.pack(fill='x', padx=10, pady=10)
    dropdown.bind("<<ComboboxSelected>>", lambda e: update_layout())  # Bind selection event to update_layout function

    # Setup for text input frame
    text_frame = tk.Frame(root)
    label_tb = tk.Label(text_frame, text="Enter DNA sequence:")  # Label for text input
    label_tb.pack(side='top', fill='x', padx=10, pady=(10, 0))
    text_box = tk.Text(text_frame, height=5, width=50)  # Text box for input
    text_box.pack(side='top', fill='both', expand=True, padx=10)
    text_box.bind("<KeyRelease>", validate_input)  # Bind key release event to validate_input function
    label_invalid = tk.Label(text_frame, text="Invalid input", fg="red")  # Label for invalid input warning

    # Setup for file input frame
    file_frame = tk.Frame(root)
    label_browse = tk.Label(file_frame, text="Select FASTA file:")  # Label for file selection
    label_browse.pack(side='left', padx=10, pady=10)
    browse_button = ttk.Button(file_frame, text="Browse", command=browse_file)  # Button to open file dialog
    browse_button.pack(side='left', pady=10)
    path_label = tk.Label(file_frame, width=30)  # Label to display selected file path

    # Submit button setup
    sub_button = ttk.Button(root, text="Submit", command=submit)

    # Function to update layout based on dropdown selection
    def update_layout():
        if dropdown.get() == "raw text":
            file_frame.pack_forget()  # Hide file input frame
            text_frame.pack(fill='both', expand=True)  # Show text input frame
            sub_button.pack_forget()  # Hide submit button initially
            validate_input()  # Validate existing text input
        elif dropdown.get() == "FASTA file":
            text_frame.pack_forget()  # Hide text input frame
            file_frame.pack(fill='both', expand=True)  # Show file input frame
            label_invalid.pack_forget()  # Hide invalid input warning
            sub_button.pack(side='left', padx=(10, 0))  # Show submit button

    root.mainloop()  # Start the GUI event loop
    return root.return_value  # Return the user input after closing the GUI


def para_gui(length_default):
    # Validates numeric input, ensuring it's within a specified range or empty
    def on_validate(P):
        if P.strip() == "":  # Allow empty input to be cleared or edited
            return True
        try:
            value = int(P)
            return 1 <= value <= length_default  # Check if within allowable range
        except ValueError:  # Non-numeric input is invalid
            return False

    # Adjust the start entry based on its current value and the end entry's value
    def adjust_start(*args):
        try:
            start_value = int(start_entry.get())
            end_value = int(end_entry.get() or length_default)  # Use default if end entry is empty
            # Correct start value if out of range
            if start_value < 1:
                start_entry.delete(0, tk.END)
                start_entry.insert(0, "1")
            elif start_value > end_value:
                start_entry.delete(0, tk.END)
                start_entry.insert(0, str(end_value))
        except ValueError:  # Reset to default if invalid
            start_entry.delete(0, tk.END)
            start_entry.insert(0, "1")

    # Adjust the end entry based on its current value and the start entry's value
    def adjust_end(*args):
        try:
            start_value = int(start_entry.get() or 1)  # Use 1 if start entry is empty
            end_value = int(end_entry.get())
            # Correct end value if out of range
            if end_value > length_default:
                end_entry.delete(0, tk.END)
                end_entry.insert(0, str(length_default))
            elif end_value < start_value:
                end_entry.delete(0, tk.END)
                end_entry.insert(0, str(start_value))
        except ValueError:  # Reset to default if invalid
            end_entry.delete(0, tk.END)
            end_entry.insert(0, str(length_default))

    # Updates the state of input fields based on the selected amplicon placement option
    def update_input_fields():
        # Disable all fields initially
        length_entry.config(state='disabled')
        start_entry.config(state='disabled')
        end_entry.config(state='disabled')
        # Enable fields based on the selected option
        if radio_value.get() == 1:
            length_entry.config(state='normal')
        elif radio_value.get() == 2:
            start_entry.config(state='normal')
            end_entry.config(state='normal')

    # Enable or disable GC% spinbox based on its checkbox state
    def toggle_gc():
        gc_spinbox.config(state='normal' if activate_gc_var.get() else 'disabled')

    # Enable or disable Tm spinbox based on its checkbox state
    def toggle_tm():
        tm_spinbox.config(state='normal' if activate_tm_var.get() else 'disabled')

    # Enable or disable ΔG spinbox based on its checkbox state
    def toggle_dg():
        dg_spinbox.config(state='normal' if activate_dg_var.get() else 'disabled')
    
    # Enable or disable primer length spinbox based on its checkbox state
    def toggle_pl():
        pl_spinbox.config(state='normal' if activate_pl_var.get() else 'disabled')

    # Collect input values, store them in a dictionary, and close the GUI
    def submit_and_close():
        result['length'] = int(length_entry.get()) if radio_value.get() == 1 else 0
        result['start'] = int(start_entry.get()) if radio_value.get() == 2 else 0
        result['end'] = int(end_entry.get()) if radio_value.get() == 2 else 0
        result['amplicon_option'] = radio_value.get()
        result['pl'] = int(pl_spinbox.get()) if activate_pl_var.get() else 20
        result['gc'] = int(gc_spinbox.get()) if activate_gc_var.get() else 40
        result['tm'] = float(tm_spinbox.get()) if activate_tm_var.get() else 55
        result['dg'] = float(dg_spinbox.get()) if activate_dg_var.get() else -10
        root.destroy()

    root = tk.Tk()
    root.title("Nucleotide Options")

    result = {}  # Dictionary to store the results

    vcmd = (root.register(on_validate), '%P')

    font = ('Arial', 10)
    entry_width = 10
    spinbox_width = 5
    padding = {'padx': 5, 'pady': 5}

    radio_value = tk.IntVar(value=1)
    activate_gc_var = tk.BooleanVar(value=False)
    activate_tm_var = tk.BooleanVar(value=False)
    activate_dg_var = tk.BooleanVar(value=False)
    activate_pl_var = tk.BooleanVar(value=False)

    # UI setup for selecting amplicon placement and input fields
    amplicon_placement_label = tk.Label(root, text="Amplicon Placement", font=('Arial', 12, 'bold'))
    amplicon_placement_label.pack(fill='x', **padding)

    # Length input option
    length_frame = tk.Frame(root)
    length_frame.pack(fill='x', **padding)
    tk.Radiobutton(length_frame, text="Length", variable=radio_value, value=1, command=update_input_fields, font=font).pack(side=tk.LEFT, **padding)
    length_entry = tk.Entry(length_frame, width=entry_width, font=font, validate='key', validatecommand=vcmd)
    length_entry.pack(side=tk.LEFT, **padding)
    length_entry.insert(0, str(length_default))

    # Range input option
    range_frame = tk.Frame(root)
    range_frame.pack(fill='x', **padding)
    tk.Radiobutton(range_frame, text="Range", variable=radio_value, value=2, command=update_input_fields, font=font).pack(side=tk.LEFT, **padding)
    start_entry = tk.Entry(range_frame, width=entry_width, font=font)
    dash_label = tk.Label(range_frame, text="-", font=font)
    end_entry = tk.Entry(range_frame, width=entry_width, font=font)
    start_entry.pack(side=tk.LEFT, **padding)
    dash_label.pack(side=tk.LEFT)
    end_entry.pack(side=tk.LEFT, **padding)
    start_entry.insert(0, "1")
    end_entry.insert(0, str(length_default))

    start_entry.bind('<FocusOut>', adjust_start)
    end_entry.bind('<FocusOut>', adjust_end)

    # Option for applying parameters to all sequences
    all_frame = tk.Frame(root)
    all_frame.pack(fill='x', **padding)
    tk.Radiobutton(all_frame, text="All", variable=radio_value, value=3, command=update_input_fields, font=font).pack(side=tk.LEFT, **padding)

    # Other parameters setup
    other_parameters_label = tk.Label(root, text="Other Parameters", font=('Arial', 12, 'bold'))
    other_parameters_label.pack(fill='x', **padding)
    
    # Primer length option with checkbox
    pl_frame = tk.Frame(root)
    pl_frame.pack(fill='x', **padding)
    activate_pl_checkbtn = tk.Checkbutton(pl_frame, text="", variable=activate_pl_var, onvalue=True, offvalue=False, command=toggle_pl, font=font)
    pl_label = tk.Label(pl_frame, text="Primer\nLength:", font=font)
    pl_spinbox = tk.Spinbox(pl_frame, from_=15, to=50, textvariable=tk.DoubleVar(value=20), width=spinbox_width, font=font, state='disabled')
    activate_pl_checkbtn.pack(side=tk.LEFT, **padding)
    pl_label.pack(side=tk.LEFT, **padding)
    pl_spinbox.pack(side=tk.LEFT, **padding)

    # GC% option with checkbox
    gc_frame = tk.Frame(root)
    gc_frame.pack(fill='x', **padding)
    activate_gc_checkbtn = tk.Checkbutton(gc_frame, text="", variable=activate_gc_var, onvalue=True, offvalue=False, command=toggle_gc, font=font)
    gc_label = tk.Label(gc_frame, text="GC%:", font=font)
    gc_spinbox = tk.Spinbox(gc_frame, from_=0, to=100, wrap=True, textvariable=tk.IntVar(value=40), width=spinbox_width, font=font, state='disabled')
    activate_gc_checkbtn.pack(side=tk.LEFT, **padding)
    gc_label.pack(side=tk.LEFT, **padding)
    gc_spinbox.pack(side=tk.LEFT, **padding)

    # Tm option with checkbox
    tm_frame = tk.Frame(root)
    tm_frame.pack(fill='x', **padding)
    activate_tm_checkbtn = tk.Checkbutton(tm_frame, text="", variable=activate_tm_var, onvalue=True, offvalue=False, command=toggle_tm, font=font)
    tm_label = tk.Label(tm_frame, text="Tm:", font=font)
    tm_spinbox = tk.Spinbox(tm_frame, from_=0, to=100, increment=0.1, format="%.1f", wrap=True, textvariable=tk.DoubleVar(value=55), width=spinbox_width, font=font, state='disabled')
    activate_tm_checkbtn.pack(side=tk.LEFT, **padding)
    tm_label.pack(side=tk.LEFT, **padding)
    tm_spinbox.pack(side=tk.LEFT, **padding)

    # ΔG option with checkbox
    dg_frame = tk.Frame(root)
    dg_frame.pack(fill='x', **padding)
    activate_dg_checkbtn = tk.Checkbutton(dg_frame, text="", variable=activate_dg_var, onvalue=True, offvalue=False, command=toggle_dg, font=font)
    dg_label = tk.Label(dg_frame, text="Minimal\nΔG:", font=font)
    dg_spinbox = tk.Spinbox(dg_frame, from_=-100, to=0, increment=0.1, textvariable=tk.DoubleVar(value=-10), width=spinbox_width, font=font, state='disabled')
    activate_dg_checkbtn.pack(side=tk.LEFT, **padding)
    dg_label.pack(side=tk.LEFT, **padding)
    dg_spinbox.pack(side=tk.LEFT, **padding)

    # Submit button to collect and return the parameters
    submit_button = tk.Button(root, text="Submit", command=submit_and_close, font=font)
    submit_button.pack(**padding)

    # Initially disable start and end entries to match the default amplicon placement option
    start_entry.config(state='disabled')
    end_entry.config(state='disabled')

    root.mainloop()
    return result  # Return the dictionary of collected input values



def displayAndSave_gui(df: pd.DataFrame):
    # Function to prompt user for a file name and save the DataFrame to a .txt file
    def save_file():
        file_name = filedialog.asksaveasfilename(defaultextension=".txt",
                                                filetypes=[("Text documents", "*.txt"), ("All files", "*.*")])
        if file_name:  # Check if a file name was given
            df.to_csv(file_name, index=False)  # Save DataFrame to the specified file
            status_label.config(text="Saved!", fg="green")  # Update status label
            status_label.pack(pady=5)  # Display the status label

    # Function to display the DataFrame in a Treeview widget
    def display_dataframe(df):
        for i in tree.get_children():
            tree.delete(i)  # Clear existing rows in the tree
        tree["columns"] = list(df.columns)  # Set the columns in the Treeview
        tree.column("#0", width=0, stretch=tk.NO)  # Hide the first column
        # Configure column headings and widths
        for col in df.columns:
            col_width = max(df[col].astype(str).map(len).max(), len(col)) * 10  # Calculate column width
            tree.column(col, anchor=tk.CENTER, width=col_width)  # Set column width and alignment
            tree.heading(col, text=col, anchor=tk.CENTER)  # Set column headings
        # Insert data rows into the Treeview
        for row in df.itertuples(index=False, name=None):
            tree.insert("", tk.END, values=row)  # Insert each row of the DataFrame

    root = tk.Tk()
    root.title("Results")  # Title for the window

    # Configure style for the Treeview widget
    style = ttk.Style()
    style.configure("Treeview", rowheight=25, background="#D3D3D3",
                    fieldbackground="#D3D3D3", foreground="black")
    style.configure("Treeview.Heading", font=('Calibri', 13, 'bold'), background="lightblue")
    style.map('Treeview', background=[('selected', '#A9A9A9')])  # Change selection color

    tree = ttk.Treeview(root)
    display_dataframe(df)  # Display the DataFrame in the Treeview
    tree.pack(pady=20, fill='x', expand=True)  # Pack the Treeview widget

    save_button = tk.Button(root, text="Save", command=save_file)  # Button to save the DataFrame
    save_button.pack(pady=5)

    status_label = tk.Label(root, text="", fg="green")  # Label to show save status

    root.mainloop()  # Start the GUI event loop