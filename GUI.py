import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
def validate_input(event):
    input_text = text_box.get("1.0", "end-1c")
    if all(char in "gtacGTAC" for char in input_text):
        sub_button.grid(row=2, column=1,sticky='W')
        label_invalid.grid_forget()
    else:
        label_invalid.grid(row=2, column=1, sticky="W")
        sub_button.grid_forget()


def on_select(event):
    selected_option = dropdown.get()
    if selected_option == "raw text":
        label_tb.grid(row=1, column=0, sticky="N")
        root.grid_propagate(False)
        text_box.grid(row=1, column=1, sticky="E", columnspan=3)
        browse_button.grid_forget()
        label_browse.grid_forget()
        path_label.grid_forget()
        sub_button.grid(row=2, column=1,sticky='W')
        root.geometry("500x200")
        
    elif selected_option == "FASTA or Genebank file":
        label_browse.grid(row=1, column=0)
        browse_button.grid(row=2, column=0)
        text_box.grid_forget()
        browse_button.grid()
        sub_button.grid_forget()
        root.geometry("350x100")

def browse_file():
    file_path = filedialog.askopenfilename()
    print("Selected file:", file_path)
    path_label.config(text=file_path)
    path_label.grid(row=1, column=1 ,sticky="W")
    sub_button.grid(row=2, column=1,sticky='W')

root = tk.Tk()
root.title("GUI with Drop-down Menu")
# Create features
dropdown = ttk.Combobox(root, values=["-select input type-", "raw text", "FASTA or Genebank file"])
browse_button = ttk.Button(root, text="Browse", command=browse_file)
text_box = tk.Text(root, height=10, width=50)
label_tb = tk.Label(root, text="Enter DNA sequence:")
label_browse = tk.Label(root, text="Select:")
path_label = tk.Label(root, width=30)  # Set the width of the label
sub_button = ttk.Button(root, text="Submit")
text_box.bind("<KeyRelease>", validate_input)
label_invalid = tk.Label(root, text="Invalid input", fg="red")


label_select = tk.Label(root, text="Select input:")
label_select.grid(row=0, column=0, sticky="e")

dropdown.current(0)  # Set the default selected option
dropdown.bind("<<ComboboxSelected>>", on_select)
dropdown.grid(row=0, column=1, sticky="w")

root.mainloop()