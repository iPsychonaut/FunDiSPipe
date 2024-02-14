# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 12:23:15 2023

@author: ian.michael.bollinger@gmail.com/researchconsultants@critical.consulting

This is the main Graphical User Interface for a modified protocol by

Stephen Douglas Russell and paid for by the Fungal Diversity Survey (FunDiS).

Protocol Link: https://www.protocols.io/view/primary-data-analysis-basecalling-demultiplexing-a-dm6gpbm88lzp/v3?step=3

Current Default Settings:
    - Intended for Fungal ITS data

Current Allowed Input File:
    - A single ".fastq.gz" file that contains quality control checked Oxford Nanopore Guppy Basecalled sequences

Current Assumptions:
    - There is a text file in the same folder as the Input File that is named "primers.text" for minibar.py
    - There is a text file in the same folder as the Input File that is named "Index.text" for NGSpeciesID
"""

# Base Python Imports
import os, multiprocessing, math

# Required Python Imports
import tkinter as tk
from tkinter import filedialog, PhotoImage, messagebox

# Custom Python Imports
from MycoMap_Summarize import run_summary_prep
from FunDiS_NGSpeciesID import run_ngsid_prep
from FunDiS_Minbar import run_minibar_prep
from FunDiS_Tools import log_print, generate_log_file, initialize_logging_environment, run_subprocess_cmd, get_resource_values

# Global output_area variable
CPU_THREADS = 1
PERCENT_RESOURCES = 0.75
DEFAULT_LOG_FILE = None
ENVIRONMENT_TYPE = None
SAMPLE_SIZE = "500"
MIN_LENGTH_BP = "730"
MAX_STD_DEV_BP = "400"
HAP_PHASE_BOOL = True
LOGGING_ON = False

# GUI function to trigger minibar preparation
def gui_minibar_prep():
    """
    GUI function to trigger minibar preparation.

    This function is designed to be called from a graphical user interface (GUI) to start the minibar preparation. 
    It retrieves necessary paths from the GUI elements, logs the start of the process, and initiates the minibar 
    preparation sequence.

    Global Variables:
        output_area (str): A global variable used for logging output messages.

    Notes:
        - Requires certain GUI elements like file_path_entry to be defined and accessible.
    """
    minibar_fastq_gz_path= minibar_file_path_entry.get()

    # Populate this dictionary with entries from the main menu
    chopper_command_dict = {"-q": str(minibar_percent_match_entry.get()),
                            "--maxlength": str(minibar_barcode_edit_dist_entry.get()),
                            "--minlength": str(minibar_primer_edit_dist_entry.get())}

    # Populate this dictionary with entries from the main menu
    minibar_command_dict = {"-p": str(minibar_percent_match_entry.get()),
                            "-e": str(minibar_barcode_edit_dist_entry.get()),
                            "-E": str(minibar_primer_edit_dist_entry.get()),
                            "-l": str(minibar_search_length_entry.get())}

    if input_check(minibar_fastq_gz_path) == False:
        pass
    else:
        log_print( "MiniBar preparation started...\n")
        ngsid_output_dir = minibar_fastq_gz_path.replace(".fastq.gz","")
        minibar_path = "/mnt/d/FunDiS/minibar.py"
        minibar_index_path = f"{os.path.dirname(minibar_fastq_gz_path)}/Index.txt"
        
        run_minibar_prep(minibar_path, minibar_index_path, minibar_fastq_gz_path, ngsid_output_dir, chopper_command_dict, minibar_command_dict)

# GUI function to trigger NGSpeciesID preparation
def gui_ngsid_prep():
    """
    GUI function to trigger NGSpeciesID preparation.

    This function is designed to be invoked from a GUI to initiate the NGSpeciesID preparation process. It collects 
    necessary parameters from the GUI elements, logs the start of the process, and calls the function to run 
    NGSpeciesID prep.

    Global Variables:
        output_area (str): A global variable used for logging output messages.

    Notes:
        - Retrieves parameters such as sample size, minimum length, and maximum standard deviation from GUI elements.
        - Requires certain GUI elements like sample_size_entry, min_length_bp_entry, etc., to be defined and accessible.
    """
    ngsid_folder_path = ngsid_folder_path_entry.get()
    if input_check(ngsid_folder_path) == False:
        pass
    else:
        log_print( "NGSpeciesID preparation started...\n")
        
        sample_size = sample_size_entry.get()
        min_length_bp = min_length_bp_entry.get()
        max_std_dev_bp = max_std_dev_bp_entry.get()
        hap_phase_bool = hap_phase_var.get()
        
        input_fastq_file = minibar_file_path_entry.get().replace(".gz", "")
        ngsid_output_dir = ngsid_folder_path.replace(".fastq.gz", "")
        ngsid_primers_path = f"{os.path.dirname(ngsid_folder_path)}/primers.txt"
        
        run_ngsid_prep(ngsid_primers_path, input_fastq_file, sample_size, min_length_bp, max_std_dev_bp, hap_phase_bool, ngsid_output_dir, CPU_THREADS)

# GUI function to trigger NGSpeciesID preparation  
def gui_summary_prep():
    """
    GUI function to trigger summary preparation.

    This function is designed for use within a GUI to start the summary preparation process. It initializes and starts 
    a new thread dedicated to running the summary preparation.

    Global Variables:
        output_area (str): A global variable used for logging output messages.

    Notes:
        - Initiates the summary preparation in a separate thread to avoid blocking the GUI.
    """
    summary_folder_path = summary_folder_path_entry.get()
    if input_check(summary_folder_path) == False:
        pass
    else:
        log_print( "MycoMap Summarize preparation started...\n")
        
        hap_phase_bool = hap_phase_var.get()
        summary_output_dir = minibar_file_path_entry.get().replace(".fastq.gz", "")
        
        run_summary_prep(summary_output_dir, hap_phase_bool)

def choose_file_or_folder(file_path_entry):
    """
    Opens a dialog for selecting either a file or a folder and initializes the logging environment.

    This function first asks the user if they want to select a file or a folder. Based on the user's choice, it 
    triggers a GUI file dialog or a folder dialog to allow the user to make a selection. Once a selection is made, 
    it updates a GUI element with the file or folder path and calls `initialize_logging_environment` to set up the 
    logging based on the chosen file or folder.

    Notes:
        - Assumes the existence of specific GUI elements like file_path_entry for updating the file or folder path.
        - Directly interacts with the user interface to facilitate file or folder selection.
    """
    # Ask the user if they want to select a file or a folder
    choice = messagebox.askquestion("Choose File or Folder", "Do you want to select a file?")
    if choice == 'yes':
        selection = filedialog.askopenfilename(initialdir="/", title="Select file", filetypes=(("gzip files", "*.gz"), ("all files", "*.*")))
    else:
        selection = filedialog.askdirectory(initialdir="/", title="Select folder")
    
    # Update the file_path_entry with the selected file or folder path
    file_path_entry.delete(0, tk.END)
    file_path_entry.insert(0, selection)
    
    if LOGGING_ON == False:
        # Initialize the logging environment with the path to the input file
        initialize_logging_environment(selection)
        LOGGING_ON == True

def choose_minibar_file():
    """
    Opens a file dialog for selecting a .fastq.gz file for minibar and updates related fields.
    """
    filename = filedialog.askopenfilename(initialdir="/", title="Select .fastq.gz file", filetypes=(("gzip files", "*.gz"), ("all files", "*.*")))
    minibar_file_path_entry.delete(0, tk.END)
    minibar_file_path_entry.insert(0, filename)
    
    # Assuming filename is a path to a .fastq.gz file, prepare the base directory for other uses
    base_output_dir = filename.replace(".gz", "").replace(".fastq","")
    
    # Update ngsid_folder_path_entry and summary_folder_path_entry
    ngsid_folder_path_entry.delete(0, tk.END)
    ngsid_folder_path_entry.insert(0, base_output_dir)
    
    summary_folder_path_entry.delete(0, tk.END)
    summary_folder_path_entry.insert(0, base_output_dir)

    if LOGGING_ON == False:
        initialize_logging_environment(filename)
        LOGGING_ON == True

def about():
    """
    Displays an 'About' message box with information about the application.
    
    This function shows a message box containing information about the application, including the author's email, 
    references to the development protocol, and funding sources.
    
    Notes:
        - Utilizes the Tkinter messagebox for displaying the information.
        - Contains formatted text with multiple lines, including email addresses and a URL.
    """
    about_message = ("Code written by:\n"
                     "ian.michael.bollinger@gmail.com/\n"
                     "researchconsultants@critical.consulting\n"
                     "based on a modified protocol developed by:\n"
                     "Stephen Douglas Russell\n"
                     "paid for by the Fungal Diversity Survey (FunDiS)\n"
                     "Protocol link: https://www.protocols.io/view/primary-data-analysis-basecalling-demultiplexing-a-dm6gpbm88lzp/v3?step=3")
    messagebox.showinfo("About", about_message)

def exit_app():
    """
    Closes the main application window.
    
    This function is responsible for terminating the application by destroying the main application window.
    
    Notes:
        - Typically linked to an 'Exit' option in the application menu.
    """
    root.destroy()

def full_pipe():
    """
    

    Returns
    -------
    None.

    """
    gui_minibar_prep()
    gui_ngsid_prep()
    gui_summary_prep()

def input_check(file_path):
    if os.path.exists(file_path):
        return True
    else:
        log_print(f"NOTE:\tSelected file is invalid, please select a valid file: {file_path}")
        return False

# Test files
# "D:/FunDiS/ONT04/combined.fastq.gz"
# "/mnt/d/FunDiS/ONT04/combined.fastq.gz"

# Calculate the number of threads as specified percentage of available CPUs & RAM
CPU_THREADS, _ = get_resource_values(PERCENT_RESOURCES)

# Create root window
root = tk.Tk()
root.title("FunDiS Nanopore Barcoding Pipeline")

# Set the window icon
icon_image = PhotoImage(file="./fundis_icon.png")
root.iconphoto(False, icon_image)

# Initialize your Entry variables with default values in the main script body (not in a function)
sample_size_entry = tk.Entry(root, width=50)
sample_size_entry.insert(0, SAMPLE_SIZE)
min_length_bp_entry = tk.Entry(root, width=50)
min_length_bp_entry.insert(0, MIN_LENGTH_BP)
max_std_dev_bp_entry = tk.Entry(root, width=50)
max_std_dev_bp_entry.insert(0, MAX_STD_DEV_BP)
hap_phase_var = tk.BooleanVar(value=HAP_PHASE_BOOL)

# Set initial window size and allow resizing
root.geometry("500x900")  # Width x Height
root.resizable(True, True)  # Allow resizing in both directions

# Create a menu bar
menu_bar = tk.Menu(root)

# Create a File menu and add commands
file_menu = tk.Menu(menu_bar, tearoff=0)
file_menu.add_separator()
file_menu.add_command(label="Exit", command=exit_app)
menu_bar.add_cascade(label="File", menu=file_menu)

# Create a Help menu and add commands
help_menu = tk.Menu(menu_bar, tearoff=0)
help_menu.add_command(label="About", command=about)
menu_bar.add_cascade(label="Help", menu=help_menu)

# Display the menu
root.config(menu=menu_bar)

# Display Main Window Frame
main_frame = tk.Frame(root)
main_frame.pack(side=tk.LEFT, padx=10, pady=10, fill=tk.BOTH, expand=True)

# Load the logo image and display it
logo_image = PhotoImage(file="./fundis_logo.png")
logo_label = tk.Label(main_frame, image=logo_image)
logo_label.pack(side=tk.TOP, padx=10, pady=10)

# Set Percent Resources to use
percent_resources_label = tk.Label(main_frame, text="Percent CPU Resources (0.00-1.00):")
percent_resources_label.pack()
percent_resources_entry = tk.Entry(main_frame, width=50)
percent_resources_entry.insert(0, str(PERCENT_RESOURCES))
percent_resources_entry.pack()

# MiniBar File Path Entry
minibar_file_path_label = tk.Label(main_frame, text="MiniBar Input fastq.gz Path:")
minibar_file_path_label.pack()
minibar_file_path_entry = tk.Entry(main_frame, width=50)
minibar_file_path_entry.pack()

# File Browser Button for Minibar file input
file_browse_button = tk.Button(main_frame, text="Browse", command=choose_minibar_file)
file_browse_button.pack()

# Chopper Read Quality
chopper_read_quality_label = tk.Label(main_frame, text="Chopper Lowest Quality Read Allowed:")
chopper_read_quality_label.pack()
chopper_read_quality_entry = tk.Entry(main_frame, width=50)
chopper_read_quality_entry.insert(0, "10")  # Default value
chopper_read_quality_entry.pack()

# Chopper Max Length
chopper_max_length_label = tk.Label(main_frame, text="Chopper Maximum Seuqence Length (bp):")
chopper_max_length_label.pack()
chopper_max_length_entry = tk.Entry(main_frame, width=50)
chopper_max_length_entry.insert(0, "2500")  # Default value
chopper_max_length_entry.pack()

# Chopper Min Length
chopper_min_length_label = tk.Label(main_frame, text="Chopper Minimum Sequence Length (bp):")
chopper_min_length_label.pack()
chopper_min_length_entry = tk.Entry(main_frame, width=50)
chopper_min_length_entry.insert(0, "450")  # Default value
chopper_min_length_entry.pack()

# Minibar Percent Match
minibar_percent_match_label = tk.Label(main_frame, text="Minibar Percent Match (0.00-1.00):")
minibar_percent_match_label.pack()
minibar_percent_match_entry = tk.Entry(main_frame, width=50)
minibar_percent_match_entry.insert(0, "")  # Default value
minibar_percent_match_entry.pack()

# Minibar Barcode Edit Distance
minibar_barcode_edit_dist_label = tk.Label(main_frame, text="Minibar Barcode Edit Distance (bp):")
minibar_barcode_edit_dist_label.pack()
minibar_barcode_edit_dist_entry = tk.Entry(main_frame, width=50)
minibar_barcode_edit_dist_entry.insert(0, "")  # Default value
minibar_barcode_edit_dist_entry.pack()

# Minibar Primer Edit Distance
minibar_primer_edit_dist_label = tk.Label(main_frame, text="Minibar Primer Edit Distance (bp):")
minibar_primer_edit_dist_label.pack()
minibar_primer_edit_dist_entry = tk.Entry(main_frame, width=50)
minibar_primer_edit_dist_entry.insert(0, "")  # Default value
minibar_primer_edit_dist_entry.pack()

# Minibar Search Length
minibar_search_length_label = tk.Label(main_frame, text="Minibar Search Length (bp):")
minibar_search_length_label.pack()
minibar_search_length_entry = tk.Entry(main_frame, width=50)
minibar_search_length_entry.insert(0, "")  # Default value
minibar_search_length_entry.pack()

# MiniBar Button
minibar_prep_button = tk.Button(main_frame, text="Run MiniBar", command=gui_minibar_prep)
minibar_prep_button.pack(fill=tk.BOTH, expand=True)

# Advanced Settings Directly on Main Menu
ngsid_folder_path_label = tk.Label(main_frame, text="NGSID Input Folder Path:")
ngsid_folder_path_label.pack()
ngsid_folder_path_entry = tk.Entry(main_frame, width=50)
ngsid_folder_path_entry.pack()

# File Browser Button for NGSpeciesID folder input
file_browse_button = tk.Button(main_frame, text="Browse", command=lambda: choose_file_or_folder(ngsid_folder_path_entry))
file_browse_button.pack()

sample_size_label = tk.Label(main_frame, text="NGSID Sample Size (n):")
sample_size_label.pack()
sample_size_entry = tk.Entry(main_frame, width=50)
sample_size_entry.insert(0, SAMPLE_SIZE)
sample_size_entry.pack()

min_length_bp_label = tk.Label(main_frame, text="NGSID Min Length of Sequence (bp):")
min_length_bp_label.pack()
min_length_bp_entry = tk.Entry(main_frame, width=50)
min_length_bp_entry.insert(0, MIN_LENGTH_BP)
min_length_bp_entry.pack()

max_std_dev_bp_label = tk.Label(main_frame, text="NGSID Max Standard Deviation (bp):")
max_std_dev_bp_label.pack()
max_std_dev_bp_entry = tk.Entry(main_frame, width=50)
max_std_dev_bp_entry.insert(0, MAX_STD_DEV_BP)
max_std_dev_bp_entry.pack()

hap_phase_checkbox = tk.Checkbutton(main_frame, text="NGSID Haplotype Phasing", variable=hap_phase_var)
hap_phase_checkbox.pack()

# NGSpecies ID Button
ngsid_prep_button = tk.Button(main_frame, text="Run NGSpeciesID", command=gui_ngsid_prep)
ngsid_prep_button.pack(fill=tk.BOTH, expand=True)

# Summarize Input
summary_folder_path_label = tk.Label(main_frame, text="Summarize Input Folder Path:")
summary_folder_path_label.pack()
summary_folder_path_entry = tk.Entry(main_frame, width=50)
summary_folder_path_entry.pack()

# File Browser Button for Summarize folder input
file_browse_button = tk.Button(main_frame, text="Browse", command=lambda: choose_file_or_folder(summary_folder_path_entry))
file_browse_button.pack()

# Summary Button
summary_prep_button = tk.Button(main_frame, text="Run MycoMap Summary", command = gui_summary_prep)
summary_prep_button.pack(fill=tk.BOTH, expand=True)

# Full Pipeline Button
full_pipe_button = tk.Button(main_frame, text="Run Full Pipeline", command = full_pipe)
full_pipe_button.pack(fill=tk.BOTH, expand=True)

# Run main tkinter GUI loop
root.mainloop()
