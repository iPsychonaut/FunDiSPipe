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
import pathlib, os, platform, multiprocessing, math

# Required Python Imports
from termcolor import cprint
from datetime import datetime
import tkinter as tk
from tkinter import filedialog, PhotoImage, messagebox

# Custom Python Imports
from MycoMap_Summarize import run_summary_prep
from FunDiS_NGSpeciesID import run_ngsid_prep
from FunDiS_Minbar import run_minibar_prep

# Global output_area variable
PERCENT_RESOURCES = 0.2
DEFAULT_LOG_FILE = None

# Function to Generate a Log File
def generate_log_file(log_file_path, use_numerical_suffix=False):
    """
    Generates or clears a log file based on the given parameters.
    
    This function either creates a new log file or clears an existing one, depending on the specified parameters. 
    If the 'use_numerical_suffix' parameter is True, and the file already exists, a new file with a numerical suffix 
    will be created. Otherwise, the existing file will be cleared.
    
    Parameters:
        log_file_path (str): Path to the log file.
        use_numerical_suffix (bool): If True, creates new files with numerical suffixes if the file exists; otherwise, clears the existing file.
    
    Returns:
        str: Path to the log file.
    
    Notes:
    - Useful for managing log file versions without overwriting existing logs.
    - The numerical suffix increments for each new file created in the same location.
    - When 'use_numerical_suffix' is False, it refreshes the log file by clearing existing content.
    """
    if os.path.exists(log_file_path) and use_numerical_suffix:
        # If using numerical suffixes, increment until a new filename is found
        counter = 1
        new_log_file_path = f"{log_file_path.rsplit('.', 1)[0]}_{counter}.txt"
        while os.path.exists(new_log_file_path):
            counter += 1
            new_log_file_path = f"{log_file_path.rsplit('.', 1)[0]}_{counter}.txt"
        log_file_path = new_log_file_path
    else:
        # Clear the existing log file or create a new one
        open(log_file_path, 'w').close()
    
    return log_file_path

# Function to color coded print to console and save to log_file information
def log_print(input_message, log_file=None):
    """
    Logs a message to a file and prints it to the console with appropriate coloring.
    
    This function takes a message and logs it to the specified file. Additionally, the message is printed to the 
    console, potentially with specific coloring depending on the context.
    
    Parameters:
        input_message (str): Message to be logged and printed.
        log_file (str): Path to the log file.

    Notes:
        - The function uses a global default log file if none is specified.
        - Timestamps each log entry for easy tracking.
        - Utilizes color coding in the console to distinguish between different types of messages (e.g., errors, warnings).
        - Supports color coding for specific message types: NOTE, CMD, ERROR, WARN, and PASS.
        - Falls back to default (white) color if the message type is unrecognized.
    """
    # Access the global variable
    global DEFAULT_LOG_FILE
    
    # Use the default log file if none specified
    if log_file is None:
        log_file = DEFAULT_LOG_FILE  
    
    # Establish current date-time
    now = datetime.now()
    message = f'[{now:%Y-%m-%d %H:%M:%S}]\t{input_message}'

    # Determine the print color based on the input_message content
    message_type_dict = {'NOTE': ['blue'],
                         'CMD': ['cyan'],
                         'ERROR': ['red'],
                         'WARN': ['yellow'],
                         'PASS': ['green'],}
    print_color = ['white']  # Default color
    for key, value in message_type_dict.items():
        if key.lower() in input_message.lower():
            print_color = value
            break

    try:
        # Writing the message to the log file
        with open(log_file, 'a') as file:
            print(message, file=file)
    except TypeError:
        print(f"UNLOGGED ERROR:\tUnable to load the log file provided: {log_file}")

    # Handling different message types for colored printing
    try:
        cprint(message, print_color[0])
    except (KeyError, IndexError):
        cprint(message, print_color[1] if len(print_color) > 1 else 'white')

# Function to initialize the logging environment
def initialize_logging_environment(input_file_path):
    """
    Initializes the logging environment based on the given input file path.

    This function sets up the logging environment by adjusting file paths according to the operating system in use, 
    ensuring file existence, and then generating a log file. It sets the global DEFAULT_LOG_FILE variable to the path 
    of the generated log file.

    Parameters:
        input_file_path (str): Path to the input file which influences the log file generation.

    Global Variables:
        DEFAULT_LOG_FILE (str): The default path for the log file used throughout the logging process.

    Notes:
        - Supports Windows, Linux/WSL, and Darwin (macOS) environments.
        - Prints unlogged messages to the console regarding environment detection and file existence.
        - Modifies the global DEFAULT_LOG_FILE variable.
    """
    global DEFAULT_LOG_FILE
    
    # Determine the operating system
    os_name = platform.system()
    
    # Depending on the operating system, handle the input file path differently
    if os_name == "Windows":
        # On Windows, ensure the file exists
        if not os.path.exists(input_file_path):
            print(f'UNLOGGED ERROR:\tThe specified file does not exist: {input_file_path}')
            return
        print('UNLOGGED:\tWINDOWS ENVIRONMENT')
    elif os_name in ["Linux", "Darwin"]:  # Darwin is the system name for macOS
        # For Linux/WSL/Mac, convert the Windows drive letter to the appropriate mount path
        drive, path_without_drive = os.path.splitdrive(input_file_path)
        if drive:
            # Replace the drive letter with the WSL-style path (e.g., '/mnt/c')
            drive_letter = drive.strip(":\\/")
            path_without_drive_mod = path_without_drive.replace("\\", "/")
            input_file_path = f'/mnt/{drive_letter.lower()}{path_without_drive_mod}'
        else:
            # If there's no drive letter, it's already a POSIX-style path
            pass
        print('UNLOGGED:\tLINUX/WSL/MAC ENVIRONMENT')
    else:
        print(f'UNLOGGED ERROR:\tUnsupported OS: {os_name}')
        return

    # Check if the file exists on the adjusted path
    if not os.path.exists(input_file_path):
        print(f'UNLOGGED ERROR:\tThe provided log file does not exist: {input_file_path}')
        return
    
    # Generate the log file based on the input file
    file_extension = pathlib.Path(input_file_path).suffix
    run_log = input_file_path.replace(f".fastq{file_extension}", '.tsv')
    run_log = generate_log_file(run_log, use_numerical_suffix=False)
    
    # Set the default log file to the generated run_log
    DEFAULT_LOG_FILE = run_log

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
    ngsid_fastq_gz_path = file_path_entry.get()
    if input_check(ngsid_fastq_gz_path) == False:
        pass
    else:
        log_print( "MiniBar preparation started...\n")
    
        ngsid_output_dir = ngsid_fastq_gz_path.replace(".fastq.gz","")
        minibar_path = "/mnt/d/FunDiS/minibar.py"
        minibar_index_path = f"{os.path.dirname(ngsid_fastq_gz_path)}/Index.txt"
        
        run_minibar_prep(minibar_path, minibar_index_path, ngsid_fastq_gz_path, ngsid_output_dir)

# GUI function to trigger NGSpeciesID preparation
def gui_ngsid_prep():
    """
    GUI function to trigger NGSpeciesID preparation.

    This function is designed to be invoked from a GUI to initiate the NGSpeciesID preparation process. It collects 
    necessary parameters from the GUI elements, logs the start of the process, and calls the function to run 
    NGSpeciesID prep.

    Global Variables:
        output_area (str): A global variable used for logging output messages.
        cpu_threads (int): Number of CPU threads to be used in processing.

    Notes:
        - Retrieves parameters such as sample size, minimum length, and maximum standard deviation from GUI elements.
        - Requires certain GUI elements like sample_size_entry, min_length_bp_entry, etc., to be defined and accessible.
    """
    global cpu_threads
    ngsid_fastq_gz_path = file_path_entry.get()
    if input_check(ngsid_fastq_gz_path) == False:
        pass
    else:
        log_print( "NGSpeciesID preparation started...\n")
    
        sample_size = sample_size_entry.get()
        min_length_bp = min_length_bp_entry.get()
        max_std_dev_bp = max_std_dev_bp_entry.get()
        hap_phase_bool = hap_phase_var.get()
        input_fastq_file = ngsid_fastq_gz_path.replace(".gz", "")
        ngsid_output_dir = ngsid_fastq_gz_path.replace(".fastq.gz", "")
        ngsid_primers_path = f"{os.path.dirname(ngsid_fastq_gz_path)}/primers.txt"
        
        run_ngsid_prep(ngsid_primers_path, input_fastq_file, sample_size, min_length_bp, max_std_dev_bp, hap_phase_bool, ngsid_output_dir)

# GUI function to trigger NGSpeciesID preparation  
def gui_summary_prep(ngsid_fastq_gz_path):
    """
    GUI function to trigger summary preparation.

    This function is designed for use within a GUI to start the summary preparation process. It initializes and starts 
    a new thread dedicated to running the summary preparation.

    Global Variables:
        output_area (str): A global variable used for logging output messages.

    Notes:
        - Initiates the summary preparation in a separate thread to avoid blocking the GUI.
    """
    if input_check(ngsid_fastq_gz_path) == False:
        pass
    else:
        log_print( "MycoMap Summarize preparation started...\n")
        
        hap_phase_bool = hap_phase_var.get()
        ngsid_output_dir = ngsid_fastq_gz_path.replace(".fastq.gz", "")
    
        run_summary_prep(ngsid_output_dir, hap_phase_bool)

def choose_file():
    """
    Opens a file dialog for selecting a file and initializes the logging environment.

    This function triggers a GUI file dialog to allow the user to select a file. Once a file is selected, it updates 
    a GUI element with the file path and calls `initialize_logging_environment` to set up the logging based on the 
    chosen file.

    Notes:
        - Assumes the existence of specific GUI elements like file_path_entry for updating the file path.
        - Directly interacts with the user interface to facilitate file selection.
    """
    filename = filedialog.askopenfilename(initialdir="/", title="Select .fastq.gz file", filetypes=(("gzip files", "*.gz"),("all files", "*.*")))
    file_path_entry.delete(0, tk.END)
    file_path_entry.insert(0, filename)
    # Initialize the logging environment with the path to the input file
    initialize_logging_environment(filename)

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

def open_settings_window():
    """
    Opens a new window for adjusting advanced settings of the application.
    
    This function creates a new top-level window containing various input fields and settings related to the application. 
    It includes settings for CPU resource allocation, NGSpeciesID parameters, and options for haplotype phasing.
    
    Notes:
        - The settings window includes fields for percent CPU resources, sample size, minimum sequence length, 
          maximum standard deviation, and an option for haplotype phasing.
        - Each setting field is initialized with a default value.
        - The function declares several global variables that are used to store the settings.
        - The settings window has a fixed size and does not allow resizing.
    """
    settings_window = tk.Toplevel(root)
    settings_window.title("Advanced Settings")

    # Set the window icon
    icon_image = PhotoImage(file="./fundis_icon.png")
    settings_window.iconphoto(False, icon_image)

    # Set initial window size and allow resizing
    settings_window.geometry("250x250")  # Width x Height
    settings_window.resizable(False, False)

    # General FunDiS Pipeline Settings Label
    tk.Label(settings_window, text="FunDiS Pipeline Settings").pack()
    
    # PERCENT_RESOURCES Entry
    tk.Label(settings_window, text="Percent CPU Resources (0.00-1.00):").pack()
    global PERCENT_RESOURCES  # Declare as global if you need to access outside the function
    PERCENT_RESOURCES = tk.Entry(settings_window, width=50)
    PERCENT_RESOURCES.pack()
    PERCENT_RESOURCES.insert(0, "0.2")  # Default value for PERCENT_RESOURCES

    # TODO: Minibar Settings Label
    tk.Label(settings_window, text="Minibar Advanced Settings").pack()

    # NGSpeciesID Settings Label
    tk.Label(settings_window, text="NGSpeciesID Advanced Settings").pack()

    # Sample Size Entry
    tk.Label(settings_window, text="Sample Size (n):").pack()
    global sample_size_entry  # Declare as global if you need to access outside the function
    sample_size_entry = tk.Entry(settings_window, width=50)
    sample_size_entry.pack()
    sample_size_entry.insert(0, "500")  # Default value for Sample Size

    # Minimum Length of Sequence Entry
    tk.Label(settings_window, text="Minimum Length of Sequence (bp):").pack()
    global min_length_bp_entry
    min_length_bp_entry = tk.Entry(settings_window, width=50)
    min_length_bp_entry.pack()
    min_length_bp_entry.insert(0, "730")  # Default value for Minimum Length
    
    # Maximum Standard Deviation Entry
    tk.Label(settings_window, text="Maximum Standard Deviation (bp):").pack()
    global max_std_dev_bp_entry
    max_std_dev_bp_entry = tk.Entry(settings_window, width=50)
    max_std_dev_bp_entry.pack()
    max_std_dev_bp_entry.insert(0, "400")  # Default value for Maximum Standard Deviation
    
    # Checkbox for Haplotype Phasing
    global hap_phase_var
    hap_phase_var = tk.IntVar(value=1)
    hap_phase_checkbox = tk.Checkbutton(settings_window, text="Haplotype Phasing", variable=hap_phase_var)
    hap_phase_checkbox.pack()

def input_check(file_path):
    if os.path.exists(file_path):
        return True
    else:
        log_print(f"NOTE:\tSelected file is invalid, please select a valid file: {file_path}")
        return False

# Get the number of CPUs available on the system
num_cpus = multiprocessing.cpu_count()

# Calculate the number of threads as 80% of available CPUs & RAM
cpu_threads = int(math.floor(num_cpus * PERCENT_RESOURCES))

# Create root window
root = tk.Tk()
root.title("FunDiS Nanopore Barcoding Pipeline")

# Set the window icon
icon_image = PhotoImage(file="./fundis_icon.png")
root.iconphoto(False, icon_image)

# Set initial window size and allow resizing
root.geometry("250x500")  # Width x Height
root.resizable(True, True)  # Allow resizing in both directions

# Create a menu bar
menu_bar = tk.Menu(root)

# Create a File menu and add commands
file_menu = tk.Menu(menu_bar, tearoff=0)
file_menu.add_command(label="Browse", command=choose_file)
file_menu.add_command(label="Settings", command=open_settings_window)  # Add this line
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

# File Path Entry
file_path_label = tk.Label(main_frame, text="NGSID gz Fastq Path:")
file_path_label.pack()
file_path_entry = tk.Entry(main_frame, width=50)
file_path_entry.pack()

# File Browser Button
file_browse_button = tk.Button(main_frame, text="Browse", command=choose_file)
file_browse_button.pack()

# MiniBar Button
minibar_prep_button = tk.Button(main_frame, text="Run MiniBar", command=gui_minibar_prep)
minibar_prep_button.pack(fill=tk.BOTH, expand=True)

# NGSpecies ID Button
ngsid_prep_button = tk.Button(main_frame, text="Run NGSpeciesID", command=gui_ngsid_prep)
ngsid_prep_button.pack(fill=tk.BOTH, expand=True)

# Summary Button
summary_prep_button = tk.Button(main_frame, text="Run MycoMap Summary", command=gui_summary_prep)
summary_prep_button.pack(fill=tk.BOTH, expand=True)

# Run main tkinter GUI loop
root.mainloop()
