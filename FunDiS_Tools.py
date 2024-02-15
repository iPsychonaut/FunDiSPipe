# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 15:04:51 2024

@author: theda
"""
# Base Python Imports
import subprocess, os, platform, multiprocessing, psutil, math

# Required Python Imports
from datetime import datetime

# Global output_area variable
CPU_THREADS = 1
PERCENT_RESOURCES = 0.75
DEFAULT_LOG_FILE = None
ENVIRONMENT_TYPE = None

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
        - This function is essential for managing log file versions, especially in long-running applications or in situations where log file preservation is crucial.
        - The numerical suffix increments for each new file created in the same location.
        - When 'use_numerical_suffix' is False, it refreshes the log file by clearing existing content.
    
    Considerations:
        - Ensure the directory for the log file exists, or handle directory creation within the function or externally.
    
    Examples:
        log_file_path = "logs/my_log.txt"
        generate_log_file(log_file_path, use_numerical_suffix=True)
    """
    # Extract directory and base from the log file path
    directory, base_filename = os.path.split(log_file_path)
    filename, ext = os.path.splitext(base_filename)

    # Construct the initial log file name
    if use_numerical_suffix and os.path.exists(log_file_path):
        counter = 1
        new_log_file_name = f"{filename}_{counter}{ext}"
        new_log_file_path = os.path.join(directory, new_log_file_name)

        # Increment the counter until a new file name is found
        while os.path.exists(new_log_file_path):
            counter += 1
            new_log_file_name = f"{filename}_{counter}{ext}"
            new_log_file_path = os.path.join(directory, new_log_file_name)

        log_file_path = new_log_file_path
    else:
        # Clear the existing log file or create a new one
        log_file_path = os.path.join(directory, f"{filename}{ext}")
        open(log_file_path, 'w').close()

    return log_file_path

def initialize_logging_environment(INPUT_FOLDER):
    """
    Initializes the logging environment based on the given input file path.

    This function sets up the logging environment by adjusting file paths according to the operating system in use, 
    ensuring file existence, and then generating a log file. It sets the global DEFAULT_LOG_FILE variable to the path 
    of the generated log file.

    Parameters:
        INPUT_FOLDER (str): Path to the folder which influences the log file generation.

    Global Variables:
        DEFAULT_LOG_FILE (str): The default path for the log file used throughout the logging process.

    Notes:
        - Adapts to different operating systems, making the logging system more flexible.
        - Prints unlogged messages to the console regarding environment detection and file existence.
        - Modifies the global DEFAULT_LOG_FILE variable.
        
    Considerations:
        - Verify the input folder's existence and accessibility before attempting to create a log file.

    Examples:
        input_folder = "data/input_data"
        initialize_logging_environment(input_folder)
    """
    global DEFAULT_LOG_FILE, ENVIRONMENT_TYPE

    # Extract the directory from INPUT_FOLDER and form the log file path
    input_dir = os.path.dirname(INPUT_FOLDER)
    log_file_name = os.path.basename(INPUT_FOLDER).split('.')[0] + "_log.txt"
    input_file_path = os.path.join(input_dir, log_file_name)

    # Determine the operating system
    os_name = platform.system()
    
    # Handle path adjustments for different operating systems
    if os_name == "Windows":
        print('UNLOGGED:\tWINDOWS ENVIRONMENT')
        ENVIRONMENT_TYPE = "WIN"
    elif os_name in ["Linux", "Darwin"]:  # Darwin is the system name for macOS
        print('UNLOGGED:\tLINUX/WSL/MAC ENVIRONMENT')
        ENVIRONMENT_TYPE = "LINUX/WSL/MAC"
    else:
        print(f'UNLOGGED ERROR:\tUnsupported OS: {os_name}')
        return

    # Generate the log file
    run_log = generate_log_file(input_file_path, use_numerical_suffix=False)
    DEFAULT_LOG_FILE = run_log

def log_print(input_message, log_file=None):
    """
    Logs a message to a file and prints it to the console with appropriate coloring.
    
    This function takes a message and logs it to the specified file. Additionally, the message is printed to the 
    console, potentially with specific coloring depending on the context.
    
    Parameters:
        input_message (str): Message to be logged and printed.
        log_file (str): Path to the log file.

    Notes:
        - This function serves as a centralized way to manage both logging and console output, ensuring consistency across the application.
        - The function uses a global default log file if none is specified.
        - Timestamps each log entry for easy tracking.
        - Utilizes color coding in the console to distinguish between different types of messages (e.g., errors, warnings).
        - Supports color coding for specific message types: NOTE, CMD, ERROR, WARN, and PASS.
        - Falls back to default (white) color if the message type is unrecognized.
    
    Considerations:
        - Consider the security implications of logging sensitive information.
        
    Examples:
        log_print("NOTE: Starting process")
        log_print("ERROR: An error occurred", log_file="error_log.txt")
    """
    # Access the global variable
    global DEFAULT_LOG_FILE
    # ANSI escape sequences for colors
    COLORS = {"grey": "\033[90m",
              "red": "\033[91m",
              "green": "\033[92m",
              "yellow": "\033[93m",
              "blue": "\033[94m",
              "magenta": "\033[95m",
              "cyan": "\033[96m",
              "white": "\033[97m",
              "reset": "\033[0m"}
    
    if log_file is None:
        log_file = DEFAULT_LOG_FILE

    now = datetime.now()
    message = f'[{now:%Y-%m-%d %H:%M:%S}]\t{input_message}'

    # Determine the print color
    message_type_dict = {
        'NOTE': 'blue',
        'CMD': 'cyan',
        'ERROR': 'red',
        'WARN': 'yellow',
        'PASS': 'green',
    }
    print_color = 'white'  # Default color
    for key, value in message_type_dict.items():
        if key.lower() in input_message.lower():
            print_color = value
            break

    # Writing the message to the log file
    try:
        with open(log_file, 'a') as file:
            print(message, file=file)
    except TypeError:
        print(f"UNLOGGED ERROR:\tUnable to load the log file provided: {log_file}")

    # Print the message with color
    color_code = COLORS.get(print_color, COLORS['white'])
    print(f"{color_code}{message}{COLORS['reset']}")

def run_subprocess_cmd(cmd_list, shell_check):
    """
    Executes a command using the subprocess.Popen and displays its output in real-time.

    This function is designed to execute shell commands from within a Python script. It uses subprocess.Popen to
    provide real-time output of the command to the command line window. It also logs the command execution details.

    Parameters:
        cmd_list (str or list of str): The command to be executed. Can be a single string or a list of strings
                                       representing the command and its arguments.
        shell_check (bool): If True, the command is executed through the shell. This is necessary for some 
                            commands, especially those that are built into the shell or involve shell features
                            like wildcard expansion.

    Features:
        - Uses subprocess.Popen for more control over command execution.
        - Captures the command's standard output and errors in real-time and displays them as they occur.
        - Waits for the command to complete and checks the return code to determine success or failure.
        - Logs the command, its real-time output, and any execution errors.

    Usage and Considerations:
        - Useful for executing commands where live feedback is important, especially for long-running commands.
        - Requires careful use of 'shell_check' due to potential security risks with shell commands.

    Example:
        run_subprocess_cmd(["ls", "-l"], shell_check=False)
        # This would execute 'ls -l' and display its output in real-time, while handling logging.
    """
    if isinstance(cmd_list, str):
        log_print(f"CMD:\t{cmd_list}")    
        process = subprocess.Popen(cmd_list, shell=shell_check, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    else:
        log_print(f"CMD:\t{' '.join(cmd_list)}")    
        process = subprocess.Popen(cmd_list, shell=shell_check, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

    # Iterating over the output
    for line in process.stdout:
        print(line, end='')

    # Wait for the process to complete and fetch the return code
    process.wait()

    if process.returncode != 0:
        log_print(f"ERROR:\tCommand failed with return code {process.returncode}")
    else:
        log_print(f"PASS:\tSuccessfully processed command: {' '.join(cmd_list)}" if isinstance(cmd_list, list) else cmd_list)
        
# Function to convert the user input PERCENT_RESOURCES into usuable cpu_threads and ram_gb values
def get_resource_values(PERCENT_RESOURCES):
    """
    Converts user input PERCENT_RESOURCES into usuable cpu_threads and ram_gb values.

    Arg:
        PERCENT_RESOURCES (float): Percentage of resources to use.

    Returns:
        cpu_threads (str): A count of the CPUs available for processing.
        ram_gb (str): A count of the RAM (in GB) available for processing.

    Notes:
        - Allows dynamic allocation of resources based on the system's current state.

    Considerations:
        - Ensure that the PERCENT_RESOURCES value is within an acceptable range to avoid over-utilization of system resources.

    Examples:
        cpu_threads, ram_gb = get_resource_values(0.5)  # Use 50% of system resources
    """
    # Get the number of CPUs available on the system
    num_cpus = multiprocessing.cpu_count()
    
    # Get the amount of RAM (GB) currently available
    mem_info = psutil.virtual_memory()
    
    # Calculate the number of threads as 80% of available CPUs & RAM
    cpu_threads = int(math.floor(num_cpus * PERCENT_RESOURCES))
    ram_gb = int(mem_info.total / (1024.0 ** 3) * PERCENT_RESOURCES)
    
    return cpu_threads, ram_gb 
