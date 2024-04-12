"""
USpekPy Application Module

This module defines the main window class for the USpekPy application. It provides functionality to interact with
various widgets, such as selecting files and folders, and running the application.

Classes:
    MainWindow: Main window class for the USpekPy application.

Functions:
    main: Main function to create and run the USpekPy application.

Usage:
    To run the USpekPy application, execute the script directly.

Example:
    python uspekpy.py
"""
import tkinter as tk
from tkinter import filedialog
from tkinter import ttk

from ttkthemes import ThemedTk


class MainWindow:
    """
    Main window class for the USpekPy application.

    This class represents the main window of the USpekPy application, providing functionality to interact
    with various widgets and perform actions such as selecting files and folders, and running the application.

    Attributes:
        root (tk.Tk or ThemedTk): The root window of the application.
        frame (ttk.Frame): The frame containing all the widgets.
        header_font (tuple): Font style for header text.
        normal_font (tuple): Font style for normal text.
        frame_pad_y (int): Padding for the frame.
        widget_pad_x (int): Horizontal padding for widgets.
        widget_pad_y (int): Vertical padding for widgets.
        widget_sticky (str): Alignment property for widgets.
        entry_width (int): Width of entry widgets.
        style1 (ttk.Style): Custom style for buttons.
        style2 (ttk.Style): Custom style for buttons.
        button_sticky (str): Alignment property for buttons.
        frame_header (ttk.Frame): Frame to contain header widgets.
        frame_input (ttk.Frame): Frame to contain input widgets.
        frame_buttons (ttk.Frame): Frame to contain button widgets.
        header (ttk.Label): Label for displaying the application header.
        label1 (ttk.Label): Label for describing the first input field.
        label2 (ttk.Label): Label for describing the second input field.
        label3 (ttk.Label): Label for describing the third input field.
        label4 (ttk.Label): Label for describing the fourth input field.
        label5 (ttk.Label): Label for describing the fifth input field.
        entry1 (ttk.Entry): Entry widget for the first input field.
        entry2 (ttk.Entry): Entry widget for the second input field.
        entry3 (ttk.Entry): Entry widget for the third input field.
        entry4 (ttk.Entry): Entry widget for the fourth input field.
        entry5 (ttk.Entry): Entry widget for the fifth input field.
        button1 (ttk.Button): Button for selecting a file path for the first input field.
        button2 (ttk.Button): Button for selecting a file path for the second input field.
        button5 (ttk.Button): Button for selecting a folder path for the fifth input field.
        button_run (ttk.Button) - Button for running the USpekPy application.

    Methods:
        __init__(root): Initialize the main window.
        styling_widgets(): Define styles and fonts for the widgets.
        create_widgets(): Create all the widgets for the main window.
        grid_widgets(): Grid all the widgets within the main window.
        select_file(entry): Open a file dialog to select a file path and insert it into an entry.
        select_folder(entry): Open a folder dialog to select a folder path and insert it into an entry.
        run(): Run the USpekPy application.
    """

    def __init__(self, root):
        """
        Initialize the main window.

        Parameters
        ----------
        root : tk.Tk or ThemedTk
            The root window of the application.
        """
        # Initialize the main root window of the Tkinter application, setting its title and its icon
        self.root = root
        self.root.title("USpekPy")
        self.root.iconbitmap('assets/letter-u.ico')

        # Create a frame to hold all the widgets with padding
        self.frame = ttk.Frame(self.root, padding="50 0 50 50")
        self.frame.grid(row=0, column=0)

        # Set up styles and fonts for the widgets
        self.styling_widgets()

        # Create all the widgets
        self.create_widgets()

        # Organize the layout of the widgets
        self.grid_widgets()

    def styling_widgets(self):
        """
        Define styles and fonts for the widgets.

        This method sets various styles and fonts used by the widgets within the main window.
        """
        # Define font styles for header and normal text
        self.header_font = ('', 20, 'bold')
        self.normal_font = ('', 14, '')

        # Define padding for the frame and widgets
        self.frame_pad_y = 20
        self.widget_pad_x = 20
        self.widget_pad_y = 10

        # Define sticky property for widget alignment
        self.widget_sticky = 'w'

        # Define width for entry widgets
        self.entry_width = 35

        # Configure custom button styles
        self.style1 = ttk.Style()
        self.style1.configure('Custom1.TButton', font=('', 14, 'bold'), foreground='RoyalBlue2')
        self.style2 = ttk.Style()
        self.style2.configure('Custom2.TButton', font=self.normal_font)

        # Define sticky property for button alignment
        self.button_sticky = 'e'

    def create_widgets(self):
        """
        Create all the widgets for the main window.

        This method initializes all the widgets including labels, entries, buttons, etc.
        """
        # Create frames to organize the layout
        self.frame_header = ttk.Frame(self.frame)
        self.frame_input = ttk.Frame(self.frame)
        self.frame_buttons = ttk.Frame(self.frame)

        # Create labels for various inputs
        self.header = ttk.Label(self.frame_header, text="USpekPy: ISO 4037:2019 Magnitudes & uncertainties",
                                font=self.header_font)
        self.label1 = ttk.Label(self.frame_input, text="Mono-energetic conversion coefficients", font=self.normal_font)
        self.label2 = ttk.Label(self.frame_input, text="Mass transmission coefficients", font=self.normal_font)
        self.label3 = ttk.Label(self.frame_input, text="Uncertainty of mass transmission coefficients",
                                font=self.normal_font)
        self.label4 = ttk.Label(self.frame_input, text="Number of simulations", font=self.normal_font)
        self.label5 = ttk.Label(self.frame_input, text="Select output folder", font=self.normal_font)

        # Create entry widgets for input fields
        self.entry1 = ttk.Entry(self.frame_input, width=self.entry_width, font=self.normal_font, state='disabled')
        self.entry2 = ttk.Entry(self.frame_input, width=self.entry_width, font=self.normal_font, state='disabled')
        self.entry3 = ttk.Entry(self.frame_input, width=self.entry_width, font=self.normal_font)
        self.entry4 = ttk.Entry(self.frame_input, width=self.entry_width, font=self.normal_font)
        self.entry5 = ttk.Entry(self.frame_input, width=self.entry_width, font=self.normal_font, state='disabled')

        # Create buttons for file/folder selection and running the application
        self.button1 = ttk.Button(self.frame_input, text="Select", command=lambda: self.select_file(self.entry1),
                                  style='Custom2.TButton')
        self.button2 = ttk.Button(self.frame_input, text="Select", command=lambda: self.select_file(self.entry2),
                                  style='Custom2.TButton')
        self.button5 = ttk.Button(self.frame_input, text="Select", command=lambda: self.select_folder(self.entry5),
                                  style='Custom2.TButton')
        self.button_run = ttk.Button(self.frame_buttons, text='Run', command=self.run, style='Custom1.TButton')

    def grid_widgets(self):
        """
        Grid all the widgets within the main window.

        This method organizes the layout of all the widgets within the main window.
        """
        # Grid frames
        self.frame_header.grid(row=0, column=0, pady=self.frame_pad_y)
        self.frame_input.grid(row=1, column=0, pady=self.frame_pad_y)
        self.frame_buttons.grid(row=2, column=0, pady=self.frame_pad_y, sticky=self.button_sticky)

        # Grid labels
        self.header.grid(row=0, column=0)
        self.label1.grid(row=1, column=0, pady=self.widget_pad_y, padx=self.widget_pad_x, sticky=self.widget_sticky)
        self.label2.grid(row=2, column=0, pady=self.widget_pad_y, padx=self.widget_pad_x, sticky=self.widget_sticky)
        self.label3.grid(row=3, column=0, pady=self.widget_pad_y, padx=self.widget_pad_x, sticky=self.widget_sticky)
        self.label4.grid(row=4, column=0, pady=self.widget_pad_y, padx=self.widget_pad_x, sticky=self.widget_sticky)
        self.label5.grid(row=5, column=0, pady=self.widget_pad_y, padx=self.widget_pad_x, sticky=self.widget_sticky)

        # Grid entries
        self.entry1.grid(row=1, column=1, padx=self.widget_pad_x, sticky=self.widget_sticky)
        self.entry2.grid(row=2, column=1, padx=self.widget_pad_x, sticky=self.widget_sticky)
        self.entry3.grid(row=3, column=1, padx=self.widget_pad_x, sticky=self.widget_sticky)
        self.entry4.grid(row=4, column=1, padx=self.widget_pad_x, sticky=self.widget_sticky)
        self.entry5.grid(row=5, column=1, padx=self.widget_pad_x, sticky=self.widget_sticky)

        # Grid buttons
        self.button1.grid(row=1, column=2, padx=self.widget_pad_x, sticky=self.widget_sticky)
        self.button2.grid(row=2, column=2, padx=self.widget_pad_x, sticky=self.widget_sticky)
        self.button5.grid(row=5, column=2, padx=self.widget_pad_x, sticky=self.widget_sticky)
        self.button_run.grid(row=0, column=2, padx=self.widget_pad_x, sticky=self.button_sticky)

    @staticmethod
    def select_file(entry):
        """
        Open a file dialog to select a file path and insert it into an entry.

        Parameters
        ----------
        entry : ttk.Entry
            The entry widget to insert the selected file path.
        """
        # Open a file dialog to select a file path
        file_path = filedialog.askopenfilename()

        # Enable the entry widget to allow editing
        entry.config(state='normal')

        # Delete any existing text in the entry widget
        entry.delete(0, tk.END)

        # Insert the selected file path into the entry widget
        entry.insert(0, file_path)

        # Disable the entry widget to prevent further editing
        entry.config(state='disabled')

    @staticmethod
    def select_folder(entry):
        """
        Open a folder dialog to select a folder path and insert it into an entry.

        Parameters
        ----------
        entry : ttk.Entry
            The entry widget to insert the selected folder path.
        """
        # Open a folder dialog to select a folder path
        folder_path = filedialog.askdirectory()

        # Enable the entry widget to allow editing
        entry.config(state='normal')

        # Delete any existing text in the entry widget
        entry.delete(0, tk.END)

        # Insert the selected folder path into the entry widget
        entry.insert(0, folder_path)

        # Disable the entry widget to prevent further editing
        entry.config(state='disabled')

    def run(self):
        """
        Run the USpekPy application.

        This method retrieves input values from the entry widgets and displays them in a new window.
        """
        # Get input values
        conversion_coefficients_file = self.entry1.get()
        transmission_coefficients_file = self.entry2.get()
        transmission_coefficients_uncertainty = self.entry3.get()
        simulations_number = self.entry4.get()
        output_folder = self.entry5.get()

        # Create a new window to display input values
        new_window = tk.Toplevel(self.root)
        new_window.title("USpekPy")
        new_window.iconbitmap('assets/letter-u.ico')

        # Prepare items to display
        items = [
            f"Mono-energetic conversion coefficients: {conversion_coefficients_file}",
            f"Mass transmission coefficients: {transmission_coefficients_file}",
            f"Uncertainty of mass transmission coefficients: {transmission_coefficients_uncertainty}",
            f"Number of simulations: {simulations_number}",
            f"Output folder: {output_folder}"]

        # Create a listbox to display items
        listbox = tk.Listbox(new_window, width=50, height=20)
        for item in items:
            listbox.insert(tk.END, item)
        listbox.grid(row=0, column=0)


# Main function to create and run the USpekPy application
def main():
    """
    Main function to create and run the USpekPy application.

    This function creates the main Tkinter root window, initializes the main application window, and starts the event loop.

    Returns
    -------
    None
        This function does not return anything.
    """
    # Create a themed Tkinter root window with the "arc" theme
    root = ThemedTk(theme="arc")

    # Initialize the main window of the application
    window = MainWindow(root)

    # Start the Tkinter event loop
    root.mainloop()


# Entry point of the script
if __name__ == "__main__":
    # If the script is executed directly, call the main function
    main()
