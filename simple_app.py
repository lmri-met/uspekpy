"""
# USpekPy Application

This script defines a graphical user interface (GUI) application using Tkinter for USpekPy, a tool for handling ISO 4037:2019 magnitudes and uncertainties.

## MainWindow Class

### Attributes
- `root`: ThemedTk root widget.
- `frame`: ttk.Frame containing UI widgets.
- `header`: ttk.Label displaying the application title.
- `label1`, `label2`, `label3`, `label4`, `label5`: ttk.Labels for various information or instructions.
- `entry1`, `entry2`, `entry3`, `entry4`, `entry5`: ttk.Entries for user input.
- `button`: ttk.Button to trigger a specific action (in this case, to run the `get_iso_magnitudes` method).

### Methods

#### `__init__(self, root)`
Initialize the main window of the application.

#### `grid_widgets(self)`
Position widgets on the grid layout within the main window.

#### `get_iso_magnitudes(self)`
Open a new window with "Hello World" text.

## main() Function

### Purpose
The main entry point function for the application.

### Steps
1. Creates a themed Tkinter root window with the "arc" theme.
2. Sets the title of the root window to "USpekPy".
3. Sets the application icon to 'assets/letter-u.ico'.
4. Creates an instance of the `MainWindow` class, passing the root window as an argument.
5. Positions the widgets within the main window using the `grid_widgets()` method of the `MainWindow` instance.
6. Starts the main event loop, which waits for events such as user input, using `root.mainloop()`.
"""
import tkinter as tk
from tkinter import ttk

from ttkthemes import ThemedTk


class MainWindow:
    """
    Class representing the main window of the application.

    This class manages the main window of the application, including the layout and interaction with UI widgets.

    Attributes
    ----------
    root : object
        ThemedTk root widget.
    frame : ttk.Frame
        Frame containing UI widgets.
    header : ttk.Label
        Label displaying the application title.
    label1, label2, label3, label4, label5 : ttk.Label
        Labels for various information or instructions.
    entry1, entry2, entry3, entry4, entry5 : ttk.Entry
        Entry fields for user input.
    button : ttk.Button
        Button to trigger a specific action (in this case, to run the `get_iso_magnitudes` method).
    """

    def __init__(self, root):
        """
        Initialize the main window.

        This method initializes the main window of the application. It performs the following steps:

        1. Initializes the main window with a frame to hold widgets and positions it within the root window.
        2. Sets up styling for widgets such as entry fields and labels.
        3. Creates labels, entry fields, and a button inside the frame for user interaction.

        Parameters
        ----------
        root : object
            ThemedTk root widget.
        """
        # Initializes the main window with a frame to hold widgets and position it within the root window
        self.root = root
        self.frame = ttk.Frame(self.root, padding="50 0 50 50")
        self.frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # Styling the widgets
        entry_width = 35
        header_font = ('', 20, 'bold')
        normal_font = ('', 14, '')

        # Apply style to button
        style = ttk.Style()
        style.configure('Custom.TButton', font=('', 14, 'bold'))

        # Create labels, entries, and a button
        self.header = ttk.Label(self.frame, text="USpekPy: ISO 4037:2019 Magnitudes & uncertainties", font=header_font)
        self.label1 = ttk.Label(self.frame, text="Mono-energetic conversion coefficients", font=normal_font)
        self.label2 = ttk.Label(self.frame, text="Mass transmission coefficients", font=normal_font)
        self.label3 = ttk.Label(self.frame, text="Uncertainty of mass transmission coefficients", font=normal_font)
        self.label4 = ttk.Label(self.frame, text="Number of simulations", font=normal_font)
        self.label5 = ttk.Label(self.frame, text="Select output folder", font=normal_font)
        self.entry1 = ttk.Entry(self.frame, width=entry_width, font=normal_font)
        self.entry2 = ttk.Entry(self.frame, width=entry_width, font=normal_font)
        self.entry3 = ttk.Entry(self.frame, width=entry_width, font=normal_font)
        self.entry4 = ttk.Entry(self.frame, width=entry_width, font=normal_font)
        self.entry5 = ttk.Entry(self.frame, width=entry_width, font=normal_font)
        self.button = ttk.Button(self.frame, text='Run', command=self.get_iso_magnitudes, style='Custom.TButton')

    def grid_widgets(self):
        """
        Position widgets on the grid.

        This method is responsible for positioning the widgets (labels, entries, button) on the grid layout within the main window.

        Parameters
        ----------
        self : object
            Instance of the MainWindow class.

        Returns
        -------
        None

        Notes
        -----
        This method utilizes the grid geometry manager to arrange widgets in rows and columns within the main window's frame.

        - The `pad_x` and `pad_y` variables define the padding between widgets.
        - The `sticky` variable determines how widgets should stick to the sides of the cells they occupy.
        - The `header_pad_y` variable provides additional padding for the header label.
        - The `button_sticky` variable specifies the alignment of the button within its cell.

        The grid layout is as follows:
        - Header label spans two columns at the top.
        - Labels, entries, and the button are placed in different rows and columns with specified padding and alignment.
        """
        pad_x = 20
        pad_y = 10
        sticky = 'w'
        header_pad_y = 50
        button_sticky = 'e'

        self.header.grid(row=0, column=0, columnspan=2, pady=header_pad_y)
        self.label1.grid(row=1, column=0, pady=pad_y, padx=pad_x, sticky=sticky)
        self.label2.grid(row=1, column=1, pady=pad_y, padx=pad_x, sticky=sticky)
        self.entry1.grid(row=2, column=0, pady=pad_y, padx=pad_x, sticky=sticky)
        self.entry2.grid(row=2, column=1, pady=pad_y, padx=pad_x, sticky=sticky)
        self.label3.grid(row=3, column=0, pady=pad_y, padx=pad_x, sticky=sticky)
        self.label4.grid(row=3, column=1, pady=pad_y, padx=pad_x, sticky=sticky)
        self.entry3.grid(row=4, column=0, pady=pad_y, padx=pad_x, sticky=sticky)
        self.entry4.grid(row=4, column=1, pady=pad_y, padx=pad_x, sticky=sticky)
        self.label5.grid(row=5, column=0, pady=pad_y, padx=pad_x, sticky=sticky)
        self.entry5.grid(row=6, column=0, pady=pad_y, padx=pad_x, sticky=sticky)
        self.button.grid(row=6, column=1, pady=pad_y, padx=pad_x, sticky=button_sticky)

    def get_iso_magnitudes(self):
        """
        Open a new window with "Hello World" text.

        This method is responsible for opening a new window when called. It performs the following steps:

        1. Creates a new Toplevel window (`new_window`) which is a separate window on top of the main application window.
        2. Creates a label (`label`) widget inside the new window with the text "Hello World".
        3. Packs the label widget within the new window, which makes it visible to the user.

        This method is typically used as a placeholder or example for opening new windows and displaying information.
        """
        new_window = tk.Toplevel(self.root)
        label = ttk.Label(new_window, text="Hello World")
        label.pack()


def main():
    """
    The main entry point function for the application.

    This function serves as the entry point for the application. It performs the following steps:

    1. Creates a themed Tkinter root window with the "arc" theme.
    2. Sets the title of the root window to "USpekPy".
    3. Sets the application icon to 'assets/letter-u.ico'.
    4. Creates an instance of the `MainWindow` class, passing the root window as an argument.
    5. Positions the widgets within the main window using the `grid_widgets()` method of the `MainWindow` instance.
    6. Starts the main event loop, which waits for events such as user input, using `root.mainloop()`.

    Returns
    -------
    None
    """
    root = ThemedTk(theme="arc")
    root.title("USpekPy")
    root.iconbitmap('assets/letter-u.ico')
    window = MainWindow(root)
    window.grid_widgets()
    root.mainloop()


if __name__ == "__main__":
    main()
