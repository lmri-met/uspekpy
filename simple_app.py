import tkinter as tk
from tkinter import ttk
from ttkthemes import ThemedTk


class MainWindow:
    """
    Main window class for handling UI widgets.

    Parameters
    ----------
    root : object
        ThemedTk root widget.
    """

    def __init__(self, root):
        self.root = root
        self.frame = ttk.Frame(self.root, padding="50 0 50 50")  # Create a frame with padding for margins
        self.frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))  # Make it stick to all sides
        # self.root.columnconfigure(0, weight=1)
        # self.root.rowconfigure(0, weight=1)

        entry_width = 35
        header_font = ('', 20, 'bold')
        normal_font = ('', 14, '')

        # Apply style to button
        style = ttk.Style()
        style.configure('Custom.TButton', font=('', 14, 'bold'))

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
        self.button = ttk.Button(self.frame, text='Run', command=self.update_label, style='Custom.TButton')

    def grid_widgets(self):
        """
        Method to position widgets on the grid.

        Returns
        -------
        None.
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

    def update_label(self):
        """
        Open new window with "Hello World" text

        Returns
        -------
        None.
        """
        new_window = tk.Toplevel(self.root)
        label = ttk.Label(new_window, text="Hello World")
        label.pack()


def main():
    """
    The main entry point function for the standard GUI application.

    Returns
    -------
    None.
    """
    root = ThemedTk(theme="arc")
    root.title("USpekPy")
    root.iconbitmap('assets/letter-u.ico')
    window = MainWindow(root)
    window.grid_widgets()
    root.mainloop()


if __name__ == "__main__":
    main()
