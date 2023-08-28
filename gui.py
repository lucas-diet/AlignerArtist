import tkinter as tk
from tkinter import ttk
from tkinter.messagebox import showinfo

from fasta import Fasta
from blosum62 import blosum62
from globalAl import NeedlemannWunsch as nw
from localAl import SmithWaterman as sw


class GuiApp(tk.Tk):
    def __init__(self):
        super().__init__()
        
        self.title('My Awesome App')
        self.geometry('300x300')

        self.label = ttk.Label(self, text='Hello, Tkinter!')
        self.label.pack()

        self.button = ttk.Button(self, text='Click Me')
        #self.button['command'] = self.button_clicked
        self.button.pack()

if __name__ == '__main__':
    app = GuiApp()
    app.mainloop()