import tkinter as tk
from tkinter import ttk

from fasta import Fasta
from blosum62 import blosum62
from globalAl import NeedlemannWunsch as nw
from localAl import SmithWaterman as sw

class GuiApp(tk.Tk):
    
    def __init__(self):
        super().__init__()
        self.title('Menu')
        self.geometry('300x200')

        self.label = ttk.Label(self, text='Choose a tool')
        self.label.pack()
        
        self.options = [
            'Needlemann-Wunsch Algorithm',
            'Smith-Waterman Algorithm',
            'Best Cost Matrix']
        
        self.tools = tk.StringVar()
        self.tools.set(self.options[0])
        self.drop = tk.OptionMenu(self, self.tools, *self.options)
        self.drop.config(fg='black')
        self.drop.place(x=30,y=30)

        self.button = tk.Button(self, text='Apply', command=self.openToolWindow)
        self.button.place(x=10,y=150)
    
    def openToolWindow(self):
        toolWindow = tk.Tk()
        
        sequenceOptions = ['Nucleotidsequence', 'Aminoacidsequence']
        sequenceTypes = tk.StringVar(toolWindow)
        sequenceTypes.set(sequenceOptions[0])

        if self.tools.get() == self.options[0]:
            toolWindow.title('Needlemann-Wunsch Algorithm')
            toolWindow.geometry('700x600')
            
            toolWindow.sequenceType = tk.OptionMenu(toolWindow, sequenceTypes, *sequenceOptions)
            toolWindow.sequenceType.config(fg='black')
            toolWindow.sequenceType.place(x=20,y=20)

            toolWindow.seq1 = tk.Entry(toolWindow,width=50)
            toolWindow.seq1.config(bg='white',fg='black')
            toolWindow.seq1.place(x=90,y=80)

            toolWindow.seq2 = tk.Entry(toolWindow,width=50)
            toolWindow.seq2.config(bg='white',fg='black')
            toolWindow.seq2.place(x=90,y=120)
        
        if self.tools.get() == self.options[1]:
            toolWindow.title('Smith-Waterman Algorithm')
            toolWindow.geometry('700x600')

            toolWindow.sequenceType = tk.OptionMenu(toolWindow, sequenceTypes, *sequenceOptions)
            toolWindow.sequenceType.config(fg='black')
            toolWindow.sequenceType.place(x=20,y=20)

            toolWindow.seq1 = tk.Entry(toolWindow,width=50)
            toolWindow.seq1.config(bg='white',fg='black')
            toolWindow.seq1.place(x=60,y=80)

            toolWindow.seq2 = tk.Entry(toolWindow,width=50)
            toolWindow.seq2.config(bg='white',fg='black')
            toolWindow.seq2.place(x=60,y=120)

if __name__ == '__main__':
    app = GuiApp()
    app.mainloop()

