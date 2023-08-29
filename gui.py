import tkinter as tk

from fasta import Fasta
from blosum62 import blosum62
from globalAl import NeedlemannWunsch as nw
from localAl import SmithWaterman as sw

class GuiApp(tk.Tk):
    def __init__(self):
        super().__init__()
        
        self.title('Menu')
        self.geometry('300x200')
        #self.configure(background='grey')

        self.label = tk.Label(self, text='Choose a tool')
        self.label.pack()
        
        options = [
            'Needlemann-Wunsch Algorithm',
            'Smith-Waterman Algorithm',
            'Best Cost Matrix']
        
        tools = tk.StringVar()
        tools.set('Neeleman-Wunsch Algorithm')
        drop = tk.OptionMenu(self, tools, *options)
        #drop.configure(bg='grey')
        drop.place(x=50,y=30)

        self.button = tk.Button(self, text='Apply')
        #self.button.configure(bg='grey')
        self.button.place(x=10,y=150)

if __name__ == '__main__':
    app = GuiApp()
    app.mainloop()

