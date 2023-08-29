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

        self.label = tk.Label(self, text='Choose a tool')
        self.label.pack()
        
        self.options = [
            'Needlemann-Wunsch Algorithm',
            'Smith-Waterman Algorithm',
            'Best Cost Matrix']
        
        self.tools = tk.StringVar()
        self.tools.set(self.options[0])
        self.drop = tk.OptionMenu(self, self.tools, *self.options)
        self.drop.place(x=30,y=30)

        self.button = tk.Button(self, text='Apply', command=self.openToolWindow)
        self.button.place(x=10,y=150)
    
    def openToolWindow(self):
        toolWindow = tk.Tk()
        if self.tools.get() == self.options[0]:
            toolWindow.title('Needlemann-Wunsch Algorithm')
            toolWindow.geometry('500x500')
        
        if self.tools.get() == self.options[1]:
            toolWindow.title('Smith-Waterman Algorithm')
            toolWindow.geometry('500x500')

if __name__ == '__main__':
    app = GuiApp()
    app.mainloop()

