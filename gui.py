import tkinter as tk
from tkinter import ttk

from fasta import Fasta
from blosum62 import blosum62
from globalAl import NeedlemannWunsch as NW
from localAl import SmithWaterman as SW

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
        
        self.toolWindow = tk.Tk()

        sequenceOptions = ['Nucleotidsequence', 'Aminoacidsequence']
        sequenceTypes = tk.StringVar(self.toolWindow)
        sequenceTypes.set(sequenceOptions[0])

        if self.tools.get() == self.options[0]:
            self.toolWindow.title('Needlemann-Wunsch Algorithm')
            self.toolWindow.geometry('700x600')
            
            self.sequenceType = tk.OptionMenu(self.toolWindow, sequenceTypes, *sequenceOptions)
            self.sequenceType.config(fg='black')
            self.sequenceType.place(x=20,y=20)

            self.seq1 = tk.Entry(self.toolWindow, width=50)
            self.seq1.config(bg='white', fg='black')
            self.seq1.place(x=90,y=80)

            self.seq2 = tk.Entry(self.toolWindow,width=50)
            self.seq2.config(bg='white', fg='black')
            self.seq2.place(x=90,y=120)

            self.al_button = tk.Button(self.toolWindow, text='Align', command=self.alignInput)
            self.al_button.place(x=10,y=300)

            self.result = tk.StringVar()
            self.output = tk.Label(self.toolWindow, height=15, width=65)
            self.output.config(state='disabled') #bg='white', fg='black', 
            self.output.place(x=90,y=340)
        
        if self.tools.get() == self.options[1]:
            self.toolWindow.title('Smith-Waterman Algorithm')
            self.toolWindow.geometry('700x600')

            self.toolWindow.sequenceType = tk.OptionMenu(self.toolWindow, sequenceTypes, *sequenceOptions)
            self.toolWindow.sequenceType.config(fg='black')
            self.toolWindow.sequenceType.place(x=20,y=20)

            self.toolWindow.seq1 = tk.Entry(self.toolWindow,width=50)
            self.toolWindow.seq1.config(bg='white',fg='black')
            self.toolWindow.seq1.place(x=60,y=80)

            self.toolWindow.seq2 = tk.Entry(self.toolWindow,width=50)
            self.toolWindow.seq2.config(bg='white',fg='black')
            self.toolWindow.seq2.place(x=60,y=120)

            self.toolWindow.button = tk.Button(self.toolWindow, text='Align')
            self.toolWindow.button.place(x=10,y=300)

            self.toolWindow.output = tk.Text(self.toolWindow, height=15, width=65)
            self.toolWindow.output.config(bg='white', fg='black', state='disabled')
            self.toolWindow.output.place(x=90,y=340)

    def alignInput(self):
        s1 = self.seq1.get()
        s2 = self.seq2.get()
        t = ''
        
        '''
        if self.sequenceType == 'Nucleotidsequence':
            t = 'nt'
        elif self.sequenceType == 'Aminoacidsequence':
            t = 'aa'
        
        self.output.config(bg='white', fg='black', state='normal')
        nw = NW()
        dp = nw.calcualteDP('nt',s1,s2)
        als = nw.trackbackGlobalAlignments(dp,'nt',s1,s2,len(s1),len(s2))
        print(als)
        self.output.config(text=nw.printGlobalAlignments(als))
        self.output.config(bg='white', fg='black', state='disabled')
        '''

if __name__ == '__main__':
    app = GuiApp()
    app.mainloop()

