import tkinter as tk
from tkinter import ttk

from fasta import Fasta
from blosum62 import blosum62
from globalAl import NeedlemannWunsch as NW
from localAl import SmithWaterman as SW

'''
menu_window = Tk()
menu_window.title('Menu')
menu_window.geometry('300x200')

tool_choose = Label(menu_window, text='Choose a tool').pack()

tool_options = [
            'Needlemann-Wunsch Algorithm',
            'Smith-Waterman Algorithm',
            'Best Cost Matrix']

tools = StringVar()
tools.set(tool_options[0])
drop_menu = OptionMenu(menu_window, tools, *tool_options)
drop_menu.config(fg='black')
drop_menu.place(x=30,y=30)

def openToolWindow():
        
        tool_window = Tk()

        sequence_options = ['Nucleotidsequence', 'Aminoacidsequence']
        sequence_types = StringVar(tool_window)
        sequence_types.set(sequence_options[0])

        if tools.get() == tool_options[0]:
            tool_window.title('Needlemann-Wunsch Algorithm')
            tool_window.geometry('700x600')
            
            sequence_type = OptionMenu(tool_window, sequence_types, *sequence_options)
            sequence_type.config(fg='black')
            sequence_type.place(x=20,y=20)

            seq1 = Entry(tool_window, width=50)
            seq1.config(bg='white', fg='black')
            seq1.place(x=90,y=80)

            seq2 = Entry(tool_window,width=50)
            seq2.config(bg='white', fg='black')
            seq2.place(x=90,y=120)

            al_button = Button(tool_window, text='Align', command=None)
            al_button.place(x=10,y=300)

            result = StringVar()
            output = Label(tool_window, height=15, width=65)
            output.config(state='disabled') #bg='white', fg='black', 
            output.place(x=90,y=340)

menu_btn = Button(menu_window, text='Apply', command=openToolWindow)
menu_btn.place(x=10,y=150)

menu_window.mainloop()
'''



class GuiApp(tk.Tk):
    
    def __init__(self):
        super().__init__()
        self.title('Menu')
        self.geometry('300x200')

        self.label = ttk.Label(self, text='Choose a tool')
        self.label.pack()
        
        self.tool_options = [
            'Needlemann-Wunsch Algorithm',
            'Smith-Waterman Algorithm',
            'Multiples-Sequence-Alignment',
            'Best Cost Matrix']
        
        self.tools = tk.StringVar()
        self.tools.set(self.tool_options[0])
        self.tool_drop_menu = tk.OptionMenu(self, self.tools, *self.tool_options)
        self.tool_drop_menu.config(fg='black')
        self.tool_drop_menu.place(x=30,y=30)

        self.tool_choice_btn = tk.Button(self, text='Apply', command=self.openToolWindow)
        self.tool_choice_btn.place(x=10,y=150)
    
    def openToolWindow(self):
        
        self.tool_window = tk.Tk()

        sequence_options = ['Nucleotidsequence', 'Aminoacidsequence']
        sequence_types = tk.StringVar(self.tool_window)
        sequence_types.set(sequence_options[0])

        if self.tools.get() == self.tool_options[0]:
            self.tool_window.title('Needlemann-Wunsch Algorithm')
            self.tool_window.geometry('700x600')
            
            self.sequence_type = tk.OptionMenu(self.tool_window, sequence_types, *sequence_options)
            self.sequence_type.config(fg='black')
            self.sequence_type.place(x=20,y=20)

            self.seq1 = tk.Entry(self.tool_window, width=50)
            self.seq1.config(bg='white', fg='black')
            self.seq1.place(x=90,y=80)

            self.seq2 = tk.Entry(self.tool_window,width=50)
            self.seq2.config(bg='white', fg='black')
            self.seq2.place(x=90,y=120)

            self.al_button = tk.Button(self.tool_window, text='Align', command=self.alignInput)
            self.al_button.place(x=10,y=300)

            #self.result = tk.StringVar()
            self.output_als = tk.Text(self.tool_window, wrap='none', height=15, width=40)
            self.output_als.pack(fill="both", expand=True)
            self.output_als.place(x=90,y=340)

            
            self.output_dp = tk.Text(self.tool_window, wrap='none', height=15, width=40)
            self.output_dp.pack(fill="both", expand=True)
            self.output_dp.place(x=400,y=340)

            self.output_als_scrol_v = tk.Scrollbar(self.tool_window, command=self.output_als.yview, orient='vertical')
            self.output_als_scrol_h = tk.Scrollbar(self.tool_window, command=self.output_als.xview, orient='horizontal')
            #self.output_als_scrol_v.pack(fill='y')
            #self.output_als_scrol_h.pack(fill="y")

            self.output_als.config(xscrollcommand=self.output_als_scrol_h.set, yscrollcommand=self.output_als_scrol_v.set)

            self.output_dp_scrol_v = tk.Scrollbar(self.tool_window, command=self.output_dp.yview, orient='vertical')
            self.output_dp_scrol_h = tk.Scrollbar(self.tool_window, command=self.output_dp.xview, orient='horizontal')
            #self.output_dp_scrol_v.pack(fill='y')
            #self.output_dp_scrol_h.pack(fill="y")

            self.output_dp.config(xscrollcommand=self.output_dp_scrol_h.set,yscrollcommand=self.output_dp_scrol_v.set)
            
            self.output_als.config(bg='white', fg='black', state='disabled')
            self.output_dp.config(bg='white', fg='black', state='disabled')
      
    
    def alignInput(self):
        
        s1 = self.seq1.get()
        s2 = self.seq2.get()
        t = ''

        
        if self.sequence_type == 'Nucleotidsequence':
            t = 'nt'
        elif self.sequence_type == 'Aminoacidsequence':
            t = 'aa'
        
        self.output_als.config(state='normal')
        self.output_dp.config(state='normal')
        
        nw = NW()
        dp = nw.calcualteDP('nt',s1,s2)
        als = nw.trackbackGlobalAlignments(dp,'nt',s1,s2,len(s1),len(s2))

        formatted_als = ''
        self.output_als.delete('1.0',tk.END)
        for al in als:
            formatted_als += '\n'.join(map(str, al)) + '\n\n'
        
        als.clear()
        self.output_als.delete('1.0', tk.END)
        self.output_als.insert('end', formatted_als)
        
        formatted_dp = ''
        for line in dp:
            formatted_dp += '\t'.join(map(str, line)) +'\n'
        print(formatted_dp)
        self.output_dp.delete('1.0', 'end')
        self.output_dp.update()
        self.output_dp.insert('end', formatted_dp)
        self.output_dp.insert('end', '\n')

        self.output_als.config(state='disabled')
        self.output_dp.config(state='disabled')
      

if __name__ == '__main__':
    app = GuiApp()
    app.mainloop()

