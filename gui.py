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

        self.label = tk.Label(self, text='Choose a tool')
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

        self.sequence_options = ['Nucleotidsequence', 'Aminoacidsequence']
        self.sequence_types = tk.StringVar(self.tool_window)
        self.sequence_types.set(self.sequence_options[0])

        if self.tools.get() == self.tool_options[0]:
            self.tool_window.title('Needlemann-Wunsch Algorithm')
            self.tool_window.geometry('700x600')
            
            self.sequence_type = tk.OptionMenu(self.tool_window, self.sequence_types, *self.sequence_options)
            self.sequence_type.config(fg='black')
            self.sequence_type.place(x=20,y=20)

            self.seq1 = tk.Entry(self.tool_window, width=50)
            self.seq1.config(bg='white', fg='black')
            self.seq1.place(x=90,y=80)

            self.seq2 = tk.Entry(self.tool_window,width=50)
            self.seq2.config(bg='white', fg='black')
            self.seq2.place(x=90,y=120)

            self.al_button = tk.Button(self.tool_window, text='Align', command=self.alignInputNW)
            self.al_button.place(x=10,y=300)

            self.output_als = tk.Text(self.tool_window, wrap='none', height=15, width=40)
            self.output_als.pack(fill="both", expand=True)
            self.output_als.place(x=90,y=340)

            self.output_dp = tk.Text(self.tool_window, wrap='none', height=15, width=40)
            self.output_dp.pack(fill="both", expand=True)
            self.output_dp.place(x=400,y=340)

            self.output_als_scrol_v = tk.Scrollbar(self.tool_window, command=self.output_als.yview, orient='vertical')
            self.output_als_scrol_h = tk.Scrollbar(self.tool_window, command=self.output_als.xview, orient='horizontal')

            self.output_als.config(xscrollcommand=self.output_als_scrol_h.set, yscrollcommand=self.output_als_scrol_v.set)

            self.output_dp_scrol_v = tk.Scrollbar(self.tool_window, command=self.output_dp.yview, orient='vertical')
            self.output_dp_scrol_h = tk.Scrollbar(self.tool_window, command=self.output_dp.xview, orient='horizontal')

            self.output_dp.config(xscrollcommand=self.output_dp_scrol_h.set,yscrollcommand=self.output_dp_scrol_v.set)
            
            self.output_als.config(bg='white', fg='black', state='disabled')
            self.output_dp.config(bg='white', fg='black', state='disabled')
        
        if self.tools.get() == self.tool_options[1]:
            self.tool_window.title('Smith-Waterman Algorithm')
            self.tool_window.geometry('700x600')
            
            self.sequence_type = tk.OptionMenu(self.tool_window, self.sequence_types, *self.sequence_options)
            self.sequence_type.config(fg='black')
            self.sequence_type.place(x=20,y=20)

            self.seq1 = tk.Entry(self.tool_window, width=50)
            self.seq1.config(bg='white', fg='black')
            self.seq1.place(x=90,y=80)

            self.seq2 = tk.Entry(self.tool_window,width=50)
            self.seq2.config(bg='white', fg='black')
            self.seq2.place(x=90,y=120)

            self.al_button = tk.Button(self.tool_window, text='Align', command=self.alignInputSW)
            self.al_button.place(x=10,y=300)

            self.output_als = tk.Text(self.tool_window, wrap='none', height=15, width=40)
            self.output_als.pack(fill="both", expand=True)
            self.output_als.place(x=90,y=340)

            self.output_dp = tk.Text(self.tool_window, wrap='none', height=15, width=40)
            self.output_dp.pack(fill="both", expand=True)
            self.output_dp.place(x=400,y=340)

            self.output_als_scrol_v = tk.Scrollbar(self.tool_window, command=self.output_als.yview, orient='vertical')
            self.output_als_scrol_h = tk.Scrollbar(self.tool_window, command=self.output_als.xview, orient='horizontal')

            self.output_als.config(xscrollcommand=self.output_als_scrol_h.set, yscrollcommand=self.output_als_scrol_v.set)

            self.output_dp_scrol_v = tk.Scrollbar(self.tool_window, command=self.output_dp.yview, orient='vertical')
            self.output_dp_scrol_h = tk.Scrollbar(self.tool_window, command=self.output_dp.xview, orient='horizontal')

            self.output_dp.config(xscrollcommand=self.output_dp_scrol_h.set,yscrollcommand=self.output_dp_scrol_v.set)
            
            self.output_als.config(bg='white', fg='black', state='disabled')
            self.output_dp.config(bg='white', fg='black', state='disabled')
         
    def alignInputNW(self):
        
        self.s1 = self.seq1.get().lower()
        self.s2 = self.seq2.get().lower()

        self.type = ''

        self.f1 = Fasta()
        self.f2 = Fasta()

        self.error = tk.Label(self.tool_window)
        self.error.place(x=90,y=150)
        self.msg = tk.StringVar()
        
        if len(self.s1) == 0 or len(self.s2) == 0:
            self.msg.set('Please enter a sequence')
            self.error.config(text=self.msg.get(), fg='yellow', font=('20'))
        
        else:
            self.msg.set('\t \t ')
            self.error.config(text=self.msg.get())
                
            self.output_als.config(state='normal')
            self.output_dp.config(state='normal')
                
            nw = NW()
            
            if self.sequence_types.get() == self.sequence_options[0]:
                self.type = 'nt'
                self.f1.setSequenceType(self.type)
                self.f2.setSequenceType(self.type)
                self.f1.setSequence(self.type, self.s1)
                self.f2.setSequence(self.type, self.s2)
                
            elif self.sequence_types.get() == self.sequence_options[1]:
                self.type = 'aa'
                self.f1.setSequenceType(self.type)
                self.f2.setSequenceType(self.type)
                self.f1.setSequence(self.type, self.s1)
                self.f2.setSequence(self.type, self.s2)
            
            if self.f1.getSequence() == 'ERROR' or self.f2.getSequence() == 'ERROR':
                #print(self.f1.getSequence())
                self.msg.set('Illigal sequence')
                self.error.config(text=self.msg.get(), fg='yellow', font=('20'))
            
            else:
                #print(self.f1.getSequence())
                dp = nw.calcualteDP(self.f1.getSequenceType(), self.s1, self.s2)
                als = nw.trackbackGlobalAlignments(dp, self.f1.getSequenceType(), self.s1, self.s2, len(self.s1), len(self.s2))
                
                formatted_als = ''
                self.output_als.delete('1.0',tk.END)
                for al in als:
                    formatted_als += '\n'.join(map(str, al)) + '\n\n'
                    
                als.clear()
                self.output_als.delete('1.0', tk.END)
                self.output_als.insert('end', formatted_als.upper())
                        
                formatted_dp = ''
                for line in dp:
                    formatted_dp += '\t'.join(map(str, line)) +'\n'
                        
                self.output_dp.delete('1.0', 'end')
                self.output_dp.update()
                self.output_dp.insert('end', formatted_dp)
                self.output_dp.insert('end', '\n')

                self.output_als.config(state='disabled')
                self.output_dp.config(state='disabled')

                self.seq1.delete(0,tk.END)
                self.seq2.delete(0,tk.END)
              
    def alignInputSW(self):
        
        self.s1 = self.seq1.get().lower()
        self.s2 = self.seq2.get().lower()

        self.type = ''

        self.f1 = Fasta()
        self.f2 = Fasta()

        self.error = tk.Label(self.tool_window)
        self.error.place(x=90,y=150)
        self.msg = tk.StringVar()
        
        if len(self.s1) == 0 or len(self.s2) == 0:
            self.msg.set('Please enter a sequence')
            self.error.config(text=self.msg.get(), fg='yellow', font=('20'))
        
        else:
            self.msg.set('\t \t ')
            self.error.config(text=self.msg.get())
                
            self.output_als.config(state='normal')
            self.output_dp.config(state='normal')
                
            sw = SW()
            
            if self.sequence_types.get() == self.sequence_options[0]:
                self.type = 'nt'
                self.f1.setSequenceType(self.type)
                self.f2.setSequenceType(self.type)
                self.f1.setSequence(self.type, self.s1)
                self.f2.setSequence(self.type, self.s2)
                
            elif self.sequence_types.get() == self.sequence_options[1]:
                self.type = 'aa'
                self.f1.setSequenceType(self.type)
                self.f2.setSequenceType(self.type)
                self.f1.setSequence(self.type, self.s1)
                self.f2.setSequence(self.type, self.s2)
            
            if self.f1.getSequence() == 'ERROR' or self.f2.getSequence() == 'ERROR':
                #print(self.f1.getSequence())
                self.msg.set('Illigal sequence')
                self.error.config(text=self.msg.get(), fg='yellow', font=('20'))
            
            else:
                #print(self.f1.getSequence())
                dp = sw.calcualteDP(self.f1.getSequenceType(), self.s1, self.s2)
                als = sw.trackbackLocalAlignment(dp, self.f1.getSequenceType(), self.s1, self.s2)
                
                formatted_als = ''
                self.output_als.delete('1.0',tk.END)
                for al in als:
                    formatted_als += '\n'.join(map(str, al)) + '\n\n'
                    
                als.clear()
                self.output_als.delete('1.0', tk.END)
                self.output_als.insert('end', formatted_als.upper())
                        
                formatted_dp = ''
                for line in dp:
                    formatted_dp += '\t'.join(map(str, line)) +'\n'
                        
                self.output_dp.delete('1.0', 'end')
                self.output_dp.update()
                self.output_dp.insert('end', formatted_dp)
                self.output_dp.insert('end', '\n')

                self.output_als.config(state='disabled')
                self.output_dp.config(state='disabled')

                self.seq1.delete(0,tk.END)
                self.seq2.delete(0,tk.END)
                 

if __name__ == '__main__':
    app = GuiApp()
    app.mainloop()

