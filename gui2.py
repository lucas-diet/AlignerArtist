import tkinter as tk
from tkinter import ttk

from fasta import Fasta
from blosum62 import blosum62
from globalAl import NeedlemannWunsch as NW
from localAl import SmithWaterman as SW
from msa import MultipleSequenzalignment as MSA
from bcm import BestCostMatrix as BCM

class App(tk.Tk):

    def __init__(self):
        super().__init__()
        
        #######################################
        ###            Start Window         ###
        #######################################
        
        self.title('Menu')
        self.geometry('300x200')
        self.label = tk.Label(self, text='Choose a tool')
        self.label.pack()
        
        self.tool_options = [
            'Needlemann-Wunsch Algorithm',
            'Smith-Waterman Algorithm',
            'Multiples-Sequence-Alignment',
            'Best Cost Matrix',
            'Carillo-Lipman-Barrier']
        
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

            self.dp_button = tk.Button(self.tool_window, text='Show DP', command=self.showDPNW)
            self.dp_button.place(x=100,y=300)
            self.dp_button.config(state='disabled')

            self.output_als = tk.Text(self.tool_window, wrap='none', height=15, width=80)
            #self.output_als.pack(fill="both", expand=True)
            self.output_als.place(x=90,y=340)

            self.output_als_scrol_v = tk.Scrollbar(self.tool_window, command=self.output_als.yview, orient='vertical')
            self.output_als_scrol_h = tk.Scrollbar(self.tool_window, command=self.output_als.xview, orient='horizontal')
            self.output_als.config(xscrollcommand=self.output_als_scrol_h.set, yscrollcommand=self.output_als_scrol_v.set)
            self.output_als.config(bg='white', fg='black', state='disabled')

        
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

            self.dp_button = tk.Button(self.tool_window, text='Show DP', command=self.showDPSW)
            self.dp_button.place(x=100,y=300)
            self.dp_button.config(state='disabled')

            self.output_als = tk.Text(self.tool_window, wrap='none', height=15, width=80)
            #self.output_als.pack(fill="both", expand=True)
            self.output_als.place(x=90,y=340)

            self.output_als_scrol_v = tk.Scrollbar(self.tool_window, command=self.output_als.yview, orient='vertical')
            self.output_als_scrol_h = tk.Scrollbar(self.tool_window, command=self.output_als.xview, orient='horizontal')
            self.output_als.config(xscrollcommand=self.output_als_scrol_h.set, yscrollcommand=self.output_als_scrol_v.set)
            self.output_als.config(bg='white', fg='black', state='disabled')

        
        if self.tools.get() == self.tool_options[2]: 
        
            self.tool_window.title('Multiple-Sequenzalignment')
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

            self.seq3 = tk.Entry(self.tool_window, width=50)
            self.seq3.config(bg='white', fg='black')
            self.seq3.place(x=90,y=160)

            self.al_button = tk.Button(self.tool_window, text='Align', command=self.alignInputMSA)
            self.al_button.place(x=10,y=300)

            self.output_als = tk.Text(self.tool_window, wrap='none', height=15, width=80)
            self.output_als.pack(fill="both", expand=True)
            self.output_als.place(x=90,y=340)

            self.output_als_scrol_v = tk.Scrollbar(self.tool_window, command=self.output_als.yview, orient='vertical')
            self.output_als_scrol_h = tk.Scrollbar(self.tool_window, command=self.output_als.xview, orient='horizontal')

            self.output_als.config(xscrollcommand=self.output_als_scrol_h.set, yscrollcommand=self.output_als_scrol_v.set)
            self.output_als.config(bg='white', fg='black', state='disabled')


        if self.tools.get() == self.tool_options[3]:
            self.tool_window.title('Best-Cost-Matrix')
            self.tool_window.geometry('700x300')
            
            self.sequence_type = tk.OptionMenu(self.tool_window, self.sequence_types, *self.sequence_options)
            self.sequence_type.config(fg='black')
            self.sequence_type.place(x=20,y=20)

            self.seq1 = tk.Entry(self.tool_window, width=50)
            self.seq1.config(bg='white', fg='black')
            self.seq1.place(x=90,y=80)

            self.seq2 = tk.Entry(self.tool_window,width=50)
            self.seq2.config(bg='white', fg='black')
            self.seq2.place(x=90,y=120)

            self.al_button = tk.Button(self.tool_window, text='Align', command=self.alignInputBCM)
            self.al_button.place(x=10,y=200)

            self.dp_button = tk.Button(self.tool_window, text='DP', command=self.showDPNW)
            self.dp_button.place(x=100,y=200)
            self.dp_button.config(state='disabled')

            self.dprev_button = tk.Button(self.tool_window, text='DP-Rev', command=self.showDP_rev)
            self.dprev_button.place(x=170,y=200)
            self.dprev_button.config(state='disabled')

            self.m_button = tk.Button(self.tool_window, text='M', command=self.showM)
            self.m_button.place(x=270,y=200)
            self.m_button.config(state='disabled')


        if self.tools.get() == self.tool_options[4]:
            self.tool_window.title('Needlemann-Wunsch Algorithm')
            self.tool_window.geometry('700x300')
            
            self.sequence_type = tk.OptionMenu(self.tool_window, self.sequence_types, *self.sequence_options)
            self.sequence_type.config(fg='black')
            self.sequence_type.place(x=20,y=20)

            self.seq1 = tk.Entry(self.tool_window, width=50)
            self.seq1.config(bg='white', fg='black')
            self.seq1.place(x=90,y=80)

            self.seq2 = tk.Entry(self.tool_window,width=50)
            self.seq2.config(bg='white', fg='black')
            self.seq2.place(x=90,y=120)

            self.seq3 = tk.Entry(self.tool_window,width=50)
            self.seq3.config(bg='white', fg='black')
            self.seq3.place(x=90,y=160)

            self.al_score = tk.Entry(self.tool_window,width=5)
            self.al_score.config(bg='white', fg='black')
            self.al_score.place(x=90,y=200)

            self.al_button = tk.Button(self.tool_window, text='Align', command=None)
            self.al_button.place(x=10,y=250)

            self.dp_button = tk.Button(self.tool_window, text='Show M', command=None)
            self.dp_button.place(x=100,y=250)
            self.dp_button.config(state='disabled')
            
            '''
            self.output_als = tk.Text(self.tool_window, wrap='none', height=15, width=80)
            #self.output_als.pack(fill="both", expand=True)
            self.output_als.place(x=90,y=340)
            
            self.output_als_scrol_v = tk.Scrollbar(self.tool_window, command=self.output_als.yview, orient='vertical')
            self.output_als_scrol_h = tk.Scrollbar(self.tool_window, command=self.output_als.xview, orient='horizontal')
            self.output_als.config(xscrollcommand=self.output_als_scrol_h.set, yscrollcommand=self.output_als_scrol_v.set)
            self.output_als.config(bg='white', fg='black', state='disabled')
            '''
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
                self.msg.set('Illigal sequence \t \t')
                self.error.config(text=self.msg.get(), fg='yellow', font=('20'))

            else:
                #print(self.f1.getSequence())
                self.dp_button.config(state='active')

                self.dp = nw.calcualteDP(self.f1.getSequenceType(), self.s1, self.s2)
                als = nw.trackbackGlobalAlignments(self.dp, self.f1.getSequenceType(), self.s1, self.s2, len(self.s1), len(self.s2))
                
                formatted_als = ''
                self.output_als.delete('1.0',tk.END)
                for al in als:
                    formatted_als += '\n'.join(map(str, al)) + '\n\n'
                    
                als.clear()
                self.output_als.delete('1.0', tk.END)
                self.output_als.insert('end', formatted_als.upper())

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
                
            self.sw = SW()

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
                self.msg.set('Illigal sequence \t \t')
                self.error.config(text=self.msg.get(), fg='yellow', font=('20'))

            else:
                #print(self.f1.getSequence())
                self.dp_button.config(state='active')
                self.dp = self.sw.calcualteDP(self.f1.getSequenceType(), self.s1, self.s2)
                als = self.sw.trackbackLocalAlignment(self.dp, self.f1.getSequenceType(), self.s1, self.s2)
                
                formatted_als = ''
                self.output_als.delete('1.0',tk.END)
                for al in als:
                    formatted_als += '\n'.join(map(str, al)) + '\n\n'
                    
                als.clear()
                self.output_als.delete('1.0', tk.END)
                self.output_als.insert('end', formatted_als.upper())

    def alignInputMSA(self):
        
        self.s1 = self.seq1.get()
        self.s2 = self.seq2.get()
        self.s3 = self.seq3.get()

        self.f1 = Fasta()
        self.f2 = Fasta()
        self.f3 = Fasta()

        self.error = tk.Label(self.tool_window)
        self.error.place(x=90,y=190)
        self.msg = tk.StringVar()
        
        if len(self.s1) == 0 or len(self.s2) == 0 or len(self.s3) == 0:
            self.msg.set('Please enter a sequence')
            self.error.config(text=self.msg.get(), fg='yellow', font=('20'))
        
        else:
            self.msg.set('\t \t ')
            self.error.config(text=self.msg.get())
                
            self.output_als.config(state='normal')
                
            msa = MSA()

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
                self.msg.set('Illigal sequence \t \t')
                self.error.config(text=self.msg.get(), fg='yellow', font=('20'))
                
            else:
                #print(self.f1.getSequence())
                dp = msa.calcualteDP(self.s1, self.s2, self.s3)
                als = msa.trackbackMSA(dp, self.s1, self.s2, self.s3, len(self.s1), len(self.s2), len(self.s3)) 
                
                formatted_als = ''
                self.output_als.delete('1.0',tk.END)
                for al in als:
                    formatted_als += '\n'.join(map(str, al)) + '\n\n'
        
                als.clear()
                self.output_als.delete('1.0', tk.END)
                self.output_als.insert('end', formatted_als.upper())

                self.output_als.config(state='disabled')

                self.seq1.delete(0,tk.END)
                self.seq2.delete(0,tk.END)
                self.seq3.delete(0, tk.END)

    def alignInputBCM(self):
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

            bcm = BCM()

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
                self.msg.set('Illigal sequence \t \t')
                self.error.config(text=self.msg.get(), fg='yellow', font=('20'))
            
            else:
                self.dp_button.config(state='active')
                self.dprev_button.config(state='active')
                self.m_button.config(state='active')

                self.dp = NW().calcualteDP(self.f1.getSequenceType(), self.s1, self.s2)
                self.dprev = bcm.calculateDPRev(self.s1, self.s2)
                self.m = bcm.calculateM(self.dp, self.dprev)
                
    def showDPNW(self):

        self.dp_win = tk.Tk()
        self.dp_win.title('DP-Matrix')
        self.dp_win.geometry('500x500')

        canvas = tk.Canvas(self.dp_win)
        canvas.pack(fill=tk.BOTH, expand=True)

        content_frame = tk.Frame(canvas)
        canvas.create_window((0, 0), window=content_frame, anchor='nw')

        v_scrollbar = tk.Scrollbar(canvas, orient=tk.VERTICAL)
        v_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        canvas.config(yscrollcommand=v_scrollbar.set)
        v_scrollbar.config(command=canvas.yview)

        h_scrollbar = tk.Scrollbar(self.dp_win, orient=tk.HORIZONTAL, command=canvas.xview)
        h_scrollbar.pack(side=tk.BOTTOM, fill=tk.X)
        canvas.config(xscrollcommand=h_scrollbar.set)

        for x in range(0, len(self.dp)):
            for y in range(0, len(self.dp[0])):
                entry = tk.Text(content_frame, width=5, height=2)
                entry.grid(row=x, column=y)
                entry.insert(tk.END, self.dp[x][y])
                entry.config(state='disabled', bg='white', fg='black')
                if x == len(self.dp)-1 and y == len(self.dp[0])-1:
                    entry.config(state='disabled', bg='red', fg='black')

                entry.tag_configure("center", justify="center", font=('Arial', '20'), )
                entry.tag_add("center", "1.0", "end")

        content_frame.update_idletasks()
        canvas.config(scrollregion=canvas.bbox("all"))

    def showDPSW(self):

        self.dp_win = tk.Tk()
        self.dp_win.title('DP-Matrix')
        self.dp_win.geometry('500x500')

        canvas = tk.Canvas(self.dp_win)
        canvas.pack(fill=tk.BOTH, expand=True)

        content_frame = tk.Frame(canvas)
        canvas.create_window((0, 0), window=content_frame, anchor='nw')

        v_scrollbar = tk.Scrollbar(canvas, orient=tk.VERTICAL)
        v_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        canvas.config(yscrollcommand=v_scrollbar.set)
        v_scrollbar.config(command=canvas.yview)

        h_scrollbar = tk.Scrollbar(self.dp_win, orient=tk.HORIZONTAL, command=canvas.xview)
        h_scrollbar.pack(side=tk.BOTTOM, fill=tk.X)
        canvas.config(xscrollcommand=h_scrollbar.set)

        for x in range(0, len(self.dp)):
            for y in range(0, len(self.dp[0])):
                entry = tk.Text(content_frame, width=5, height=2)
                entry.grid(row=x, column=y)
                entry.insert(tk.END, self.dp[x][y])
                entry.config(state='disabled', bg='white', fg='black')

                if self.dp[x][y] == self.sw.getMaximalSimilarities(self.dp)[0][0]: 
                    entry.config(state='disabled', bg='red', fg='black')
                
                entry.tag_configure("center", justify="center", font=('Arial', '20'), )
                entry.tag_add("center", "1.0", "end")

        # Update the scroll region of the canvas
        content_frame.update_idletasks()
        canvas.config(scrollregion=canvas.bbox("all"))

    def showDP_rev(self):
        self.dp_win = tk.Tk()
        self.dp_win.title('DP-Reverse-Matrix')
        self.dp_win.geometry('500x500')

        canvas = tk.Canvas(self.dp_win)
        canvas.pack(fill=tk.BOTH, expand=True)

        content_frame = tk.Frame(canvas)
        canvas.create_window((0, 0), window=content_frame, anchor='nw')

        v_scrollbar = tk.Scrollbar(canvas, orient=tk.VERTICAL)
        v_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        canvas.config(yscrollcommand=v_scrollbar.set)
        v_scrollbar.config(command=canvas.yview)

        h_scrollbar = tk.Scrollbar(self.dp_win, orient=tk.HORIZONTAL, command=canvas.xview)
        h_scrollbar.pack(side=tk.BOTTOM, fill=tk.X)
        canvas.config(xscrollcommand=h_scrollbar.set)

        for x in range(0, len(self.dprev)):
            for y in range(0, len(self.dprev[0])):
                entry = tk.Text(content_frame, width=5, height=2)
                entry.grid(row=x, column=y)
                entry.insert(tk.END, self.dprev[x][y])
                entry.config(state='disabled', bg='white', fg='black')

                entry.tag_configure("center", justify="center", font=('Arial', '20'))
                entry.tag_add("center", "1.0", "end")

        # Update the scroll region of the canvas
        content_frame.update_idletasks()
        canvas.config(scrollregion=canvas.bbox("all"))

    def showM(self):
        self.dp_win = tk.Tk()
        self.dp_win.title('DP-Matrix')
        self.dp_win.geometry('500x500')

        canvas = tk.Canvas(self.dp_win)
        canvas.pack(fill=tk.BOTH, expand=True)

        content_frame = tk.Frame(canvas)
        canvas.create_window((0, 0), window=content_frame, anchor='nw')

        v_scrollbar = tk.Scrollbar(canvas, orient=tk.VERTICAL)
        v_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        canvas.config(yscrollcommand=v_scrollbar.set)
        v_scrollbar.config(command=canvas.yview)

        h_scrollbar = tk.Scrollbar(self.dp_win, orient=tk.HORIZONTAL, command=canvas.xview)
        h_scrollbar.pack(side=tk.BOTTOM, fill=tk.X)
        canvas.config(xscrollcommand=h_scrollbar.set)

        for x in range(0, len(self.m)):
            for y in range(0, len(self.m[0])):
                entry = tk.Text(content_frame, width=5, height=2)
                entry.grid(row=x, column=y)
                entry.insert(tk.END, self.m[x][y])
                entry.config(state='disabled', bg='white', fg='black')
                
                entry.tag_configure("center", justify="center", font=('Arial', '20'), )
                entry.tag_add("center", "1.0", "end")

        # Update the scroll region of the canvas
        content_frame.update_idletasks()
        canvas.config(scrollregion=canvas.bbox("all"))

if __name__ == '__main__':
    app = App()
    app.mainloop()

