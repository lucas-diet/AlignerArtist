from sequence import Sequence as sq

class Fasta():

    def __init__(self, header, sequence):
        self.header = header
        self.sequence = sequence

    def checkHeader(self, header):
        if header[0] != '>':
            correctHeader = '>'+ header
            return correctHeader
        else:
            return header
        
    def setHeader(self, header):
        self.header = header
    
    def getHeader(self):
        return self.header
    

h = '>jj'
s = 'ACG'
a = 'nt'
f = Fasta(h,s).checkHeader(h)
print(f) 