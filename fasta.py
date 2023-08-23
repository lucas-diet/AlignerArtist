from sequence import Sequence as sq

class Fasta():
    def __init__(self, header, sequence):
        self.header = header
        self.sequence = sequence

    def setHeader(self, header):
        self.header = header
    
    def getHeader(self):
        return self.header
    
    