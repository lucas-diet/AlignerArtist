from fasta import Fasta

class NeedlemannWunsch():

    def __init__(self, match=0,mismatch=1,gap=1):
        self.match = match
        self.mismatch = mismatch
        self.gap = gap