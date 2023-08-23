
class Sequence():
    def __init__(self, type, sequence):
        self.type = type
        self.sequence = sequence

    def checkSequence(self, type, sequence):
        nts = ['a','c','g','t']
        aas = ['a','c','h','m','t','r','q','i','f','w','n','e','l','p','y','d','g','k','s','v']
        
        if type == 'nt':
            for nt in nts:
                if nt not in sequence:
                    print('not in')
                    break

    def setSequenceType(self, type):
        if type == 'nt':
            self.type = type
        elif type == 'aa':
            self.type = type
    
    def getSequencetype(self):
        return self.type

    def setSequence(self, sequence):
        self.sequence = sequence
    
    def getSequence(self):
        return self.sequence

    


s = 'atg'
seq = Sequence('nt',s)
print(seq.checkSequence('nt',s))
    