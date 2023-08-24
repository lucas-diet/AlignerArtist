
class Sequence():
    
    def __init__(self, type,sequence):
        self.type = type
        self.sequence = sequence

    def checkSequence(self, type,sequence):
        nts = ['a','c','g','t']
        aas = ['a','c','h','m','t','r','q','i','f','w','n','e','l','p','y','d','g','k','s','v']

        if type == 'nt':
            for nt in sequence:
                if nt not in nts:
                    return 'error'
                    break
                else:
                    self.sequence = sequence
                    return sequence
        elif type == 'aa':
            for nt in sequence:
                if nt not in aas:
                    return 'error'
                else:
                    self.sequence = sequence
                    return sequence

    def setSequenceType(self, type):
        if type == 'nt':
            self.type = type
        elif type == 'aa':
            self.type = type
    
    def getSequencetype(self):
        return self.type

    def setSequence(self, sequence):
        type = self.type
        chekedSequence = self.checkSequence(type,sequence)
        if chekedSequence == 'error':
            return print('error')
        else:
            self.sequence = sequence
    
    def getSequence(self):
        return self.sequence
    