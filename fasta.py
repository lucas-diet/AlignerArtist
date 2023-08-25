
class Fasta():

    def __init__(self, type,header,sequence):
        self.type = type
        self.header = header
        self.sequence = sequence

    def checkHeader(self, header):
        if header[0] != '>':
            correctHeader = '>'+ header
            return correctHeader
        else:
            return header

    def checkSequence(self, type,sequence):
        __nts = ['a','c','g','t']
        __aas = ['a','c','h','m','t','r','q','i','f','w','n','e','l','p','y','d','g','k','s','v']

        if type == 'nt':
            for nt in sequence:
                if nt not in __nts:
                    return 'error'
                    break
                else:
                    self.sequence = sequence
                    return sequence
        elif type == 'aa':
            for nt in sequence:
                if nt not in __aas:
                    return 'error'
                else:
                    self.sequence = sequence
                    return sequence
       
    def setHeader(self, header):
        checkedHeader = self.checkHeader(header)
        self.header = checkedHeader
    
    def getHeader(self):
        return self.header

    def setSequenceType(self, type):
        if type == 'nt':
            self.type = type
        elif type == 'aa':
            self.type = type
    
    def getSequencetype(self):
        return self.type
    
    def setSequence(self, type,sequence):
        if type == 'nt':
            chekedSequence = self.checkSequence(type,sequence)
            if chekedSequence == 'error':
                return print('error')
            else:
                self.sequence = sequence
        elif type == 'aa':
            chekedSequence = self.checkSequence(type,sequence)
            if chekedSequence == 'error':
                return print('error')
            else:
                self.sequence = sequence
    
    def getSequence(self):
        return self.sequence

    def printFasta(self, fasta):
        f_h = fasta.getHeader()
        f_s = fasta.getSequence()

        print(f_h)
        count = 0
        for i in range(0,len(f_s)):
            if count <= 80:
                print(f_s[i], end='')
                count += 1
            else:
                print()
                count = 0
    

h = '>jj'
s = 'acgaaaaaagaggagagagaggagagaagagagagaggagagggagagagagagaggagagagagaggagaggagagagagaggagagagagaggagagagagaggagagagaggaag'
a = 'nt'

fasta = Fasta(a,h,s)
fasta.setHeader(h)
fasta.setSequence(a,s)

fasta.printFasta(fasta)