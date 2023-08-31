
class Fasta():

    def __init__(self, type='',header='',sequence=''):
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
                    return 'ERROR' #print('ERROR')
                    break
                else:
                   self.sequence = sequence
        elif type == 'aa':
            for aa in sequence:
                if aa not in __aas:
                    return 'ERROR' #print('ERROR')
                    break
                else:
                   self.sequence = sequence

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
    
    def getSequenceType(self):
        return self.type
    
    def setSequence(self, type,sequence):
        if type == 'nt':
            chekedSequence = self.checkSequence(type,sequence)
            if chekedSequence == 'ERROR':
                #return 'ERROR'
                self.sequence = 'ERROR'
            elif self.sequence == sequence:
                self.sequence = sequence
        elif type == 'aa':
            chekedSequence = self.checkSequence(type,sequence)
            if chekedSequence == 'ERROR':
                #return 
                self.sequence = 'ERROR'
            elif self.sequence == sequence:
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
        print()
        print()

'''
s1 = 'acgt'
t = 'nt'
h1 = 'h1'

f1 = Fasta()
f1.setHeader(h1)
f1.setSequenceType(t)
f1.setSequence(t,s1)
f1.printFasta(f1)


s2 = 'acggg'
h2 = 'h1'

f2 = Fasta()
f2.setHeader(h2)
f2.setSequenceType(t)
f2.setSequence(t,s2)
f1.printFasta(f2)
'''