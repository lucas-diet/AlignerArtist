class Fasta():

    def __init__(self, type='',header='',sequence=''):
        self.type = type
        self.header = header
        self.sequence = sequence

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

    def setSequence(self, type,sequence):
        if type == 'nt':
            chekedSequence = self.checkSequence(type,sequence)
            if chekedSequence == 'ERROR':
                print('ERROR')
                self.sequence = ''
            elif self.sequence == sequence:
                self.sequence = sequence
        elif type == 'aa':
            chekedSequence = self.checkSequence(type,sequence)
            if chekedSequence == 'ERROR':
                print('ERROR')
                self.sequence = ''
            elif self.sequence == sequence:
                self.sequence = sequence

    def getSequence(self):
        return self.sequence
        
    
'''
s1 = 'acg'
t = 'aa'
h1 = 'h1'

f1 = Fasta()
f1.setSequence(t,s1)
print(f1.getSequence())
'''