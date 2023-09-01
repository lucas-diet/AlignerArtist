
class Fasta():

    def __init__(self, type='',header='',sequence=''):
        self.type = type
        self.header = header
        self.sequence = sequence

    def checkHeader(self, header):
        """_summary_
            Nimmt einen String (header) und überprüft, ob erstes Zeichen ein '>' ist.
            Wenn nicht dann wird der String so angepasst, dass die Bedingung für ein
            Header erfüllt ist.
        Args:
            header (String): Eine Art Name für eine Sequenz

        Returns:
            String: Den Header
        """
        if header[0] != '>':
            correctHeader = '>'+ header
            return correctHeader
        else:
            return header

    def checkSequence(self, type,sequence):
        """_summary_
            Überprüft für die zwei Arten von Sequenzen (Nukleotidsequenz (nt), Aminosäuresequenz (aa)), 
            ob es keine verbotenen Buchstben in der Sequenz gibt. Falls es verbotene Buchstaben gibt, dann
            wird 'ERROR' zurückgegeben.
        Args:
            type (String):Beschreibt die Art der Sequenz, ob es sich um eine Nukleotidsequenz 'nt' oder eine
            Aminosäuresequenz 'aa' handelt.
            sequence (String): _description_
        Returns:
            String: Error-String
        """
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
        """_summary_
            Legt den Header fest. Bevor der Header aber festgelegt wird, wird dieser auf Korrektheit überprüft.
        Args:
            header (String): Name einer Sequenz
        """
        checkedHeader = self.checkHeader(header)
        self.header = checkedHeader
    
    def getHeader(self):
        return self.header

    def setSequenceType(self, type):
        """_summary_
            Soll die Art der Sequenz festlegen.
        Args:
            type (String): Beschreibt die Art der Sequenz (Nukleotide oder Aminosäuren)
        """
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