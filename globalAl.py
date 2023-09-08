
from fasta import Fasta as Fasta
from blosum62 import blosum62

class NeedlemannWunsch():

    def __init__(self, s1='', s2=''):
        self.s1 = s1
        self.s2 = s2
        self.penalty = {'match':0, 'mismatch':1, 'gap':1}

    def setPenalty(self, match, mismatch, gap):
        """_summary_
            Die Kosten für das Einfügen eines Gaps, Match oder eines Mismatches werden hier festgelegt.
        Args:
            match (int): Beschreibt die Kosten bei einem Match
            mismatch (int): Beschreibt die Kosten bei einem Mismatch
            gap (int): Beschreibt die Kosten bei einem Gap
        """
        self.penalty = {'match':match, 'mismatch':mismatch, 'gap':gap}

    def getPenalty(self):
        return self.penalty
    
    def costs(self, c1, c2):
        """_summary_
            Hier werden zwei Buchstaben verglichen, um damit die entsprechenden Kosten zu bestimmen.
        Args:
            c1 (char): Ist ein gewählter Buchstabe
            c2 (char): Ist eine gewählter Buchstabe

        Returns:
            _type_: Die entsprechdnen Kosten, je nach dem welcher Wert des Dict zurückgegeben wird.
        """
        if c1 != c2:
            if c1 == '-' or c2 == '-':
                return self.penalty['gap']
            elif c1 != '-' or c2 != '-':
                return self.penalty['mismatch']
        else:
            return self.penalty['match']
                
    def initDP(self, s1, s2):
        """_summary_
            Es werden zwei Sequenzen übergeben, mit denen dann eine Matrix erzeugt wird.
            Die Nullte Zeile bzw. Spalte wird dann entsprechend der Kostenfunktion initialisiert.
        Args:
            s1 (str): Erste Sequenz
            s2 (str): Zweite Sequenz

        Returns:
            _type_: Eine Matrix
        """
        dp = [[0 for _ in range(len(s2)+1)] for _ in range(len(s1)+1)]
        for i in range(0,len(s1)+1):
            dp[i][0] = self.penalty['gap'] * i
        for j in range(0,len(s2)+1):
            dp[0][j] = self.penalty['gap'] * j
        
        return dp
    
    def calcualteDP(self, type, s1, s2):
        """_summary_
            Hier wird die initialisierte Matrix dann gefüllt, in dem das Prinzip der dynamischen Programmierung 
            angewendet wird. Bzw. die Rekurrenz für den Needlemann-Wunsch-Algorithmus ausgeführt wird.
            Der gewählte Wert für eine Zelle der Matrix ist dann immer das Minimum aus den drei Fällen.
        Args:
            type (str): Art der Sequent (nt/aa)
            s1 (str): Sequenz 1
            s2 (str): Sequenz 2

        Returns:
            _type_: Eine gefüllte Matrix
        """
        dp = self.initDP(s1,s2)
        if type == 'nt':
            for i in range(1,len(s1)+1):
                for j in range(1,len(s2)+1):
                    diag = dp[i-1][j-1] + self.costs(s1[i-1],s2[j-1])
                    hori = dp[i][j-1] + self.costs('-',s2[j-1])
                    vert = dp[i-1][j] + self.costs(s1[i-1],'-')
                    dp[i][j] = min(diag,hori,vert)

        elif type == 'aa':
            for i in range(1,len(s1)+1):
                for j in range(1,len(s2)+1):
                    diag = dp[i-1][j-1] + blosum62[s1[i-1]][s2[j-1]]
                    hori = dp[i][j-1] + self.costs('-',s2[j-1])
                    vert = dp[i-1][j] + self.costs(s1[i-1],'-')
                    dp[i][j] = min(diag,hori,vert)
            
        return dp
    
    def getMinimalCosts(self, dp_mat):
        """_summary_
            Hier soll der letzte Wert einer Matrix aisgegebn werden. Welcher dann den Kosten zweier Sequenzen entspricht.
        Args:
            dp_mat (_type_): Matrix

        Returns:
            int: gibt den letzten Wert einer übergeben Matrix zurück 
        """
        return dp_mat[-1][-1]
    
    def trackbackGlobalAlignments(self, dp_mat, type, s1, s2, i, j, al1='', al2='', alignments=[]):
        """_summary_
            Es wird rekursiv durch die Matrix gegangen und geschaut ob Bedingungen der Rukkrenz werfüllt sind.
            Wenn ja dann wird die Funktion erneut aufgerufen und das bisher aufgebaute Alignment um eine Positioen 
            erweitert.
            Das geschieht so lange, bis die Laufvariablen gleich Null sind, wo die ermittelten Alignments dann in die 
            Liste aller ALignments aufgenommen werden.
        Args:
            dp_mat (_type_): Eine gefüllte Matrix
            type (str): Beschreibt die Art der Sequenz
            s1 (str): Sequenz 1
            s2 (str): Sequenz 2
            i (int): Eine Laufvariable 
            j (int): Eine Laufvariable
            al1 (str, optional): _description_. Defaults to ''.
            al2 (str, optional): _description_. Defaults to ''.
            alignments (list, optional): _description_. Defaults to [].

        Returns:
            list: Beinhaltet alle möglichen optimalen Alignments
        """
        if i == 0 and j == 0:
            alignments.append((al2, al1))

        if type == 'nt':
            if i > 0 and dp_mat[i][j] == dp_mat[i-1][j] + self.costs(s1[i-1],'-'):
                self.trackbackGlobalAlignments(dp_mat, type, s1, s2, i-1, j, s1[i-1]+al1, '-'+al2, alignments)
            
            if j > 0 and dp_mat[i][j] == dp_mat[i][j-1] + self.costs('-',s2[j-1]):
                self.trackbackGlobalAlignments(dp_mat, type, s1, s2, i, j-1, '-'+al1, s2[j-1]+al2, alignments)
            
            if i > 0 and j > 0 and dp_mat[i][j] == dp_mat[i-1][j-1] + self.costs(s1[i-1],s2[j-1]):
                    self.trackbackGlobalAlignments(dp_mat, type, s1, s2, i-1, j-1, s1[i-1]+al1, s2[j-1]+al2, alignments)
        
        elif type == 'aa':
            if i > 0 and dp_mat[i][j] == dp_mat[i-1][j] + self.costs(s1[i-1],'-'):
                self.trackbackGlobalAlignments(dp_mat, type, s1, s2, i-1, j, s1[i-1]+al1, '-'+al2, alignments)
            
            if j > 0 and dp_mat[i][j] == dp_mat[i][j-1] + self.costs('-',s2[j-1]):
                self.trackbackGlobalAlignments(dp_mat, type, s1, s2, i, j-1, '-' + al1, s2[j-1]+al2, alignments)

            if i > 0 and j > 0 and dp_mat[i][j] == dp_mat[i-1][j-1] + blosum62[s1[i-1]][s2[j-1]]:
                self.trackbackGlobalAlignments(dp_mat, type, s1, s2, i-1, j-1, s1[i-1]+al1, s2[j-1]+al2, alignments)

        return alignments
    
    def printGlobalAlignments(self,als):
        for al in als:
            print()
            for s in al:
                print(s)
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


s2 = 'at'
h2 = 'h2'

f2 = Fasta()
f2.setHeader(h2)
f2.setSequenceType(t)
f2.setSequence(t,s2)
f1.printFasta(f2)

nw = NeedlemannWunsch()
dp = nw.calcualteDP(t, f1.getSequence(),f2.getSequence())


nw = NeedlemannWunsch()
dp = nw.calcualteDP(t,f1.getSequence(),f2.getSequence())

for i in range(0,len(dp)):
   print(dp[i])
print()

c = nw.getMinimalCosts(dp)
print(c)

als = nw.trackbackGlobalAlignments(dp,t,s1,s2,len(s1),len(s2))

nw.printGlobalAlignments(als)
'''