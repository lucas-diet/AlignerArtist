from fasta import Fasta as Fasta
from blosum62 import blosum62

class MultipleSequenzalignment():

    def __init__(self, s1='', s2='', s3=''):
        self.s1 = s1
        self.s2 = s2
        self.s3 = s3
    
    def sumOfPair(self, c1, c2, c3):
        """_summary_
            Berehcnet den Sum of Pairs Score für drei Buchsten und gibt diesen zurück.
        Args:
            c1 (char): Buchstabe aus Sequenz 1
            c2 (char): Buchstabe aus Sequenz 2
            c3 (char): Buchstabe aus Sequenz 3

        Returns:
            int: Sum of Pairs Score
        """
        sc = 0
        
        if c1 != c2:
            sc += 1
        if c1 != c3:
            sc += 1
        if c2 != c3:
            sc += 1
        else:
            sc += 0
        
        return sc

    def initDP(self, s1, s2, s3):
        """_summary_
            Es werden drei Sequenzen übergeben und anhand der länge die 3D-Matrix initialisiert.
        Args:
            s1 (str): Sequenz 1
            s2 (str): Sequenz 2
            s3 (str): Sequenz 3

        Returns:
            _type_: Matrix mit den Initalisierten der Nullten Spalten und Zeilen.
        """
        dp = [[[0 for _ in range(len(s3)+1)] for _ in range(len(s2)+1)] for _ in range(len(s1)+1)]
        for i in range(0,len(s1)+1):
            dp[i][0][0] = self.sumOfPair(s1[i-1], '-', '-') * i
        for j in range(0,len(s2)+1):
            dp[0][j][0] = self.sumOfPair('-', s2[j-1], '-') * j
        for k in range(1,len(s3)+1):
            dp[0][0][k] = self.sumOfPair('-', '-', s3[k-1]) * k
        
        for i in range(1,len(s1)+1):
            for j in range(1,len(s2)+1):
                f1 = dp[i-1][j-1][0] + self.sumOfPair(s1[i-1],s2[j-1],'-')
                f2 = dp[i-1][j][0] + self.sumOfPair(s1[i-1],'-','-')
                f3 =  dp[i][j-1][0] + self.sumOfPair('-',s2[j-1],'-')
                
                dp[i][j][0] = min(f1,f2,f3)

        for i in range(1,len(s1)+1):
            for k in range(1,len(s3)+1):
                f1 = dp[i-1][0][k-1] + self.sumOfPair(s1[i-1],'-',s3[k-1])
                f2 = dp[i-1][0][k] + self.sumOfPair(s1[i-1],'-','-')
                f2 = dp[i][0][k-1] + self.sumOfPair('-','-',s3[k-1])
                
                dp[i][0][k] = min(f1,f2,f3)

        for j in range(1,len(s2)+1):
            for k in range(1,len(s3)+1):
                f1 = dp[0][j-1][k-1] + self.sumOfPair('-',s2[j-1],s3[k-1])
                f2 = dp[0][j-1][k] + self.sumOfPair('-',s2[j-1],'-')
                f3 = dp[0][j][k-1] + self.sumOfPair('-','-',s3[k-1])
                
                dp[0][j][k] = min(f1,f2,f3)

        return dp
    
    def calcualteDP(self, s1, s2, s3):
        """_summary_
            Um den Wert für eine Zelle zu berechnen, werden durch 3 for-loops durch die Sequenzen gelaufen und dann 
            für jede Zelle der miniamle Score berechnt.
        Args:
            s1 (str): Sequenz 1
            s2 (str): Sequenz 2
            s3 (str): Sequenz 3

        Returns:
            _type_: 3D-Matrix
        """
        dp = self.initDP(s1, s2, s3)
        for i in range(1, len(s1)+1):
             for j in range(1, len(s2)+1):
                    for k in range(1, len(s3)+1):
                        f1 = dp[i-1][j-1][k-1] + self.sumOfPair(s1[i-1], s2[j-1], s3[k-1])
                        f2 = dp[i-1][j-1][k] + self.sumOfPair(s1[i-1], s2[j-1], '-')
                        f3 = dp[i-1][j][k-1] + self.sumOfPair(s1[i-1], '-', s3[k-1])
                        f4 = dp[i][j-1][k-1] + self.sumOfPair('-', s2[j-1], s3[k-1])
                        f5 = dp[i][j][k-1] + self.sumOfPair('-', '-', s3[k-1])
                        f6 = dp[i][j-1][k] + self.sumOfPair('-', s2[j-1], '-')
                        f7 = dp[i-1][j][k] + self.sumOfPair(s1[i-1], '-', '-')
                        
                        dp[i][j][k] = min(f1,f2,f3,f4,f5,f6,f7)
                      
        return dp
    
    def getMinimalCosts(self, dp_mat):
        return dp_mat[-1][-1][-1]

    def trackbackMSA(self, dp_mat, s1, s2, s3, i, j, k, al1='', al2='', al3='', alignments=[]):
        """_summary_

        Args:
            dp_mat (_type_): 3D-Matrix
            s1 (str): Sequenz 1
            s2 (str): Sequenz 2
            s3 (str): Sequenz 3
            i (int): Laufvariable für Sequenz 1
            j (int): Laufvariable für Sequenz 2
            k (int): Laufvariable für Sequenz 3
            al1 (str, optional): Alignments-Sequenz für Sequenz 1. Defaults to ''.
            al2 (str, optional): Alignment-Sequnez für Sequenz 2. Defaults to ''.
            al3 (str, optional): Alignment-Sequenz für Sequenz 3. Defaults to ''.
            alignments (list, optional): _description_. Defaults to [].

        Returns:
            _type_: _description_
        """
        if i == 0 and j == 0 and k == 0:
            alignments.append((al1, al2, al3))
        
        if i > 0 and j > 0 and dp_mat[i][j][k] == dp_mat[i-1][j-1][k] + self.sumOfPair(s1[i-1],s2[j-1],'-'):
            self.trackbackMSA(dp_mat, s1, s2, s3, i-1, j-1, k, s1[i-1]+al1, s2[j-1], '-'+al3, alignments)
        
        if i > 0 and k > 0 and dp_mat[i][j][k] == dp_mat[i-1][j][k-1] + self.sumOfPair(s1[i-1],'-',s3[k-1]):
            self.trackbackMSA(dp_mat, s1, s2, s3, i-1, j, k-1, s1[i-1]+al1, '-'+al3, s3[k-1]+al3, alignments)
        
        if j > 0 and k > 0 and dp_mat[i][j][k] == dp_mat[i][j-1][k-1] + self.sumOfPair('-',s2[j-1],s3[k-1]):
            self.trackbackMSA(dp_mat, s1, s2, s3, i, j-1, k-1, '-'+al1, s2[j-1]+al2, s3[k-1]+al3, alignments)
        
        if i > 0 and dp_mat[i][j][k] == dp_mat[i-1][j][k] + self.sumOfPair(s1[i-1],'-','-'):
            self.trackbackMSA(dp_mat, s1, s2, s3, i-1, j, k, s1[i-1]+al1, '-'+al2, '-'+al3, alignments)
        
        if j > 0 and dp_mat[i][j][k] == dp_mat[i][j-1][k] + self.sumOfPair('-',s2[j-1],'-'):
            self.trackbackMSA(dp_mat, s1, s2, s3, i, j-1, k, '-'+al1, s2[j-1]+al2, '-'+al3, alignments)

        if k > 0 and dp_mat[i][j][k] == dp_mat[i][j][k-1] + self.sumOfPair('-','-',s3[k-1]):
            self.trackbackMSA(dp_mat, s1, s2, s3, i, j, k-1, '-'+al1, '-'+al2, s3[j-1]+al3, alignments)

        if i > 0 and j > 0 and k > 0 and dp_mat[i][j][k] == dp_mat[i-1][j-1][k-1] + self.sumOfPair(s1[i-1],s2[j-1],s3[k-1]):
            self.trackbackMSA(dp_mat,  s1, s2, s3, i-1, j-1, k-1, s1[i-1]+al1, s2[j-1]+al2, s3[k-1]+al3, alignments)

        return alignments

    def printMSA(self, als):
        for al in als:
            for s in al:
                print(s)
            print()

'''
s1 = 'taca'
s2 = 'ctac'
s3 = 'gtag'
msa = MultipleSequenzalignment()

print()

dp = msa.calcualteDP(s1, s2, s3)

c = msa.getMinimalCosts(dp)
print(c)

als = msa.trackbackMSA(dp, s1, s2, s3, len(s1), len(s2), len(s3))

msa.printMSA(als)
'''