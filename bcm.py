from globalAl import NeedlemannWunsch as NW

class BestCostMatrix():

    def __init__(self, s1='',s2=''):
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
        if c1 != c2:
            if c1 == '-' or c2 == '-':
                return self.penalty['gap']
            elif c1 != '-' or c2 != '-':
                return self.penalty['mismatch']
        else:
            return self.penalty['match']
    
    def initDP(self, s1, s2):
        dp = [[0 for _ in range(len(s2)+1)] for _ in range(len(s1)+1)]

        for i in range(1,len(s1)+1):
            dp[i][0] = self.costs(s1[i-1],'-') * i
        
        for j in range(1,len(s2)+1):
            dp[0][j] = self.costs('-',s2[j-1]) * j
        
        return dp
    
    def calcualteDP(self, s1, s2):
        dp = self.initDP(s1,s2)
        for i in range(1,len(s1)+1):
            for j in range(1,len(s2)+1):
                diag = dp[i-1][j-1] + self.costs(s1[i-1],s2[j-1])
                hori = dp[i][j-1] + self.costs('-',s2[j-1])
                vert = dp[i-1][j] + self.costs(s1[i-1],'-')
                dp[i][j] = min(diag,hori,vert)
        return dp
    
    def initRevDP(self, s1, s2):
        """_summary_
            Initailaisiert eine Matrix, wobei die letzt Spalte und Zeile werden entsprechend der 
            Kostenfunktion initialisiert.
        Args:
            s1 (str): Sequenz 1
            s2 (str): Sequenz 2
        Returns:
            _type_: Matrix
        """
        dp = [[0 for _ in range(len(s2)+1)] for _ in range(len(s1)+1)]

        for i in range(len(s1)-1, -1, -1):
            dp[i][len(s2)] = dp[i+1][len(s2)] + self.costs(s1[i-1],'-')

        for j in range(len(s2)-1, -1, -1):
            dp[len(s1)][j] = dp[len(s1)][j+1] + self.costs('-',s2[j-1])
        
        return dp
    
    def calculateDPRev(self, s1, s2):
        """_summary_
            Berechnet die Werte der Matrix, wobei nun die Matrix von (m,n) -> (0,0) gefüllt werden.
        Args:
            s1 (str): Sequenz 1
            s2 (str): Sequenz 1

        Returns:
            _type_: Gefüllte Matrix
        """
        dp = self.initRevDP(s1, s2)

        for i in range(len(s1)-1, -1, -1):
            for j in range(len(s2)-1, -1, -1):
                f1 = dp[i+1][j+1] + self.costs(s1[i], s2[j])
                f2 = dp[i+1][j] + self.costs(s1[i], '-')
                f3 = dp[i][j+1] + self.costs('-', s2[j])
                dp[i][j] = min(f1, f2, f3)
        return dp
    
    def calculateM(self, dp, dpRev):
        """_summary_
            Berechnet die Beste-Kosten-Matrix, indem die DP-Matrix mit der DP-Matrix (Reverse) addiert wird
        Args:
            dp (_type_): DP-Matrix der Sequenzen  
            dpRev (_type_): DP-Matrix der reversen Sequenzen (String+'epsilon')

        Returns:
            _type_: Best-Kosten-Matrix
        """
        m = [[0 for _ in range(len(dp[0]))] for _ in range(len(dp))]
        
        for i in range(0,len(dp)):
            for j in range(0,len(dp[0])):
                m[i][j] = dp[i][j] + dpRev[i][j]
        return m
    
    def printMatrix(self, text, mat):
        print(text)
        for row in mat:
            for col in row:
                print(col, end=' ')
            print()

'''
s1 = 'AGATC'
s2 = 'TACATA'

bcm = BestCostMatrix()

bcm.setPenalty(10,1,1)
print(bcm.getPenalty())

d = bcm.calcualteDP(s1, s2)
drev = bcm.calculateDPRev(s1,s2)
m = bcm.calculateM(d,drev)

bcm.printMatrix('D:', d)
print()
bcm.printMatrix('D_rev:', drev)
print()
bcm.printMatrix('M:', m)
'''