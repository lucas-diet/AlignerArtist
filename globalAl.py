
from fasta import Fasta as Fasta
from blosum62 import blosum62

class NeedlemannWunsch():

    def __init__(self, s1='', s2=''):
        self.s1 = s1
        self.s2 = s2
        self.penalty = {'match':0, 'mismatch':1, 'gap':1}

    def setPenalty(self, match,mismatch,gap):
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
        for i in range(0,len(s1)+1):
            dp[i][0] = self.penalty['gap'] * i
        for j in range(0,len(s2)+1):
            dp[0][j] = self.penalty['gap'] * j
        
        return dp
    
    def calcualteDP(self, type, s1, s2):
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
        return dp_mat[-1][-1]
    
    def trackbackGlobalAlignments(self, dp_mat, type, s1, s2, i, j, al1='', al2='', alignments=[]):
        
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