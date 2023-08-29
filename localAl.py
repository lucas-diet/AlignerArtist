
from fasta import Fasta as Fasta
from blosum62 import blosum62

class SmithWaterman():

    def __init__(self, s1='',s2=''):
        self.s1 = s1
        self.s2 = s2
        self.penalty = {'match':1, 'mismatch':-1, 'gap':-1}

    def setPenalty(self, match,mismatch,gap):
        self.penalty = {'match':match, 'mismatch':mismatch, 'gap':gap}

    def getPenalty(self):
        return self.penalty
    
    def similarities(self, c1,c2):
        if c1 == c2:
            return self.penalty['match']
        elif c1 != c2:
            if c1 != '-' and c2 != '-':
                return self.penalty['mismatch']
            else:
                return self.penalty['gap']
    
    def initDP(self, s1,s2):
        dp = [[0 for _ in range(len(s2)+1)] for _ in range(len(s1)+1)]
    
        for i in range(0,len(s1)+1):
            dp[i][0] = 0
        for j in range(0,len(s2)+1):
            dp[0][j] = 0
        
        return dp
    
    def calcualteDP(self, type,s1,s2):
        dp = self.initDP(s1,s2)
        if type == 'nt':
            for i in range(1,len(s1)+1):
                for j in range(1,len(s2)+1):
                    diag = dp[i-1][j-1] + self.similarities(s1[i-1],s2[j-1])
                    hori = dp[i][j-1] + self.similarities('-',s2[j-1])
                    vert = dp[i-1][j] + self.similarities(s1[i-1],'-')
                    dp[i][j] = max(0,diag,hori,vert)
              
        elif type == 'aa':
            for i in range(1,len(s1)+1):
                for j in range(1,len(s2)+1):
                    diag = dp[i-1][j-1] + blosum62[s1[i-1]][s2[j-1]]
                    hori = dp[i][j-1] + self.similarities('-',s2[j-1])
                    vert = dp[i-1][j] + self.similarities(s1[i-1],'-')
                    dp[i][j] = max(0,diag,hori,vert)
 
        return dp

    def getMaximalSimilarities(self, dp_matrix):
        score = row = col = 0
        tmp_max = [[score, row, col]]

        for i in range(0,len(dp_matrix)):
            for j in range(0,len(dp_matrix[0])):
                if dp_matrix[i][j] >= score:
                    score = dp_matrix[i][j]
                    row = i
                    col = j
                    tmp_max.append([score,row,col])
        
        max_similarities = []
        for i in range(0,len(tmp_max)):
            if score == tmp_max[i][0]:
                max_similarities.append(tmp_max[i])
                             
        return max_similarities
    
    def trackbackLocalAlignment(self, dp_mat, type, s1, s2):
        maxScore = 0
        maxIndices = []

        for i in range(0, len(s1)+1):
            for j in range(0, len(s2)+1):
                if dp_mat[i][j] > maxScore:
                    maxScore = dp_mat[i][j]
                    maxIndices = [(i,j)]
                elif dp_mat[i][j] == maxScore:
                    maxIndices.append((i,j))

        alignments = []
        for indices in maxIndices:
            i,j = indices
            al1, al2 = '',''
            if type == 'nt':
                while i > 0 and j > 0 and dp_mat[i][j] > 0:
                    if dp_mat[i][j] == dp_mat[i-1][j] + self.similarities(s1[i-1], '-'):
                        al1 = s1[i-1] + al1
                        al2 = '-' + al2
                        i -= 1

                    elif dp_mat[i][j] == dp_mat[i][j-1] + self.similarities('-', s2[j-1]):
                        al1 = '-' + al1
                        al2 = s2[j-1] + al2
                        j -= 1

                    elif dp_mat[i][j] == dp_mat[i-1][j-1] + self.similarities(s1[i-1],s2[j-1]):
                        al1 = s1[i-1] + al1
                        al2 = s2[j-1] + al2
                        i -= 1
                        j -= 1

                alignments.append((al1, al2))

            elif type == 'aa':
                while i > 0 and j > 0 and dp_mat[i][j] > 0:
                    if dp_mat[i][j] == dp_mat[i-1][j] + self.similarities(s1[i-1], '-'):
                        al1 = s1[i-1] + al1
                        al2 = '-' + al2
                        i -= 1

                    elif dp_mat[i][j] == dp_mat[i][j-1] + self.similarities('-', s2[j-1]):
                        al1 = '-' + al1
                        al2 = s2[j-1] + al2
                        j -= 1

                    elif dp_mat[i][j] == dp_mat[i-1][j-1] + blosum62[s1[i-1]][s2[j-1]]:
                        al1 = s1[i-1] + al1
                        al2 = s2[j-1] + al2
                        i -= 1
                        j -= 1
                alignments.append((al1, al2))
                
        return alignments 
            
    def printAlignmnts(self, als):
        for al in als:
            print()
            for s in al:
                print(s)


'''
s1 = 'agtagt'
t = 'nt'
h1 = 'h1'

f1 = Fasta()
f1.setHeader(h1)
f1.setSequenceType(t)
f1.setSequence(t,s1)
f1.printFasta(f1)


s2 = 'agt'
h2 = 'h2'

f2 = Fasta()
f2.setHeader(h2)
f2.setSequenceType(t)
f2.setSequence(t,s2)
f1.printFasta(f2)

sw = SmithWaterman()
dp = sw.calcualteDP(t,f1.getSequence(),f2.getSequence())

for i in range(0,len(dp)):
    print(dp[i])

idx = sw.getMaximalSimilarities(dp)
#al = sw.buildAlignments(idx)

als = sw.trackbackLocalAlignment(dp, t,s1,s2)

sw.printAlignmnts(als)
'''