from fasta import Fasta as fasta
from blosum62 import blosum62
from collections import deque

class NeedlemannWunsch():

    def __init__(self, s1='', s2=''):
        self.s1 = s1
        self.s2 = s2
        self.penalty = {'match':0, 'mismatch':1, 'gap':1}

    def setPenalty(self, match,mismatch,gap):
        self.penalty = {'match':match, 'mismatch':mismatch, 'gap':gap}

    def getPenalty(self):
        return self.penalty
    
    def listOfCosts(self, c1,c2):
        if c1 != c2:
            if c1 == '-' or c2 == '-':
                return self.penalty['gap']
            elif c1 != '-' or c2 != '-':
                return self.penalty['mismatch']
        else:
            return self.penalty['match']
                
    def initDP(self, s1,s2):
        dp = [[0 for _ in range(len(s2)+1)] for _ in range(len(s1)+1)]
        for i in range(0,len(s1)+1):
            dp[i][0] = self.penalty['gap'] * i
        for j in range(0,len(s2)+1):
            dp[0][j] = self.penalty['gap'] * j
        
        return dp
    
    def calcualteDP(self, type,s1,s2):
        dp = self.initDP(s1,s2)
        if type == 'nt':
            for i in range(1,len(s1)+1):
                for j in range(1,len(s2)+1):
                    diag = dp[i-1][j-1] + self.listOfCosts(s1[i-1],s2[j-1])
                    hori = dp[i][j-1] + self.listOfCosts('-',s2[j-1])
                    vert = dp[i-1][j] + self.listOfCosts(s1[i-1],'-')
                    dp[i][j] = min(diag,hori,vert)

        elif type == 'aa':
            for i in range(1,len(s1)+1):
                for j in range(1,len(s2)+1):
                    diag = dp[i-1][j-1] + blosum62[s1[i-1]][s2[j-1]]
                    hori = dp[i][j-1] + self.listOfCosts('-',s2[j-1])
                    vert = dp[i-1][j] + self.listOfCosts(s1[i-1],'-')
                    dp[i][j] = min(diag,hori,vert)
            
        return dp
    
    def getMinimallistOfCosts(self, dp_matrix):
        return dp_matrix[-1][-1]
    
    def allAlignments(self, s1,s2):

        def helper(s1,s2):
            if len(s1) == 0 and len(s2) == 0:
                yield deque()

            scenarios = []
            count = 0
            if len(s1) > 0 and len(s2) > 0:
                scenarios.append((s1[0],s1[1:],s2[0],s2[1:]))
            if len(s1) > 0:
                scenarios.append((s1[0],s1[1:],None,s2))
            if len(s2) > 0:
                scenarios.append((None,s1,s2[0],s2[1:]))

            for s1h,s1t,s2h,s2t in scenarios:
                for alignment in helper(s1t,s2t):
                    alignment.appendleft((s1h,s2h))
                    yield alignment
		
        alignments = helper(range(len(s1)),range(len(s2)))
        return map(list,alignments)
    
    def buildAlignments(self,sequence1,sequence2,alignment):
        al = []
        nt1 = ['-' if i is None else sequence1[i] for i, _ in alignment]
        nt2 = ['-' if j is None else sequence2[j] for _, j in alignment]
        al.append([nt1,nt2])
        alignments = [val for sublist in al for val in sublist]
        return alignments

    def optAlignments(self,min_cost,listOfAlignments):
        '''
		Durchläuft die Liste aller möglichen Alignments und berechnet die Kosten der Alignmets.
		Dann wird eine Liste erzeugt mit Kosten und den zugehörigen Sequenzen.
		Filtere die optimalen Alignments -> Return.
		'''
        alignment = listOfAlignments
        cost = 0
        listOfCosts = []
        optAls = []
        for nt in range(1,len(alignment),2):
            s1 = alignment[nt-1]
            s2 = alignment[nt]
            for i in range(0,len(s1)):
                for j in range(i,len(s2)):
                    if s1[i] == s2[j]:
                        cost += self.penalty['match']
                        break
                    elif s1[i] != s2[j]:
                        if s1[i] == '-':
                            cost += self.penalty['gap']
                            break
                        elif s2[j] == '-':
                            cost += self.penalty['gap']
                            break
                        else:
                            cost += self.penalty['mismatch']
                            break
        listOfCosts.append([cost,s1,s2])
        listOfCosts_al = [val for sublist in listOfCosts for val in sublist]
        for i in range(0,len(listOfCosts_al)):
            if listOfCosts_al[i] == min_cost:
                optAls.append(listOfCosts_al[i+1])
                optAls.append(listOfCosts_al[i+2])
        
        return optAls
    
    def printAlignments(self,opts):
        if len(opts) > 0:
            for seq in opts:
                print()
                for nt in seq:
                    print(nt,end='')
            print()
                
s1 = 'agtt'
s2 = 'agta'

nw = NeedlemannWunsch()
dp = nw.calcualteDP('nt',s1,s2)

for i in range(0,len(dp)):
    print(dp[i])

listOfCosts = nw.getMinimallistOfCosts(dp)
print(listOfCosts)

als = nw.allAlignments(s1,s2)

for al in als:
    alignments = nw.buildAlignments(s1,s2,al)
    opt = nw.optAlignments(listOfCosts,alignments)
    nw.printAlignments(opt)