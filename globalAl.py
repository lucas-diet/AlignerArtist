from fasta import Fasta

class NeedlemannWunsch():

    def __init__(self, s1='', s2=''):
        self.s1 = s1
        self.s2 = s2
        self.penalty = {'match':0, 'mismatch':1, 'gap':1}

    def setPenalty(self, match,mismatch,gap):
        self.penalty = {'match':match, 'mismatch':mismatch, 'gap':gap}

    def getPenalty(self):
        return self.penalty
    
    def costs(self, c1,c2):
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
    
    def calcualteDP(self, s1,s2):
        dp = self.initDP(s1,s2)

        for i in range(1,len(s1)+1):
            for j in range(1,len(s2)+1):
                diag = dp[i-1][j-1] + self.costs(s1[i-1],s2[j-1])
                hori = dp[i][j-1] + self.costs('-',s2[j-1])
                vert = dp[i-1][j] + self.costs(s1[i-1],'-')
                dp[i][j] = min(diag,hori,vert)
        return dp
    

s1 = 'agt'
s2 = 'agta'

dp = NeedlemannWunsch().calcualteDP(s1,s2)

for i in range(0,len(dp)):
    print(dp[i])

    