
class MultipleSequenzalignment():

    def __init__(self, s1='', s2='', s3=''):
        self.s1 = s1
        self.s2 = s2
        self.s3 = s3
        self.penalty = {'match':0, 'mismatch':1, 'gap':1}

    def setPenalty(self, match,mismatch,gap):
        self.penalty = {'match':match, 'mismatch':mismatch, 'gap':gap}

    def getPenalty(self):
        return self.penalty
    
    def sumOfPair(self, c1, c2, c3):
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
        dp = [[[0 for _ in range(len(s3)+1)] for _ in range(len(s2)+1)] for _ in range(len(s1)+1)]
        for i in range(0,len(s1)+1):
            dp[i][0][0] = self.sumOfPair(s1[i-1], '-', '-') * i
        for j in range(0,len(s2)+1):
            dp[0][j][0] = self.sumOfPair('-', s2[j-1], '-') * j
        for k in range(1,len(s3)+1):
            dp[0][0][k] = self.sumOfPair('-', '-', s3[k-1]) * k
        
        for i in range(1, len(s1)+1):
            for j in range(1, len(s2)+1):
                f1 = dp[0][j-1][k-1] + self.sumOfPair(s1[i-1], '-', '-')
                f2 = dp[0][j-1][k] + self.sumOfPair('-', s2[j-1], '-')
                f3 = dp[0][j][k-1] + self.sumOfPair('-', '-', s3[k-1])
                dp[0][j][k] = min (f1,f2,f3)
        
            for k in range(1, len(s3)+1):
                f1 = dp[i-1][0][k-1] + self.sumOfPair(s1[i-1], '-', s3[k-1])
                f2 = dp[i-1][0][k] + self.sumOfPair(s1[i-1], '-', '-')
                f3 = dp[i][0][k-1] + self.sumOfPair('-', '-', s3[k-1]) 
                dp[i][0][k] = min (f1,f2,f3) 
        
        for j in range(1, len(s2)+1):
             for k in range(1, len(s3)+1):
                f1 = dp[i-1][j-1][0] + self.sumOfPair('-', s2[j-1], s3[k-1])
                f2 = dp[i-1][j][0] + self.sumOfPair(s1[i-1], '-', '-')
                f3 = dp[i][j-1][0] + self.sumOfPair('-', s2[j-1], '-')
                dp[i][j][0] = min(f1,f2,f3)
        
        return dp
    
    def calcualteDP(self, s1, s2, s3):
        dp = self.initDP(s1, s2, s3)
        for i in range(1, len(s1)+1):
             for j in range(1, len(s2)+1):
                    for k in range(1, len(s3)+1):
                        f1 = dp[i-1][j-1][k-1] + self.sumOfPair(s1[i-1], s2[j-1], s3[k-1])
                        f2 = dp[i-1][j-1][k] + self.sumOfPair(s1[i-1], s2[j-1], '-')
                        f3 = dp[i-1][j][k-1] + self.sumOfPair(s1[i-1], '-', s2[k-1])
                        f4 = dp[i][j-1][k-1] + self.sumOfPair('-', s2[j-1], s3[k-1])
                        f5 = dp[i][j][k-1] + self.sumOfPair('-', '-', s3[k-1])
                        f6 = dp[i][j-1][k] + self.sumOfPair('-', s2[j-1], '-')
                        f7 = dp[i-1][j][k] + self.sumOfPair(s1[i-1], '-', '-')
                        
                        dp[i][j][k] = min(f1,f2,f3,f4,f5,f6,f7)
                      
        return dp
s1 = 'TACA'
s2 = 'CTAC'
s3 = 'GTAG'
msa = MultipleSequenzalignment()

dp = msa.initDP(s1, s2, s3)
for i in range(0,len(dp)):
    print(dp[i])

print()

dp = msa.calcualteDP(s1, s2, s3)

for i in range(0,len(dp)):
    print(dp[i])