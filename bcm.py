class BestCostMatrix():

    def __init__(self, s1='',s2=''):
        self.s1 = s1
        self.s2 = s2
        self.penalty = {'match':0, 'mismatch':1, 'gap':1}

    def setPenalty(self, ma, mi, ga):
        self.penalty = {'match':ma, 'mismatch':mi, 'gap':ga}

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
        dp = [[0 for _ in range(len(s2)+1)] for _ in range(len(s1)+1)]

        for i in range(len(s1)-1, -1, -1):
            dp[i][len(s2)] = dp[i+1][len(s2)] + self.costs(s1[i-1],'-')

        for j in range(len(s2)-1, -1, -1):
            dp[len(s1)][j] = dp[len(s1)][j+1] + self.costs('-',s2[j-1])
        
        return dp
    
    def calculateDPRev(self, s1, s2):
        dp = self.initRevDP(s1, s2)

        for i in range(len(s1)-1, -1, -1):
            for j in range(len(s2)-1, -1, -1):
                f1 = dp[i+1][j+1] + self.costs(s1[i], s2[j])
                f2 = dp[i+1][j] + self.costs(s1[i], '-')
                f3 = dp[i][j+1] + self.costs('-', s2[j])
                dp[i][j] = min(f1, f2, f3)
        return dp
    
    def calculateM(self, dp, dpRev):
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


s1 = 'AGATC'
s2 = 'TACATA'

bcm = BestCostMatrix()

d = bcm.calcualteDP(s1,s2)
drev = bcm.calculateDPRev(s1,s2)
m = bcm.calculateM(d,drev)

bcm.printMatrix('D:', d)
print()
bcm.printMatrix('D_rev:', drev)
print()
bcm.printMatrix('M:', m)
    