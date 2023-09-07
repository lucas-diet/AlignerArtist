from bcm import BestCostMatrix as BCM

class CarilloLipmanBarrier():

    def __init__(self, alScore = 0, s1='', s2=''):
        self.s1 = s1
        self.s2 = s2
        self.alScore = alScore

    def setAlScore(self, alScore):
        self.alScore = alScore

    def getAlScore(self):
       return self.alScore
    
    def printSmallMatrix(self, m, u):
        for i in range(0,len(m)):
            for j in range(0,len(m[0])):
                if m[i][j] <= u:
                    print(m[i][j], end=' ')
                else:
                    print('  ', end='')
            print()

'''
s1 = 'AGATC'
s2 = 'TACATA'
s3 = 'GAGAT'

bcm = BCM()

dp12 = bcm.calcualteDP(s1,s2)
dp13 = bcm.calcualteDP(s1,s3)
dp23 = bcm.calcualteDP(s2,s3)

dprev12 = bcm.calculateDPRev(s1,s2)
dprev13 = bcm.calculateDPRev(s1,s3)
dprev23 = bcm.calculateDPRev(s2,s3)

m12 = bcm.calculateM(dp12, dprev12)
m13 = bcm.calculateM(dp13, dprev13)
m23 = bcm.calculateM(dp23, dprev23)

print()
clb = CarilloLipmanBarrier()
clb.setAlScore(10)

alScore = clb.getAlScore()

u12 = alScore - (dp13[-1][-1] + dp23[-1][-1])
u13 = alScore - (dp12[-1][-1] + dp23[-1][-1])
u23 = alScore - (dp12[-1][-1] + dp13[-1][-1])

clb.printSmallMatrix(m12, u12)
clb.printSmallMatrix(m13, u13)
clb.printSmallMatrix(m23, u23)

'''