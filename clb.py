from bcm import BestCostMatrix as BCM

class CarilloLipmanBarrier():

    def __init__(self, dp, dpRev, m, s1='', s2=''):
        self.s1 = s1
        self.s2 = s2
        self.dp = dp
        self.dpRev = dpRev
        self.m = m

    def calculateBarrier(self, alScore, *dps):
       upperBarrier = alScore - (dps)
       return upperBarrier

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

bcm.printMatrix('D12:', dp12)
print()
bcm.printMatrix('DRev12:', dprev12)

print()
bcm.printMatrix('D13:', dp13)
print()
bcm.printMatrix('DRev13:', dprev13)

print()
bcm.printMatrix('D23:', dp23)
print()
bcm.printMatrix('DRev23:', dprev23)
