
from globalAl import NeedlemannWunsch as NW

class CenterStar():

    def __init__(self, type='', s1='', s2='', s3='', s4=''):
        self.s1 = s1
        self.s2 = s2
        self.s3 = s3
        self.s4 = s4
        self.type = type

    def getCenterSequence(self, type, s1, s2, s3, s4):
        nw = NW()
        cost12 = nw.calcualteDP(type, s1,s2)[-1][-1]
        cost13 = nw.calcualteDP(type, s1,s3)[-1][-1]
        cost14 = nw.calcualteDP(type, s1,s4)[-1][-1]

        cost23 = nw.calcualteDP(type, s2,s3)[-1][-1]
        cost24 = nw.calcualteDP(type, s2,s4)[-1][-1]
        
        cost34 = nw.calcualteDP(type, s3,s4)[-1][-1]
        
        return min(cost12, cost13, cost14, cost23, cost24, cost34)
    
s1 = 'TACA'
s2 = 'CTAC'
s3 = 'GTAG'
s4 = 'ATGC'

dca = CenterStar()

cs = dca.getCenterSequence('nt', s1,s2,s3,s4)

print(cs)

