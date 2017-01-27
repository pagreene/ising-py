#!/usr/bin/python

import ising

if __name__ == "__main__":
    I = ising.IsingLattice(1,10,10,0.1,1)
    
    m = I.cSettle(100000, trackFuncList = [I.cGetm], nTrack = 1000)[0]
