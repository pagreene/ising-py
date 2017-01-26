#!/usr/bin/python

import numpy as np
import matplotlib.pyplot  as mplot
import ctypes
import re
import os
from py_queue import Queuer
from datetime import datetime
from copy import copy

def getProcInfo():
    '''
    Read and convert the /proc/meminfo file
    '''
    patt = re.compile("([\w \(\)]+):\s+(\d+) kB")
    f = open("/proc/meminfo", "r")
    fstr = f.read()
    f.close()
    
    l = patt.findall(fstr)
    infoDict = dict()
    for n, (name, val) in enumerate(l):
        infoDict["%s" % (name)] = int(val)/float(l[0][1])
    
    return infoDict

class IsingLattice(object):
    def __init__(self, J, h, L, b, R):
        self.J = J
        self.h = h
        self.L = L
        self.b = b
        self.R = R
        
        self.s = np.ones([L,L], dtype = np.int)
        self.lib = ctypes.cdll.LoadLibrary('./libIsing.so')
        
        self.measureList = []
        return
    
    def __getitem__(self, i):
        return self.s[tuple(i)]
    
    def __setitem__(self, i, val):
        self.s[tuple(i)] = val
        return
    
    def fixBounds(self, j):
        '''
        This function enforces the boundary conditions. Given an index j, it will
        adjust the index such that it remains in the boundary. Currently,
        torroidal cyclic boundary conditions are implemented.
        '''
        if j >= self.L:
            j -= self.L
        elif j < 0:
            j += self.L
        return j
    
    def fixInts(self, cArr):
        '''
        When we get back the c-arrays from c ints, negative values are mapped to
        the end of the integers. This fixes that problem. Any array that was
        passed to c code as a pointer should be cleansed by this function.
        
        Inputs:
        cArr -- the array that was modified by c-code
        
        Returns:
        None -- the array is modified in-place.
        '''
        # This is the max value of the returned integers...there's a bit
        # of a glitch going from c ints to python ints. It maps negatives
        # such that an int i, i<0 becomes N + i = N - |i|. I have to
        # correct that when the lattice is returned.
        N = 4294967296
        
        # Here is the unfortunate correction.
        cArr[cArr>N/2] = cArr[cArr>N/2] - N
        return
    
    def calcEPlusMinus(self, j):
        '''
        Calculates the energy of a spin with regard to the field (h) and it's
        neighbors within range R for the up and down configurations.
        '''
        # I do most of the calculation without the spin, because I can just
        # multiply by +/- at the end.
        
        # Account for the field.
        E_divby_sj = -self.h
        
        # And now the surrounding spins...if J is set to zero, just skip this.
        if self.J is not 0:
            
            # I use two loops to iterate over the surrounding region.
            for dx in range(-self.R,self.R+1):
                jx = self.fixBounds(j[0] + dx)
                
                for dy in range(-self.R,self.R+1):
                    # The spin does not give itself energy.
                    if dx is 0 and dy is 0:
                        continue
                    
                    jy = self.fixBounds(j[1] + dy)
                    
                    E_divby_sj -= self.J*self[(jx,jy)]/float((2*self.R + 1)**2 - 1)
        
        # It is here that I multiply by +/-1 to account for the sign of the
        # spin at index j.
        return np.array([E_divby_sj, -E_divby_sj])
    
    def lowlevel_cSettle(self, N):
        '''
        This function should perform the exact same function as settle,
        however it calls a method written in c, which should make it MUCH
        faster. There is some translation that needs to be done first though.
        '''
        # Pass the stuff to the c code.
        csettle = self.lib.settle
        csettle(ctypes.c_void_p(self.s.ctypes.data),
                ctypes.c_int(self.L),
                ctypes.c_int(self.R),
                ctypes.c_double(self.h),
                ctypes.c_double(self.J),
                ctypes.c_double(self.b),
                ctypes.c_int(np.random.randint(0,100000)),
                ctypes.c_int(N) )
        
        # This is an unfortunate workaround that is needed. There is some
        # type incompatability issue, I think.
        #self.s[abs(self.s) > 1] = -np.sign(self.s[abs(self.s) > 1])
        
        self.fixInts(self.s)
        return
    
    def plotMeasures(self):
        '''
        Plots everything that was measured during the last run.
        '''
        i = 0
        for thing in self.measureList:
            if not isinstance(thing[0], np.ndarray):
                mplot.figure(i)
                mplot.clf()
                mplot.plot(thing, '.')
                i += 1
        return
    
    def lookForDroplet(self, nBatch, nTrackPerBatch, nIterPerBatch):
        patt = re.compile(".+?(\d+) kB")
        mFullArr = np.empty(nBatch*nTrackPerBatch)
        sArr = None
        try:
            for n in range(nBatch):
                mArr, sArr = self.cSettle(nIterPerBatch*self.L**2, 
                                 trackFuncList = [self.cGetm, self.getLatticeCopy],
                                 nTrack = nTrackPerBatch)
                mFullArr[n*nTrackPerBatch:(n+1)*nTrackPerBatch] = mArr[:]
                
                # Check for nucleation
                if np.sign(mArr[0]) is -np.sign(mArr[-1]):
                    break
                
                # Check Available Memory
                f = open("/proc/meminfo", 'r')
                fstr = f.read()
                l = patt.findall(fstr)
                p = int(l[1])/float(l[0])
                if p < 0.1:
                    print "Not enough RAM", p
                    break
                else:
                    print "Still OK on RAM", p
                print "On batch", n
        finally:
            return mFullArr, sArr
    
    def cSettle(self, N, trackFuncList = None, nTrack = None, showy = False, 
                endOnFlip = False):
        '''
        Settles faster using c
        
        Args:
        N -- the number of iterations to be performed
        
        Kwargs:
        trackFuncList -- (default None) a list of functions that are
            members of this object, that return single values.
        nTrack  -- an int which is the number of measurements and the
            number of shows.
        showy -- (T/F) If True, then show plots as you go along. The
            number of plots plots is given by nTrack
        endOnFlip -- (T/F) if True, end when the magnetization becomes
            negative.
        '''
        # Make a list to track any values that need to be tracked.
        if trackFuncList is not None:
            self.measureList = []
            for trackFunc in trackFuncList:
                dtype = None
                if isinstance(trackFunc(),object):
                    dtype = object
                self.measureList.append(np.empty(nTrack, dtype=dtype))
            if nTrack is None:
                nTrack = N
        else:
            nTrack = 1
        
        # If I end the simulation midway, I'd still like to get the
        # results I had so far.
        try:
            # Run through the settling, while measuring.
            for n in range(nTrack):
                print n
                 
                # Make the measurement(s)
                if trackFuncList is not None:
                    for i, trackFunc in enumerate(trackFuncList):
                        self.measureList[i][n] = trackFunc()
                
                # This is done whether I'm tracking or not. Thus nTrack has
                # two meanings, and may be meaningful without track functions.
                if showy:
                    self.plotMeasures()
                    self.show()
                    mplot.draw()
                    mplot.show(block = False)
                
                # Check to see if the spins have flipped, and end if that's what
                # the user wanted.
                if endOnFlip and self.cGetm() < 0:
                    break
                
                # Settle some more
                self.lowlevel_cSettle(N/nTrack)
        except:
            print "We encountered a problem. Attempting to return the measures."
            raise
        finally:
            # If nothing was tracked, return None
            if trackFuncList is not None:
                ret = self.measureList
            else:
                ret = None
        
        return ret
    
    def show(self, *arg, **kwarg):
        '''
        This function plots a colorplot of the ising lattice. I will always
        plot to the figure labeled "Lattice Plot". If no such figure exists,
        it will create one. However, the active figure will be reset once the
        function has concluded.
        '''
        # Get the old figure (if there are any figures)
        replaceFig = False
        if mplot.get_fignums():
            oldFig = mplot.gcf()
            replaceFig = True
        
        # Set/create the new active figure.
        fig = mplot.figure(10000)
        
        # Make the plot
        ret = mplot.imshow(self.s, *arg, interpolation='none', vmin = -1, vmax = 1,                             **kwarg)
        
        # Go back to the old figure (If I ever changed)
        if replaceFig:
            mplot.figure(oldFig.number)
        
        return ret
    
    def cGetM(self):
        cgetM = self.lib.getM 
        return cgetM(ctypes.c_void_p(self.s.ctypes.data), 
                     ctypes.c_int(self.L))
    
    def cGetInteVal(self):
        cgetInterVal = self.lib.getInterVal
        return cgetInterVal(ctypes.c_void_p(self.s.ctypes.data),
                            ctypes.c_int(self.L),
                            ctypes.c_int(self.R))
    
    def cGetInterLattice(self, R=None):
        '''
        Return a lattice with the interaction potential calculated for each
        spin lattice site, modulo J. In other words,
        
            ret[i,j] = s[i,j]*sum(s[k,l])
        
        where the sum is performed over neighbors within a block with sides
        2R+1. If R = None, then the range of interaction for this lattice is
        used. Note that this function passes things to the libIsing c function,
        and is thus quite fast.
        
        Inputs:
        R -- (default: None) The radius of interaction. If None, it will be
                set to this lattice's radius of interaction.
        
        Returns:
        iLatt -- the lattice of interaction potentials divided by J.
        
        '''
        # By default, we will just calculat the interaction part of the
        # ising energy.
        if R is None:
            R = self.R
        
        # Initialize the lattice. Make sure the type is set to int, else
        # this process will return nonsense.
        iLatt = np.zeros(self.s.shape, dtype = np.int)
        
        # Get the c-function wrapper.
        cgetInterLattice = self.lib.getInterLattice
        cgetInterLattice(ctypes.c_void_p(self.s.ctypes.data),
                         ctypes.c_void_p(iLatt.ctypes.data),
                         ctypes.c_int(self.L),
                         ctypes.c_int(R),
                         ctypes.c_double(self.h))
        
        self.fixInts(iLatt)
        return iLatt
    
    def cGetInterSpin(self, i, j):
        cGetInterSpin = self.lib.getInterSpin
        return cGetInterSpin(ctypes.c_void_p(self.s.ctypes.data),
                             ctypes.c_int(self.L),
                             ctypes.c_int(self.R),
                             ctypes.c_int(i),
                             ctypes.c_int(j))
    
    def getM(self):
        return sum(self.s.flatten())
    
    def cGetm(self):
        return float(self.cGetM())/self.L**2
    
    def cGetmLattice(self, R = None):
        if R is None:
            R = self.R
        
        cGetMLatt = self.lib.getMLattice
        
        MLatt = np.zeros(self.s.shape, dtype=np.int)
        
        cGetMLatt(ctypes.c_void_p(self.s.ctypes.data),
                  ctypes.c_void_p(MLatt.ctypes.data),
                  ctypes.c_int(self.L),
                  ctypes.c_int(R))
        
        self.fixInts(MLatt)
        return MLatt/(2.*R + 1)**2
                            
    
    def getm(self):
        return float(self.getM())/self.L**2
    
    def getLatticeCopy(self):
        return self.s.copy()
    
    def getLCSize(self):
        clusters = self.getClusters()
        maxSize = 0
        for cluster in clusters:
            if len(cluster) > maxSize:
                maxSize = len(cluster)
        return maxSize
    
    def getClusters(self):
        '''
        Get cluster statistics form the lattice. Partly implemented in C.
        '''
        # We are interested in the clustering of lattice sites that are
        # in the same direction as the magnetic field. These will be our
        # "sites".
        iArr, jArr = np.where(self.h*self.s > 0)
        ijArr = np.empty([len(iArr), 2])
        ijArr[:,0] = iArr[:]
        ijArr[:,1] = jArr[:]
        
        print iArr.shape
        siteDict = {}
        for (i,j) in ijArr:
            siteDict[(int(i),int(j))] = Site(int(i),int(j))
        
        # We want to connect up the sites in a random predetermined order.
        order = np.arange(len(siteDict))
        np.random.shuffle(order)
        
        # This is the bond probability.
        P = 1 - np.exp(-2*self.b*self.J)
        print P
        
        def connect():
            return np.random.choice([True, False], p = [P, 1-P])
        
        # Go through the sites
        clusters = []
        for n in order:
            (i, j), thisSite = siteDict.items()[n]
            for i_ in range(i-self.R, i + self.R + 1):
                i_ = self.__fixBounds(i_)
                for j_ in range(j-self.R, j+self.R +1):
                    j_ = self.__fixBounds(j_)
                    
                    # A site can't connect to itself.
                    if i_ == i and j_ == j:
                        continue
                    
                    # If the site is "Occupied"
                    if self.s[i_,j_]*self.h > 0 and connect():
                        thatSite = siteDict[(i_, j_)]
                        if not thisSite.cluster and not thatSite.cluster:
                            clusters.append([thisSite, thatSite])
                            thisSite.cluster = clusters[-1]
                            thatSite.cluster = clusters[-1]
                        elif not thatSite.cluster:
                            thisSite.cluster.append(thatSite)
                            thatSite.cluster = thisSite.cluster
                        elif not thisSite.cluster:
                            thatSite.cluster.append(thisSite)
                            thisSite.cluster = thatSite.cluster
                        elif thisSite.cluster is thatSite.cluster:
                            pass
                        else:
                            clusters.remove(thisSite.cluster)
                            clusters.remove(thatSite.cluster)
                            newCluster = thisSite.cluster + thatSite.cluster
                            clusters.append(newCluster)
                            for aSite in newCluster:
                                aSite.cluster = newCluster
        
        return clusters

class Site(object):
    def __init__(self, i, j, cluster = None):
        self.loc = (i,j)
        self.cluster = cluster
        return
    
    def __repr__(self):
        return "Site(%d, %d)" % self.loc
    
    def __str__(self):
        return self.__repr__()
    
class Cluster(object):
    '''
    This object allows of the easy handling of clusters.
    '''
    def __init__(self, sites):
        self.sites = copy(sites)
        return
    
    def __add__(self, c):
        if not isinstance(c, Cluster):
            raise Exception("Cannot add type %s to Cluster" % type(c))
        
        return Cluster(self.sites + c.sites)
    
    def __len__(self):
        return len(self.sites)
    
    def __contains__(self, site):
        return (site in self.sites)
    
    def addSite(self, site):
        self.sites.append(site)

def settleMany(I, nRuns, nIter, nProcs, verbose=False):
    '''
    Given a scenario where there may be a nucleating droplet, run from
    `nRuns' settles from that same point using different random seeds
    to determine if we really are at the nucleating droplet.
    
    I -- Ising object with lattice at supposed critical droplet.
    nRuns -- the number of different runs
    nIter -- the number of iterations for each run. Note that an "iteration"
             here actionally L^2 iterations, in other words in one iteration,
             every spin should have been touched an average of once.
    '''
    Q = Queuer(nProcs)
    Q.startMonitor()
    
    print "Running the simulations..."
    myDir = datetime.now().strftime("Outputs/batch_%Y%m%d-%H%M%S")
    os.mkdir(myDir)
    for i in range(nRuns):
        Q.addToQueue(I.cSettle, 
                     nIter*I.L**2, 
                     trackFuncList = [I.cGetm, I.getLatticeCopy], 
                     nTrack = nIter,
                     baseName = "%s/Run_%d" % (myDir,i),
                     seed = np.random.randint(0,100000),
                     verbose = verbose)
    Q.waitForAll()
    
    print "Compiling the data..."
    ret = []
    for i in range(nRuns):
        m = np.load("%s/Run_%d_cGetm.npy" % (myDir, i))
        ret.append(m.copy())
    
    return ret

def makeGif(sArr, name, start = 0, end = -1, loc = '.', 
            vmin = -1, vmax=1, **kwargs):
    '''
    Function to convert all or part of an array of 2D arrays into
    a gif. The 2D arrays are plotted using imshow.
    
    Args:
    sArr    -- the array of arrays
    name    -- the name of the output file (no filetype suffix please)
    
    Kwargs:
    start   -- (default 0) the index within sArr to begin
    end     -- (default -1) the index within sArr at which to end.
    loc     -- (default '.') the directory in which to do all this.
    (Also any keyword arguments appropriate for matplotlib.pyplot.imshow)
    '''
    tmpDir = "%s_pngs" % (name,)
    loc += "/"
    
    # Make a unique file directory 
    n = 1
    newTmpDir = tmpDir
    while newTmpDir in os.listdir(loc):
        newTmpDir = "%s-%d" % (tmpDir, n)
        n+=1
    os.mkdir(loc + newTmpDir)
    
    # Plot the 2D arrays using imshow.
    fig = mplot.figure()
    for i, s in enumerate(sArr[start:end]):
        fig.clf()
        mplot.imshow(s, interpolation = 'none', vmin = vmin, vmax = vmax, **kwargs)
        mplot.title("Timestep: %d" % (start + i))
        fig.savefig(loc + newTmpDir + "/%04d" % i)
        print "Saved at iteration", i
    
    # Convert the set of images to a gif.
    os.system("convert -delay 50 -loop 0 %s%s/*.png %s.gif" % (loc,newTmpDir,name))
    
    return
