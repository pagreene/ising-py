#!/usr/bin/python

import numpy as np
import matplotlib.pyplot  as mplot
import ctypes
import re
import os
from py_queue import Queuer
from datetime import datetime
from copy import copy
from time import sleep

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
    '''
    This object is the primary workhorse of this code. The methods of this
    object control the actions that can be taken on the lattice, by changing
    the boltzman factor `b` (for `beta`), the magnetic field `h`, the
    interaction strength `J`, and the interaction radius `R`. The dimensions
    of the lattice are set at instantiation, but cannot be changed.
    '''
    def __init__(self, J, h, L, b, R):
        self.J = J
        self.h = h
        self.__L = L
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
        if j >= self.__L:
            j -= self.__L
        elif j < 0:
            j += self.__L
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
    
    def lowlevel_cSettle(self, N, seed=None):
        '''
        This function should perform the exact same function as settle,
        however it calls a method written in c, which should make it MUCH
        faster. There is some translation that needs to be done first though.
        '''
        # Pass the stuff to the c code.
        if seed is None:
            seed = np.random.randint(0,100000)
        
        csettle = self.lib.settle
        csettle(ctypes.c_void_p(self.s.ctypes.data),
                ctypes.c_int(self.__L),
                ctypes.c_int(self.R),
                ctypes.c_double(self.h),
                ctypes.c_double(self.J),
                ctypes.c_double(self.b),
                ctypes.c_int(seed),
                ctypes.c_int(N) )
        
        # This is an unfortunate workaround that is needed. There is some
        # type incompatability issue, I think.
        self.s[abs(self.s) > 1] = -np.sign(self.s[abs(self.s) > 1])
        
        #self.fixInts(self.s)
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
                mArr, sArr = self.cSettle(nIterPerBatch*self.__L**2, 
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
                baseName=None, verbose = True, seed = None, eoFlip = False):
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
        eoFlip -- (T/F) if True, end when the magnetization becomes
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
                if verbose:
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
                if eoFlip and self.cGetm() < 0:
                    break
                
                # Settle some more
                self.lowlevel_cSettle(N/nTrack, seed=seed)
        except:
            print "We encountered a problem. Attempting to return the measures."
            raise
        finally:
            # If nothing was tracked, return None
            if trackFuncList is not None:
                if baseName:
                    for i, measureArr in enumerate(self.measureList):
                        np.save("%s_%s" % (baseName,trackFuncList[i].func_name),
                                measureArr)
                
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
        ret = mplot.imshow(self.s, *arg, interpolation='none', vmin = -1, vmax = 1,
                           **kwarg)
        
        # Go back to the old figure (If I ever changed)
        if replaceFig:
            mplot.figure(oldFig.number)
        
        return ret
    
    def cGetM(self):
        cgetM = self.lib.getM 
        return cgetM(ctypes.c_void_p(self.s.ctypes.data), 
                     ctypes.c_int(self.__L))
    
    def cGetInterVal(self):
        cgetInterVal = self.lib.getInterVal
        return cgetInterVal(ctypes.c_void_p(self.s.ctypes.data),
                            ctypes.c_int(self.__L),
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
                         ctypes.c_int(self.__L),
                         ctypes.c_int(R),
                         ctypes.c_double(self.h))
        
        self.fixInts(iLatt)
        return iLatt/(float(2*R + 1)**2)
    
    def cGetInterSpin(self, i, j):
        cGetInterSpin = self.lib.getInterSpin
        return cGetInterSpin(ctypes.c_void_p(self.s.ctypes.data),
                             ctypes.c_int(self.__L),
                             ctypes.c_int(self.R),
                             ctypes.c_int(i),
                             ctypes.c_int(j))
    
    def getM(self):
        return sum(self.s.flatten())
    
    def cGetm(self):
        return float(self.cGetM())/self.__L**2
    
    def cGetmLattice(self, R = None):
        if R is None:
            R = self.R
        
        cGetMLatt = self.lib.getMLattice
        
        MLatt = np.zeros(self.s.shape, dtype=np.int)
        
        cGetMLatt(ctypes.c_void_p(self.s.ctypes.data),
                  ctypes.c_void_p(MLatt.ctypes.data),
                  ctypes.c_int(self.__L),
                  ctypes.c_int(R))
        
        self.fixInts(MLatt)
        return MLatt/(2.*R + 1)**2
                            
    
    def getm(self):
        return float(self.getM())/self.__L**2
    
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

fnamePatt = "Run_%%03d"
def settleMany(I, nRuns, nIter, nTrack, nProcs, verbose = False, eoFlip = False,
                clean = False, trackFuncList = None, retFuncList = None):
    '''
    Given a scenario where there may be a nucleating droplet, run `nRuns' from
    that same point using different random seeds to determine if we really are 
    at the nucleating droplet.
    
    I       --  Ising object with lattice at supposed critical droplet.
    nRuns   --  the number of different runs
    nIter   --  the number of iterations for each run. Note that an "iteration"
                here actionally L^2 iterations, in other words in one iteration,
                every spin should have been touched an average of once.
    nTrack  --  Number of measurements taken on the system.
    nProcs  --  the number of processes that may be operating at any given time.
    verbose --  whether each individual process should be allowed to spew its
                output. In general, this just makes for very messy output.
    eoFlip  --  whether each simulation should end once the spins have 
                flipped, specifically when the magnetization changes sign to match
                the magnetic field.
    clean   --  If true, remove all the files.
    '''
    # By default, we will not return anything.
    if trackFuncList is None:
        trackFuncList = []
    
    # By default, if we don't get a retFuncList, we assume we return everything
    # in the trackFuncList
    if retFuncList is None:
        retFuncList = trackFuncList[:]
    
    # Start the queing aparatus.
    Q = Queuer(nProcs)
    Q.startMonitor()
    
    # This is the name of the directory and the base of the file names for
    # where files will be saved. If clean=True, these will all be erased
    # after the run is done.
    myDir = datetime.now().strftime("Outputs/batch_%Y%m%d-%H%M%S")
    basename = ("%s/" + fnamePatt) % myDir
    
    # Add things to the loop.
    print "Running the simulations..."
    os.mkdir(myDir)
    for i in range(nRuns):
        Q.addToQueue(I.cSettle, 
                     nIter*I.getL()**2, 
                     trackFuncList = trackFuncList, 
                     nTrack = nTrack,
                     baseName = basename % i,
                     seed = np.random.randint(0,100000),
                     verbose = verbose,
                     eoFlip = eoFlip)
    Q.waitForAll()
    
    print "Compiling the data..."
    ret = {}
    for i in range(nRuns):
        
        for trackFunc in trackFuncList:
            funcName = trackFunc.func_name
            fname = (basename + "_%s.npy") % (i,funcName)
            
            # If we are returning this value, then add it to the ret dict.
            if trackFunc in retFuncList:
                # Get the data that was saved.
                m = np.load(fname)
                
                # Make sure there is an entry in the ret dict for this type
                # function output.
                if funcName not in ret.keys():
                    ret[funcName] = []
                ret[funcName].append(m.copy())
            
            # If we are cleaning up, then clean up.
            if clean:
                os.remove(fname)
    
    if clean:
        os.rmdir(myDir) 
    
    return (ret, myDir)

def didFlip(mArr):
    '''
    Check if a trajectory flipped, given the list of ave. magnetizations.
    '''
    i = -1
    while mArr[i] is None:
        i -= 1
    
    return mArr[0]*mArr[i] < 0

def getPercentFlipped(mArrArr):
    '''
    A function to calculate what percent of trajectories ended in nucleation.
    '''
    nFlipped = 0
    for mArr in mArrArr:
        if didFlip(mArr):
            nFlipped += 1
    return float(nFlipped)/len(mArrArr)

class IsingData(IsingLattice):
    '''
    Object to handle the data generated from many Ising lattices, with the goal
    of using that data in a machine learning algorithm.
    '''
    def __init__(self, L, defaultR, psList = None):
        if psList is None:
            psList = []
        self.psList = psList
        
        self.__L = L
        self.R = defaultR
        self.s = None
        self.J = 0
        return
    
    def addLattice(self, s, p):
        self.psList.append((p,s))
        return
    
    def get_mLatts(self):
        for i, ps in enumerate(self.psList):
            self.s = ps[1]
            
        return
    
    def get_iLatts(self):
        return
    

def generateTrainingData(nStart, nIter, nTrack, h, R, maxBack=50, minBack=0):
    '''
    A function to generate the training data.
    '''
    # Create the original lattice
    I = IsingLattice(4, h, 100, 9./16, R)
    
    # Run a few different simulations with different random seeds.
    retDict, dataDir = settleMany(I, nStart, nIter, nTrack, 5, eoFlip = True,
                                 trackFuncList = [I.cGetm, I.getLatticeCopy],
                                 retFuncList = [I.cGetm])
    mArrArr = retDict[I.cGetm.func_name]
    
    ret = []
    
    # Now we will cycle through each of the simulations, generating a set of
    # classified data sets.
    baseName = ("%s/" + fnamePatt) % dataDir
    for i, mArr in enumerate(mArrArr):
        # If the spins flipped...
        if didFlip(mArr):
            
            # Load the s array
            sArr = np.load((baseName + "_getLatticeCopy.npy") % i)
            mArr_ = mArr[mArr > -2] # filter out the nones
            sArr_ = sArr[mArr > -2] # again, filter out Nones.
            
            # Find where they flipped (to flip, they must go through zero)
            w = np.where(abs(mArr_) == min(abs(mArr_)))[0][0]
            
            # Now we run intervention from all the lattices within some range of
            # the flipped lattice.
            dw_i = min(maxBack, w)
            dw_f = min(minBack, w)
            for dw in range(-dw_i, -dw_f):
                I.s = sArr_[w + dw].copy()
                testRetDict, _ = settleMany(I, 100, 500, 50, 5, eoFlip = True,
                                            clean=True, trackFuncList = [I.cGetm])
                test_mArrArr = testRetDict[I.cGetm.func_name]
                mplot.plot(test_mArrArr)
                mplot.draw()
                mplot.show()
                sleep(0.1)
                p = getPercentFlipped(test_mArrArr)
                ret.append((p, I.s.copy()))
                print "\n%d/%d for run %d | p = %.03f\n" % (-dw_i-dw, dw_i-dw_f, i, p)
                if p == 1:
                    print "They're all collapsing now, no need to keep going"
                    break
            
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
        mplot.colorbar()
        mplot.title("Timestep: %d" % (start + i))
        fig.savefig(loc + newTmpDir + "/%04d" % i)
        print "Saved at iteration", i
    
    # Convert the set of images to a gif.
    os.system("convert -delay 50 -loop 0 %s%s/*.png %s.gif" % (loc,newTmpDir,name))
    
    return

def xData(data, a, R):
    '''
    This function turns the raw lattice data into the input into the
    feature data.
    
    Arguments:
    data -- Iterable of tuples (p, s), where p is the probability of flipping
            and s is the lattice configuration.
    a    -- The threshold for `close enough` to being the nucleating droplet.
    R    -- The radius of interaction to be considered
    
    Returns:
    X   -- the feature data
    y   -- the classification
    '''
    N = data[0][1].shape[0]
    X = np.empty([len(data), N*N*2])
    y = np.empty([len(data)])
    I = IsingLattice(4, -0.1, N, 9./16, R)
    for i, (p,s) in enumerate(data):
        I.s = s.copy()
        X[i,:N*N] = I.cGetmLattice().flatten()
        X[i,N*N:] = I.cGetInterLattice().flatten()
        y[i] = round(abs(2*(p-0.5)) + a)
    return (X, y)

def renorm(X, L):
    Xnew = np.empty([L,L])
    N = X.shape[0]/L
    for i in range(L):
        for j in range(L):
            Xnew[i,j] = mean(X[N*i:N*(i+1),N*j:N*(j+1)])
    return Xnew
