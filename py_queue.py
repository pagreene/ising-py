from multiprocessing import Process
from threading import Thread, Lock
from time import sleep

#================================================================================
# Queuing object
#================================================================================

class QueueError(Exception):
    pass

class Queuer(object):
    '''
    This object acts as a pool of a set number of processors to which
    to send tasks. Each new task is in fact a newly spawned process,
    however there is a limit to the number of processes that may be
    created. The purpose of this object is to maintain that limit.

    All access to the multiprocessing module should go through this
    function, as it implement thread safety procedures that are lacking
    in multiprocessing. If you want to create a single process, simply
    create a 1-process pool and add your task to the queue.
    
    Along with the above, and somewhat in explanation, this object
    should NEVER BE SHARED BY MULTIPLE PROCESSES. EVER. Just don't
    do it. It is designed to be shared by multiple threads.
    '''
    
    # It's really a pain, but multipocessing is NOT thread
    # safe (kind of stupid, right?) so I need to lock any
    # time anything with multiprocessing is done.
    __GLOBAL_LOCK = Lock()
    
    def __init__(self, np):
        # Give this a name
        self.name = "Queue"
        
        # The total number of processes that may be spawned
        # at any given time.
        self.__total = np
        
        # A lock for this object
        self.__lock = Lock()

        # Dictionary of all tasks ever created
        self.__proc_dict = {}

        # List of ids of tasks waiting to run
        self.__queue = []
        
        # List of ids of tasks currently running
        self.__running = []
        
        # List of ids of tasks that have already run
        self.__done = []

        # The most recent task id
        self.__n = 0

        # Signals for the monitor
        self.__ending = False
        self.__closing = False
        self.__open = False

        return
    
    def changeNum(self, np):
        '''
        API to change the number of processes allowed to spawn.
        '''
        self.__total = np
        return
    
    def getQueue(self):
        return self.__queue[:]

    def getRunning(self):
        return self.__running[:]
        
    def getDone(self):
        return self.__done[:]

    def isOpen(self):
        return self.__open
    
    def getNumProcs(self):
        return self.__total
    
    def addToQueue(self, target, *args, **kwargs):
        '''
        API to add a task to the queue. Arguments are the arguments for process 
        spawning. The monitor MUST be running, or this method will except.
        '''
        fmt = "'%s' - Sorry, we're %s. No tasks are being accepted."
        if not self.__open:
            raise QueueError(fmt % ("not open yet", self.name))

        if self.__ending or self.__closing:
            raise QueueError(fmt % ("closing", self.name))
        
        # Create the process object. These locks may not be strictly necessary,
        # however, as a precaution I am locking anytime I interract with the
        # multiprocessing module, because it is not thread safe.
        self.__GLOBAL_LOCK.acquire()
        try:
            p = Process(name = str(self.__n + 1), target = target, args = args, 
                        kwargs = kwargs)
        finally:
            self.__GLOBAL_LOCK.release()
        
        # Locking here because different threads could be attempting to modify
        # __n, and I need the proc ids to be unique.
        self.__lock.acquire()
        try:
            self.__n += 1
            qid = self.__n
            self.__proc_dict[qid] = p
            self.__queue.append(qid)
        finally:
            self.__lock.release()

        return qid
        
    def startMonitor(self):
        '''
        API to spawn the threads that monitor and handle the queue
        '''
        if self.__ending or self.__closing:
            print "WARNING: '%s' has already been openned and closed." % self.name
            print "Are you sure you want to open again?"
            self.__ending = False
            self.__closing = False
        
        # Create the thread objects
        self.__queue_thread = Thread(target = self.__queueMonitor, 
                                     name = self.name + "_queued")
        self.__running_thread = Thread(target = self.__runningMonitor, 
                                       name = self.name + "_running")
        self.__queue_thread.daemon = True
        self.__running_thread.daemon = True

        # Start the threads
        self.__queue_thread.start()
        self.__running_thread.start()
        
        # We are now open!
        self.__open = True
        return

    def __runningMonitor(self):
        '''
        Internal API to monitor the processes that are running. This
        is run as a thread (see 'startMonitors').
        '''
        # Define a method for checking '__running' (list of tasks
        # currently running)
        def checkRunning():
            for qid in self.__running[:]:
                proc = self.__proc_dict[qid]
                
                # Check if process is done and join if it is. 
                self.__GLOBAL_LOCK.acquire()
                try:
                    proc_is_dead = False
                    if not proc.is_alive():
                        proc.join()
                        proc_is_dead = True
                finally:
                    self.__GLOBAL_LOCK.release()
                
                if proc_is_dead:
                    # Update the lists
                    self.__running.remove(qid)
                    self.__done.append(qid)
            
            print "Checked the running processes", len(self.__running), len(self.__queue)
            sleep(1)
            return
        
        # While we're open, check __running no matter what, when we're done or
        # we except, still keep checking.
        try:
            while not (self.__closing or self.__ending):
                checkRunning()
            
            while len(self.__queue):
                checkRunning()
        finally:
            while len(self.__running):
                checkRunning()
        print "Yay i'm done"
        return
        
    def __queueMonitor(self):
        '''
        Internal API to monitor the number of processes and handle queing. 
        This is run as a thread (see 'startMonitors').
        '''
        # Define a method for checking '__queue' (list of tasks waiting to run)
        def checkQueue():
            if len(self.__running) < self.__total:
                for qid in self.__queue[:]:
                    # Start the process
                    self.__GLOBAL_LOCK.acquire()
                    try:
                        self.__proc_dict[qid].start()
                    finally:
                        self.__GLOBAL_LOCK.release()
                    
                    # Update the lists
                    self.__queue.remove(qid)
                    self.__running.append(qid)
                    
                    # Check if we're at (or over) our limit.
                    if len(self.__running) >= self.__total:
                        return
            sleep(1)
            return
        
        # Until the end or closing value is set, keep monitoring
        while not self.__ending and not self.__closing:
            checkQueue()
        else:
            if self.__closing:
                print "We're closing! Waiting for %d task(s)." % len(self.__queue)
                while len(self.__queue):
                    checkQueue()
        return
    
    def waitForAll(self, timeout = None):
        '''
        API to block until the processes in the queue run and complete. 
        A timeout may be specified in seconds; default is None.
        '''
        self.__closing = True
        self.__queue_thread.join(timeout)
        self.__running_thread.join(timeout)
        self.__open = False
        return len(self.__queue)
    
    def waitForRunning(self, timeout = None):
        '''
        API to block until the processes currently running complete. A 
        timeout may be specified in seconds; default is None.
        '''
        self.__ending = True
        self.__queue_thread.join(timeout)
        self.__running_thread.join(timeout)
        self.__open = False
        return len(self.__running)
    
    def waitFor(self, qid, timeout = None):
        '''
        API to block until the specified pid is completed. Timeout may be
        specified in seconds if desired. If timeout is exceded, this function
        will raise.
        '''
        class TimeoutError(Exception):
            pass
        
        def wait(t):
            sleep(1)
            t += 1
            if timeout != None:
                if t > timeout:
                    msg = "Timeout exceded while waiting for process '%s'." % self.__proc_dict[qid].name
                    raise TimeoutError(msg)
            return t
        
        t = 0
        # Wait for the process to start
        while qid in self.__queue:
            if self.__ending:
                print "Monitor is ending, process '%s' will not be run." % self.__proc_dict[qid].name
                return
            t = wait(t)

        # Wait for the process to complete
        while qid in self.__running:
            t = wait(t)
        
        return
