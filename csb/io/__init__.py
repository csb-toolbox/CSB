"""
Basic IO Function
"""
def Dump(this, filename, gzip = 0, lock = None):
    """
    Writes a python object to a file, using python cPickle
    Supports also '~' or '~user'.

    @param this: The object, which will be written to disk
    @type param: Any python object
    @param filename: Filename of the new file
    @type filename: String
    @param gzip: Use gzip to compress the file
    @type gzip: Boolean
    @param lock: Use a lockfile to restrict access
    """

    import os
    
    ## check whether file is locked
    if lock is not None:
        import os

        path, file = os.path.split(filename)
        lockfile = '%s/%s.lock' % (__lockPath, file)

        if os.path.exists(lockfile):
            raise '%s is locked (%s)' % (filename, lockfile)

        ## create lock file
        open(lockfile, 'w').close()

    filename = os.path.expanduser(filename)
        
    if gzip:
        import gzip
        f = gzip.GzipFile(filename, 'w')
    else:
        f = open(filename, 'w')
        
    if type(this).__name__ == 'array':
        import Numeric
        p = Numeric.Pickler(f)
        p.dump(this)

    else:
        import cPickle
        cPickle.dump(this, f, 1)
        
    f.close()

    if lock is not None:
        ## release lock
        import os
        try:
            os.remove(lockfile)
        except:
            raise 'missing lockfile %s' % lockfile


def Load(filename, gzip = 0, lock = None):
    """
    Unpickle an object from filename
    
    @param filename: Filename pickled object
    @param gzip: Use gzip to uncompress the file
    @type gzip: Boolean    
    @param lock: Use a lockfile to restrict access
    @return: Python object unpickled from file
    """
    import cPickle, os

    ## file locked?
    if lock is not None:
        import os

        path, file = os.path.split(filename)
        lockfile = '%s/%s.lock' % (__lockPath, file)
        
        if os.path.exists(lockfile):
            raise '%s is locked (%s)' % (filename, lockfile)

        ## lock file
        open(lockfile, 'w').close()

    filename = os.path.expanduser(filename)
        
    if gzip:
        import gzip
        try:
            f = gzip.GzipFile(filename)
        except:
            return
        
    f = open(filename)

    try:
        this = cPickle.load(f)
    except:
        import Numeric

        f.close()
        try:
            f = gzip.GzipFile(filename)
        except:
            f = open(filename)
        
        u = Numeric.Unpickler(f)
        this = u.load()

    f.close()

    if lock is not None:
        ## release lock
        import os
        try:
            os.remove(lockfile)
        except:
            raise 'missing lockfile %s' % lockfile

    return this
