class UnbufferedStreamWrapper(object):
    """
    Wrapper around a buffered stream which automatically calls flush()
    after each write(). 
    
    @param stream: the stream object to wrap
    @type stream: file
    """
    
    def __init__(self, stream):
        self.stream = stream
    
    def write(self, data):
        self.stream.write(data)    
        self.stream.flush()
        
    def flush(self):
        self.stream.flush()

class EntryReader(object):
    """
    Generic flatfile reader. Provides efficient iterable interface over the entries
    in the specified stream.
    
    @param stream: the source data stream to read
    @type stream: file
    @param start_marker: a string marker which denotes the beginning of a new entry
    @type start_marker: str
    @param end_marker: a string marker which signals the end of the file
    @type end_marker: str, None
    """
    def __init__(self, stream, start_marker, end_marker):
        
        if not (hasattr(stream, 'seek') and hasattr(stream, 'readline')):
            raise TypeError('The entry reader requires an opened stream.')
        
        stream.seek(0)
        self._stream = stream
        
        self.start_marker = start_marker
        self.end_marker = end_marker 
    
    def entries(self):
        """
        Return an iterator over all entries from the stream/flat file.
        
        @return: iterable over all entries read from the stream
        @rtype: generator
        """

        entry = ''
        in_entry = False
        
        while True:            
            line = self._stream.readline()
            
            if not in_entry:
                if line.startswith(self.start_marker):
                    in_entry = True
                    entry = line                    
            else:
                if line.startswith(self.start_marker):
                    yield entry
                    entry = line
                elif not line or line.strip() == self.end_marker:
                    yield entry
                    break
                else:
                    entry += line                    
                    
            if not line:
                break

        self._stream.close()     
        
    def readall(self):
        """
        Return a list of all entries in the stream.
        
        @rtype: list
        """
        return list(self.entries())

def Dump(this, filename, gzip=False, lock=None):
    """
    Writes a python object to a file, using python cPickle
    Supports also '~' or '~user'.

    @param this: The object, which will be written to disk
    @type this: Any python object
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


def Load(filename, gzip=False, lock=None):
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

