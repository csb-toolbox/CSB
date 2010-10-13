"""
Common IO-related utility functions and classes. 
"""

NEWLINE = '\n'
    
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
        
    def __del__(self):
        
        try:
            self._stream.close()
        except:
            pass

    def entries(self):
        """
        Return an iterator over all entries from the stream/flat file.
        
        @return: iterable over all entries read from the stream
        @rtype: generator
        """

        self._stream.seek(0)
        
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

    def readall(self):
        """
        Return a list of all entries in the stream.
        
        @rtype: list
        """
        return list(self.entries())

class EntryWriter(object):
    """
    A simple stream writer. The best way to use it is::
    
        with EntryWriter(output_file, close=True) as out:
            out.write(object)
    
    In this way the stream is automatically closed at the end of the block.
    
    @param destination: output file name or opened stream
    @type destination: str or stream
    @param newline: new line string (the default is L{csb.io.NEWLINE})
    @type newline: str
    @param close: if True (the default), the stream is automatically
                  closed when the object is destroyed
    @type close: bool  
    """
    
    def __init__(self, destination, close=True, newline=NEWLINE):
        
        self._stream = None
        self.newline = newline
        self.autoclose = close
        
        if isinstance(destination, file) or hasattr(destination, 'write'):
            self._stream = destination
        elif isinstance(destination, basestring):
            self._stream = open(destination, 'w')
            self.autoclose = True
            
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self.autoclose:
            self.close()
    
    def __del__(self):
        if self.autoclose:
            self.close()
    
    @property
    def stream(self):
        return self._stream
    
    def close(self):
        """
        Close the destination stream.
        """
        try:
            self._stream.close()
        except:
            pass        
            
    def write(self, data):
        """
        Write a chunk of sting data to the destination stream. 
        """
        self._stream.write(data)
    
    def writeline(self, data):
        """
        Same as C{write}, but appends a newline character at the end.
        """        
        self._stream.write(data)
        self._stream.write(self.newline)
        
    def writeall(self, entries, delimiter=NEWLINE):
        """
        Write all C{entries} to the destination stream, separating them with
        C{delimiter}
        
        @param entries: a collection of objects
        @type entries: iterable
        @param delimiter: append this string after each entry (the default is a 
                          C{self.newline} character)
        @type delimiter: str 
        """        
        if delimiter == NEWLINE:
            delimiter = self.newline
        for entry in entries:
            self.write(entry)
            self.write(delimiter)

def dump(this, filename, gzip=False, lock=None, lock_path='~/'):
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

        _path, file = os.path.split(filename)
        lockfile = '%s/%s.lock' % (lock_path, file)

        if os.path.exists(lockfile):
            raise '%s is locked (%s)' % (filename, lockfile)

        ## create lock file
        open(lockfile, 'w').close()

    filename = os.path.expanduser(filename)

    if gzip:
        import gzip
        f_handle = gzip.GzipFile(filename, 'w')
    else:
        f_handle = open(filename, 'w')

    if type(this).__name__ == 'array':
        import Numeric
        p = Numeric.Pickler(f_handle)
        p.dump(this)

    else:
        import cPickle
        cPickle.dump(this, f_handle, 1)

    f_handle.close()

    if lock is not None:
        ## release lock
        try:
            os.remove(lockfile)
        except:
            raise IOError('missing lockfile %s' % lockfile)

def load(filename, gzip=False, lock=None, lock_path='~/'):
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
        _path, file = os.path.split(filename)
        lockfile = '%s/%s.lock' % (lock_path, file)

        if os.path.exists(lockfile):
            raise IOError('%s is locked (%s)' % (filename, lockfile))

        ## lock file
        open(lockfile, 'w').close()

    filename = os.path.expanduser(filename)

    if gzip:
        import gzip
        try:
            f_handle = gzip.GzipFile(filename)                          #@UnusedVariable
        except:
            return

    f_handle = open(filename)

    try:
        this = cPickle.load(f_handle)
    except:
        import Numeric
        f_handle.close()
        try:
            f_handle = gzip.GzipFile(filename)
        except:
            f_handle = open(filename)

        unpickler = Numeric.Unpickler(f_handle)
        this = unpickler.load()

    f_handle.close()

    if lock is not None:
        ## release lock
        try:
            os.remove(lockfile)
        except:
            raise IOError('missing lockfile %s' % lockfile)

    return this

