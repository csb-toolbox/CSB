"""
Common IO-related utility functions and classes. 
"""

import os
import time
import errno
import shutil
import tempfile
import csb.pyutils

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
    
try:
    import cPickle as Pickle
except ImportError:
    import pickle as Pickle

try:
    import urllib.request as urllib
except ImportError:
    import urllib
        

NEWLINE = '\n'


class MemoryStream(StringIO):
    """
    In-memory stream object.
    """
    def __enter__(self):
        return self
    
    def __exit__(self, *a, **k):
        try:
            self.close()
        except:
            pass
    
class AutoFlushStream(csb.pyutils.Proxy):
    """
    Wrapper around a buffered stream which automatically calls flush()
    after each write(). This is essentially a proxy/decorator. 
    
    @param stream: the stream object to wrap
    @type stream: file
    """

    def __init__(self, stream):
        super(AutoFlushStream, self).__init__(stream)

    def write(self, data):
        self._subject.write(data)
        self._subject.flush()
        
class TempFile(csb.pyutils.Proxy):
    """
    Create a temporary file and take care of deleting it upon object
    destruction. The file can be opened multiple times on any platform, unlike
    the case with tempfile.NamedTemporaryFile (does not work on Windows).
    
    @param dispose: automatically delete the file
    @type dispose: bool
    @param mode: file open mode (text, binary), default=t
    @type text: str
    """

    def __init__(self, dispose=True, mode='t'):
        
        fd, file = tempfile.mkstemp()
        
        self.__file = file
        self.__fd = fd
        self.__fh = open(self.__file, 'w' + mode)
        self.__mode = mode
        self.__dispose = bool(dispose)
        
        super(TempFile, self).__init__(self.__fh)
        
    def __del__(self):
        
        if self.__dispose:
            try:
                self.close()
            except:
                pass
        
    def __enter__(self):
        return self
    
    def __exit__(self, *args):
        self.close() 
        
    def close(self):
        """
        Flush, close and delete the file.
        """
        
        if not self.__fh.closed:
            self.__fh.flush()
            self.__fh.close()
            os.close(self.__fd)
            
        if os.path.exists(self.__file):
            os.remove(self.__file)
            
    def content(self):
        """
        @return: the current content of the file.
        @rtype: str or bytes
        """
        self.flush()
        with open(self.name, 'r' + self.__mode) as f:
            return f.read()
            
    @property
    def name(self):
        return self.__file
    
class TempFolder(object):
    """
    Create a temporary directory which is automatically wiped when the object
    is closed.
    
    @param dispose: automaticlaly delete the folder and its contents
    @type dispose: bool                
    """
    
    def __init__(self, dispose=True):
         
        name = tempfile.mkdtemp()

        self.__name = os.path.abspath(name)
        self.__dispose = bool(dispose)        
        
    def __del__(self):
        
        if self.__dispose:
            try:
                self.close()
            except:
                pass
        
    def __enter__(self):
        return self
    
    def __exit__(self, *args):
        self.close() 
        
    def close(self):
        """
        Delete the entire directory and its contents.
        """
        if os.path.exists(self.name):
            shutil.rmtree(self.name)
            
    @property
    def name(self):
        return self.__name    

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
        
        if isinstance(destination, csb.pyutils.string):
            self._stream = open(destination, 'w')
            self.autoclose = True
            
        elif hasattr(destination, 'write'):
            self._stream = destination
        
        else:
            raise TypeError(destination)         
            
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

def dump(this, filename, gzip=False, lock=None, timeout=None):
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

    ## check whether file is locked
    ## file locked?
    filename = os.path.expanduser(filename)

    if lock is not None:
        lockdir = filename + '.lock'

        if timeout is not None and timeout > 0:
            end_time = timeout + time.time()

        while True:
            try:
                os.mkdir(lockdir)
            except OSError as e:
                # File is already locked
                if e.errno == errno.EEXIST:
                    if timeout is not None and time.time() > end_time:
                        raise IOError("Failed to acquire Lock")
                else:
                    raise IOError("Failed to acquire Lock")
            else:
                break

    if gzip:
        import gzip
        stream = gzip.GzipFile(filename, 'wb')
    else:
        stream = open(filename, 'wb')

    try:
        if type(this).__name__ == 'array':
            import Numeric                                                  #@UnresolvedImport
            p = Numeric.Pickler(stream)
            p.dump(this)
        else:
            Pickle.dump(this, stream, 1)
    finally:
        stream.close()

    if lock is not None:
        ## release lock
        try:
            os.rmdir(lockdir)
        except:
            raise IOError('missing lockfile %s' % lockdir)

def load(filename, gzip=False, lock=None, timeout=None):
    """
    Unpickle an object from filename
    
    @param filename: Filename pickled object
    @param gzip: Use gzip to uncompress the file
    @type gzip: Boolean    
    @param lock: Use a lockfile to restrict access
    
    @return: Python object unpickled from file
    """
    ## file locked?
    filename = os.path.expanduser(filename)

    if lock is not None:
        lockdir = filename + '.lock'

        if timeout is not None and timeout > 0:
            end_time = timeout + time.time()

        while True:
            try:
                os.mkdir(lockdir)
            except OSError as e:
                # File is already locked
                if e.errno == errno.EEXIST:
                    if timeout is not None and time.time() > end_time:
                        raise IOError("Failed to acquire Lock")
                else:
                    raise IOError("Failed to acquire Lock")
            else:
                break

    if gzip:
        import gzip
        stream = gzip.GzipFile(filename)
        try:
            stream.readline() 
            stream.seek(0)
        except:
            stream.close()
            raise

    else:
        stream = open(filename, 'rb')

    try:
        this = Pickle.load(stream)
    except:                                
        stream.close()
        import Numeric                                                      #@UnresolvedImport        
        try:
            stream = gzip.GzipFile(filename)
        except:
            stream = open(filename, 'rb')

        try:
            unpickler = Numeric.Unpickler(stream)
            this = unpickler.load()
        except:
            stream.close()
            raise

    stream.close()

    if lock is not None:
        ## release lock
        try:
            os.rmdir(lockdir)
        except:
            raise IOError('missing lockfile %s' % lockdir)

    return this


