"""
Generic I/O utility objects.

Here is a list of the most essential classes in this module:

    1. temporary file system objects: L{TempFile}, L{TempFolder} 
    2. special/decorated streams: L{MemoryStream}, L{AutoFlushStream}
    3. reusable stream readers and writers: L{EntryReader}, L{EntryWriter}
    4. convenient communication with the shell: L{Shell}

In addition, csb.io is also part of the CSB compatibility layer. In order to
ensure cross-interpreter compatibility, always use the following csb.io objects:

    - L{MemoryStream} instead of (c)StringIO
    - csb.io.Pickle instead of pickle or cPickle
    - csb.io.urllib instead of urllib or urllib.request
    
See also L{csb.core} for additional notes on compatibility.    
"""

import os
import time
import errno
import shlex
import shutil
import tempfile
import subprocess

import csb.core


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
        

NEWLINE = "\n"


class Shell(object):
    
    POLL = 1.0

    @staticmethod        
    def run(cmd, timeout=None):
        """
        Run a shell command and return the output.
        
        @param cmd: shell command with its arguments
        @param cmd: tuple or str
        @param timeout: maximum duration in seconds
        @type timeout: float or None
        
        @rtype: L{ShellInfo}
        @raise InvalidCommandError: on invalid executable
        @raise TimeoutError: when the timeout is expired
        """
        
        if isinstance(cmd, csb.core.string):
            cmd = shlex.split(cmd)
        
        try:
            cmd = tuple(cmd)
            start = time.time()    
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            
            if timeout is not None:
                while True:
                    if process.poll() == 0:
                        break
                    elif time.time() >= (start + timeout):
                        try:
                            process.kill()
                        except:
                            pass
                        raise TimeoutError(cmd, timeout)
                    else:
                        time.sleep(Shell.POLL)
    
            stdout, stderr = process.communicate()                
            code = process.returncode
            
        except OSError as oe:
            if oe.errno == 2:
                raise InvalidCommandError(oe.strerror, cmd)
            else:
                raise
                
        return ShellInfo(code, stdout.decode() or '', stderr.decode() or '', cmd)

    @staticmethod
    def runstrict(cmd, timeout=None):
        """
        Same as L{Shell.run()}, but raises L{ProcessError} on bad exit code.
        
        @param cmd: shell command with its arguments
        @param cmd: tuple or str
        @param timeout: maximum duration in seconds
        @type timeout: float or None        
        
        @rtype: L{ShellInfo}
        @raise ProcessError: on bad exit code
        @raise TimeoutError: when the timeout is expired        
        """
        si = Shell.run(cmd, timeout=timeout)
        
        if si.code == 0:
            return si
        else:
            raise ProcessError(si)            
   
class ProcessError(Exception):
    """
    Raised on L{Shell.run()} failures.
    @type context: L{ShellInfo} 
    """    
    def __init__(self, context, *args):
        self.context = context
        super(ProcessError, self).__init__(context, [])
        
    def __str__(self):
        return 'Bad exit code: #{0.code}'.format(self.context)

class TimeoutError(ProcessError):
    """
    Raised on L{Shell.run()} timeouts. 
    """    
    def __init__(self, cmd, timeout):
        
        self.timeout = timeout
        context = ShellInfo(None, '', '', cmd)
        
        super(TimeoutError, self).__init__(context)
        
    def __str__(self):
        return 'The process "{0.context.cmd}" did not finish in {0.timeout}s'.format(self)
        
class InvalidCommandError(ValueError):
    """
    Raised when L{Shell.run()} encounters an OSError.
    """
    def __init__(self, message, cmd):
        
        self.program = cmd[0]
        if csb.core.iterable(cmd):
            cmd = ' '.join(cmd)
        self.cmd = cmd
        self.msg = message
        
        super(InvalidCommandError, self).__init__(message, cmd)
        
    def __str__(self):
        return self.msg

class ShellInfo(object):
    """
    Shell command execution info
    """
    
    def __init__(self, code, stdout, stderr, cmd):
        
        self.code = code
        self.stdout = stdout or ''
        self.stderr = stderr or ''
        self.cmd = ' '.join(cmd)
        

class MemoryStream(StringIO):
    """
    In-memory stream object. Can be used in a context manager.
    """
    def __enter__(self):
        return self
    
    def __exit__(self, *a, **k):
        try:
            self.close()
        except:
            pass
    
class AutoFlushStream(csb.core.Proxy):
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
        
class TempFile(csb.core.Proxy):
    """
    Create a temporary file and take care of deleting it upon object
    destruction. The file can be opened multiple times on any platform, unlike
    the case with tempfile.NamedTemporaryFile (does not work on Windows).
    
        >>> with TempFile() as tmp:
                tmp.write(...)
                open(tmp.name)...
    
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
        """
        Full path and file name
        @rtype: str
        """
        return self.__file
    
class TempFolder(object):
    """
    Create a temporary directory which is automatically wiped when the object
    is closed.

        >>> with TempFolder() as tmp:
                # put some files in tmp.name...
    
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
        """
        Full directory name
        @rtype: str
        """        
        return self.__name    

class EntryReader(object):
    """
    Generic flatfile reader. Provides efficient iterable interface over the entries
    in the specified stream.
    
    @param stream: the source data stream to read
    @type stream: file
    @param start_marker: a string marker which marks the beginning of a new entry
    @type start_marker: str
    @param end_marker: a string marker which signals the end of the file
    @type end_marker: str, None
    """
    def __init__(self, stream, start_marker, end_marker=None):

        if not (hasattr(stream, 'seek') and hasattr(stream, 'readline')):
            raise TypeError('The entry reader requires an opened stream.')

        stream.seek(0)
        self._stream = stream
        self._start_marker = None
        self._end_marker = None

        self.start_marker = start_marker
        self.end_marker = end_marker        
        
    @property
    def start_marker(self):
        return self._start_marker
    @start_marker.setter
    def start_marker(self, value):
        if value is not None:
            value = str(value)        
        self._start_marker = value

    @property
    def end_marker(self):
        return self._end_marker
    @end_marker.setter
    def end_marker(self, value):
        if value is not None:
            value = str(value)
        self._end_marker = value
                
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
        self._newline = NEWLINE
        self._autoclose = True

        self.newline = newline
        self.autoclose = close        
        
        if isinstance(destination, csb.core.string):
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
        """
        Destination stream
        @rtype: stream
        """
        return self._stream
    
    @property
    def newline(self):
        return self._newline
    @newline.setter
    def newline(self, value):
        self._newline = str(value)

    @property
    def autoclose(self):
        return self._autoclose
    @autoclose.setter
    def autoclose(self, value):
        self._autoclose = bool(value)        
    
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
            Pickle.dump(this, stream, 2)
    finally:
        stream.close()

    if lock is not None:
        ## release lock
        try:
            os.rmdir(lockdir)
        except:
            raise IOError('missing lockfile {0}'.format(lockdir))

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
            raise IOError('missing lockfile {0}'.format(lockdir))

    return this


