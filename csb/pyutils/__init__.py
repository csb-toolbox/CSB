"""
Common high-level Python code utilities.
 
@warning: This module will try to use an OrderedDict from Python 2.7+'s 
          collections module. If such class is missing, the importer will
          silently fallback to ActiveState's OrderedDict implementation.
          
@todo: Remove the backup OrderedDict once CSB is officially migrated to Python
       version 2.7+ 
"""

import re
import sys
import time
import shlex
import subprocess
import threading

from abc import ABCMeta, abstractproperty, abstractmethod


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
        
        if isinstance(cmd, basestring):
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
                
        return ShellInfo(code, stdout or '', stderr or '', cmd)

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
        if hasattr(cmd, '__iter__'):
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

class InterruptibleThread(threading.Thread):
    
    def __init__(self, target, name=None, args=[], kwargs={}):
        
        super(InterruptibleThread, self).__init__(target=target, name=name, args=args, kwargs=kwargs)
        self.setDaemon(True)
        self.__result = None
        self.__target = target
        self.__name = name
        self.__args = args
        self.__kwargs = kwargs
        
    @property
    def result(self):
        return self.__result
    
    def run(self):
        try:
            self.__result = self.__target(*self.__args, **self.__kwargs)
        except Exception as ex:
            sys.stderr.write(ex)
            self.__result = None

def singleton(klass):
    """
    Singleton class decorator.
    """
    instances = {}
    
    def get():
        if klass not in instances:
            instances[klass] = klass()
        return instances[klass]
    
    return get
            
class Proxy(object):
    """
    Base class implementing the proxy pattern. Subclasses that define
    a customized constructor must call super().__init__(subject) for proper
    initialization.
    
    @param subject: an arbitrary object, hidden behind the proxy
    """

    def __init__(self, subject):
        self._subject = subject
        
    def __getattr__(self, name):
        try:
            return object.__getattribute__(self, name)
        except AttributeError:
            subject = object.__getattribute__(self, '_subject')
            return getattr(subject, name)                        
        
class REMatchProxy(object):
    
    def __init__(self, match):

        self._start = [match.start()]
        self._end = [match.end()]
        self._groups = match.groups()        
                
        for i, dummy in enumerate(self._groups, start=1):
            self._start.append(match.start(i))
            self._end.append(match.end(i))            
    
    def start(self, group=0):
        try:
            if not group >= 0:
                raise IndexError
            return self._start[group]
        except IndexError:
            raise IndexError('no such group') 
    
    def end(self, group=0):
        try:
            if not group >= 0:
                raise IndexError
            return self._end[group]
        except IndexError:
            raise IndexError('no such group') 
    
    def groups(self):
        return self._groups

class Stack(list):
    
    def push(self, item):
        """
        Push an item onto the top of the stack
        """
        self.append(item)
        
    def peek(self):
        """
        Return the object at the top of the stack without removing it
        """
        if self.empty():
            raise IndexError('peek in empty list')
        return self[-1]
    
    def empty(self):
        """
        Return True if the stack is empty
        """
        return len(self) == 0
        
def metaclass(meta, base=object):
    """
    Return a new class with parent class C{base} and metaclass C{meta}.
    This works in both python 2 and 3.
    """
    return meta("NewBase", (base,), {})

def deepcopy(obj, recursion=100000):
    """
    Perform a deep copy of obj using cPickle. Faster than copy.deepcopy()
    for large objects.
    
    @param obj: the object to copy
    @return: a deep copy of obj
    @param recursion: maximum recursion limit
    @type recursion: int
    """
    import cPickle
    
    current = sys.getrecursionlimit()
    sys.setrecursionlimit(recursion)
    
    tmp = cPickle.dumps(obj, cPickle.HIGHEST_PROTOCOL)
    copy = cPickle.loads(tmp)
    
    sys.setrecursionlimit(current)
    return copy    

def _deserialize_enum(enum, name):
    return getattr(enum, name)

class EnumValueError(ValueError):
    pass

class EnumMemberError(AttributeError):
    pass
    
class EnumItem(object):
    
    def __init__(self, name, value, enum):
        self.__name = name
        self.__value = value
        self.__container = enum
    
    def __reduce__(self):
        return (_deserialize_enum, (self.enum, self.name))
    def __deepcopy__(self, memo):
        return self
    def __copy__(self):
        return self         
    def __repr__(self):
        return '{0}'.format(self.__name, self.__value)
    def __str__(self):
        return str(self.__value)
    def __int__(self):
        return int(self.__value)
    def __float__(self):
        return float(self.__value)
    def __hash__(self):
        return hash(self.__value)
    def __cmp__(self, other):
        if isinstance(other, EnumItem):
            return cmp(self.__value, other.value)
        else:
            return cmp(self.__value, other)

    @property
    def name(self):
        return self.__name

    @property
    def value(self):
        return self.__value

    @property
    def enum(self):
        return self.__container
    
    def _attach(self, enum):
        assert getattr(enum, self.__name) is self
        self.__container = enum
        
class EnumMeta(type):
    """
    Metaclass for enum types.
    """

    def __new__(cls, name, bases, items):
        
        enumdict = {}
        
        enumdict['_name'] = name
        enumdict['_names'] = {}
        enumdict['_namesci'] = {}
        enumdict['_values'] = {}
        enumdict['_valuesci'] = {}    

        for attribute in items:
            value = items[attribute]
                                        
            if attribute.startswith('_'):
                enumdict[attribute] = value
            else:            
                if value in enumdict['_values']:
                    raise EnumValueError('Duplicate value {0} in enum {1}.'.format(value, name))
    
                enumdict['_names'][attribute] = attribute
                enumdict['_namesci'][attribute.lower()] = attribute
                enumdict['_values'][value] = attribute
                if isinstance(value, basestring):
                    enumdict['_valuesci'][value.lower()] = attribute   
                
                enumdict[attribute] = EnumItem(attribute, value, None)      
            
        enum = super(EnumMeta, cls).__new__(cls, name, bases, enumdict)     

        for attribute in items:
            if not attribute.startswith('_'):            
                value = items[attribute]
                enum.__dict__[attribute]._attach(enum)
            
        return enum
                
    def __setattr__(self, attribute, value):
        raise AttributeError("enum types are immutable")
    
    def __repr__(self):
        items = list(self._names)[:2]
        preview = [ '{0}={1}'.format(i, getattr(self, i)) for i in items ]
        if len(preview) < len(self._names):
            preview.append('...')         
        return '<{0}: {1}>'.format(self._name, ', '.join(preview))
    
    def __str__(self):
        return repr(self)

EnumBase = metaclass(EnumMeta)
 
class enum(EnumBase):
    """
    Base class for all enumeration types. Supports both string and integer
    enumeration values. Examples:

        >>> class MolTypes(enum): DNA, RNA = range(2)
        <MolTypes: RNA=1, DNA=0>        
        >>> class MolTypes(enum): DNA=1; RNA=2
        <MolTypes: RNA=2, DNA=1>        
        >>> MolTypes.DNA
        1
        >>> MolTypes.DNA == 1
        True
        >>> int(MolTypes.DNA)
        1
        >>> repr(MolTypes.DNA)
        'DNA'
        >>> type(MolTypes.DNA)
        L{EnumItem} instance
        
    @note: The recommended way to create an enum is to define a public
           subclass of L{enum} in the global namespace of your module. Nesting
           the enum in another class for example will break pickling.
    """
    def __init__(self):
        raise TypeError("Can't instantiate static enum type {0}".format(self.__class__))

class Enum(object):
    """
    A collection of efficient static methods for working with L{enum}
    classes.
    """

    @staticmethod
    def create(classname, **members):
        """
        Dynamically create a new enum from a list of key:value pairs.
        Note that each key must be a valid python identifier, and the
        values must be unique.
        
        @param classname: class name for the new enum
        @type classname: str
        
        @note: The recommended way to create an enum is to define a public
               subclass of L{enum} in the global namespace of your module.
               You should avoid creating enums dynamically if static
               construction is possible, because dynamically created enums
               cannot be pickled.
        """
        
        return type(classname, (enum,), members)        
                
    @staticmethod
    def members(enumclass):
        """
        Return all member items of the C{enumclass}.
        
        @param enumclass: the enumeration class to traverse
        @type enumclass: type
        
        @return: a set of all enumclass members
        @rtype: frozenset
        """
        return frozenset([enumclass.__dict__[i] for i in enumclass._names])

    @staticmethod
    def names(enumclass):
        """
        Return the names of all items in the C{enumclass}.
        
        @param enumclass: the enumeration class to traverse
        @type enumclass: type
        
        @return: a set of all enumclass member names
        @rtype: frozenset        
        """
        return frozenset(enumclass._names)

    @staticmethod
    def values(enumclass):
        """
        Return all values of the C{enumclass}.
        
        @param enumclass: the enumeration class to traverse
        @type enumclass: type
        
        @return: a set of all enum values
        @rtype: frozenset        
        """
        return frozenset(enumclass._values)

    @staticmethod
    def parse(enumclass, value, ignore_case=True):
        """
        Parse C{value} as an C{enumclass} member value.
        
        @param enumclass: target L{enum} subclass
        @type enumclass: type        
        @type value: str, int, float
        @param ignore_case: if set to True, triggers case insensitive parsing
        @type ignore_case: bool
        
        @return: a member of enumclass, having such value
        @rtype: L{EnumItem}
        
        @raise EnumValueError: when such value is not found in enumclass
        """

        if isinstance(value, basestring) and ignore_case:
            values = enumclass._valuesci
            value = value.lower()
        else:
            values = enumclass._values

        if value in values:
            return enumclass.__dict__[ values[value] ]
        else:
            raise EnumValueError('No such value {0} in {1}'.format(value, enumclass))

    @staticmethod
    def parsename(enumclass, name, ignore_case=True):
        """
        Parse C{name} as a member of C{enumclass}.
        
        @param enumclass: target L{enum} subclass
        @type enumclass: type
        @param name: enumclass member name to parse
        @type name: str
        @param ignore_case: if set to True, triggers case insensitive parsing
        @type ignore_case: bool
        
        @return: a member of enumclass, having such name
        @rtype: L{EnumItem}
        
        @raise EnumValueError: when such name is not found in enumclass's members  
        """

        if isinstance(name, basestring) and ignore_case:
            names = enumclass._namesci
            name = name.lower()
        else:
            names = enumclass._names

        if name in names:
            return enumclass.__dict__[ names[name] ]
        else:
            raise EnumMemberError('No such item {0} in {1}'.format(name, enumclass))

    @staticmethod
    def tostring(item):
        """
        Return a string representation of the enum item.
        
        @param item: an enum item to be converted to a string
        @type item: L{EnumItem}
        
        @return: the value of the enum member
        @rtype: str
        """
        return item.name

    @staticmethod
    def ismember(item, enumclass):
        """
        Return True if item is a member of enumclass.
        
        @param enumclass: target enumeration type
        @type enumclass: type
        @param item: the enum item to be tested
        @type item: L{EnumItem}
        @rtype: bool        
        """
        if not issubclass(enumclass, enum):
            raise TypeError(enumclass)
        return item.enum is enumclass
    
class ItemNotFoundError(KeyError):
    pass

class InvalidKeyError(KeyError):
    pass

class DuplicateKeyError(InvalidKeyError):
    pass

class AbstractIndexer(object):
    
    @abstractmethod
    def __getitem__(self, i):
        pass
    
    @abstractmethod
    def _getinternal(self, i):
        """
        Implementing classes are expected to provide direct access to the
        requested element in the internal data structure via this method.
        """
        pass   

class BaseDictionaryContainer(AbstractIndexer):
    """ 
    Base class which defines the behavior of a read only key-value collection 
    container. 
    
    @note: Methods for editing an existing dictionary are also defined in the 
           base class, but left internal on purpose. Subclasses which are
           supposed to be write-enabled containers must provide their own 
           public methods for editing which might do some custom work and then
           finally call any of the internal methods in the base class to do the
           real data editing.

    @param items: an initialization dictionary
    @param restrict: a list of keys allowed for this dictionary
    """

    def __init__(self, items=None, restrict=None):

        self._keys = None
        self._items = OrderedDict({})

        if restrict:
            self._keys = frozenset(restrict)
        if items is not None:
            self._set_items(items)
    
    def __getitem__(self, key):
        try:
            return self._items[key]
        except KeyError:
            raise ItemNotFoundError(key)

    def _getinternal(self, key):
        return self._items[key]

    def __contains__(self, item):
        return item in self._items

    def __len__(self):
        return len(self._items)

    def __nonzero__(self):
        return len(self) > 0
    
    @property
    def length(self):
        return len(self)

    def keys(self):
        return self._items.keys()

    def __iter__(self):
        return iter(self._items)

    def __repr__(self):
        return '<{0.__class__.__name__}: {0.length} items>'.format(self)

    def _set_items(self, new_items):
        new_items = OrderedDict(new_items)

        if self._keys and not self._keys.issuperset(new_items):
            raise InvalidKeyError("One or more of the keys provided are not allowed for this collection.")

        self._items = new_items

    def _update_items(self, new_items={ }, **named_items):
        new_items = OrderedDict(new_items)

        if self._keys:
            if not set(self).issuperset(new_items) or not set(self).issuperset(named_items):
                raise ItemNotFoundError("One or more of the keys provided were not found in this collection.")

        self._items.update(new_items)
        self._items.update(named_items)
        
    def _append_item(self, key, item):
        
        if self._keys and key not in self._keys:
            raise InvalidKeyError("Key {0} is not allowed for this collection.".format(key))
        if key in self:
            raise DuplicateKeyError("Key {0} already exists in the collection.".format(key))
        self._items[key] = item        
        
    def _remove_item(self, key):
        
        if key not in self._items:
            raise ItemNotFoundError(key)
        del self._items[key]

class ReadOnlyDictionaryContainer(BaseDictionaryContainer):
    """
    This is a write-once container, which is filled with items only at object
    construction.
    """
    pass

class DictionaryContainer(BaseDictionaryContainer):
    """ 
    Write-enabled Dictionary Container. New items can be added using a public
    C{append} method. Subclasses may also override the internal C{_update}
    and C{_set} to expose them to the clients.
    """
    def __init__(self, items=None, restrict=None):

        super(DictionaryContainer, self).__init__(items, restrict)

    def append(self, key, item):
        """
        Append a new key-value to the collection.
        
        @param key: key
        @param item: value

        @raise InvalidKeyError: when the key is not allowed for this container
        @raise DuplicateKeyError: when such a key already exists 
        """
        self._append_item(key, item)

    def _set(self, new_items):
        """ 
        Set the collection to the dictionary provided in the new_items parameter.
        
        @param new_items: the new value of the dictionary container
        @type new_items: dict
        """
        self._set_items(new_items)

    def _update(self, new_items={ }, **named_items):
        """ 
        Update the collection with the dictionary provided in the C{new_items} parameter.
        Update also specific items with the values provided as keyword arguments.
        
        @param new_items: a dictionary that provides new values for certain keys
        @type new_items: dict
        """
        self._update_items(new_items, **named_items)
        
    def _remove(self, item):
        """ 
        Delete C{item}.
        
        @param item: key to remove
        """
        self._remove_item(item)        

class CollectionIndexError(IndexError):
    pass

class BaseCollectionContainer(AbstractIndexer):
    """ 
    Base class which defines the behavior of a read-only collection container.

    @note: Methods for editing an existing collection are also defined in the 
           base class, but left internal on purpose. Subclasses which are
           supposed to be write-enabled containers must provide their own 
           public methods for editing which might do some custom work and then
           finally call any of the internal methods in the base class to do the
           real data editing.
               
    @param items: initialization list
    @type items: list
    @param type: if defined, makes the container type-safe
    @type type: type
    @param start_index: the index of the zero element
    @type start_index: int
    """

    def __init__(self, items=None, type=None, start_index=0):

        self._items = [ ]
        self._type = type

        if not (isinstance(start_index, int) or start_index >= 0):
            raise ValueError('start_index must be a positive integer.')

        self._start = start_index

        if items is not None:
            self._update_items(items)

    def __fix(self, i):
        if i >= 0:
            if i < self._start:
                return None
            return i - self._start
        else:
            return i

    def __getitem__(self, i):
        try:
            if isinstance(i, slice):
                sl = slice(self.__fix(i.start), self.__fix(i.stop), i.step)
                return self._items[sl]
            else:
                if 0 <= i < self._start:
                    raise IndexError
                return self._items[self.__fix(i)]

        except IndexError:
            if len(self) > 0:
                raise CollectionIndexError('List index {0} out of range: {1}..{2}'.format(i, self.start_index, self.last_index))
            else:
                raise CollectionIndexError('List index out of range. The collection is empty.')
    
    def _getinternal(self, i):
        return self._items[i]

    def __contains__(self, item):
        return item in self._items

    def __len__(self):
        return len(self._items)

    def __nonzero__(self):
        return len(self) > 0

    @property
    def length(self):
        return len(self)

    @property
    def start_index(self):
        return self._start

    @property
    def last_index(self):
        length = len(self._items)
        if length > 0:
            return length + self._start - 1
        else:
            return None

    def __iter__(self):
        return iter(self._items)

    def __repr__(self):
        return '<{0.__class__.__name__}: {0.length} items>'.format(self)

    def _append_item(self, item):

        if self._type:
            if not isinstance(item, self._type):
                raise TypeError("Item {0} is not of the required {1} type.".format(item, self._type.__name__))
        self._items.append(item)

        return len(self) + self._start - 1

    def _update_items(self, new_items):

        if self._type:
            for a in new_items:
                if not isinstance(a, self._type):
                    raise TypeError("Item {0} is not of the required {1} type.".format(a, self._type.__name__))
        self._items = list(new_items)
        
    def _sort(self):
        self._items.sort()

class ReadOnlyCollectionContainer(BaseCollectionContainer):
    """
    This is a write-once container, which is filled with items only at object
    construction.
    """
    pass

class CollectionContainer(BaseCollectionContainer):
    """ 
    Write-enabled Collection Container. New items can be added using a public
    C{append} method. Subclasses may also override the internal C{_update}
    to expose it to the clients. 
    """

    def __init__(self, items=None, type=None, start_index=0):

        super(CollectionContainer, self).__init__(items, type, start_index)

    def append(self, item):
        """
        Append a new item to the collection.
        
        @param item: the new item to append

        @return: the index of the new element        
        @rtype: int
        
        @raise TypeError: when the container is type-safe and item has an 
                          incorrect type
        """
        return self._append_item(item)

    def _update(self, new_items):
        """ 
        Set the collection to the list provided in the new_items parameter.
        
        @param new_items: a list providing the new items for the container
        @type new_items: list
        """
        self._update_items(new_items)

class AbstractContainer(object):
    """
    Defines the behavior of a high-level object, which can hold an array of
    elements. Implementing classes automatically provide iterable and index/key
    based access to those objects in a read-only encapsulated manner.
    
    This is an abstract class with an abstract property C{_children}. Subclasses
    must override this property. The overridden implementation is usually 
    extremely simple - you just need to return a reference to an iterable and
    subscriptable object, containing the elements.
    """
    
    __metaclass__ = ABCMeta
    
    def __getitem__(self, i):
        return self._children[i]
    
    def __iter__(self):
        for i in self._children:
            yield i
    
    @abstractproperty
    def _children(self):
        pass
    
class Container(AbstractContainer):
    """
    Generic implementation of L{AbstractContainer}. Provides an easy way to
    encapsulate class properties that behave like read-only collections or
    dictionaries. This container is super lightweight and completely dynamic:
    it serves as a proxy to an internal list or dict and stores no data in its
    own instance. Owning classes are therefore free to modify those internal
    data structures, which provides advantages over using 
    L{ReadOnlyCollectionContainer}s. 
    """
    
    def __init__(self, data):
        self._data = data

    @property
    def _children(self):
        return self._data        
    
class AbstractNIContainer(AbstractContainer):
    """
    Same as the L{AbstractContainer}, but provides access to the child elements
    through C{AbstractNativeIndexer._getinternal} instead of the standard 
    __getitem__.
    
    Therefore, the C{self._children} property must return an object which 
    implements L{AbstractIndexer}. 
    """
    def __getitem__(self, i):
        return self._children._getinternal(i)
                
try:        
    from collections import OrderedDict as _dict
    class OrderedDict(_dict):
        pass
except ImportError:
    import UserDict
    
    class OrderedDict(dict, UserDict.DictMixin):
    
        def __init__(self, *args, **kwds):
            if len(args) > 1:
                raise TypeError('expected at most 1 arguments, got %d' % len(args))
            try:
                self.__end
            except AttributeError:
                self.clear()
            self.update(*args, **kwds)
    
        def clear(self):
            self.__end = end = []
            end += [None, end, end]         # sentinel node for doubly linked list
            self.__map = {}                 # key --> [key, prev, next]
            dict.clear(self)
    
        def __setitem__(self, key, value):
            if key not in self:
                end = self.__end
                curr = end[1]
                curr[2] = end[1] = self.__map[key] = [key, curr, end]
            dict.__setitem__(self, key, value)
    
        def __delitem__(self, key):
            dict.__delitem__(self, key)
            key, prev, next = self.__map.pop(key)
            prev[2] = next
            next[1] = prev
    
        def __iter__(self):
            end = self.__end
            curr = end[2]
            while curr is not end:
                yield curr[0]
                curr = curr[2]
    
        def __reversed__(self):
            end = self.__end
            curr = end[1]
            while curr is not end:
                yield curr[0]
                curr = curr[1]
    
        def popitem(self, last=True):
            if not self:
                raise KeyError('dictionary is empty')
            if last:
                key = reversed(self).next()
            else:
                key = iter(self).next()
            value = self.pop(key)
            return key, value
    
        def __reduce__(self):
            items = [[k, self[k]] for k in self]
            tmp = self.__map, self.__end
            del self.__map, self.__end
            inst_dict = vars(self).copy()
            self.__map, self.__end = tmp
            if inst_dict:
                return (self.__class__, (items,), inst_dict)
            return self.__class__, (items,)
    
        def keys(self):
            return list(self)
    
        setdefault = UserDict.DictMixin.setdefault
        update = UserDict.DictMixin.update
        pop = UserDict.DictMixin.pop
        values = UserDict.DictMixin.values
        items = UserDict.DictMixin.items
        iterkeys = UserDict.DictMixin.iterkeys
        itervalues = UserDict.DictMixin.itervalues
        iteritems = UserDict.DictMixin.iteritems
    
        def __repr__(self):
            if not self:
                return '%s()' % (self.__class__.__name__,)
            return '%s(%r)' % (self.__class__.__name__, self.items())
    
        def copy(self):
            return self.__class__(self)
    
        @classmethod
        def fromkeys(cls, iterable, value=None):
            d = cls()
            for key in iterable:
                d[key] = value
            return d
    
        def __eq__(self, other):
            if isinstance(other, OrderedDict):
                return len(self)==len(other) and self.items() == other.items()
            return dict.__eq__(self, other)
    
        def __ne__(self, other):
            return not self == other

