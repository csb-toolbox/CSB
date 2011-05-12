"""
Read, query and update textual tables via flexible SQL interface.

L{Table}s can be created and populated with data from scratch, built from TSV
files, 2D lists or other tables. Once the data is loaded in memory, each
storage operation on the table object is delegated via bridge to an SQL
storage backend (by default this is SQLite). However the table uses the backend
only as a temp storage to ensure maximum portability of the data. Tables can be
stored persistently as text (TSV) files and then loaded back in memory when
needed.

These Tables can be queried and updated in a vast number of ways; each query
returns a new L{Table}:  

    1. Using slice expressions. The general form of a slice expression is
    C{[rows, columns]}, where C{rows} can be:
    
    
        - a row index, 0-based, e.g. C{5}
        - a tuple of row indices, e.g. C{(1, 3, 6)}
        - a standard Python slice, e.g. C{1:3} or C{:5} or C{:}
        - omitted (means: all rows)
        
    and C{columns} can be:
    
        - a column index, 0-based, e.g. C{5}    
        - a tuple of columns indices, 0-based
        - a column name, e.g. C{'TmScore'}        
        - a tuple of column names, e.g. C{('ID', 'TmScore')}
        - a standard Python slice using column indices
        - a slice using column names, e.g. C{'ID':'TM'} or C{:'TM'} or C{:}
        - omitted (means: all columns)
    
    2. Using query expressions, for example:
    
    
    >>> table.where('ID').between(1, 5).select('TmScore', 'RMSD')
    Table ('TmScore', 'RMSD')
    
    >>> table.where('ID').between(1, 5).update('RMSD', 0.2)
    Table (the same table)
        
    3. With SQL queries:
    
    
    >>> t.query(r'''SELECT  ColumnB * ColumnA AS ComputedValue
                    FROM    {0.name}
                    WHERE   ColumnC IN ({1}, {1})'''.format(t, Predicate.PH),
                [12, 55])
    iterable
    
The data contained in a Table can be extracted in several ways:

    - if you need a single (scalar) value -- with the C{table[row, column]}
    indexing expression or with the dedicated C{table.scalar(row, column)} method.
    - by treating the table as an iterator; each cycle will then yield a L{DataRow}
    object
    - with text (TSV) serialization: simply call C{table.dump(file)}.
"""

import os
import csb.io
import cStringIO
import __builtin__

from abc import ABCMeta, abstractmethod, abstractproperty


class RepositoryImp(object):
    """
    Abstract SQL backend interface. Defines a number of platform-specific
    operations, that each concrete backend implementor must provide. 
    """
    
    __metaclass__ = ABCMeta

    def __init__(self, tablename):
        self._table = tablename
    
    @abstractproperty
    def pk(self):
        pass
    
    @property
    def table(self):
        return self._table    
    
    def query(self, sql, params=None):
        """
        Execute a native SQL query against the backend, as-is.
        
        @param sql: SQL query
        @type sql: str
        @param params: query bound parameters, if any
        @type params: tuple
        
        @return: data reader (2D iterable)
        """
        raise NotImplementedError()
    
    @abstractmethod
    def execute(self, expression):
        """
        Perform a select operation given L{expression}.
        
        @type expression: L{Expression}    
        @return: data reader (2D iterable)        
        """
        pass
    
    @abstractmethod
    def update(self, expression):
        """
        Perform an update operation given L{expression}.
        
        @type expression: L{Expression}    
        @return: void        
        """
        pass
    
    @abstractmethod
    def insert(self, row):
        """
        Insert a new tuple in the table.
        
        @type row: tuple    
        @return: void        
        """
        pass    
    
    @abstractmethod
    def create(self, metadata):
        """
        Create a table given L{metadata}.
        
        @type metadata: tuple of L{ColumnInfo}    
        @return: void        
        """        
        pass
    
    @abstractmethod
    def close(self):
        """
        Perform cleanup (e.g. close connections).
        """
        pass 

class InvalidColumnError(KeyError):
    pass

class UnsupportedTypeError(ValueError):
    pass

class SQLiteRepository(RepositoryImp):
    """
    SQLite-based concrete repository implementor.
    This is the default L{Table} backend. 
    """

    PK = 'ROWID'
    TYPES = { int: 'INT', long: 'BIGINT', float: 'REAL', str: 'VARCHAR' }
        
    def __init__(self, tablename):
        import sqlite3
        
        self._conn = sqlite3.connect(':memory:')
        self._pk = SQLiteRepository.PK
        
        super(SQLiteRepository, self).__init__(tablename)
        
    @property
    def pk(self):
        return self._pk
    
    def query(self, sql, params=None):
        
        sql = sql.replace(Predicate.PH, '?')
        if not params:
            params = []
        
        return self._conn.execute(sql, params).fetchall() 
        
    def execute(self, exp):

        query = 'SELECT  {0}\nFROM    {1}\n'.format(', '.join(exp.select), self.table)
        
        if exp.where:
            predicate = str(exp.predicate).replace(Predicate.PH, '?')
            query += 'WHERE   {0} {1}\n'.format(exp.where, predicate)

        query += 'ORDER BY  {0} ASC\n'.format(self.pk)                
        return self.query(query, exp.params)
    
    def update(self, exp):

        params = [exp.data]
        query = 'UPDATE  {0}\n  SET   {1} = ?\n'.format(self.table, exp.select[0])
        
        if exp.where:
            predicate = str(exp.predicate).replace(Predicate.PH, '?')
            query += 'WHERE   {0} {1}\n'.format(exp.where, predicate)
            if exp.params:
                params.extend(list(exp.params))            
        
        return self.query(query, params)
    
    def insert(self, row):

        row = list(row)
        params = ','.join(['?' for dummy in row])
        query = 'INSERT INTO {0} VALUES({1})'.format(self.table, params)
        self.query(query, row)
        
    def create(self, metadata):
        
        cols = []
        
        for ci in metadata:
            type = self._gettype(ci.type)
            cols.append('{0}  {1}'.format(ci.name, type))
            
        statement = 'CREATE TABLE {0} (\n    {1}\n);'.format(self.table, ',\n    '.join(cols))
        
        self._conn.execute(statement)
    
    def _gettype(self, type):
        try:
            return SQLiteRepository.TYPES[type]
        except KeyError:
            raise UnsupportedTypeError(type)
    
    def close(self):
        try:
            return self._conn.close()
        except:
            pass

class ColumnInfo(object):
    """
    Holder object for column metadata.
    
    @param name: column name
    @type name: str
    @param type: column data type (Python)
    @type type: type
    """
    
    def __init__(self, name, type):
        self._name = name
        self._type = type
    
    @property
    def name(self):
        return self._name
    
    @property
    def type(self):
        return self._type    
    
    def __str__(self):
        return '{0.name}:{0.type.__name__}'.format(self)
    
    def copy(self):
        """
        @return: a deep copy of C{self}
        """
        return ColumnInfo(self.name, self.type)
    
class DataRow(object):
    """
    Represents a table data row. This is basically what a table iterator
    yields for each row in a table. Provides both index (position) and
    column name-based access to the data.
    """
    
    def __init__(self, columns, row):
                
        self._row = tuple(row)
        self._columns = {}
        for i, c in enumerate(columns):
            self._columns[c] = i
        assert len(self._columns) == len(self._row)            

    def __iter__(self):
        return iter(self._row)
        
    def __getitem__(self, i):
        
        if isinstance(i, int):
            return self._row[i]
        else:
            return self._row[self._columns[i]]
        
    def __len__(self):
        return len(self._row)
    
    def __repr__(self):
        return '{0}: {1}'.format(self.__class__.__name__, repr(self._row))
    
    def __str__(self):
        return self.dump()
    
    def dump(self, delimiter='\t'):
        """
        Dump the row as a string.
        
        @param delimiter: column separator (defaults to tab)
        @type delimiter: str 
        """
        return delimiter.join(map(str, self._row))
    
    @property
    def columns(self):
        return tuple(self._columns)
                
class Table(object):
    """
    Build and query a TSV Table. See the documentation of L{csb.io.tsv} for
    details and examples.
    
    @param definition: column definition string: C{ColA:typeA colB:typeB ...},
                       where C{ColN} is a column name and C{typeN} is one of the
                       base Python data types: str, int, long, float.
                       Alternatively, the table definition may be specified
                       directly as a list of metadata objects.
    @type definition: str, tuple of L{ColumnInfo}
    @param name: name of the table on the SQL backend. Useful when you need to
                 execute native SQL queries against the table.
    @type name: str
    @param backend: table backend storage engine. This must be a proper
                    L{RepositoryImp} bridge implementor.
    @type backend: type (reference to a L{RepositoryImp} subclass)
    
    @raise UnsupportedTypeError: when an unsupported type is used in the table
                                 C{definition}
    @raise ValueError: if the C{definition} is not valid
    """
    
    """
    Table header string, used when saving and restoring TSV files.
    """
    HEADER = '# @TSV '
    
    def __init__(self, definition, name='TSV', backend=SQLiteRepository):
        
        if not issubclass(backend, RepositoryImp):
            raise TypeError('The Table Backend must be a Repository Implementor')

        self._name = name
        self._backend = backend
        self._imp = backend(name)
        
        try:
            if isinstance(definition[0], ColumnInfo):
                self._metadata = [ c.copy() for c in definition ]
            else:
                if isinstance(definition, basestring):
                    definition = [ (d.split(':')[0], getattr(__builtin__, d.split(':')[1])) for d in definition.split() ]
                self._metadata = [ ColumnInfo(c[0], c[1]) for c in definition ]
            if len(self._metadata) < 1:
                raise ValueError()
        except UnsupportedTypeError:
            raise            
        except (TypeError, IndexError, ValueError, NameError, AttributeError):
            raise ValueError('Invalid table definition')
        
        self._imp.create(self._metadata)
        
    @staticmethod
    def from_tsv(tsv, definition=None, delimiter='\t', skip=0, name='TSV', backend=SQLiteRepository):
        """
        Table factory: build a L{Table} from a TSV file.
        
        @param tsv: TSV path and filename. This can be either a conventional
                    TSV/CSV file, or a file created with C{table.dump(tsv)}
        @type tsv: str
        @param definition: table column definition (see L{Table}). If defined,
                           this parameter will determine the structure of the
                           table. Otherwise, the table definition will be
                           extracted from the TSV header. If the file contains
                           no TSV header, this parameter is mandatory.
        @type definition: str, tuple of L{ColumnInfo}                           
        @param delimiter: column separator used in the file
        @type delimiter: str
        @param skip: skip the first N number of rows (the header can still be
                     extracted from those however)
        @type skip: int
        
        @rtype: L{Table}
        
        @raise ValueError: if neither a table C{definition} is provided,
                           nor the C{tsv} file has a header line 
        """

        if not definition:
            with open(tsv) as tsvfile:
                for line in tsvfile:
                    if line.startswith(Table.HEADER):
                        definition = line[ len(Table.HEADER) : ]
                        
        if not definition:
            raise ValueError('No header definition found')
                    
        table = Table(definition, name=name, backend=backend) 
        
        with open(tsv) as tsvfile:
            for i, line in enumerate(tsvfile, start=1):
                if (skip and i <= skip) or line.startswith(Table.HEADER):
                    continue
                table.insert(line.rstrip(os.linesep).split(delimiter))
        
        return table
    
    @staticmethod
    def from_iterable(iterable, definition, name='TSV', backend=SQLiteRepository):
        """
        Table factory: build a L{Table} from a 2D iterable/data reader.
        
        @param iterable: data container
        @type iterable: iterable (2D)
        @param definition: table column definition (see L{Table}).
        @type definition: str, tuple of L{ColumnInfo}

        @rtype: L{Table}
        """        
        table = Table(definition, name=name, backend=backend)
         
        for row in iterable:
            table.insert(list(row))
        
        return table
    
    @staticmethod
    def from_table(table, data=False, name='TSV', backend=SQLiteRepository):
        """
        Table factory: build a L{Table} with the definition of another L{Table}.
        
        @param table: template table
        @type table: L{Table}
        @param data: if True, also copy the data from the source C{table}
        @type data: bool 

        @rtype: L{Table}
        """
        if data:
            return Table.from_iterable(table, table._metadata, name=name, backend=backend)
        else:
            return Table(table._metadata, name=name, backend=backend)            
        
    def __del__(self):
        self._imp.close()
        
    def __iter__(self):
        exp = Expression(self.columns)
        for row in self._imp.execute(exp):
            yield DataRow(self.columns, row)
            
    def __getstate__(self):
        
        temp = cStringIO.StringIO()
        self.dump(temp)
        return temp.getvalue()
        
    def __setstate__(self, state):

        with csb.io.TempFile() as temp:
            temp.write(state)
            temp.flush()
            clone = Table.from_tsv(temp.name)
        
        self.__init__(definition=clone._metadata, name=clone.name, backend=clone._backend)

        for row in clone:
            self.insert(row)
            
    def __setitem__(self, i, value):

        exp = self._interpret(i)
        
        if len(exp.select) != 1:
            raise NotImplementedError('single-column expression expected')
        if hasattr(value, '__iter__'):
            raise NotImplementedError("single-value assignment expected")     
        
        exp.data = value
        self._update(exp)
            
    def __getitem__(self, i):
                
        exp = self._interpret(i)
        
        if exp.scalar:
            return self.scalar(i[0], exp.select[0])
        else:
            return self._execute(exp)
    
    def _interpret(self, i):
        """
        Parse a table slice and convert it into an L{Expression}.
        @rtype: L{Expression}
        """
        
        if not hasattr(i, '__iter__'):
            i = [i, slice(None, None)]        
        else:
            i = list(i)
            
        if len(i) not in (1, 2):
            raise ValueError('Tables are only 2 dimensional')
        if len(i) == 1:
            i.append(slice(None, None))
        
        exp = Expression(self.columns)
        columns = self._getcols(i[1])
        if len(columns) < 1:
            raise ValueError('Column slices must return at least one column')
        exp.select = columns
        exp.where = self.pk
        
        if isinstance(i[0], int):
            self._checkrow(i[0])
            if len(columns) == 1 and isinstance(i[1], int):
                exp.scalar = True
            exp.predicate = Equals(i[0] + 1)
            
        elif hasattr(i[0], '__iter__'):
            params = list(i[0])
            self._checkrow(params)
            params = map(lambda x: x+1, params)
            exp.predicate = In(params)
            
        elif isinstance(i[0], slice):
            
            sl = i[0]
            if sl.step is not None:
                raise NotImplementedError('Row slice steps are not supported')              
            
            if sl == slice(None, None):
                exp.where = None
            elif sl.start is None:
                self._checkrow(sl.stop)                
                exp.predicate = Lower(sl.stop + 1)
            elif sl.stop is None:
                self._checkrow(sl.start)                
                exp.predicate = GreaterOrEquals(sl.start + 1)
            else:
                self._checkrow([sl.start, sl.stop])       
                exp.predicate = Between(sl.start + 1, sl.stop)
                
        else:
            raise TypeError("Can't handle row slice expression: {0}".format(i[0]))
        
        return exp
        
    def _checkrow(self, i):
        
        if isinstance(i, int):
            if i < 0:
                raise NotImplementedError('Negative row indices are not supported')
        elif hasattr(i, '__iter__'):
            for j in i:
                self._checkrow(j)
        else:
            raise TypeError(i)
             
    def _getcols(self, spec, ifnull=None):
        
        columns = list(self.columns)
        
        if spec is None and ifnull is not None:
            return [ifnull]
        
        elif isinstance(spec, int):
            try:
                return [columns[spec]]
            except:
                raise IndexError('Column {0} out of range'.format(spec))
            
        elif isinstance(spec, basestring):
            if spec in columns:
                return [spec]
            else:
                raise InvalidColumnError(spec)
        
        elif isinstance(spec, slice):
            start = self._getcols(spec.start, columns[0])
            start = columns.index(start[0])

            end = self._getcols(spec.stop, columns[-1])
            end = columns.index(end[0])
            if spec.stop is None:
                end += 1
                
            return [columns[i] for i in range(start, end, spec.step or 1)]
        
        elif hasattr(spec, '__iter__'):
            return [self._getcols(i)[0] for i in spec]
        
        else:
            raise TypeError("Can't handle column slice expression: {0}".format(spec))

    @property
    def name(self):
        return self._name

    @property
    def columns(self):
        return tuple(i.name for i in self._metadata)
    
    @property
    def pk(self):
        return self._imp.pk
    
    def dump(self, file):
        """
        Dump the table in a file.
        
        @param file: destination stream or filename
        @type file: file (stream) or str (filename)
        """

        import csb.io
                
        with csb.io.EntryWriter(file, close=False) as out:
            
            definition = map(str, self._metadata)
            out.write(Table.HEADER)
            out.writeall(definition, delimiter=' ')
            out.write(csb.io.NEWLINE)
            
            for row in self:
                out.writeline(row.dump(delimiter='\t'))
    
    def query(self, sql, params=None):
        """
        Execute a native SQL query against the storage engine.
        
        @param sql: SQL query text. May contain parameter binding placeholders
                    (see L{Predicate.PH}). The SQL dialect of the query depends
                    on the SQL C{backend} being used by the table.
        
        @return: native data reader
        @rtype: iterable (2D)        
        """
        return self._imp.query(sql, params)
    
    def insert(self, row):
        """
        Insert a new row in the table.
        
        @param row: a tuple of the appropriate length
        @type row: tuple 
        """
        self._imp.insert(row)
    
    def _project(self, columns):
        
        metadata = dict((c.name, c) for c in self._metadata)
        try:
            return [metadata[cn].copy() for cn in columns]
        except KeyError as ke:
            raise InvalidColumnError(ke.message)
        except:
            raise
        
    def _execute(self, exp):
        
        newdef = self._project(exp.select)        
        reader = self._imp.execute(exp)
        return Table.from_iterable(reader, newdef, name=self.name, backend=self._backend)
    
    def _update(self, exp):
        
        if exp.select[0] not in self.columns:
            raise InvalidColumnError(exp.select[0])
        
        self._imp.update(exp)
        return self
        
    def where(self, column):
        """
        @param column: column name        
        @type column: str
        @raise InvalidColumnError: when an invalid column is requested
        """
        exp = Expression(self.columns)
        return Where(self, exp, column)
    
    def select(self, *columns):
        """
        @return: a new L{Table}
        
        @param columns: column names; defaults to all columns
        @type columns: str, tuple of str
        @raise InvalidColumnError: when an invalid column is requested
        """        
        columns = Expression.array(columns)
        
        exp = Expression(self.columns)
        exp.select = columns
        
        return self._execute(exp)
    
    def update(self, column, value):
        """
        Update C{column} for all rows in the table.
        
        @param column: column to update (name)        
        @type column: str
        @param value: new column value
        @raise InvalidColumnError: when an invalid column is referenced        
        """          
        exp = Expression(self.columns)
        exp.select = [column]
        exp.data = value
        
        return self._update(exp)
    
    def scalar(self, row=None, column=None):
        """
        @return: a scalar value at the specified row and column.

        @param row: row index; if not specified - take the first row
        @type row: int
        @param column: column name; if not specified - take the first
        @type column: str
        
        @raise IndexError: when an invalid row is requested
        @raise InvalidColumnError: when an invalid column is requested        
        """
        
        if row is None:
            row = 0
        row += 1
        if column is None:
            column = self.columns[0]
        elif column not in self.columns:
            raise InvalidColumnError(column)
            
        exp = Expression(self.columns)
        exp.select = [column]
        exp.where = self.pk
        exp.predicate = Equals([row])
        
        reader = self._imp.execute(exp)
        if len(reader) > 0:
            return reader[0][0]
        else:
            raise IndexError()

class Expression(object):
    """
    Metadata container: represents a table select or update expression.
    """
    
    def __init__(self, columns):
        
        self._table = None
        self._columns = []
        
        self._columns = list(columns)
        self._select = []
        self._where = None
        self._predicate = None
        self.data = None
        self.scalar = False
        
        self.select = '*'
    
    @staticmethod
    def array(args):
        if len(args) == 1 and hasattr(args[0], '__iter__'):
            args = args[0]
        return list(args)
    
    @property
    def all(self):
        return tuple(self._columns)
    
    @property
    def params(self):
        if self.where and self.predicate:
            return self.predicate.params
        else:
            return None
    
    @property
    def select(self):
        return self._select
    @select.setter
    def select(self, value):
        self._select = []
        if not value:
            self._select = list(self.all)
        else:        
            if not hasattr(value, '__iter__'):
                value = [value]        
            for i in value:
                if i == '*':
                    self._select.extend(self.all)
                else:
                    if i not in self._columns:
                        raise InvalidColumnError(i)
                    self._select.append(i)
            
    @property
    def where(self):
        return self._where
    @where.setter
    def where(self, value):
        if not value:
            self._where = None
            self._predicate = None
        else:
            self._where = value
    
    @property
    def predicate(self):
        return self._predicate
    @predicate.setter
    def predicate(self, value):
        if not value:
            self._where = None
            self._predicate = None
        else:
            self._predicate = value
      
class Step(object):
    
    def __init__(self, table, expression):

        self._table = table        
        self._expression = expression
        
    @property
    def table(self):
        return self._table
    
    @property
    def expression(self):
        return self._expression      
    
class Where(Step):
    
    def __init__(self, table, expression, column):
        
        if column not in table.columns and column != table.pk:
            raise InvalidColumnError(column) 
        
        expression.where = column
        super(Where, self).__init__(table, expression)      
    
    def in_(self, *values):
        return Operator(self.table, self.expression, In(values))
    
    def between(self, start, end):
        return Operator(self.table, self.expression, Between(start, end))
    
    def equals(self, value):
        return Operator(self.table, self.expression, Equals(value))
    
    def greater(self, value):
        return Operator(self.table, self.expression, Greater(value))

    def lower(self, value):
        return Operator(self.table, self.expression, Lower(value))
        
class Operator(Step):

    def __init__(self, table, expression, predicate):
        
        expression.predicate = predicate
        super(Operator, self).__init__(table, expression)  
        
    def select(self, *columns):
        """
        @return: a new L{Table}
        
        @param columns: column names; defaults to all columns
        @type columns: str, tuple of str
        @raise InvalidColumnError: when an invalid column is requested
        """        
        exp = self.expression
        exp.select = columns

        return self.table._execute(exp)
    
    def update(self, column, value):
        """
        Update C{column} for all rows in the table.
        
        @param column: column to update (name)        
        @type column: str
        @param value: new column value
        @raise InvalidColumnError: when an invalid column is referenced
        """        
        exp = self.expression
        exp.select = [column]
        exp.data = value
        
        return self.table._update(exp)
    
class Predicate(object):
    
    __metaclass__ = ABCMeta
    
    PH = '?'
    
    def __init__(self, params):
        
        self._params = []
        
        if not hasattr(params, '__iter__'):
            params = [params]
        
        for p in params:
            if hasattr(p, '__iter__'):
                self._params.extend(p)
            else:
                self._params.append(p)
                
        self._validate()

    @property
    def params(self):
        return tuple(self._params)
    
    def _validate(self):
        
        if len(self._params) < 1:
            raise ValueError('{0} predicate with no params'.format(self.__class__.__name__))
    
    @abstractproperty
    def sql(self):
        pass
    
    def __str__(self):
        return self.sql
    
class In(Predicate):
    
    @property
    def sql(self):
        p = [Predicate.PH for dummy in self.params]
        return 'IN ({0})'.format(', '.join(p))

class Equals(Predicate):
        
    @property
    def sql(self):
        return '= {0}'.format(Predicate.PH)
    
class Greater(Predicate):
        
    @property
    def sql(self):
        return '> {0}'.format(Predicate.PH)    

class GreaterOrEquals(Predicate):
        
    @property
    def sql(self):
        return '>= {0}'.format(Predicate.PH)

class Lower(Predicate):
        
    @property
    def sql(self):
        return '< {0}'.format(Predicate.PH)    

class LowerOrEquals(Predicate):
        
    @property
    def sql(self):
        return '<= {0}'.format(Predicate.PH)
            
class Between(Predicate):

    def __init__(self, start, end):
        super(Between, self).__init__([start, end])           
    @property
    def sql(self):
        return 'BETWEEN {0} AND {0}'.format(Predicate.PH)
