
import os
import sys
import numpy
import types
import csb.io
import csb.io.tsv
import csb.test as test
import csb.bio.structure

from csb.bio.io.wwpdb import LegacyStructureParser
# Workaround I just need some class  which I can abuse to create a graph
from csb.statistics.pdf import Gamma as Node


@test.unit
class TestShell(test.Case):
    
    def setUp(self):
        
        super(TestShell, self).setUp()
        
        cmd = '{0.python} -c "print({0.text})"'
        
        self.python = sys.executable.replace("\\", "/") # "\" fails on win7
        self.text = "123"
        self.output = csb.io.Shell.run(cmd.format(self))

        if not self.python:
            self.skipTest("Can't get interpreter's path")
        
    def testSTDOUT(self):
        self.assertEqual(self.output.stdout.strip(), self.text)
        
    def testSTDERR(self):
        self.assertEqual(self.output.stderr, '')
        
    def testExitCode(self):
        self.assertEqual(self.output.code, 0)  
        
        
@test.unit
class TestTable(test.Case):

    def setUp(self):
        super(TestTable, self).setUp()
        
        self.tsv = self.config.getTestFile('csb.tsv')
        self.table = csb.io.tsv.Table.from_tsv(self.tsv)
        self.t = csb.io.tsv.Table.from_table(self.table, data=True)
        
    def testIteration(self):
        
        data = list(self.table)
        self.assertEqual(len(data), 3)
        
    def testArrayCompatibility(self):
        
        array = numpy.array(self.table)
        self.assertEqual(array.shape, (3, 3))
        
        t = csb.io.tsv.Table('a:int b:float')
        t.insert([1, 1.2])
        array = numpy.array(t)
        self.assertEqual(list(array[0]), [1, 1.2])        

    def testDataAccess(self):
        
        firstrow = list(self.table)[0]
        
        self.assertEqual(firstrow[1], 11.1)    
        self.assertEqual(firstrow['ID'], 11)
        self.assertEqual(firstrow['B'], 'Row eleven')
        
        self.assertEqual(firstrow.columns, ('ID', 'A', 'B'))
        self.assertEqual(firstrow.dump('|'), '11|11.1|Row eleven')
        
        self.assertRaises(IndexError, lambda:firstrow[4])
        self.assertRaises(KeyError, lambda:firstrow['Missing'])        
        
    def testSQLDataAccess(self):

        reader = self.table.query('SELECT A, ID FROM TSV WHERE ID > ?', [0])
        data = list(reader)
        firstrow = data[0]
                
        self.assertEqual(len(data), 3)
        self.assertEqual(firstrow[0], 11.1)
        self.assertEqual(firstrow[1], 11)
        
    def testInsert(self):
        
        self.t.insert([333, 333.0, '3 3 3']) 
        self.assertEqual(self.t[3, 2], '3 3 3')
        
        self.assertRaises(Exception, lambda:self.t.insert([1]))
        
    def testUpdate(self):

        self.assertRaises(csb.io.tsv.InvalidColumnError, lambda:self.t.where('ID').equals(11).update('Missing', 0))            
        self.t.where('ID').equals(11).update('ID', 133)
        self.assertEqual(self.t[0, 'ID'], 133)
        self.assertEqual(self.t[1, 'ID'], 12)
                
        self.assertRaises(csb.io.tsv.InvalidColumnError, lambda:self.t.update('Missing', 0))
        self.assertRaises(csb.io.tsv.InvalidColumnError, lambda:self.t[0, 'Missing'])        
        self.t.update('ID', 13)
        
        for r in self.t[:, 'ID']:
            self.assertEqual(r['ID'], 13)                    
            
        self.t[0, 'ID'] = 11
        self.t[1, 'ID'] = 12
        self.assertEqual(self.t[0, 'ID'], 11)
        self.assertEqual(self.t[2, 'ID'], 13)        
    
    def testDump(self):
        
        with self.config.getTempStream() as dump:
            
            self.table.dump(dump)
            dump.flush()
            
            self.assertEqual(open(self.tsv).read(), open(dump.name).read())
            
    def testScalar(self):
        
        self.assertEqual(self.table.scalar(0, 'ID'), 11)
        self.assertEqual(self.table[0, 'ID'], 11)
        self.assertEqual(self.table[0, 0], 11)
        
        self.assertTrue(isinstance(self.table[:, :], csb.io.tsv.Table))
        self.assertTrue(isinstance(self.table[:, ('A',)], csb.io.tsv.Table))
        self.assertTrue(isinstance(self.table[:, 'A':], csb.io.tsv.Table))
        self.assertTrue(isinstance(self.table[1:2, 0], csb.io.tsv.Table))
        self.assertTrue(isinstance(self.table[(1, 2), 0], csb.io.tsv.Table))                        
        
    def testSelect(self):
        
        # test column selection
        self.assertEqual(self.table.select().columns, ('ID', 'A', 'B'))
        self.assertEqual(self.table.select('*').columns, ('ID', 'A', 'B'))
        self.assertEqual(self.table[:, :].columns, ('ID', 'A', 'B'))
        self.assertEqual(self.table[:, :'A'].columns, ('ID',))
        self.assertEqual(self.table[:, 'ID':'B'].columns, ('ID', 'A'))                
                
        self.assertEqual(self.table.select(['B', 'A']).columns, ('B', 'A'))
        self.assertEqual(self.table.select('A', 'B').columns, ('A', 'B'))
        self.assertEqual(self.table[:, ('B', 'A')].columns, ('B', 'A'))
                
        self.assertEqual(self.table.select(['A']).columns, ('A',))
        self.assertEqual(self.table.select(['A']).columns, ('A',))
        
        def fr(t):
            return list(t)[0]
        
        # test data        
        self.assertEqual(len(list(self.table[1:2, :])), 1)
        self.assertEqual(len(list(self.table[1:, :])), 2)
        self.assertEqual(len(list(self.table[(1, 2), :])), 2)
        self.assertEqual(len(list(self.table[(0, 1, 2), :])), 3)                           
        self.assertEqual(len(list(self.table[:, :])), 3)
                
        firstrow = fr(self.table.select('B', 'A'))
        self.assertEqual(firstrow[0], 'Row eleven')
        self.assertEqual(firstrow[1], 11.1)

        self.assertEqual(fr(self.table[:, :])[0], 11)
                
        self.assertEqual(fr(self.table[:, 'A'])[0], 11.1)
        self.assertEqual(fr(self.table[:, 'B':])[0], 'Row eleven')
                
        self.assertEqual(fr(self.table[2, :])[0], 13)
        self.assertEqual(fr(self.table[(1, 2), :])[0], 12)
        self.assertEqual(fr(self.table[1:9, :])[0], 12)  
        
    def testWhere(self):
        
        self.assertEqual(self.table.where('ID').equals(11).select('A').scalar(), 11.1)
        self.assertEqual(self.table.where('ID').greater(12).select('A').scalar(), None)        
        self.assertEqual(self.table.where('ID').between(11, 12).select('A').scalar(), 11.1)        
        self.assertEqual(self.table.where('ID').in_(11).select('A').scalar(), 11.1)
        self.assertEqual(self.table.where('ID').in_(11, 12).select('A').scalar(), 11.1)                

        
@test.unit
class TestAutoFlushStream(test.Case):
    
    def setUp(self):
        
        super(TestAutoFlushStream, self).setUp()
        self.temp = self.config.getTempStream()
        
    def testWrite(self):
        
        stream = csb.io.AutoFlushStream(self.temp)
        stream.write('$')
        
        content = open(self.temp.name).read()
        self.assertEqual(content, '$')          # will fail if not flushed (this is the actual test)
        stream.flush()
        self.assertEqual(content, '$')
        
    def tearDown(self):
        self.temp.close()


@test.unit
class TestMemoryStream(test.Case):
    
    def setUp(self):
        
        super(TestMemoryStream, self).setUp()
        self.stream = csb.io.MemoryStream()
        
    def testWrite(self):
        
        self.stream.write('tesT')
        self.assertEqual(self.stream.getvalue(), 'tesT')

    def testFlush(self):
        
        self.stream.write('tesT')
        self.stream.flush()
        self.assertEqual(self.stream.getvalue(), 'tesT')
        
    def testClose(self):
        self.stream.close()
        self.assertTrue(self.stream.closed)
                
    def testContextManager(self):
        with csb.io.MemoryStream() as stream:
            stream.write('tesT')
            stream.flush()
            self.assertEqual(stream.getvalue(), 'tesT')        
            

@test.unit
class TestTempFile(test.Case):                           
            
    def testCallForwarding(self):
        
        with csb.io.TempFile() as temp:
            
            self.assertTrue(hasattr(temp, 'write'))        
            self.assertTrue(hasattr(temp, 'flush'))        
            self.assertTrue(hasattr(temp, 'name'))
            self.assertTrue(hasattr(temp, 'close'))
            self.assertTrue(hasattr(temp, 'read'))                                    
            self.assertTrue(hasattr(temp, 'readline'))
            
            self.assertFalse(hasattr(temp, 'jxYUgg8'))
            self.assertRaises(AttributeError, lambda:getattr(temp, 'jxYUgg8'))

    def testContextManager(self):
                
        with csb.io.TempFile() as temp:
            self.assertTrue(os.path.exists(temp.name))
        self.assertFalse(os.path.exists(temp.name))

    def testClose(self):
        
        temp = csb.io.TempFile()

        self.assertTrue(os.path.isfile(temp.name))
        temp.close()
        self.assertFalse(os.path.exists(temp.name)) 
            
    def testAutoCleanup(self):
        
        for dispose in [True, False]:
            temp = csb.io.TempFolder(dispose=dispose)
            name = temp.name
        
            self.assertTrue(os.path.isdir(name))
            del temp
            self.assertEqual(os.path.isdir(name), not dispose)
        
    def testMultipleHandles(self):
        
        with csb.io.TempFile() as temp:
            
            data = '_Test_'
            temp.write(data)
            temp.flush()
            
            self.assertEqual(data, open(temp.name).read())        
        
        
@test.unit
class TestTempFolder(test.Case):                           

    def testContextManager(self):
        
        with csb.io.TempFolder() as temp:
            self.assertTrue(os.path.isdir(temp.name))
            
        self.assertFalse(os.path.exists(temp.name))
        
    def testClose(self):
        
        temp = csb.io.TempFolder()

        for i in range(3):
            name = os.path.join(temp.name, 'test{0}.txt'.format(i))
            with open(name, 'w') as f:
                f.write('.')
                
        self.assertTrue(os.path.isdir(temp.name))
        temp.close()
        self.assertFalse(os.path.exists(temp.name))   

    def testAutoCleanup(self):
        
        for dispose in [True, False]:
            temp = csb.io.TempFolder(dispose=dispose)
            name = temp.name
        
            self.assertTrue(os.path.isdir(name))
            del temp
            self.assertEqual(os.path.isdir(name), not dispose)
                    
                                        
@test.unit
class TestEntryReader(test.Case):
    
    def setUp(self):
        
        super(TestEntryReader, self).setUp()
        
        self.stream = csb.io.MemoryStream()
        self.data = r"""

>>E0
text
>>E1
text1
text2

text3
>>E2
//      """        
        self.stream.write(self.data)

    def testStartMarker(self):
        
        reader = csb.io.EntryReader(self.stream, '>>', None)        
        e = reader.readall()
        
        self.assertEqual(e[0].splitlines(), ['>>E0', 'text'])
        self.assertEqual(e[1].splitlines(), ['>>E1', 'text1', 'text2', '', 'text3'])
        self.assertEqual(e[2].splitlines(), ['>>E2', '//      '])
                
    def testAllMarkers(self):
        
        reader = csb.io.EntryReader(self.stream, '>>', '//')        
        e = reader.readall()
                
        self.assertEqual(e[0].splitlines(), ['>>E0', 'text'])
        self.assertEqual(e[1].splitlines(), ['>>E1', 'text1', 'text2', '', 'text3'])
        self.assertEqual(e[2].splitlines(), ['>>E2'])
        
    def testReadAll(self):
        
        r1 = csb.io.EntryReader(self.stream, '>>', '//')        
        self.assertEqual(len(r1.readall()), 3)

        r2 = csb.io.EntryReader(self.stream, '>>', None)        
        self.assertEqual(len(r2.readall()), 3)
    
    def testRead(self):
        
        reader = csb.io.EntryReader(self.stream, '>>', '//')
        self.assertEqual(reader.readall(), list(reader.entries()))
        
        self.assertEqual(type(reader.entries()), types.GeneratorType)
        
    def testDel(self):
        
        stream = csb.io.MemoryStream()
        reader = csb.io.EntryReader(stream, None, None)

        self.assertFalse(stream.closed)        
        del reader
        self.assertTrue(stream.closed)


@test.unit
class TestEntryWriter(test.Case):
            
    def testConstructor(self):

        # test with a file name
        tempFileName = self.config.getTempStream().name   
        writer = csb.io.EntryWriter(tempFileName, close=False)  # see below:
        self.assertTrue(hasattr(writer.stream, 'write')) 
        self.assertTrue(writer.autoclose)                       # close is ignored in init when destination is a path, not a stream
        
        # with a stream, autoclose
        temp = self.config.getTempStream()   
        writer = csb.io.EntryWriter(temp, close=True)
        # - test .stream
        self.assertEqual(writer.stream, temp)
        # - test .autoclose
        self.assertTrue(writer.autoclose)
        # - test del
        del writer
        self.assertTrue(temp.closed)  

        # with a stream, no autoclose        
        temp = self.config.getTempStream()   
        writer = csb.io.EntryWriter(temp, close=False)
        self.assertEqual(writer.stream, temp)
        self.assertFalse(writer.autoclose)             
        del writer
        self.assertFalse(temp.closed)
        temp.close()
             
    def testContextManager(self):
        
        with csb.io.EntryWriter(self.config.getTempStream(), close=True) as w:
            w.write('.')
            self.assertFalse(w.stream.closed)            
        self.assertTrue(w.stream.closed)    
          
        with csb.io.EntryWriter(self.config.getTempStream(), close=False) as w:
            w.write('.')
            self.assertFalse(w.stream.closed)            
        self.assertFalse(w.stream.closed)        
    
    def testWrite(self):
        
        with self.config.getTempStream() as temp:   
            
            writer = csb.io.EntryWriter(temp, close=False)
            writer.write('D')
            writer.stream.flush()
            
            self.assertEqual(open(temp.name).read(), 'D')
        
    def testWriteLine(self):
        
        with self.config.getTempStream() as temp:   
        
            writer = csb.io.EntryWriter(temp, close=False)
            writer.writeline('D')
            writer.stream.flush()
            
            self.assertEqual(open(temp.name).read(), 'D' + csb.io.NEWLINE)
        
    def testWriteAll(self):
        
        with self.config.getTempStream() as temp:  
        
            writer = csb.io.EntryWriter(temp, close=False)
            writer.writeall(('A', 'B', 'C',), '/')
            writer.stream.flush()
            
            self.assertEqual(open(temp.name).read(), 'A/B/C/')  



@test.functional
class TestDumpLoad(test.Case):

    def setUp(self):
        super(TestDumpLoad, self).setUp()
        
        self.lists = [[],
                      list(range(1000)),
                      list("Although that way may not be" + 
                           "obvious at first" + 
                           "unless you're Dutch.")] 
        self.arrays = [numpy.array([]),
                       numpy.random.random(1000),
                       numpy.arange(1000), ]
                       
        self.strings = ["",
                        "Although that way may not be" + \
                        "obvious at first" + \
                        "unless you're Dutch.",
                        "([0-9a-zA-Z]([-.\w]*[0-9a-zA-Z])"]

        # Completly connnected graph
        
        self.big_graph = []
        for _i in range(250):
            n = Node()
            self.big_graph.append(n)
        for n in self.big_graph:
            n.connections = set(self.big_graph)
            n.connections.remove(n)

        # Protein
        pdbfile = self.config.getTestFile('ake-xray-ensemble-ca.pdb')
        self.protein = LegacyStructureParser(pdbfile).parse_models()[0]

        self.objs = [self.lists, self.arrays, self.strings,
                     self.protein]

    def testSimple(self):
        load_func = csb.io.load
        dump_func = csb.io.dump
        self._test(load_func, dump_func)

    def testGzip(self):
        load_func = lambda x: csb.io.load(x, gzip=True)
        dump_func = lambda x, y: csb.io.dump(x, y, gzip=True)
        self._test(load_func, dump_func)

    def testGzipFail(self):
        load_func = lambda x: csb.io.load(x, gzip=True)
        dump_func = lambda x, y: csb.io.dump(x, y, gzip=False)
        with self.assertRaises(IOError):
            self._test(load_func, dump_func)
        
    def testLock(self):
        load_func = lambda x: csb.io.load(x, lock=True)
        dump_func = lambda x, y: csb.io.dump(x, y, lock=True)
        self._test(load_func, dump_func)

    def testLockGzip(self):
        load_func = lambda x: csb.io.load(x, gzip=True, lock=True)
        dump_func = lambda x, y: csb.io.dump(x, y, gzip=True, lock=True)
        self._test(load_func, dump_func)

    def testRecusionLimit(self):
        with self.assertRaises(RuntimeError):
            with csb.io.TempFile(mode='b') as temp:
                csb.io.dump(self.big_graph, temp.name)

    def _test(self, load_func, dump_func):
        for ob in self.objs:
            with csb.io.TempFile(mode='b') as temp:
                dump_func(ob, temp.name)
                ob2 = load_func(temp.name)
                if isinstance(ob, list):
                    for a, b in zip(ob2, ob):
                        if isinstance(a, numpy.ndarray):
                            self.assertTrue(numpy.all(a == b))
                        else:
                            self.assertEqual(a, b)
                elif isinstance(ob, csb.bio.structure.Structure):
                    self.assertEqual(ob.to_fasta(),
                                     ob2.to_fasta())
                    self.assertTrue(numpy.all(ob.get_coordinates()
                                              == ob2.get_coordinates()))
                else:
                    self.assertEqual(ob, ob2)

            
if __name__ == '__main__':
    
    test.Console()        
