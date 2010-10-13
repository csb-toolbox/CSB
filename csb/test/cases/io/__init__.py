"""
@todo: write tests for C{dump} and C{load}
"""

import StringIO
import types
import csb.io
import csb.test as test


@test.unit
class TestUnbufferedStream(test.Case):
    
    def setUp(self):
        
        super(TestUnbufferedStream, self).setUp()
        self.temp = self.config.getTempStream()
        
    def testWrite(self):
        
        stream = csb.io.UnbufferedStreamWrapper(self.temp)
        stream.write('$')
        
        content = open(self.temp.name).read()
        self.assertEqual(content, '$')          # will fail if not flushed (this is the actual test)
        stream.flush()
        self.assertEqual(content, '$')
                           
                            
@test.unit
class TestEntryReader(test.Case):
    
    def setUp(self):
        
        super(TestEntryReader, self).setUp()
        
        self.stream = StringIO.StringIO()
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
        
        stream = StringIO.StringIO()
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
        self.assertTrue(isinstance(writer.stream, file)) 
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
        
        temp = self.config.getTempStream()   
        
        writer = csb.io.EntryWriter(temp, close=False)
        writer.write('D')
        writer.stream.flush()
        
        self.assertEqual(open(temp.name).read(), 'D')
        
    def testWriteLine(self):
        
        temp = self.config.getTempStream()   
        
        writer = csb.io.EntryWriter(temp, close=False)
        writer.writeline('D')
        writer.stream.flush()
        
        self.assertEqual(open(temp.name).read(), 'D' + csb.io.NEWLINE)
        
    def testWriteAll(self):
        
        temp = self.config.getTempStream()   
        
        writer = csb.io.EntryWriter(temp, close=False)
        writer.writeall(('A', 'B', 'C',), '/')
        writer.stream.flush()
        
        self.assertEqual(open(temp.name).read(), 'A/B/C/')  

        
        
if __name__ == '__main__':
    
    test.Console()        
