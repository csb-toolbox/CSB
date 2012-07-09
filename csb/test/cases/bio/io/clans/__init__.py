import csb.test as test

from csb.bio.io import ClansParser
from csb.bio.io.clans import Clans, ClansEntry, ClansParams, ClansSeqgroup,\
     Color


@test.unit
class TestClansColor(test.Case):

    def setUp(self):
        
        super(TestClansColor, self).setUp()

    def testColorInit(self):
        color = Color()

        self.assertEqual(color.r, 0)
        self.assertEqual(color.g, 0)
        self.assertEqual(color.b, 0)

    def testColorSetter(self):
        color = Color()
        for color_name in ('r', 'g', 'b'):
            self.assertRaises(ValueError, color.__setattr__, color_name, -1)
            self.assertRaises(ValueError, color.__setattr__, color_name, 256)

    def testParseClansColorWithCorrectInput(self):

        correct_input = '83;92;3'
        color = Color.from_string(correct_input)
        self.assertEqual(color.r, 83)
        self.assertEqual(color.g, 92)
        self.assertEqual(color.b, 3)

    def testParseClansColorWithWrongInput(self):
        color = Color()

        wrong_input_1 = (83, 92, 3)
        self.assertRaises(TypeError,
                          color.from_string, wrong_input_1)

        wrong_input_2 = '83;92;3;'
        self.assertRaises(ValueError, color.from_string, wrong_input_2)

        wrong_input_3 = '83;92'
        self.assertRaises(ValueError, color.from_string, wrong_input_3)

    def testToClansColor(self):
        color = Color()

        self.assertEqual(color.to_clans_color(), '0;0;0;255')

        testValues = (83, 92, 3, 87)
        color.r = testValues[0]
        color.g = testValues[1]
        color.b = testValues[2]
        color.a = testValues[3]

        self.assertEqual(color.to_clans_color(),
                         ';'.join(map(str, testValues)))


@test.functional
class TestClansParams(test.Case):

    def setUp(self):
        
        super(TestClansParams, self).setUp()

    def testInstatiation(self):
        cp = ClansParams()
        for attribute_name, default_value in cp._DEFAULTS.items():
            if attribute_name == 'colors':
                continue
            self.assertEqual(cp.__getattribute__(attribute_name),
                             default_value)

    def testUnknownParamFail(self):
        self.assertRaises(KeyError, ClansParams, **{'unknownParam': True})

    def testForbiddenAssignments(self):
        self.assertRaises(ValueError, ClansParams, **{'attfactor': 'a'})
        self.assertRaises(ValueError, ClansParams, **{'attvalpow': 'a'})
        self.assertRaises(ValueError, ClansParams, **{'avgfoldchange': 'a'})
        self.assertRaises(ValueError, ClansParams, **{'blastpath': 3})
        self.assertRaises(ValueError, ClansParams, **{'cluster2d': 'a'})
        self.assertRaises(ValueError, ClansParams, **{'colors': 'a'})
        self.assertRaises(ValueError, ClansParams, **{'complexatt': 'a'})
        self.assertRaises(ValueError, ClansParams, **{'cooling': 'a'})
        self.assertRaises(ValueError, ClansParams, **{'currcool': 'a'})
        self.assertRaises(ValueError, ClansParams, **{'dampening': 'a'})
        self.assertRaises(ValueError, ClansParams, **{'dotsize': 'a'})
        self.assertRaises(ValueError, ClansParams, **{'formatdbpath': 3})
        self.assertRaises(ValueError, ClansParams, **{'groupsize': 'a'})
        self.assertRaises(ValueError, ClansParams, **{'maxmove': 'a'})
        self.assertRaises(ValueError, ClansParams, **{'minattract': 'a'})
        self.assertRaises(ValueError, ClansParams, **{'ovalsize': 'a'})
        self.assertRaises(ValueError, ClansParams, **{'pval': 'a'})
        self.assertRaises(ValueError, ClansParams, **{'repfactor': 'a'})
        self.assertRaises(ValueError, ClansParams, **{'repvalpow': 'a'})
        self.assertRaises(ValueError, ClansParams, **{'showinfo': 'a'})
        self.assertRaises(ValueError, ClansParams, **{'usefoldchange': 'a'})
        self.assertRaises(ValueError, ClansParams, **{'usescval': 'a'})
        self.assertRaises(ValueError, ClansParams, **{'zoom': 'a'})


@test.functional
class TestClans(test.Case):

    def setUp(self):
        
        super(TestClans, self).setUp()
        
    def testClansInit(self):
        '''
        Test creating an empty L{Clans} instance.
        '''

        c = Clans()

        param_names = ['attfactor', 'attvalpow', 'avgfoldchange', 'blastpath',
                       'cluster2d', 'colors', 'complexatt', 'cooling',
                       'currcool', 'dampening', 'dotsize', 'formatdbpath',
                       'groupsize', 'maxmove', 'minattract', 'ovalsize',
                       'pval', 'repfactor', 'repvalpow', 'showinfo',
                       'usefoldchange', 'usescval', 'zoom']

        for param_name in param_names:
            self.assertTrue(hasattr(c.params, param_name))
                
        self.assertEqual(c.filename, None)
        self.assertEqual(c.rotmtx.shape, (3, 3))
        self.assertEqual(len(c.entries), 0)
        self.assertEqual(len(c.seqgroups), 0)

    def testClansEntryAddingAndSorting(self):
        c = Clans()

        names = ['a', 'b', 'c', 'd', 'e', 'f', 'g']
        shuffled_names = ['g', 'f', 'b', 'd', 'e', 'c', 'a']

        for name in shuffled_names:
            c.add_entry(ClansEntry(name=name))

        c.sort()

        for i, e in enumerate(c):
            self.assertEqual(e.name, names[i])

    def testGetEntry(self):
        c = Clans()

        ## get non-existant entry from empty clans instance
        self.assertRaises(ValueError, c.get_entry, 'a')

        names = ['a', 'b', 'c', 'd', 'e', 'f', 'g']
        entries = [ClansEntry(name=name) for name in names]
        [c.add_entry(e) for e in entries]

        ## check whether entries fetched by name match those created
        for i, name in enumerate(names):
            self.assertEqual(c.get_entry(name), entries[i])

        ## check pedantic flag for duplicate name='a' entries
        c.add_entry(ClansEntry(name='a'))

        self.assertTrue(c.get_entry('a', False).name == 'a')

        self.assertRaises(ValueError, c.get_entry, 'a', True)


@test.functional
class TestClansSeqgroup(test.Case):

    def setUp(self):
        
        super(TestClansSeqgroup, self).setUp()

    def testInit(self):
        sg = ClansSeqgroup()

        self.assertTrue(sg.is_empty())

    def testAddingAndRemovingSeqgroups(self):
        c = Clans()

        names = ['a', 'b', 'c', 'd', 'e', 'f', 'g']
        for i, name in enumerate(names):
            c.add_group(ClansSeqgroup(name=name))
            self.assertEqual(len(c.seqgroups), i + 1)

        removed = 0
        while len(c.seqgroups) != 0:
            c.remove_group(c.seqgroups[-1])
            removed += 1

        self.assertEqual(removed, len(names))
        self.assertEqual(len(c.seqgroups), 0)

        testGroup = ClansSeqgroup()
        self.assertRaises(TypeError, testGroup.add, 23)
        self.assertRaises(TypeError, testGroup.remove, 23)

    def testAddingClansEntries(self):
        c = Clans()

        sg = ClansSeqgroup()
        c.add_group(sg)

        e = ClansEntry()
        c.add_entry(e)

        ## add entry to seqgroup
        sg.add(e)
        self.assertEqual(len(sg), 1)
        self.assertEqual(len(e.groups), 1)

        ## adding the same entry is forbidden
        self.assertRaises(ValueError, sg.add, e)

        ## adding s.th. else than a ClansEntry
        self.assertRaises(TypeError, sg.add, 23)

    def testRemovingClansEntries(self):
        c = Clans()

        sg = ClansSeqgroup()
        c.add_group(sg)

        e = ClansEntry()
        c.add_entry(e)

        sg.add(e)
        sg.remove(e)
        self.assertEqual(len(sg), 0)
        self.assertEqual(len(e.groups), 0)
        self.assertRaises(TypeError, sg.remove, 23)
        self.assertRaises(ValueError, sg.remove, e)
        

@test.functional
class TestClansParser(test.Case):

    def setUp(self):
        
        super(TestClansParser, self).setUp()
                
        self.filename = self.config.getTestFile('out.clans')

    def testPrematureGetter(self):
        '''
        Test whether the premature (before parsing) access to clans_instance is
        properly handled.
        '''
        cp = ClansParser()
        self.assertRaises(ValueError, cp.__getattribute__, 'clans_instance')

    def testParseFile(self):
        '''
        Test parsing of a small dummy file with known values
        '''
        
        from numpy import array

        cp = ClansParser()
        
        self.clans_instance = cp.parse_file(self.filename)
        
        self.assertEqual(len(self.clans_instance), 41)
        self.assertRaises(IndexError, self.clans_instance.__getitem__, 41)

        correct_rotmtx = array([[0.75614862, 0.65439992, 0.],
                                [-0.65439992, 0.75614862, 0.],
                                [0., 0., 1.]])
        self.assertEqual(self.clans_instance.rotmtx.shape, (3, 3))
        self.assertTrue(
            (self.clans_instance.rotmtx - correct_rotmtx < 1e-6).all())

        self.assertEqual(len(self.clans_instance.seqgroups), 2)

        seqgroup_names = ('insect hypoth. protein (2 copies, C term)',
                          'allergens')
        seqgroup_sizes = (20, 17)

        for i, seqgroup in enumerate(self.clans_instance.seqgroups):
            self.assertEqual(len(seqgroup), seqgroup_sizes[i])
            self.assertEqual(seqgroup.name, seqgroup_names[i])


@test.functional
class TestClansFileWriter(test.Case):

    def setUp(self):
        
        super(TestClansFileWriter, self).setUp()
                
        self.filename = self.config.getTestFile('out.clans')
        self.temp = self.config.getTempStream()

    def testWrittenIsIdenticalToOriginal(self):
        cp = ClansParser()
        clans_instance = cp.parse_file(self.filename)

        clans_instance.write(self.temp.name)
        self.temp.flush()

        with open(self.filename) as original_file:
            original_lines = original_file.readlines()

        with open(self.temp.name) as written_file:
            written_lines = written_file.readlines()

        self.assertEqual(len(original_lines), len(written_lines))

        in_hsps = False
        start_tag_hsp = '<hsp>\n'
        end_tag_hsp = '</hsp>\n'
        colorarr_tag = 'colorarr='
        color_tag = 'color='
        
        for i, original_line in enumerate(original_lines):

            if original_line == start_tag_hsp:
                in_hsps = True
                continue
            if original_line == end_tag_hsp:
                in_hsps = False
                continue

            if original_line.startswith(colorarr_tag):
                ## remove colorarr_tag from beginning of line
                original_line = original_line[len(colorarr_tag):].strip().strip(':')
                self.assertTrue(written_lines[i].startswith(colorarr_tag))
                written_line = written_lines[i][len(colorarr_tag):].strip().strip(':')

                original_colors = original_line.replace('(', ''). replace(')', '').split(':')
                written_colors = written_line.replace('(', ''). replace(')', '').split(':')

                self.assertEqual(len(original_colors), len(written_colors))
                
                for j, original_color_string in enumerate(original_colors):
                    original_color = Color.from_string(original_color_string)
                    written_color = Color.from_string(written_colors[j])
                    self.assertEqual(original_color.r, written_color.r)
                    self.assertEqual(original_color.g, written_color.g)
                    self.assertEqual(original_color.b, written_color.b)
                    self.assertEqual(original_color.a, written_color.a)

                continue

            if original_line.startswith(color_tag):
                original_color_string = original_line[len(color_tag):].strip()
                self.assertTrue(written_lines[i].startswith(color_tag))
                written_color_string = written_lines[i][len(color_tag):].strip()

                original_color = Color.from_string(original_color_string)
                written_color = Color.from_string(written_color_string)
                self.assertEqual(original_color.r, written_color.r)
                self.assertEqual(original_color.g, written_color.g)
                self.assertEqual(original_color.b, written_color.b)
                self.assertEqual(original_color.a, written_color.a)

                continue
                
            if in_hsps:
                original_start_end, original_value \
                                    = original_line.strip().split(':')
                written_start_end, written_value \
                                   = written_lines[i].strip().split(':')
                self.assertEqual(original_start_end, written_start_end)
                self.assertTrue((float(original_value) - float(written_value)) < 1e-6)

            else:
                self.assertEqual(original_line, written_lines[i])
            

    def tearDown(self):
        self.temp.close()
    

if __name__ == '__main__':
    
    test.Console()
