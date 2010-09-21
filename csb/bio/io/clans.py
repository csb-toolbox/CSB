'''
Classes for parsing/manipulating/writing CLANS (by Tancred Frickey) files

Author: Klaus Kopec
MPI fuer Entwicklungsbiologie, Tuebingen
Copyright (C) 2010 Klaus Kopec
All rights reserved.
No warranty implied or expressed.
'''
import os, re, sys
from numpy import array, float64, eye, random


class ClansParser:
    '''
    Class for parsing CLANS files.
    '''

    def __init__(self):
        self.clans_instance = None
        self.data_block_dict = {}

    def __repr__(self):
        return 'ClansParser instance'

    def parse_file(self, filename):
        '''Parse a CLANS file.

        @param filename: name of the CLANS file.
        @type filename: string

        @rtype: L{Clans} instance
        @return: a Clans instance containing the parsed data
        '''
        self.clans_instance = Clans()
        self.clans_instance.filename = filename

        self._read_block_dict()  # read and preprocess the CLANS file

        self._parse_param()
        self._parse_rotmtx()
        seq = self._parse_seq()
        seqgroups = self._parse_seqgroups()
        pos = self._parse_pos()

        if 'hsp' in self.data_block_dict:
            hsp = self._parse_hsp()
        elif 'mtx' in self.data_block_dict:
            hsp = self._parse_mtx()
        elif 'att' in self.data_block_dict:
            hsp = self._parse_att()

        ## print unknown blocks in case further implementations are needed
        known_block_tags = set(('param', 'rotmtx', 'seq', 'seqgroups', 'pos',
                                'hsp', 'mtx', 'att'))
        unprocessed_block_tags = set(self.data_block_dict.keys()).difference(
            known_block_tags)

        for unprocessed_tag in unprocessed_block_tags:
            print '[ClansParser.parse_file:WARNING] tag "%s" unknown.' \
                  % unprocessed_tag \
                  + '    File corrupt or further implementations needed!'

        ## if no entries exist, we cannot add pos, seqgroup and hsp data
        if len(seq) > 0:

            ## add Entries
            if len(pos) > 0:
                self.clans_instance.entries = [
                    ClansEntry(seq[i][0], seq[i][1],
                               pos[i], parent=self.clans_instance)
                    for i in pos]

            ## add groups
            self.clans_instance.seqgroups = []
            if len(seqgroups) > 0:
                for group_raw_data in seqgroups:

                    group = ClansSeqgroup(name=group_raw_data['name'],
                                          type=group_raw_data['type'],
                                          size=group_raw_data['size'],
                                          hide=group_raw_data['hide'],
                                          color=group_raw_data['color'])

                    ## get members corresponding to the IDs in this group
                    members = [self.clans_instance.entries[number]
                               for number in group_raw_data['numbers']]

                    self.clans_instance.add_group(group, members)

            ## add hsp values
            if len(hsp) > 0:
                [self.clans_instance.entries[a].add_hsp(
                    self.clans_instance.entries[b], value)
                 for ((a, b), value) in hsp.items()]

        self.clans_instance._update_index()

        return self.clans_instance

    def _read_block_dict(self):
        '''Extracts all <tag>DATA</tag> blocks from file
        self.clans_instance.filename.

        @rtype: dict
        @return: data in the form: dict[tag] = DATA.
        '''
        # read file and remove the first line, i.e. sequence=SEQUENCE_COUNT
        data_blocks = file(os.path.expanduser(
            self.clans_instance.filename)).read().split('\n', 1)[1]

        ## flag re.DOTALL is necessary to make . match newlines
        data = re.findall(r'(<(\w+)>(.+)</\2>)', data_blocks,
                          flags=re.DOTALL)
        self.data_block_dict = dict([(tag, datum.strip().split('\n'))
                                     for _tag_plus_data, tag, datum in data])

    def _parse_param(self):
        '''
        Parse a list of lines in the CLANS <param> format:

        parameter1=data1\n
        parameter2=data2\n
        ...
        '''
        if 'param' not in self.data_block_dict:
            print 'WARNING: CLANS file contains no <param> block.'
            return

        block = self.data_block_dict['param']

        tmp_params = dict([block[i].split('=') for i in range(len(block))])

        ## create colors entry from colorcutoffs and colorarr
        colorcutoffs = [float(val) for val in
                        tmp_params.pop('colorcutoffs').strip(';').split(';')]
        colors = tmp_params.pop('colorarr').strip(':')
        colors = colors.replace('(', '').replace(')', '').split(':')
        colorarr = [RGB_color(*color_definition)
                    for color_definition in [color.split(';')
                                             for color in colors]]

        tmp_params['colors'] = dict(zip(colorcutoffs, colorarr))
        self.clans_instance.params = tmp_params

    def _parse_rotmtx(self):
        '''
        Parse a list of lines in the CLANS <rotmtx> format. The data is stored
        in the clans_instance as a 3x3 numpy.array.

        @raise ValueError: if the rotmtx block does not contain exactly 3 lines
        '''
        if 'rotmtx' not in self.data_block_dict:
            print 'WARNING: CLANS file contains no <rotmtx> block.'
            return

        block = self.data_block_dict['rotmtx']

        if len(block) != 3:
            raise ValueError('CLANS <rotmtx> blocks comprise exactly 3 lines.')
        self.clans_instance.rotmtx = array(
            [[float64(val) for val in line.split(';')[:3]] for line in block])

    def _parse_seq(self):
        '''
        Parse a list of lines in the CLANS <seq> format, which are in FASTA
        format.

        @rtype: dict
        @return: dict with running numbers as key and 2-tuples (id, sequence)
                 as values
        '''
        if 'seq' not in self.data_block_dict:
            print 'WARNING: CLANS file contains no <seq> block. This is OK,'
            print '         if the file does not contain any sequences.'
            return {}

        block = self.data_block_dict['seq']

        return dict([(i, (block[2 * i][1:], block[2 * i + 1].strip()))
                     for i in range(len(block) / 2)])

    def _parse_seqgroups(self):
        '''
        Parse a list of lines in the CLANS <seqgroup> format:

        name=name of the group\n
        type=0\n
        size=12\n
        hide=0\n
        color=255;204;51\n
        numbers=0;1;2;3;4;5;6;10;13\n
        ...

        @rtype: list
        @return: list of dicts (one for each group) with the tags (name, type,
                 size, hide, ...) as keys and their typecasted data as values
                 (i.e. name will be a string, size will be an integer, etc)
        '''
        if 'seqgroups' not in self.data_block_dict:
            return []

        block = self.data_block_dict['seqgroups']

        groups = []
        for line in block:
            p, v = line.split('=')
            if p == 'name':
                groups.append({'name': v})
            elif p == 'numbers':
                groups[-1][p] = [int(val) for val in v.split(';')[:-1]]
            else:
                groups[-1][p] = v
        return groups

    def _parse_pos(self):
        '''
        Parse a list of lines in the CLANS <pos> format \'INT FLOAT FLOAT
        FLOAT\'.

        @rtype: dict
        @return: a dict using the integers as keys and a (3,1)-array created
                 from the three floats as values.
        '''
        if 'pos' not in self.data_block_dict:
            print 'WARNING: CLANS file contains no <pos> block. This is OK,'
            print '         if the file does not contain any sequences.'
            return {}

        block = self.data_block_dict['pos']

        return dict([(int(l.split()[0]),
                      array([float64(val) for val in l.split()[1:]]))
                     for l in block])

    def _parse_hsp(self):
        '''
        Parse a list of lines in the CLANS <hsp> format \'INT INT: FLOAT\'.

        NOTE: some CLANS <hsp> lines contain more than one float; we omit the
        additional numbers

        @rtype: dict
        @return: a dict using 2-tuples of the two integers as keys and the
                 float as values
        '''
        if 'hsp' not in self.data_block_dict:
            print 'WARNING: CLANS file contains no <hsp> block. This is OK,'
            print '         if the file does not contain any sequences or if'
            print '         none of the contained sequences have connections.'
            return {}

        block = self.data_block_dict['hsp']

        return dict([(tuple([int(val) for val in line.split(':')[0].split()]),
                      float(line.split(':')[1].split(' ')[0]))
                     for line in block])

    def _parse_mtx(self):
        '''
        Parse a list of lines in the CLANS <mtx> format.

        @rtype: dict
        @return: a dict using 2-tuples of the two integers as keys and the
                 float as values
        '''
        if 'mtx' not in self.data_block_dict:
            print 'WARNING: CLANS file contains no <mtx> block. This is OK,'
            print '         if the file does not contain any sequences or if'
            print '         none of the contained sequences have connections.'
            return {}

        block = self.data_block_dict['mtx']

        return dict([((i, j), float(entry))
                     for i, line in enumerate(block)
                     for j, entry in enumerate(line.split(';')[:-1])
                     if i != j if float(entry) != 0])


class ClansWriter:
    '''
    Class for writing Clans instances to files in the CLANS format.
    '''

    def __init__(self, clans_instance, output_filename):
        '''Writes the content of a Clans instance to to CLANS readable file.'''
        self.clans_instance = clans_instance

        self.file = file(os.path.expanduser(output_filename), 'w')

        self.file.write('sequences=%i\n' % len(clans_instance.entries))

        ## these methods append the CLANS file blocks to self.output_string
        self.clans_param_block()
        self.clans_rotmtx_block()
        self.clans_seq_block()
        self.clans_seqgroups_block()
        self.clans_pos_block()
        self.clans_hsp_block()

        self.file.close()

    ## methods for creating a CLANS files
    def clans_param_block(self):
        '''Adds a <param>data</param> CLANS file block to self.output_string'''
        params_dict = dict(self.clans_instance.params)

        if params_dict is None:
            return

        ## translate 'colors' into 'colorcutoffs' and 'colorarr' for clans
        if 'colors' in params_dict:
            cutoffs = sorted(params_dict['colors'].keys())
            params_dict['colorcutoffs'] = ''.join(['%.2f;' % cutoff
                                                   for cutoff in cutoffs])
            params_dict['colorarr'] = ''.join(
                ['(%s):' % params_dict['colors'][cutoff].rgb
                 for cutoff in cutoffs])
            params_dict.pop('colors')

        self.file.write('<param>\n')
        self.file.write('\n'.join(
            ['%s=%s' % (param_name, params_dict[param_name])
             for param_name in sorted(params_dict.keys())]))
        self.file.write('\n</param>\n')

    def clans_rotmtx_block(self):
        '''Adds a <rotmtx>data</rotmtx> CLANS file block to
        self.output_string

        @raise ValueError: if self.clans_instance.rotmtx is no 3x3 numpy.array
        '''
        rotmtx = self.clans_instance.rotmtx

        if rotmtx is None:
            return

        if rotmtx.shape != (3, 3):
            raise ValueError('rotmtx must be a 3x3 array')

        self.file.write('<rotmtx>\n')
        self.file.write('\n'.join(
            ['%s;%s;%s;' % tuple(rotmtx[i]) for i in range(3)]))
        self.file.write('\n</rotmtx>\n')

    def clans_seq_block(self):
        '''Adds a <seq>data</seq> CLANS file block to self.output_string'''
        self.file.write('<seq>\n')
        self.file.write(''.join([e.output_string_seq()
                                 for e in self.clans_instance.entries]))
        self.file.write('</seq>\n')

    def clans_seqgroups_block(self):
        '''Adds a <seqgroupsparam>data</seqgroups> CLANS file block
        to self.output_string'''
        seqgroups = self.clans_instance.seqgroups

        if seqgroups is None or len(seqgroups) == 0:
            return

        self.file.write('<seqgroups>\n')
        self.file.write('\n'.join([s.output_string() for s in seqgroups]))
        self.file.write('\n</seqgroups>\n')

    def clans_pos_block(self):
        '''Adds a <pos>data</pos> CLANS file block to self.output_string'''
        self.file.write('<pos>\n')
        self.file.write('\n'.join([e.output_string_pos()
                                   for e in self.clans_instance.entries]))
        self.file.write('\n</pos>\n')

    def clans_hsp_block(self):
        '''Adds a <hsp>data</hsp> CLANS file block to self.output_string'''
        self.file.write('<hsp>\n')

        ## sorting is not necessary, but makes a nicer looking clans file
        for entry1 in self.clans_instance.entries:
            entry1_id = entry1.get_id()

            for (entry2, pvalue) in sorted(entry1.hsp.items()):
                entry2_id = entry2.get_id()

                if entry1_id >= entry2_id:
                    continue

                self.file.write('%i %i:%s\n' %
                                (entry1_id, entry2_id, repr(pvalue)))

        self.file.write('</hsp>\n')


def parse_gi_and_residues(name):
    '''

    @raise ValueError: if a residue range is found that contains no end residue
    '''
    start = name.find('gi|')
    if start == -1:
        return name
    real_start = start + 3
    name = name[real_start:]

    gi_number = name.split('|', 1)[0]

    next_gi_start = name[real_start:].find('gi|')

    if next_gi_start != -1:
        name = name[:next_gi_start]

    initial_residue_number = name.find('(')
    if initial_residue_number == -1:
        return gi_number

    start = name[initial_residue_number + 1:].split('-')
    ## if start isn't an integer, we assume the '(' is not the start of a range
    try:
        start = int(start[0])
    except ValueError:
        return gi_number

    residues_end = name.find(':')
    if residues_end == -1:
        ## some entries are not (x-y:z), but only (x-y)
        residues_end = name.find(')')
        if residues_end == -1:
            raise ValueError('no end residue found in name\n\t%s' % name)

    potential_start_and_end = name[:residues_end].split('-')

    if len(potential_start_and_end) != 2:
        return gi_number
    try:
        first_res, last_res = [int(val) for val in potential_start_and_end]
    except ValueError:
        return gi_number

    return (gi_number, int(first_res), int(last_res))


class gi_comparator(object):

    def __init__(self):
        self.mapping = {}  # mapping cache for faster access

    def __call__(self, entry1, entry2):
        if entry1.name in self.mapping:
            entry1_parsed = self.mapping[entry1.name]
        else:
            entry1_parsed = parse_gi_and_residues(entry1.name)
            self.mapping[entry1.name] = entry1_parsed

        if entry2.name in self.mapping:
            entry2_parsed = self.mapping[entry2.name]
        else:
            entry2_parsed = parse_gi_and_residues(entry2.name)
            self.mapping[entry2.name] = entry2_parsed

        if entry1_parsed == entry2_parsed:
            return 0

        if len(entry1_parsed) == 3 and len(entry2_parsed) == 3:
            A = dict(zip(('gi', 'start', 'end'), entry1_parsed))
            B = dict(zip(('gi', 'start', 'end'), entry2_parsed))

            if A['gi'] != B['gi']:  # different gi numbers
                return -1

            ## switch so that A is the one that starts earlier

            if A['start'] > B['start']:
                A, B = B, A

            common_residues = A['end'] - B['start']
            if common_residues < 0:
                return -1  # B starts after A ends

            if B['end'] < A['end']:
                return 0  # as A starts before B and ends after it, B is in A

            ## > 75% of length of the shorter one are shared => identical
            if common_residues > 0.75 * min(A['end'] - A['start'],
                                            B['end'] - B['start']):
                return 0
        return -1


class RGB_color(dict):
    '''
    Class for holding and manipulating one RGB color value.
    '''

    def __init__(self, r=0, g=0, b=0):
        dict.__init__(self)
        self.r = int(r)
        self.g = int(g)
        self.b = int(b)

    def __repr__(self):
        return '%s: r%i;g%i;b%i' % (self.__class__.__name__,
                                    self.r, self.g, self.b)

    def __setitem__(self, item, value):
        if item not in ['r', 'g', 'b']:
            raise ValueError('only \'r\' \'g\', and \'b\' can be stored')

        dict.__setitem__(self, item, value)

    def __setattr__(self, attr, value):
        if attr in self.keys() or attr == 'rgb':
            self.set_color(attr, value)
        else:
            self[attr] = value

    def __getattribute__(self, attr):
        if attr in dict.keys(self):
            return self[attr]
        if attr == 'rgb':
            return '%i;%i;%i' % (self.r, self.g, self.b)

        return object.__getattribute__(self, attr)

    def set_color(self, colorname, color):
        '''Sets one or all colors to a new value.

        @param colorname: the color that will be changed. Can be \'rgb\' or any
        single letter of the three.
        @type colorname: string

        @param color: the color that is assigned to <colorname>
        @type color: a dict with keys [\'r\', \'g\', \'b\'], a semicolon
        separated string of 3 values (e.g. \'33;66;123\'), or a single value

        @raises ValueError: if any value in color is outside of range(256)
        '''
        if colorname == 'rgb':
            if isinstance(color, basestring) and color.count(';') == 2:
                color = dict(zip(['r', 'g', 'b'],
                                     [int(val) for val in color.split(';')]))

            elif not isinstance(color, dict):
                raise ValueError('unsupported color')

        elif colorname in self.keys():
            color = {colorname: int(color)}

        else:
            raise ValueError('only \'rgb\', \'r\', \'g\', and \'b\' are ' +
                             'valid colornames')

        for v in color.values():
            if v > 255 or v < 0:
                raise ValueError('valid color values are in range 0-255.')

        self.update(color)


class Clans(object):
    '''
    Class for holding and manipulating data from one CLANS file.
    '''

    def __init__(self):
        self.filename = None

        self.params = {}
        self.set_default_params()

        self.rotmtx = None
        self.set_default_rotmtx()

        self.entries = []
        self.seqgroups = []
        self._has_good_index = False

    def __repr__(self):
        return 'Clans object: %i sequences; %i seqgroups' % \
               (len(self), len(self.seqgroups))

    def __len__(self):
        return len(self.entries)

    def __getitem__(self, index):
        return self.entries[index]

    def __setitem__(self, index, data):
        self.entries[index] = data
        self._has_good_index = False

    def _update_index(self):
        '''Creates a mapping of entry names to entry indices in the Clans
        instance, speeding up entry.get_id() calls. The Index was introduced
        to get a better Clans.write() performance, which suffered from
        excessive entry.get_id() calls during HSP block generation (see clans
        _hsp_block()).

        NOTE: the index needs unique entry names, therefore remove_duplicates
        is called first and can decrease the number of entries!!!
        '''
        self.remove_duplicates()

        self._idx = dict([(e._get_unique_id(), i)
                          for i, e in enumerate(self.entries)])
        self._has_good_index = True

    def set_default_params(self):
        '''Sets the parameters to CLANS\' default values.'''
        self.params = {
            'maxmove': .1,
            'pval': 1.,
            'usescval': 'false',
            'complexatt': 'true',
            'cooling': 1.0,
            'currcool': 1.0,
            'attfactor': 10.0,
            'attvalpow': 1,
            'repfactor': 5.0,
            'repvalpow': 1,
            'dampening': 0.2,
            'minattract': 1.0,
            'cluster2d': 'false',
            'blastpath': 'blastall -p blastp',
            'formatdbpath': 'formatdb',
            'showinfo': 'true',
            'zoom': 1.0,
            'dotsize': 2,
            'ovalsize': 10,
            'groupsize': 4,
            'usefoldchange': 'false',
            'avgfoldchange': 'false',
            'colors': {.0: RGB_color(230, 230, 230),
                       .1: RGB_color(207, 207, 207),
                       .2: RGB_color(184, 184, 184),
                       .3: RGB_color(161, 161, 161),
                       .4: RGB_color(138, 138, 138),
                       .5: RGB_color(115, 115, 115),
                       .6: RGB_color(92, 92, 92),
                       .7: RGB_color(69, 69, 69),
                       .8: RGB_color(46, 46, 46),
                       .9: RGB_color(23, 23, 23)}
            }

    def set_default_rotmtx(self):
        '''Resets the rotation matrix (rotmtx) to no rotation.'''
        self.rotmtx = eye(3)

    def add_group(self, group, members=None):
        '''Adds a new group.

        @param group: the new group
        @type group: a ClansSeqgroup instance

        @param members: an iterable of members that are added to the new group
        @type members: list of ClansEntry instances

        @raise ValueError: if <group> is no ClansSeqgroup instance
        '''
        if not isinstance(group, ClansSeqgroup):
            raise ValueError('groups need to be ClansSeqgroup instances')

        self.seqgroups.append(group)
        if members is not None:
            [group.add(member) for member in members]

    def remove_group(self, group):
        '''Removes a group.

        @param group: the new group
        @type group: a ClansSeqgroup instance

        @raise ValueError: if <group> is no ClansSeqgroup instance
        '''
        if not isinstance(group, ClansSeqgroup):
            raise ValueError('groups need to be ClansSeqgroup instances')

        self.seqgroups.remove(group)
        [group.remove(member) for member in group.members]

    def add_entry(self, entry):
        '''Adds an new entry.

        @param entry: the new entry
        @type entry: a ClansEntry instance

        @raise ValueError: if <entry> is no ClansEntry instance
        '''
        if not isinstance(entry, ClansEntry):
            raise ValueError('entries need to be ClansEntry instances')

        self.entries.append(entry)
        entry.parent = self

        self._has_good_index = False

    def remove_entry(self, entry):
        '''Removes an entry.

        @param entry: the entry that shall be removed
        @type entry: string or ClansEntry instance

        @raise ValueError: if entry is neither string nor ClansEntry instance
        '''
        if isinstance(entry, basestring):
            entry = self.get_entry(entry, True)

        if not isinstance(entry, ClansEntry):
            raise ValueError('only string and ClansEntry are supported values')

        for other_entry in entry.hsp.keys():
            other_entry.remove_hsp(entry)

        for g in entry.groups:
            g.remove(entry)

        remove_groups = [g for g in self.seqgroups if g.is_empty()]
        [self.seqgroups.remove(g) for g in remove_groups]

        self.entries.remove(entry)
        self._has_good_index = False

    def get_entry(self, name, pedantic=True):
        '''Checks if an entry with name <name> exists and returns it. If
        multiple ones exists, one is returned.

        @param name: name of the sought entry
        @type name: string

        @param pedantic: If True, a ValueError is raised if multiple entries
                         with name <name> are found. If False, returns the list
                         containing all matching entries.
        @type pedantic: bool

        @raise ValueError: if no entry with name <name> is found
        @raise ValueError: if multiple entries with name <name> are found and
                           pedantic==True

        @rtype: L{ClansEntry}
        @return: entry with name <name>
        '''

        hits = [e for e in self.entries if e.name == name]

        if len(hits) >= 1:
            if pedantic:
                raise ValueError('multiple entries have name \'{}\''.format(
                    name))
            return hits

        elif len(hits) == 0:
            raise ValueError('ClansEntry %s does not exist.' % name)

        return hits

    def remove_duplicates(self, return_removed=False, cmp_function=None):
        '''Removes entries with identical unique_ids (see
        ClansEntry._get_unique_id)

        @param return_removed: if True returns a list with all removes entries.
        @type return_removed: bool

        @param cmp_function:
        @type cmp_function: function
        '''
        if cmp_function is None:
            cmp_function = gi_comparator()

        remove_us = set([e2 for i, e in enumerate(self.entries)
                         for e2 in self.entries[i + 1:]
                         if cmp_function(e, e2) == 0])

        [self.remove_entry(e) for e in remove_us]
        if return_removed:
            return remove_us

    def restrict_to_max_pvalue(self, cutoff, return_removed=False):
        ## loop to hit entries that have no HSPs left after the previous round
        removed_entries = []  # all removed entries go here
        remove_us = ['first_loop_round_starter']
        while len(remove_us) > 0:

            remove_us = []  # entries removed this round
            for entry in self.entries:
                hsp_values = entry.hsp.values()
                if len(hsp_values) == 0 or min(hsp_values) >= cutoff:
                    remove_us.append(entry)
                    removed_entries.append(entry)

            [self.remove_entry(e) for e in remove_us]

        if return_removed:
            return removed_entries

    def restrict(self, keep_names):
        [self.remove_entry(entry) for entry in
         [e for e in self.entries if e.name not in keep_names]]

    def write(self, output_filename):
        self._update_index()
        ClansWriter(self, output_filename)


class ClansEntry(object):
    '''
    Class holding the data of one CLANS sequence entry.
    '''

    def __init__(self, name=None, seq='', coords=None, hsp=None, groups=None,
                 parent=None):
        self.name = name
        self.seq = seq

        if coords is None:
            coords = random.random(3) * 2 - 1  # each CLANS coord is -1.<x<1.
        self.coords = coords

        if groups is None:
            groups = []
        self.groups = groups

        self.parent = parent

        if hsp is None:
            hsp = {}
        self.hsp = hsp

    def __repr__(self):
        if self.coords is None:
            coords_string = 'NoCoordsSet'
        else:
            coords_string = '(%.2f, %.2f, %.2f)' % tuple(self.coords)

        groups = 'not in a group'
        if len(self.groups) > 0:
            groups = 'groups: %s' % (', '.join([g.name for g in self.groups]))

        seq = 'NoSequenceSet'
        if self.seq is not '':
            seq = 'seq: %s' % self.seq

        return 'ClansEntry "%s": %s ' % \
               (self.name, '; '.join((seq, coords_string, groups)))

    def get_id(self):
        '''Returns the id of the current entry.

        @rtype: str
        @return: the entrys\' id is returned unless it has no parent in which
        case -1 is returned
        '''
        if self.parent is None:
            return -1

        if self.parent._has_good_index:
            return self.parent._idx[self._get_unique_id()]

        return self.parent.entries.index(self)

    def _get_unique_id(self):
        '''Returns a >>more<< unique ID (however this is not guaranteed to be
        really unique) than get_id. This ID determines which entries are deemed
        duplets by remove_duplicates of class Clans.

        @rtype: str
        @return: a more or less unique id
        '''
        return self.name + '<###>' + self.seq

    def add_hsp(self, other, value):
        '''Creates an HSP from self to other with the given value.'''
        self.hsp[other] = value
        other.hsp[self] = value

    def remove_hsp(self, other):
        '''Removes the HSP between self and other; if none exists, does
        nothing.
        '''
        if other in self.hsp:
            self.hsp.pop(other)

        if self in other.hsp:
            other.hsp.pop(self)

    def add_group(self, group):
        '''Adds the group to the list of groups the entry belongs to.'''
        if self.groups.count(group) == 0:
            self.groups.append(group)

    def remove_group(self, group):
        '''Removes the group from the list of groups the entry belongs to.'''
        if self.groups.count(group) == 0:
            return

        self.groups.remove(group)

    def output_string_seq(self):
        '''Creates the CLANS <seq> block format representation of the entry.

        @rtype: str
        @return: entrys\' representation in CLANS <seq> block format
        '''

        return '>%s\n%s\n' % (self.name, self.seq)

    def output_string_pos(self):
        '''Create the CLANS <pos> block format representation of the entry.

        @rtype: str
        @return: entrys\' representation in CLANS <pos> block format
        '''
        return '%i %.8f %.8f %.8f' % tuple([self.get_id()] + list(self.coords))

    def output_string_hsp(self):
        '''Creates the CLANS <hsp> block format representation of the entry.


        @rtype: str
        @return: entrys\' representation in CLANS <hsp> block format
        '''
        return '\n'.join(['%i %i:%.8f' %
                          (self.get_id(), other.get_id(), value)
                          for (other, value) in self.hsp.items()])


class ClansSeqgroup(object):
    '''
    Class holding the data of one CLANS group (seqgroup).
    '''

    def __init__(self, name, **kwargs):
        self.name = name
        self.type = kwargs.pop('type', 0)
        self.size = kwargs.pop('size', 4)
        self.hide = kwargs.pop('hide', 0)

        self.color = RGB_color()
        self.color.set_color('rgb', kwargs.pop('color', '255;255;255'))

        self.members = kwargs.pop('members', [])

    def __repr__(self):
        return 'ClansSeqgroup %s: type: %s; size: %s; ' % \
               (self.name, self.type, self.size) + \
               'hide: %s; color: %s; #members: %i' % \
               (self.hide, self.color.rgb, len(self.members))

    def is_empty(self):
        '''Checks if the group contains entries.

        @rtype: bool
        @return: True if the group contains no entries, False else.
        '''
        return len(self.members) == 0

    def add(self, new_member):
        '''Adds entry <new_member> to this group.

        @param new_member: the member that shall be added to this group
        @type new_member: a ClansEntry instance

        @raise ValueError: if <new_member> is no ClansEntry instance
        '''
        if not isinstance(new_member, ClansEntry):
            raise ValueError('only ClansEntry instances can be added as ' +
                             'group members')

        if self.members.count(new_member) == 1:
            return

        self.members.append(new_member)
        new_member.groups.append(self)

    def set_color(self, colorname, color):
        '''Set the color of this group to a new value.

        @param colorname: identifier of the color that will be changed
        @type colorname: a valid color identifier in terms of the RGB_color
        class

        @param color: the new color value
        @type color: a valid color in terms of the RGB_color class
        '''
        self.color.set_color(colorname, color)

    def remove(self, member):
        '''Removes ClansEntry \'member\' from this group. Does nothing if
        \'member\' is not contained in group self.

        @param member: the member to be removed
        @type member: a ClansEntry instance

        @raise ValueError: if member is no ClansEntry instance
        '''
        if not isinstance(member, ClansEntry):
            raise ValueError('argument must be a ClansEntry instance')

        if self.members.count(member) == 0:
            return

        self.members.remove(member)
        member.groups.remove(self)

    def output_string(self):
        '''Creates the CLANS <seqgroup> block format representation of the
        group.

        @rtype: str
        @return: entrys\' representation in CLANS <seqgroup> block format
        '''
        sorted_members = sorted([m.get_id() for m in self.members])
        return 'name=%s\ntype=%s\nsize=%s\nhide=%s\ncolor=%s\nnumbers=%s' % \
               (self.name, self.type, self.size, self.hide, self.color.rgb,
                ';'.join([str(val) for val in sorted_members]) + ';')


def transfer_groups(origin, target):
    '''Transfers the CLANS group definitions from origin to target by comparing
    entrynames. If a group from origin contains no entries in target, it is not
    created.

    @param origin: group definitions of this are added to target
    @type origin: Clans instance

    @param target: group definitions are added to this
    @type target: Clans instance
    '''
    for group in origin.seqgroups:

        new_group = ClansSeqgroup(name=group.name,
                                  type=group.type,
                                  size=group.size,
                                  hide=group.hide,
                                  color=group.color)

        for member in group.members:
            try:
                new_member = target.get_entry(member.name)
                if isinstance(new_member, tuple):
                    raise ValueError('')
            except ValueError:
                continue

            new_group.add(new_member)

        if len(new_group.members) > 0:
            target.add_group(new_group)

if __name__ == '__main__':
    if len(sys.argv) not in [4, 5]:
        print 'transfer of group definitions from source to target CLANS map.'
        print 'usage: python clans.py source-file target-file output-file ' + \
              ' [rmd]'
        print 'rmd is optional and will use the first gi number (with'
        print 'annotated residues, or the whole name if no gi is found)'
        print 'to determine duplicates and remove them.'

    else:

        cp = ClansParser()

        print 'reading source...'
        source_fn = sys.argv[1]
        source = cp.parse_file(source_fn)
        print 'source is', source

        print 'reading target...'
        target_fn = sys.argv[2]
        target = cp.parse_file(target_fn)
        print 'target is', target

        print 'transferring groups...'
        transfer_groups(source, target)
        print 'target with groups is', target

        if len(sys.argv) > 4 and sys.argv[4] == 'rmd':
            rm = target.remove_duplicates(True)
            print 'removed %i sequences as duplicates' % len(rm)
            print 'target after duplicate removal is', target

        print 'writing file...'
        output_fn = sys.argv[3]
        target.write(output_fn)

    if False:

        fn = '/tmp/cl/2be1_psiblast.clans'
        fn =  '~/tmp/20000.clans'
        cp = ClansParser()
        c = cp.parse_file(fn)
