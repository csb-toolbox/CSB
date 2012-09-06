"""
Classes for parsing/manipulating/writing CLANS (by Tancred Frickey) files

This module defines L{ClansParser} and L{ClansFileWriter} for parsing and writing CLANS format files, respectively.
Further, class L{Clans} and several helper classes are used to hold and handle the parsed data.

The most commenly used CLANS data can be accessed in an L{Clans} instance via
 - .entries <L{ClansEntryCollection} containing L{ClansEntry} instances>
  - .name
  - .seq <the amino acid sequence>
  - .hsps <connections of this L{ClansEntry} to others>
  - .groups <L{ClansSeqgroup}s the entry belongs to>
 - .seqgroups <L{ClansSeqgroupCollection} containing L{ClansSeqgroup} instances>
 - .params <L{ClansParams}>

Parse a file into L{Clans} instance C{clans_instance} by
    >>> clans_instance = ClansParser().parse_file('input.clans')

Create a new entry C{e} with name \"C{my entry}\", sequence C{AAAA} and coordinates C{(x=1, y=1, z=1)}
    >>> e = ClansEntry(name='my entry', seq='AAAA', coords=(1, 1, 1))

and add it to an existing L{Clans} instance C{clans_instance}
    >>> clans_instance.add_entry(e)

Entries can be accessed using indices of C{Clans} instances
    >>> clans_instance[0]  # access to first entry

and deleted by
    >>> clans_instance.remove_entry(e)

Equivalent functions exist for ClansSeqgroups.

Author: Klaus Kopec
MPI fuer Entwicklungsbiologie, Tuebingen
"""

import os
import re
import operator
import csb.core

from abc import ABCMeta, abstractmethod
from numpy import array, float64, eye, random


class DuplicateEntryNameError(Exception):
    """
    Raised during L{Clans.get_entry} if two entries have identical names.
    """
    pass


class DuplicateEntryError(Exception):
    """
    Raised during L{Clans._update_index} if two entries are identical in name, sequence, and coordinates.
    """
    pass
    

class MissingBlockError(Exception):
    """
    Raised if an expected tag is not found during parsing of a CLANS file.
    """
    pass


class UnknownTagError(ValueError):
    """
    Raised if an unknown tag is encountered while parsing a CLANS file.
    """
    pass


class Color(object):
    """
    RGB color handling class.
    Color is stored as r, g, b, and a (i.e. alpha) attributes.
    Default color is C{r}=C{g}=C{b}=0 (i.e. black) with a=255

    @param r: the red value
    @type r: int

    @param g: the green value
    @type g: int

    @param b: the blue value
    @type b: int

    @param a: the alpha value
    @type a: int
    """

    def __init__(self, r=0, g=0, b=0, a=255):
        self._r = None
        self.r = r
        self._g = None
        self.g = g
        self._b = None
        self.b = b
        self._a = None
        self.a = a

    def __repr__(self):
        return 'Color {0}'.format(self.to_clans_color())

    __str__ = __repr__

    @staticmethod
    def from_string(color_string, separator=';'):
        """
        Factory for a Color instance created from a string formatted as r{separator}g{separator}b{separator}a, where the final \'{separator}a\' is optional.

        @param color_string: the color string
        @type color_string: str

        @raises TypeError: if C{color_string} is not a string
        @raises ValueError: if any value in color is outside of range(256)
        """
        if not isinstance(color_string, csb.core.string):
            raise TypeError('{0} is no string'.format(color_string))
        
        if color_string.count(separator) == 2:
            r, g, b = map(int, color_string.split(';'))
            a = 255
        elif color_string.count(separator) == 3:
            r, g, b, a  = map(int, color_string.split(';'))

        else:
            raise ValueError(
                ('format needs to be \'r{0}g{0}b\' but color_string was ' +
                 '{1} [optionally with alpha value: \'r{0}g{0}b{0}a\']').format(separator, color_string))
            
        return Color(r, g, b, a)

    @property
    def r(self):
        """
        the red value of the RGB color.
        
        raises ValueError if C{value} is outside of range(256)

        @rtype: int
        """
        return self._r

    @r.setter
    def r(self, value):
        """
        Set the red value of the RGB color.
        """
        if value < 0 or value > 255:
            raise ValueError(
                'valid color values are in range(256), was \'{0}\''.format(
                    value))

        self._r = value

    @property
    def g(self):
        """
        the green value of the RGB color.

        raises ValueError if C{value} is outside of range(256)

        @rtype: int
        """
        return self._g

    @g.setter
    def g(self, value):

        if value < 0 or value > 255:
            raise ValueError('valid color values are in range(256).')

        self._g = value

    @property
    def b(self):
        """
        the blue value of the RGB color.

        raises ValueError if C{value} is outside of range(256)

        @rtype: int
        """
        return self._b

    @b.setter
    def b(self, value):

        if value < 0 or value > 255:
            raise ValueError('valid color values are in range(256).')

        self._b = value

    @property
    def a(self):
        """
        the alpha value of the RGB color.

        raises ValueError if C{value} is outside of range(256)

        @rtype: int
        """
        return self._a

    @a.setter
    def a(self, value):

        if value < 0 or value > 255:
            raise ValueError('valid color values are in range(256).')

        self._a = value

    def to_clans_color(self):
        """
        Formats the color for use in CLANS files.

        @return: the color formatted for use in CLANS files; format: r;g;b;a
        @rtype: str
        """
        return '{0.r};{0.g};{0.b};{0.a}'.format(self)


class ClansParser(object):
    """
    CLANS file format aware parser.
    """

    def __init__(self):
        self._clans_instance = None
        self._data_block_dict = {}

    def __repr__(self):
        return 'ClansParser instance'

    __str__ = __repr__

    @property
    def clans_instance(self):
        """
        the L{Clans} instance that resulted from parsing a CLANS file.

        raises a ValueError if no CLANS file has been parsed yet

        @rtype: L{Clans} instance
        """
        if self._clans_instance is None:
            raise ValueError('you need to parse a CLANS file first')

        return self._clans_instance

    def parse_file(self, filename, permissive=True):
        """
        Create a L{Clans} instance by parsing the CLANS format file C{filename}

        @param filename: name of the CLANS file.
        @type filename: str

        @param permissive: if True, tolerate missing non-essential or unknown
        blocks.
        @type permissive: bool

        @rtype: L{Clans} instance
        @return: a L{Clans} instance containing the parsed data

        @raise MissingBlockError: if C{permissive == True} and any essential
        block is missing. if C{permissive == False} and any block is missing
        @raise UnknownTagError: if C{permissive == False} and an unknown tag/
        data block is encountered
        """
        self._clans_instance = Clans()
        self._clans_instance._filename = filename

        self._read_block_dict()  # read and preprocess the CLANS file

        try:  # param and rotmtx are non-essential blocks
            self._parse_param()
            self._parse_rotmtx()
        except MissingBlockError as error:
            if not permissive:
                raise MissingBlockError(error)

        seq = {}
        try:
            seq = self._parse_seq()
        except MissingBlockError as error:
            if not permissive:
                raise MissingBlockError(error)

        seqgroups = self._parse_seqgroups()

        pos = {}
        try:
            pos = self._parse_pos()
        except MissingBlockError as error:
            if not permissive:
                raise MissingBlockError(error)

        hsp_att_mode = "hsp"
        hsp = {}
        try:
            if 'hsp' in self._data_block_dict:
                hsp = self._parse_hsp_att('hsp')

            elif 'att' in self._data_block_dict:
                hsp_att_mode = "att"
                hsp = self._parse_hsp_att('att')

            elif 'mtx' in self._data_block_dict:
                hsp = self._parse_mtx()

        except MissingBlockError as error:
            if not permissive:
                raise MissingBlockError(error)

        ## raise UnknownTagError for unknown blocks
        known_block_tags = set(('param', 'rotmtx', 'seq', 'seqgroups', 'pos',
                                'hsp', 'mtx', 'att'))
        unprocessed_block_tags = set(self._data_block_dict.keys()).difference(
            known_block_tags)

        if len(unprocessed_block_tags) > 0 and not permissive:
            raise UnknownTagError(
                ('tags unknown: {0}. File corrupt or further implementations '
                 + 'needed!').format(', '.join(unprocessed_block_tags)))

        ## if no entries exist, we cannot add pos, seqgroup and hsp data
        if len(seq) > 0:

            ## add Entries
            if len(pos) > 0:
                for i in pos:
                    self._clans_instance.add_entry(ClansEntry(seq[i][0], seq[i][1],
                               pos[i], parent=self._clans_instance))

            ## add groups
            self._clans_instance._seqgroups = ClansSeqgroupCollection()
            if len(seqgroups) > 0:
                for group_raw_data in seqgroups:

                    group = ClansSeqgroup(name=group_raw_data['name'],
                                          type=group_raw_data['type'],
                                          size=group_raw_data['size'],
                                          hide=group_raw_data['hide'] == '1',
                                          color=group_raw_data['color'])

                    ## get members corresponding to the IDs in this group
                    ## NOTE: this silently corrects files where a seqgroup
                    ##       contains the same entry multiple times
                    members = [self._clans_instance.entries[number]
                               for number in set(group_raw_data['numbers'])]

                    self._clans_instance.add_group(group, members)

            ## add hsp values
            if len(hsp) > 0:
                [self._clans_instance.entries[a].add_hsp(
                    self._clans_instance.entries[b], value)
                 for ((a, b), value) in hsp.items()]

        self._clans_instance._hsp_att_mode = hsp_att_mode

        return self._clans_instance

    def _read_block_dict(self):
        """
        Extracts all <tag>DATA</tag> blocks from file
        self.clans_instance.filename.

        @rtype: dict
        @return: data in the form: dict[tag] = DATA.
        """
        # read file and remove the first line, i.e. sequence=SEQUENCE_COUNT
        data_blocks = open(os.path.expanduser(
            self._clans_instance.filename)).read().split('\n', 1)[1]

        ## flag re.DOTALL is necessary to make . match newlines
        data = re.findall(r'(<(\w+)>(.+)</\2>)', data_blocks,
                          flags=re.DOTALL)
        self._data_block_dict = dict([(tag, datum.strip().split('\n'))
                                     for _tag_plus_data, tag, datum in data])

    def _parse_param(self):
        """
        Parse a list of lines in the CLANS <param> format:

        parameter1=data1\n
        parameter2=data2\n
        ...
        """
        if 'param' not in self._data_block_dict:
            raise MissingBlockError('file contains no <param> block.')

        block = self._data_block_dict['param']

        tmp_params = dict([block[i].split('=') for i in range(len(block))])

        ## create colors entry from colorcutoffs and colorarr
        colorcutoffs = [float(val) for val in
                        tmp_params.pop('colorcutoffs').strip(';').split(';')]
        colors = tmp_params.pop('colorarr').strip(':')
        colors = colors.replace('(', '').replace(')', '').split(':')
        colorarr = [Color(*map(int, color_definition)) for color_definition in
                    [color.split(';') for color in colors]]

        tmp_params['colors'] = tuple(zip(colorcutoffs, colorarr))

        ## convert 'true' and 'false' into Python bools
        for k, v in tmp_params.items():
            if v == 'true':
                tmp_params[k] = True
            elif v == 'false':
                tmp_params[k] = False
                
        self._clans_instance._params = ClansParams(strict=False, **tmp_params)

    def _parse_rotmtx(self):
        """
        Parse a list of lines in the CLANS <rotmtx> format. The data is stored
        in the clans_instance as a 3x3 numpy.array.

        @raise ValueError: if the rotmtx block does not contain exactly 3 lines
        """
        if 'rotmtx' not in self._data_block_dict:
            raise MissingBlockError('file contains no <rotmtx> block.')

        block = self._data_block_dict['rotmtx']

        if len(block) != 3:
            raise ValueError('CLANS <rotmtx> blocks comprise exactly 3 lines.')
        self._clans_instance.rotmtx = array(
            [[float64(val) for val in line.split(';')[:3]] for line in block])

    def _parse_seq(self):
        """
        Parse a list of lines in the CLANS <seq> format, which are in FASTA
        format.

        @rtype: dict
        @return: dict with running numbers as key and 2-tuples (id, sequence)
                 as values
        """
        if 'seq' not in self._data_block_dict:
            raise MissingBlockError(
                'file contains no <seq> block. This is OK if the file does '
                + 'not contain any sequences.')

        block = self._data_block_dict['seq']
        if len(block) % 2 == 1:
            block += ['']

        return dict([(i, (block[2 * i][1:], block[2 * i + 1].strip()))
                     for i in range(int(len(block) / 2))])

    def _parse_seqgroups(self):
        """
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
        """
        if 'seqgroups' not in self._data_block_dict:
            return ClansSeqgroupCollection()

        block = self._data_block_dict['seqgroups']

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
        """
        Parse a list of lines in the CLANS <pos> format \'INT FLOAT FLOAT
        FLOAT\'.

        @rtype: dict
        @return: a dict using the integers as keys and a (3,1)-array created
                 from the three floats as values.
        """
        if 'pos' not in self._data_block_dict:
            raise MissingBlockError(
                'file contains no <pos> block. This is OK if the file does '
                + 'not contain any sequences.')

        block = self._data_block_dict['pos']

        return dict([(int(l.split()[0]),
                      array([float64(val) for val in l.split()[1:]]))
                     for l in block])

    def _parse_hsp_att(self, mode):
        """
        Parse a list of lines in the CLANS <hsp> format \'INT INT: FLOAT\'.

        NOTE: some CLANS <hsp> lines contain more than one float; we omit the
        additional numbers

        @param mode: either "hsp" or "att" depending on the type of tag to be
        parsed
        @type mode: str

        @rtype: dict
        @return: a dict using 2-tuples of the two integers as keys and the
                 float as values
        """
        if mode not in ("hsp", "att"):
            raise ValueError('mode must be either "hsp" or "att"')

        if mode not in self._data_block_dict:
            raise MissingBlockError(
                ('file contains no <{0}> block. This is OK if the file does '
                 + 'not contain any sequences or if none of the contained '
                 + 'sequences have any connections.').format(mode))

        block = self._data_block_dict[mode]

        if mode == "hsp":
            return dict([(tuple([int(val)
                                 for val in line.split(':')[0].split()]),
                          float(line.split(':')[1].split(' ')[0]))
                         for line in block])

        else:
            return dict([(tuple([int(val) for val in line.split(' ')[:2]]),
                          float(line.split(' ')[2]))
                         for line in block])

    def _parse_mtx(self):
        """
        Parse a list of lines in the CLANS <mtx> format.

        @rtype: dict
        @return: a dict using 2-tuples of the two integers as keys and the
                 float as values
        """
        if 'mtx' not in self._data_block_dict:
            raise MissingBlockError(
                'file contains no <mtx> block. This is OK if the file does '
                + 'not contain any sequences or if none of the contained '
                + 'sequences have any connections.')

        block = self._data_block_dict['mtx']

        return dict([((i, j), float(entry))
                     for i, line in enumerate(block)
                     for j, entry in enumerate(line.split(';')[:-1])
                     if float(entry) != 0])

class ClansFileBuilder(object):
    """
    Base abstract files for building a file in CLANS format.
    Defines a common step-wise interface according to the Builder pattern.
    
    @param output: output stream (this is where the product is constructed)
    @type output: stream
    """
    
    __metaclass__ = ABCMeta

    def __init__(self, output):

        if not hasattr(output, 'write'):
            raise TypeError(output)
        
        self._out = output

    @property
    def output(self):
        """
        Destination stream
        @rtype: stream
        """
        return self._out
    
    def write(self, text):
        """
        Write a chunk of text
        """
        self._out.write(text)   

    def writeline(self, text):
        """
        Write a chunk of text and append a new line terminator
        """        
        self._out.write(text)
        self._out.write('\n')
                    
    @abstractmethod
    def add_param_block(self, block_data):
        pass

    @abstractmethod
    def add_rotmtx_block(self, block_data):
        pass

    @abstractmethod
    def add_seq_block(self, block_data):
        pass

    @abstractmethod
    def add_seqgroups_block(self, block_data):
        pass

    @abstractmethod
    def add_pos_block(self, block_data):
        pass

    @abstractmethod
    def add_hsp_block(self, block_data):
        pass
    

class ClansFileWriter(ClansFileBuilder):
    """
    Class for serializing a L{Clans} instance to a file in CLANS format.

    @param output: the output stream
    @type output: stream
    """


    def __init__(self, output):
        super(ClansFileWriter, self).__init__(output)

    def serialize(self, clans_instance):
        """
        Creates a CLANS file containing all data from C{clans_instance}

        @param clans_instance: the source of the data to be serialized
        @type clans_instance: a L{Clans} instance
        """

        self.add_sequences_line(clans_instance)
        self.add_param_block(clans_instance)
        self.add_rotmtx_block(clans_instance)
        self.add_seq_block(clans_instance)
        self.add_seqgroups_block(clans_instance)
        self.add_pos_block(clans_instance)
        self.add_hsp_block(clans_instance)
        
    def add_sequences_line(self, clans_instance):
        """
        Appends the \'sequences=<#sequences>\' line to the output.

        @param clans_instance: the source of the data to be serialized
        @type clans_instance: a L{Clans} instance
        """
        self.writeline('sequences={0}'.format(len(clans_instance.entries)))

    def add_param_block(self, clans_instance):
        """
        Appends a <param>data</param> CLANS file block to the output.

        @param clans_instance: the source of the data to be serialized
        @type clans_instance: a L{Clans} instance
        """
        param_block = clans_instance.params._to_clans_param_block()
        self.write(param_block)

    def add_rotmtx_block(self, clans_instance):
        """
        Appends a <rotmtx>data</rotmtx> CLANS file block to the output.

        @param clans_instance: the source of the data to be serialized
        @type clans_instance: a L{Clans} instance

        @raise ValueError: if self.clans_instance.rotmtx is no 3x3 numpy.array
        """
        rotmtx = clans_instance.rotmtx

        if rotmtx is None:
            return

        if rotmtx.shape != (3, 3):
            raise ValueError('rotmtx must be a 3x3 array')

        self.writeline('<rotmtx>')

        self.write('\n'.join(
            ['{0};{1};{2};'.format(*tuple(rotmtx[i])) for i in range(3)]))
        self.write('\n')

        self.writeline('</rotmtx>')

    def add_seq_block(self, clans_instance):
        """
        Appends a <seq>data</seq> CLANS file block to the output.

        @param clans_instance: the source of the data to be serialized
        @type clans_instance: a L{Clans} instance
        """
        self.writeline('<seq>')

        self.write(''.join([e.output_string_seq()
                            for e in clans_instance.entries]))

        self.writeline('</seq>')

    def add_seqgroups_block(self, clans_instance):
        """
        Appends a <seqgroupsparam>data</seqgroups> CLANS file block to the output.

        @param clans_instance: the source of the data to be serialized
        @type clans_instance: a L{Clans} instance
        """
        seqgroups = clans_instance.seqgroups

        if seqgroups is not None and len(seqgroups) > 0:

            self.writeline('<seqgroups>')

            self.write('\n'.join([s.output_string() for s in seqgroups]))
            self.write('\n')

            self.writeline('</seqgroups>')

    def add_pos_block(self, clans_instance):
        """
        Appends a <pos>data</pos> CLANS file block to the output.

        @param clans_instance: the source of the data to be serialized
        @type clans_instance: a L{Clans} instance
        """
        self.writeline('<pos>')

        self.write('\n'.join([e.output_string_pos()
                              for e in clans_instance.entries]))
        self.write('\n')

        self.writeline('</pos>')

    def add_hsp_block(self, clans_instance):
        """
        Appends a <hsp>data</hsp> CLANS file block to the output.
        If the CLANS instance has hsp_att_mode=="att" we add a <att>data<att>
        block which has the same format.

        @param clans_instance: the source of the data to be serialized
        @type clans_instance: a L{Clans} instance
        """

        self.writeline('<{0}>'.format(clans_instance._hsp_att_mode))
        
        ## sorting is not necessary, but makes a nicer looking clans file
        idToEntryMapping = [(e.get_id(), e)
                            for e in clans_instance.entries]
        idToEntryMapping.sort(key=operator.itemgetter(0))
        entryToIdMapping = dict([(entry, identifier)
                                 for (identifier, entry) in idToEntryMapping])

        for i, (entry1_id, entry1) in enumerate(idToEntryMapping):

            ## sort list of hsp targets by id
            hspTargets = [(entryToIdMapping[entry2], pvalue)
                          for (entry2, pvalue) in entry1.hsp.items()]
            hspTargets.sort(key=operator.itemgetter(0))

            for (entry2_id, pvalue) in hspTargets:
                if entry1_id >= entry2_id:
                    continue

                line_format = '{0} {1}:{2}\n'
                if clans_instance._hsp_att_mode == "att":
                    line_format = '{0} {1} {2}\n'

                self.write(
                    line_format.format(entry1_id, entry2_id, repr(pvalue)))

        self.writeline('</{0}>'.format(clans_instance._hsp_att_mode))


class ClansParams(object):
    """
    Class for handling L{Clans} parameters.
    See L{ClansParams}._DEFAULTS for accepted parameter names.

    @kwparam **kw: parameters as C{kw[parameter_name] = parameter_value}

    @raise KeyError: if a supplied parameter name is not known
    (i.e. it is not a key in _DEFAULTS)
    """

    _DEFAULTS = {'attfactor': 10.0,
                 'attvalpow': 1,
                 'avgfoldchange': False,
                 'blastpath': 'blastall -p blastp',
                 'cluster2d': False,
                 'colors': ((0.0, (230, 230, 230)),
                            (0.1, (207, 207, 207)),
                            (0.2, (184, 184, 184)),
                            (0.3, (161, 161, 161)),
                            (0.4, (138, 138, 138)),
                            (0.5, (115, 115, 115)),
                            (0.6, (92, 92, 92)),
                            (0.7, (69, 69, 69)),
                            (0.8, (46, 46, 46)),
                            (0.9, (23, 23, 23))),
                'complexatt': True,
                'cooling': 1.0,
                'currcool': 1.0,
                'dampening': 0.2,
                'dotsize': 2,
                'formatdbpath': 'formatdb',
                'groupsize': 4,
                'maxmove': 0.1,
                'minattract': 1.0,
                'ovalsize': 10,
                'pval': 1.0,
                'repfactor': 5.0,
                'repvalpow': 1,
                'showinfo': True,
                'usefoldchange': False,
                'usescval': False,
                'zoom': 1.0}

    def __init__(self, strict=True, **kw):
        self.set_default_params()

        for param_name, param_value in kw.items():
            if param_name not in self._DEFAULTS and strict:
                raise KeyError('parameter {0} (value: {1}) unknown'.format(
                    param_name, param_value))
            self.__setattr__(param_name, param_value)

    @property
    def complexatt(self):
        """
        if True, complex attraction computations are used.

        raises ValueError if set to non-boolean value

        @rtype: bool
        """
        return self._complexatt

    @complexatt.setter
    def complexatt(self, value):
        if not isinstance(value, bool):
            raise ValueError(('complexatt cannot be {0} (accepted values: True'
                              + '/False)').format(value))
        self._complexatt = value

    @property
    def attfactor(self):
        """
        factor in the attractive force

        raises ValueError if C{value} is not castable to float
        
        @rtype: float
        """
        return self._attfactor

    @attfactor.setter
    def attfactor(self, value):
        self._attfactor = float(value)
            
    @property
    def attvalpow(self):
        """
        exponent in the attractive force

        raises ValueError if C{value} is not castable to float

        @rtype: float
        """
        return self._attvalpow

    @attvalpow.setter
    def attvalpow(self, value):
        self._attvalpow = float(value)

    @property
    def repfactor(self):
        """
        factor in the repulsive force

        raises ValueError if C{value} is not castable to float
        
        @rtype: float
        """
        return self._repfactor

    @repfactor.setter
    def repfactor(self, value):
        self._repfactor = float(value)

    @property
    def repvalpow(self):
        """
        exponent in the repulsive force

        raises ValueError if C{value} is not castable to float
        
        @rtype: float
        """
        return self._repvalpow

    @repvalpow.setter
    def repvalpow(self, value):
        self._repvalpow = float(value)

    @property
    def cluster2d(self):
        """
        if True, clustering is done in 2D. Else in 3D.

        raises ValueError if set to non-boolean value
                
        @rtype: bool
        """
        return self._cluster2d


    @cluster2d.setter
    def cluster2d(self, value):
        if not isinstance(value, bool):
            raise ValueError(('cluster2d cannot be {0} (accepted values: True'
                              + '/False)').format(value))

        self._cluster2d = value

    @property
    def pval(self):
        """
        p-value cutoff that determines which connections are considered for
        the attractive force
        
        raises ValueError if C{value} is not castable to float
        
        @rtype: float
        """
        return self._pval

    @pval.setter
    def pval(self, value):
        self._pval = float(value)

    @property
    def maxmove(self):
        """
        maximal sequence (i.e. dot in the clustermap) movement per round

        raises ValueError if C{value} is not castable to float

        @rtype: float
        """
        return self._maxmove

    @maxmove.setter
    def maxmove(self, value):
        self._maxmove = float(value)

    @property
    def usescval(self):
        """
        parameter with unclear function. Check in Clans.

        raises ValueError if set to non-boolean value

        @rtype: bool
        """
        return self._usescval

    @usescval.setter
    def usescval(self, value):
        if not isinstance(value, bool):
            raise ValueError(('usescval cannot be {0} (accepted values: True'
                              + '/False)').format(value))

        self._usescval = value

    @property
    def cooling(self):
        """
        parameter  with unclear function. Check in Clans.

        raises ValueError if C{value} is not castable to float

        @rtype: float
        """
        return self._cooling

    @cooling.setter
    def cooling(self, value):
        self._cooling = float(value)

    @property
    def currcool(self):
        """
        parameter  with unclear function. Check in Clans.

        raises ValueError if C{value} is not castable to float

        @rtype: float
        """
        return self._currcool

    @currcool.setter
    def currcool(self, value):
        self._currcool = float(value)

    @property
    def dampening(self):
        """
        parameter  with unclear function. Check in Clans.

        raises ValueError if C{value} is not castable to float

        @rtype: float
        """
        return self._dampening

    @dampening.setter
    def dampening(self, value):
        self._dampening = float(value)

    @property
    def minattract(self):
        """
        parameter  with unclear function. Check in Clans.

        raises ValueError if C{value} is not castable to float

        @rtype: float
        """
        return self._minattract

    @minattract.setter
    def minattract(self, value):
        self._minattract = float(value)

    @property
    def blastpath(self):
        """
        path to the BLAST executable for protein-protein comparisons. BLAST+ is
        currently not supported by Clans.

        raises ValueError if C{value} is not a string

        @rtype: str
        """
        return self._blastpath

    @blastpath.setter
    def blastpath(self, value):
        if not isinstance(value, csb.core.string):
            raise ValueError(('blastpath cannot be {0} (accepted values: '
                              + 'strings)').format(value))
        
        self._blastpath = value

    @property
    def formatdbpath(self):
        """
        path to the formatdb executable of BLAST.

        raises ValueError if C{value} is not a string

        @rtype: str
        """
        return self._formatdbpath

    @formatdbpath.setter
    def formatdbpath(self, value):
        if not isinstance(value, csb.core.string):
            raise ValueError(('formatdbpath cannot be {0} (accepted values: '
                              + 'strings)').format(value))
              
        self._formatdbpath = value

    @property
    def showinfo(self):
        """
        if True, additional data (rotation matrix) is shown in the clustring
        window)

        raises ValueError if set to non-boolean value

        @rtype: bool
        """
        return self._showinfo

    @showinfo.setter
    def showinfo(self, value):
        if not isinstance(value, bool):
            raise ValueError(('showinfo cannot be {0} (accepted values: True'
                              + '/False)').format(value))
        
        self._showinfo = value

    @property
    def zoom(self):
        """
        zoom value (1.0 == not zoomed)

        raises ValueError if C{value} is not castable to float

        @rtype: float
        """
        return self._zoom

    @zoom.setter
    def zoom(self, value):
        self._zoom = float(value)

    @property
    def dotsize(self):
        """
        size of the central dot representing each sequence in the clustermap

        raises ValueError if C{value} is not castable to int

        @rtype: int
        """
        return self._dotsize

    @dotsize.setter
    def dotsize(self, value):
        self._dotsize = int(value)

    @property
    def ovalsize(self):
        """
        size of the circle around selected sequences
        
        raises ValueError if value not castable to int

        @rtype: int
        """
        return self._ovalsize

    @ovalsize.setter
    def ovalsize(self, value):
        self._ovalsize = int(value)

    @property
    def groupsize(self):
        """
        default for the size of circles that mark newly created groups
        
        raises ValueError if C{value} is not castable to int

        @rtype: int
        """
        return self._groupsize

    @groupsize.setter
    def groupsize(self, value):
        self._groupsize = int(value)

    @property
    def usefoldchange(self):
        """
        parameter  with unclear function. Check in Clans.

        raises ValueError if set to non-boolean value

        @rtype: bool
        """
        return self._usefoldchange

    @usefoldchange.setter
    def usefoldchange(self, value):
        if not isinstance(value, bool):
            raise ValueError(('usefoldchange cannot be {0} (accepted values: '
                              + 'True/False)').format(value))
        
        self._usefoldchange = value

    @property
    def avgfoldchange(self):
        """
        parameter  with unclear function. Check in Clans.

        raises ValueError if set to non-boolean value

        @rtype: bool
        """
        return self._avgfoldchange

    @avgfoldchange.setter
    def avgfoldchange(self, value):
        if not isinstance(value, bool):
            raise ValueError(('avgfoldchange cannot be {0} (accepted values: '
                              + 'True/False)').format(value))
        
        self._avgfoldchange = value

    @property
    def colors(self):
        """
        colors that define the coloring for different p-values/attractions

        raises ValueError if set to s.th. else than a 10-tuple of 2-tuples

        @rtype: tuple
        """
        return self._colors

    @colors.setter
    def colors(self, value):
        if not isinstance(value, tuple):
            raise ValueError('colors must be a tuple')
        if len(value) != 10:
            raise ValueError('colors must be a 10-tuple')
        lengths = [len(v) for v in value]
        if len(set(lengths)) != 1 or lengths[0] != 2:
            raise ValueError('each item of colors must be a 2-tuple')
        self._colors = value

    def set_default_params(self):
        """
        Sets the parameters to CLANS default values.
        See L{ClansParams}._DEFAULTS.
        """
        for k, v in self._DEFAULTS.items():
            if k == 'colors':
                continue
            
            self.__setattr__(k, v)

        tmp_list = []
        for i, (cutoff, color) in enumerate(ClansParams._DEFAULTS['colors']):
            tmp_list.append((cutoff, Color(*color)))
        self.colors = tuple(tmp_list)

    def _to_clans_param_block(self):
        """
        Creates a param block for a CLANS file from the L{ClansParams} values.

        @return: a CLANS file format <param>[data]</param> block
        @rtype: str
        """

        param_dict = {}

        for param_name in sorted(ClansParams._DEFAULTS):
            if param_name == 'colors':

                ## divide 'colors' into 'colorcutoffs' and 'colorarr'
                param_dict['colorcutoffs'] = ''.join(
                    ['{0:.2f};'.format(cutoff) for cutoff, color in self.colors])

                param_dict['colorarr'] = ''.join(
                    ['({0}):'.format(color.to_clans_color())
                     for cutoff, color in self.colors])

                continue

            if param_name in ('avgfoldchange', 'cluster2d', 'complexatt',
                              'showinfo', 'usefoldchange', 'usescval'):
                param_dict[param_name] = ['false', 'true'][
                    self.__getattribute__(param_name)]

                continue

            param_dict[param_name] = self.__getattribute__(param_name)

        param_block_string = '<param>\n'
        param_block_string += '\n'.join(
            ['{0}={1}'.format(param_name, param_dict[param_name])
             for param_name in sorted(param_dict)])
        param_block_string += '\n</param>\n'

        return param_block_string


class ClansEntryCollection(csb.core.ReadOnlyCollectionContainer):
    """
    Read-only container for holding (and internally managing) L{ClansEntry}
    instances.
    """

    def __init__(self):
        
        super(ClansEntryCollection, self).__init__(type=ClansEntry)

    def _remove_item(self, item):
        """
        Removes {item} from the collection.

        @param item: the item to be removed
        @type item: a L{ClansEntry} instance

        @raises TypeError: if {item} is not a L{ClansEntry} instance
        """

        if self._type:
            if not isinstance(item, self._type):
                raise TypeError("Item {0} is not of the required {1} type.".format(
                    item, self._type.__name__))
        self._items.remove(item)

    def _sort(self, key=None):
        """
        Sort entries by their {name} or by a custom key function.

        @Note: If the L{ClansEntryCollection} is part of a L{Clans} instance,
        use L{Clans.sort} instead to avoid silently corrupting the index in
        L{Clans._idx}.

        @param key: None to sort by name, a custom key function else
        @type key: function
        """
        if key is None:
            key = lambda entry: entry.name

        self._items.sort(key=key)


class ClansSeqgroupCollection(csb.core.ReadOnlyCollectionContainer):
    """
    Read-only container for holding (and internally managing) L{ClansSeqgroup}
    instances.
    """
    
    def __init__(self):
        
        super(ClansSeqgroupCollection, self).__init__(type=ClansSeqgroup)

    def _remove_item(self, item):
        """
        Removes {item} from the collection.

        @param item: the item to be removed
        @type item: a L{ClansSeqgroup} instance

        @raises TypeError: if {item} is not a L{ClansSeqgroup} instance
        """

        if self._type:
            if not isinstance(item, self._type):
                raise TypeError("Item {0} is not of the required {1} type.".format(
                    item, self._type.__name__))
        self._items.remove(item)


class Clans(object):
    """
    Class for holding and manipulating data from one CLANS file.
    Initialization is always done as empty clustermap with default parameters.
    """

    def __init__(self):
        self._filename = None

        self._params = ClansParams()

        self._rotmtx = None
        self.set_default_rotmtx()

        self._hsp_att_mode = "hsp"

        self._entries = ClansEntryCollection()
        self._seqgroups = ClansSeqgroupCollection()

        self._idx = None
        '''Index dict for fast access to entry positions'''
        self._has_good_index = False

    def __repr__(self):
        return 'Clans object: {0} sequences; {1} seqgroups'.format(
            len(self), len(self.seqgroups))

    __str__ = __repr__

    def __len__(self):
        return len(self.entries)

    def __getitem__(self, index):
        return self.entries[index]

    def __setitem__(self, index, data):
        self.entries[index] = data
        self._has_good_index = False

    @property
    def filename(self):
        """
        file from which the data was parsed

        @rtype: str or None
        """
        return self._filename

    @property
    def params(self):
        """
        L{ClansParams} that contains the parameters set for this L{Clans}
        instance.

        @rtype: L{ClansParams}
        """
        return self._params

    @property
    def rotmtx(self):
        """
        3x3 rotation matrix that indicates the rotation state of the clustermap

        raises ValueError if rotation matrix shape is not 3x3

        @rtype: numpy.array
        """
        return self._rotmtx

    @rotmtx.setter
    def rotmtx(self, value):
        if value.shape != (3, 3):
            raise ValueError('rotation matrix needs to be a 3x3 numpy array')
        self._rotmtx = value
    
    @property
    def entries(self):
        """
        list of clustermap L{ClansEntry}s.

        @rtype: list
        """
        return self._entries

    @property
    def seqgroups(self):
        """
        list of L{ClansSeqgroup}s defined in the clustermap.

        @rtype: list
        """
        return self._seqgroups

    def set_default_rotmtx(self):
        """
        Resets the rotation matrix (rotmtx) to no rotation.
        """
        self.rotmtx = eye(3)

    def _update_index(self):
        """
        Creates an index of L{ClansEntry}s to their position in the L{Clans}
        instance.

        The index is used to allow for fast access via L{ClansEntry.get_id} and
        was introduced to get a better L{Clans}.write() performance, which
        suffered from excessive entry.get_id() calls during HSP block generation
        (see L{ClansFileWriter.add_hsp_block}).

        @raises DuplicateEntryError: if two entries have the same name, sequence, and coordinates
        """
        unique_ids = [e._get_unique_id() for e in self]
        
        if len(unique_ids) != len(set(unique_ids)):
            for i, entry in enumerate(unique_ids):
                if unique_ids.count(entry) != 1:
                    raise DuplicateEntryError(
                        '{0} is contained multiple times in this Clans instance'.format(
                            repr(self.entries[i])))

        self._idx = dict(zip(unique_ids, range(len(self))))
        self._has_good_index = True

    def initialize(self):
        """Initializes the coordinates of all entries with random numbers in [-1, 1]."""
        [entry.initialize_coordinates() for entry in self]
        
    def sort(self, key=None):
        """
        Sorts the L{ClansEntry}s by their name or by a custom comparison function.

        @param key: a custom key function
        @type key: function
        """
        self._entries._sort(key)

        self._has_good_index = False

    def add_group(self, group, members=None):
        """
        Adds a new group.

        @param group: the new group
        @type group: L{ClansSeqgroup} instance

        @param members: L{ClansEntry} instances to be in the new group
        @type members: list

        @raise ValueError: if group is no ClansSeqgroup instance
        """
        self.seqgroups._append_item(group)

        if members is not None:
            [group.add(member) for member in members]

    def remove_group(self, group):
        """
        Removes a group.

        @param group: the new group
        @type group: L{ClansSeqgroup} instance
        """
        self.seqgroups._remove_item(group)

        [group.remove(member) for member in group.members]

    def append_groups_from(self, other):
        '''
        Append the L{ClansSeqgroup}-s of C{other} that contain at least one entry
        that corresponds to an entry in this instance. Entries are compared by
        their name only! Groups without any matching members in this instance are
        not created in the local instance.

        @param other: the source of the new group definitions
        @type other: L{Clans} instance
        '''

        for group in other.seqgroups:
            new_group = ClansSeqgroup(name=group.name,
                                      type=group.type,
                                      size=group.size,
                                      hide=group.hide,
                                      color=group.color)

            for member in group.members:
                try:
                    new_member = self.get_entry(member.name, pedantic=True)

                except ValueError:  # no entry with this name found
                    continue

                except DuplicateEntryNameError:
                    raise DuplicateEntryNameError('multiple entries with identical name: {0}'.format(member.name))

                new_group.add(new_member)

            if len(new_group.members) > 0:
                self.add_group(new_group)
 
    def add_entry(self, entry):
        """
        Adds an new entry.

        @param entry: the new entry
        @type entry: L{ClansEntry} instance

        @raise ValueError: if C{entry} is no L{ClansEntry} instance
        """
        if not isinstance(entry, ClansEntry):
            raise ValueError('entries need to be L{ClansEntry} instances')

        self.entries._append_item(entry)
        entry._parent = self

        self._has_good_index = False

    def remove_entry_by_name(self, entry_name):
        """
        Removes an entry fetched by its name.

        @param entry_name: name of the entry that shall be removed
        @type entry_name: string
        """
        entry = self.get_entry(entry_name, True)

        self.remove_entry(entry)

    def remove_entry(self, entry):
        """
        Removes an entry.

        @param entry: the entry that shall be removed
        @type entry: L{ClansEntry} instance
        """
        for other_entry in entry.hsp.keys():
            other_entry.remove_hsp(entry)

        groups = [g for g in entry.groups]
        for g in groups:
            g.remove(entry)

        remove_groups = [g for g in self.seqgroups if g.is_empty()]
        [self.seqgroups._remove_item(g) for g in remove_groups]

        self.entries._remove_item(entry)
        self._has_good_index = False

    def get_entry(self, name, pedantic=True):
        """
        Checks if an entry with name C{name} exists and returns it.
        
        @param name: name of the sought entry
        @type name: str

        @param pedantic: If True, a ValueError is raised if multiple entries
                         with name name are found. If False, returns the first
                         one.
        @type pedantic: bool

        @raise ValueError: if no entry with name C{name} is found
        @raise DuplicateEntryNameError: if multiple entries with name C{name}
        are found and C{pedantic == True}

        @rtype: L{ClansEntry}
        @return: entry with name C{name}
        """

        hits = [e for e in self.entries if e.name == name]

        if len(hits) == 1:
            return hits[0]

        elif len(hits) > 1:
            if pedantic:
                raise DuplicateEntryNameError(
                    'multiple entries have name \'{0}\''.format(name))
            return hits[0]

        else:
            raise ValueError('ClansEntry {0} does not exist.'.format(name))


    def restrict_to_max_pvalue(self, cutoff):
        """
        removes all L{ClansEntry}s that have no connections above the C{cutoff}

        @param cutoff: the cutoff
        @type cutoff: float
        """
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

            [self.remove_entry(e) for e in remove_us if e in self]

        return removed_entries

    def restrict(self, keep_names):
        """
        Removes all entries whose name is not in keep_names

        @param keep_names: names of entries that shall be kept
        @type keep_names: iterable
        """

        [self.remove_entry(entry) for entry in
         [e for e in self.entries if e.name not in keep_names]]

    def write(self, filename):
        """
        writes the L{Clans} instance to a file in CLANS format

        @param filename: the target file\'s name
        @type filename: str
        """

        with open(filename, 'w') as stream:
            writer = ClansFileWriter(stream)
            writer.serialize(self)


class ClansEntry(object):
    """
    Class holding the data of one CLANS sequence entry.

    @param name: the entry name
    @type name: str

    @param seq: the entry\'s amino acid sequence
    @type seq: str

    @param coords: coordinates in 3D space
    @type coords: iterable with 3 items

    @param parent: parent of this entry
    @type parent: L{Clans} instance
    """

    def __init__(self, name=None, seq='', coords=None, parent=None):
        self._name = name
        self._seq = seq

        if coords is None:
            self.initialize_coordinates()
        else:
            self._coords = coords

        self._parent = parent

        self._groups = []
        self._hsp = {}

    def __repr__(self):
        if self.coords is None:
            coords_string = 'NoCoordsSet'
        else:
            coords_string = '({0:.2f}, {1:.2f}, {2:.2f})'.format(
                *tuple(self.coords))

        groups = 'not in a group'
        if len(self.groups) > 0:
            groups = 'groups: {0}'.format(
                ', '.join([g.name for g in self.groups]))

        return 'ClansEntry "{0}": {1} '.format(
            self.name, '; '.join((coords_string, groups)))

    @property
    def name(self):
        """
        name of the entry

        raises ValueError if C{value} is not a string

        @rtype: string
        """
        return self._name

    @name.setter
    def name(self, value):
        if not isinstance(value, csb.core.string):
            raise ValueError(('name cannot be {0} (accepted values: '
                              + 'strings)').format(value))

        self._name = value

    @property
    def seq(self):
        """
        protein sequence of the entry

        raises ValueError if C{value} is not a string

        @rtype: string
        """
        return self._seq

    @seq.setter
    def seq(self, value):
        if not isinstance(value, csb.core.string):
            raise ValueError(('seq cannot be {0} (accepted values: '
                              + 'strings)').format(value))

        self._seq = value

    @property
    def coords(self):
        """
        entry coordinates in 3D space

        raises ValueError if C{value} is not an iterable with 3 items

        @rtype: string
        """
        return self._coords

    @coords.setter
    def coords(self, value):
        if len(value) != 3:
            raise ValueError(('coords cannot be {0} (accepted values: '
                              + 'iteratables with 3 items)').format(value))

        self._coords = value

    @property
    def parent(self):
        """
        L{Clans} instance that parents this L{ClansEntry}

        @rtype: L{Clans}
        """
        return self._parent

    @property
    def groups(self):
        """
        L{ClansSeqgroup}s that contain the entry

        @rtype: list
        """
        return self._groups

    @property
    def hsp(self):
        """
        connections between this and another L{ClansEntry}

        @rtype: dict
        """
        return self._hsp

    def get_id(self):
        """
        Returns the id of the current entry.

        Note: the first call to this method triggers L{Clans._update_index},
        which will make it appear slower than successive calls.

        @rtype: str
        @return: the entrys\' id is returned unless it has no parent in which
        case -1 is returned
        """

        if self.parent is None:
            return -1

        if not self.parent._has_good_index:
            self.parent._update_index()

        return self.parent._idx[self._get_unique_id()]

    def _get_unique_id(self):
        """
        Returns a >>more or less<< unique ID (however this is not guaranteed to be
        really unique) consisting of the name, sequence, and coordinates of the entry.
        If two entries have the same 'unique' id,L{Clans._update_index} will raise a
        DuplicateEntryError.

        @rtype: str
        @return: a more or less unique id
        """
        return '{0.name}<###>{0.seq}<###>{0.coords}'.format(self)

    def initialize_coordinates(self):
        """Sets the coordinates to random numbers in [-1, 1]"""
        self.coords = random.random(3) * 2 - 1
        
    def add_hsp(self, other, value):
        """
        Creates an HSP from self to other with the given value.

        @param other: the other entry
        @type other: L{ClansEntry} instance

        @param value: the value of the HSP
        @type value: float
        """
        self.hsp[other] = value
        other.hsp[self] = value

    def remove_hsp(self, other):
        """
        Removes the HSP between C{self} and C{other}; if none exists, does
        nothing.

        @param other: the other entry
        @type other: L{ClansEntry} instance
        """
        if other in self.hsp:
            self.hsp.pop(other)

        if self in other.hsp:
            other.hsp.pop(self)

    def output_string_seq(self):
        """
        Creates the CLANS <seq> block format representation of the entry.

        @rtype: str
        @return: entrys\' representation in CLANS <seq> block format
        """

        return '>{0}\n{1}\n'.format(self.name, self.seq)

    def output_string_pos(self):
        """
        Create the CLANS <pos> block format representation of the entry.

        @rtype: str
        @return: entrys\' representation in CLANS <pos> block format
        """
        return '{0} {1:.8f} {2:.8f} {3:.8f}'.format(
            *tuple([self.get_id()] + list(self.coords)))

    def output_string_hsp(self):
        """
        Creates the CLANS <hsp> block format representation of the entry.


        @rtype: str
        @return: entrys\' representation in CLANS <hsp> block format
        """
        return '\n'.join(['{0} {1}:{2:.8f}'.format(self.get_id(),
                                                   other.get_id(), value)
                          for (other, value) in self.hsp.items()])


class ClansSeqgroup(object):
    """
    Class holding the data of one CLANS group (seqgroup).

    @kwparam name: name of the seqgroup
    @type name: string

    @kwparam type: symbol used to represent the seqgroup in the graphical
    output
    @type type: int

    @kwparam size: size of the symbol used to represent the seqgroup in the
    graphical output
    @type name: int

    @kwparam hide: if True, the seqgroup\'s symbols in the graphical output are
    not drawn; default: False
    @type name: bool

    @kwparam color: color of the seqgroup
    @type color: L{Color} or string formatted like \'x;y;z\'

    @kwparam members: list of members of this seqgroup
    @type members: list
    """

    def __init__(self, **kw):
        self._name = None
        self.name = kw.pop('name', 'NO NAME')
        
        self._type = None
        self.type = kw.pop('type', 0)

        self._size = None
        self.size = kw.pop('size', 4)

        self._hide = None
        self.hide = kw.pop('hide', False)

        self._color = None
        self.color = kw.pop('color', (255, 255, 255))

        self._members = []
        if 'members' in kw:
            for member in kw['members']:
                self.add(member)

    def __repr__(self):
        return ('ClansSeqgroup {0.name}: type: {0.type}; size: {0.size}; hide:'
                + ' {0.hide}; color: {1}; #members: {2}').format(
            self, self.color.to_clans_color(), len(self.members))

    def __len__(self):
        return len(self.members)

    @property
    def name(self):
        """
        name of the seqgroup

        raises ValueError if C{value} is no string

        @rtype: string
        """
        return self._name

    @name.setter
    def name(self, value):
        if not isinstance(value, csb.core.string):
            raise ValueError('name must be a string')
        self._name = value
    
    @property
    def type(self):
        """
        symbol used to represent the seqgroup in the graphical output

        raises ValueError if C{value} is not castable to int

        @rtype: int
        """
        return self._type

    @type.setter
    def type(self, value):
        self._type = int(value)
      
    @property
    def size(self):
        """
        size of the symbol used to represent the seqgroup in the graphical
        output

        raises ValueError if C{value} is not castable to int

        @rtype: int
        """
        return self._size

    @size.setter
    def size(self, value):
        self._size = int(value)

    @property
    def hide(self):
        """
        if True, the seqgroup\'s symbols in the graphical output are not drawn

        raises ValueError if C{value} is no bool

        @rtype: int
        """
        return self._hide

    @hide.setter
    def hide(self, value):
        if not isinstance(value, bool):
            raise ValueError(('hide cannot be {0} (accepted values: '
                              + 'True/False)').format(value))

        self._hide = value
        
    @property
    def color(self):
        """
        color of the seqgroup

        raises ValueError if set to a wrongly formatted string (correct:
        \'{r};{g};{b}\')

        @rtype: L{Color}
        """
        return self._color

    @color.setter
    def color(self, value, separator=';'):
        # set values to those of existing Color instance
        if isinstance(value, Color):
            self._color = value
            return

        ## parse color from string in format 'r;g;b'
        if isinstance(value, csb.core.string):
            self._color = Color.from_string(value)
            return

        # parse 3-item iterables like (3, 5, 6)
        if len(value) == 3:
            self._color = Color(*tuple(map(int, value)))
            return

        raise ValueError('cannot parse color from \'{0}\''.format(value))


    @property
    def members(self):
        """
        the members of this seqgroup

        @rtype: list
        """
        return self._members

    def is_empty(self):
        """
        Checks if the group contains entries.

        @rtype: bool
        @return: True if the group contains no entries, else False.
        """
        return len(self) == 0

    def add(self, new_member):
        """
        Adds entry C{new_member} to this L{ClansSeqgroup}.

        @Note: L{ClansEntry}-s added using this method that are not part of the
        main L{Clans} instance need to be added to the L{Clans} instance manually.
        
        @param new_member: the member that shall be added to this
        L{ClansSeqgroup}
        @type new_member: L{ClansEntry} instance

        @raise TypeError: if C{new_member} is no L{ClansEntry} instance
        @raise ValueError: if C{new_member} is already contained in this
        L{ClansSeqgroup}
        """
        if not isinstance(new_member, ClansEntry):
            raise TypeError('only ClansEntry instances can be added as ' +
                            'group members')

        if self.members.count(new_member) > 0:
            raise ValueError(('entry {0.name} is already contained in this '
                              + 'seqgroup').format(new_member))

        self.members.append(new_member)
        new_member.groups.append(self)

    def remove(self, member):
        """
        Removes L{ClansEntry} C{member} from this group.
        
        @param member: the member to be removed
        @type member: a L{ClansEntry} instance

        @raise TypeError: if C{member} is no L{ClansEntry} instance
        @raise ValueError: if C{member} is not part of this L{ClansSeqgroup}
        """
        if not isinstance(member, ClansEntry):
            raise TypeError('argument must be a ClansEntry instance')

        if self.members.count(member) == 0:
            raise ValueError(('"{0.name}" is not a member of this '
                            + 'seqgroup').format(member))

        self.members.remove(member)
        member.groups.remove(self)

    def output_string(self):
        """
        Creates the CLANS <seqgroup> block format representation of the
        group.

        @rtype: str
        @return: entrys\' representation in CLANS <seqgroup> block format
        """
        sorted_members = sorted([m.get_id() for m in self.members])
        return ('name={0.name}\ntype={0.type}\nsize={0.size}\nhide={1}'
                + '\ncolor={2}\nnumbers={3}').format(
            self, int(self.hide), self.color.to_clans_color(),
            ';'.join([str(val) for val in sorted_members]) + ';')
