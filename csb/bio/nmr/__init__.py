"""
NMR related objects.
"""

import os
import numpy.linalg
import xml.dom.minidom

import csb.io.tsv
import csb.core as pu

from csb.statistics.pdf import GeneralizedNormal
from csb.bio.sequence import ProteinAlphabet
from csb.bio.structure import ChemElements


class InvalidResidueError(ValueError):
    pass

class EntityNotSupportedError(KeyError):
    pass


class RandomCoil(object):
    """
    Utility class containing all necessary data and methods for computing
    secondary chemical shifts.
    
    @note: You are supposed to obtain an instance of this object only via
           the dedicated factory (see L{RandomCoil.get}). The factory
           ensures a "singleton with lazy instantiation" behavior. This is
           needed since this object loads static data from the file system.
    """
    
    RESOURCES = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'resources')
    
    _instance = None    

    
    @staticmethod
    def get():
        """
        Get the current L{RandomCoil} instance (and create it, if this
        method is called for the first time).
        """
        if RandomCoil._instance is None:
            RandomCoil._instance = RandomCoil()
        
        return RandomCoil._instance  
        
    def __init__(self):
        
        if RandomCoil._instance is not None:
            raise NotImplementedError("Can't instantiate a singleton")

        RandomCoil._instance = self

        self._offsets = (-2, -1, 1, 2)
        self._reference = {}
        self._corrections = {}
        
        self._initialize()
    
    def _initialize(self):
     
        ref = os.path.join(RandomCoil.RESOURCES, 'RandomCoil.Reference.tsv')            
        cor = os.path.join(RandomCoil.RESOURCES, 'RandomCoil.Corrections.tsv')

        self._load(ref, cor)            

    def _load(self, ref, cor):

        self._reference = {}
        self._corrections = {}
                                
        header = 'Residue:str Nucleus:str Value:float'
        
        for row in csb.io.tsv.Table.from_tsv(ref, header):
            residue = pu.Enum.parsename(ProteinAlphabet, row[0])
            nucleus, value = row[1:]
            
            if residue not in self._reference:
                self._reference[residue] = {}
            
            self._reference[residue][nucleus] = value
        
        header = 'Residue:str Nucleus:str CS1:float CS2:float CS3:float CS4:float'
        
        for row in csb.io.tsv.Table.from_tsv(cor, header):   
            residue = pu.Enum.parsename(ProteinAlphabet, row[0])
            nucleus = row[1]
            values = row[2:]
            
            if residue not in self._corrections:
                self._corrections[residue] = {}
            
            self._corrections[residue][nucleus] = dict(zip(self._offsets, values))
    
    def simple_secondary_shift(self, residue, nucleus, value):
        """
        Compute a secondary shift given a raw shift C{value}.
        Residue neighborhood is not taken into account.
        
        @param residue: residue type (amino acid code)
        @type residue: str or L{EnumItem}
        @param nucleus: atom name (PDB format)
        @type nucleus: str
        @param value: raw chemical shift value
        @type value: float
        
        @return: float
        
        @raise EntityNotSupportedError: on unsupported residue or nucleus 
        """
                   
        try:
            if isinstance(residue, pu.string):
                if len(residue) == 1:
                    residue = pu.Enum.parse(ProteinAlphabet, residue)
                else:
                    residue = pu.Enum.parsename(ProteinAlphabet, residue)
            else:
                if residue.enum is not ProteinAlphabet:
                    raise TypeError(residue)
                                              
            return value - self._reference[residue][nucleus]
        
        except (pu.EnumValueError, pu.EnumMemberError):
            raise InvalidResidueError('{0} is not a protein residue'.format(residue))
        
        except KeyError as ke:
            raise EntityNotSupportedError('{0!s}, context: {1!r} {2}'.format(ke, residue, nucleus))

    def secondary_shift(self, chain, residue, nucleus, value):
        """
        Compute a secondary shift given a raw shift C{value} for a specific
        residue and its neighboring residues.
        
        @param chain: the protein chain containing the C{nucleus}
        @type chain: L{Chain}
        @param residue: the residue containing the C{nucleus}. This can be
                        a residue object, id (sequence number + insertion
                        code, string) or rank (integer, 1-based)
        @type residue: L{Residue}, str or int
        @param nucleus: atom name (PDB format)
        @type nucleus: str
        @param value: raw chemical shift value
        @type value: float            
        """
        try:
            if isinstance(residue, int):
                residue = chain.residues[residue]
            elif isinstance(residue, pu.string):
                residue = chain.find(residue)
            else:
                residue = chain.residues[residue.rank]
        except (pu.ItemNotFoundError, pu.CollectionIndexError):
            raise InvalidResidueError("Can't find residue {0} in {1}".format(residue, chain))
            
        shift = self.simple_secondary_shift(residue.type, nucleus, value)
        
        for offset in self._offsets:
            
            if 1 <= (residue.rank + offset) <= chain.length:
                try:
                    neighbor = chain.residues[residue.rank + offset]
                    shift -= self._corrections[neighbor.type][nucleus][offset * -1]
                                     
                except KeyError:
                    continue     
        
        return shift


class AtomConnectivity(object):
    
    RESOURCES = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'resources')
    
    _instance = None
    
    @staticmethod
    def get():
        """
        Get the current L{AtomConnectivity} instance (and create it if this
        method is invoked for the first time).
        @rtype: L{AtomConnectivity}
        """
        if AtomConnectivity._instance is None:
            AtomConnectivity._instance = AtomConnectivity()
        return AtomConnectivity._instance
    
    def __init__(self):
        
        self._table = {}
        self._initialize()
        
    def _initialize(self):
        
        resource = os.path.join(AtomConnectivity.RESOURCES, 'AtomConnectivity.xml')
        root = xml.dom.minidom.parse(resource)
        
        for r in root.documentElement.getElementsByTagName('residue'):
            residue = pu.Enum.parsename(ProteinAlphabet, r.getAttribute('type'))
            self._table[residue] = {}
            
            for a in r.getElementsByTagName('atom'):
                atom = a.getAttribute('name')
                self._table[residue][atom] = set()
            
            for b in r.getElementsByTagName('bond'):
                atom1 = b.getAttribute('atom1')
                atom2 = b.getAttribute('atom2')
                self._table[residue][atom1].add(atom2)
                self._table[residue][atom2].add(atom1)
                
    def connected(self, residue, atom1, atom2):
        """
        Return True if C{atom1} is covalently connected to C{atom2} in C{residue}
        
        @param residue: residue type (a member of L{ProteinAlphabet})
        @type residue: L{EnumItem} 
        @param atom1: first atom name (IUPAC)
        @type atom1: str
        @param atom2: second atom name (IUPAC)
        @type atom2: str
        
        @rtype: boolean
        """
        if residue in self._table:
            r = self._table[residue]
            if atom1 in r:
                return atom2 in r[atom1]
        
        return False
    
    def connected_atoms(self, residue, atom):
        """
        Return all atoms covalently connected to C{atom} in C{residue}.

        @param residue: residue type (a member of L{ProteinAlphabet})
        @type residue: L{EnumItem}         
        @param atom: source atom name (IUPAC)
        @type atom: str
        
        @rtype: tuple of str
        """
        if residue in self._table:
            r = self._table[residue]
            if atom in r:
                return tuple(r[atom])
        
        return tuple()
    
    def contains(self, residue, atom):
        """
        Return True if C{atom} name is contained in C{residue}.
        
        @param residue: residue type (a member of L{ProteinAlphabet})
        @type residue: L{EnumItem}  
        @param atom: atom name (IUPAC)
        @type atom: str
        
        @rtype: bool        
        """
        if residue in self._table:
            return atom in self._table[residue]
        
        return False        
    
    def get_atoms(self, residue, prefix=''):
        """
        Get all atoms contained in C{residue}.

        @param residue: residue type (a member of L{ProteinAlphabet})
        @type residue: L{EnumItem}         
        @param prefix: atom name prefix wildcard (IUPAC)
        @type prefix: str
        
        @return: set of atom names
        @rtype: frozenset of str
        """
        t = self._table[residue]
        if residue in self._table:
            return frozenset(a for a in t if a.startswith(prefix))
        
        return frozenset()
    
    
class Filters(object):
    """
    Pre-built atom filters for L{ContactMap}s. 
    """

    @staticmethod
    def ALL(a):
        return True
        
    @staticmethod
    def HYDROGENS(a):
        return a.element == ChemElements.H
    
    @staticmethod
    def CARBONS(a):
        return a.element == ChemElements.C
    
    @staticmethod
    def CALPHAS(a):
        return a.name == 'CA'

class ContactMap(object):
    """
    Describes a protein contact map. Atoms positioned at distance below
    a given cutoff are considered to be in contact.
    
    @param chain: source protein chain
    @type chain: L{csb.bio.structure.Chain} 
    @param cutoff: distance cutoff in angstroms
    @type cutoff: float
    @param filter: a callable with signature 'bool def(csb.bio.structure.Atom)',
                   invoked for every atom, which determines whether a given atom 
                   should be skipped (False) or considered (True). See L{Filters}
    @type filter: lambda 
    """
    
    DISTANCE_CUTOFF = 6.0 
    
    @staticmethod
    def load(filename):
        """
        Deserialize from a pickle.
        """
        with open(filename, 'rb') as stream:
            return csb.io.Pickle.load(stream)
    
    def __init__(self, chain, cutoff=DISTANCE_CUTOFF, filter=None):
        
        self._cutoff = float(cutoff)
        self._chain = chain
        self._atoms = []
        self._atomset = set()
        self._map = {}
        self._coords = {}
        
        if filter is None:
            filter = lambda i: True
        
        for residue in chain.residues:
            self._coords[residue.rank] = {}
            atoms = [a for a in residue.items if filter(a)]
            
            if len(atoms) == 0:
                continue
            
            step = 1.0 / len(atoms)
            n = 0
            
            for atom in atoms:
                self._atoms.append(atom)
                self._atomset.add(atom)
                self._coords[residue.rank][atom.name] = residue.rank + n * step
                n += 1        
                        
    def __iter__(self):
        return self.contacts
        
    def __contains__(self, atom):
        return atom in self._atomset
        
    @property
    def cutoff(self):
        """
        Distance cutoff in Angstroms
        @rtype: float
        """
        return self._cutoff
    
    @property
    def chain(self):
        """
        Source protein chain
        @rtype: L{Chain}        
        """        
        return self._chain
    
    @property
    def atoms(self):
        """
        All atoms involved in this map, sorted by residue number
        @rtype: tuple of L{Atom}
        """
        return tuple(self._atoms)
    
    @property
    def contacts(self):
        """
        All atom contacts: an iterator over all contacting 
        (L{Atom}, L{Atom}) pairs.
        @rtype: iterator of 2-tuples   
        """
        visited = set()
        
        for a1 in self._map:
            for a2 in self._map[a1]:
                if (a1, a2) not in visited:
                    visited.add((a1, a2))
                    visited.add((a2, a1))
                    yield (a1, a2)        
    
    def build(self):
        """
        Extract all contacts from the chain using the current distance cutoff.
        """
        
        self._map = {}
        
        for atom1 in self._atoms:
            for atom2 in self._atoms:
                if atom1 is not atom2:
                    distance = numpy.linalg.norm(atom1.vector - atom2.vector)
                    if distance <= self._cutoff:
                        self._connect(atom1, atom2)
                        
    def connect(self, atom1, atom2):
        """
        Define a contact between C{atom1} and C{atom2}.
        
        @param atom1: first atom
        @type atom1: L{Atom}
        @param atom2: second atom
        @type atom2: L{Atom}        
        """
        for atom in [atom1, atom2]:
            if atom not in self._atomset:
                raise ValueError("No such atom in contact map: {0}".format(atom))        
        
        self._connect(atom1, atom2)        
                        
    def _connect(self, atom1, atom2):
        
        if atom1 not in self._map:
            self._map[atom1] = set()
        self._map[atom1].add(atom2)
        
        if atom2 not in self._map:
            self._map[atom2] = set()
        self._map[atom2].add(atom1)
        
    def connected(self, atom1, atom2):
        """
        Return True if the specified atoms are in contact.
        
        @param atom1: first atom
        @type atom1: L{Atom}
        @param atom2: second atom
        @type atom2: L{Atom}   
        """
        if atom1 in self._map:
            return atom2 in self._map[atom1]
        return False    
        
    def atom_contacts(self, atom):
        """
        Return all atoms within C{self.cutoff} angstroms of C{atom}.
        
        @param atom: anchor atom
        @type atom: L{Atom}
        
        @rtype: frozenset of L{Atom}
        """
        
        if atom in self._map:
            return frozenset(self._map[atom])
        else:
            return frozenset()
        
    def residue_contacts(self, residue):
        """
        Return all residues, having neighboring atoms within C{self.cutoff}
        angstroms from any of the C{residue}'s atoms.
        
        @param residue: anchor residue
        @type residue: L{Residue}
        
        @rtype: frozenset of L{Residue}
        """        
        
        partners = set()
        
        for atom in residue.items:
            if atom in self._map:
                for partner in self._map[atom]:
                    partners.add(partner.residue)
                    
        return frozenset(partners)
    
    def position(self, rank, atom_name):
        """
        Compute the location of C{atom} on the contact map.
        
        @param rank: residue rank (1-based)
        @type rank: int 
        @param atom_name: atom name
        @type atom_name: str
        
        @rtype: float
        """
        residue = self._chain.residues[rank]
        atom = residue.atoms[atom_name]
                
        try:
            return self._coords[residue.rank][atom.name]
        except KeyError:
            msg = "No atom {0} at #{1} in contact map: {2}"
            raise ValueError(msg.format(atom_name, rank, self._coords[residue.rank].values()))
        
    def atom_matrix(self):
        """
        Build a 2D binary contact matrix (0=no contact, 1=contact). The order of elements
        in each dimension will match the order of atoms in the contact map
        (see L{ContactMap.atoms} and iter(L{ContactMap}). That means, the atoms in
        each dimension are sorted by residue number first.
        
        @deprecated: This method can be removed in future versions
        
        @rtype: numpy.array (2D) 
        """
        
        matrix = []
            
        for i, atom1 in enumerate(self.atoms):
            matrix.append([])
            
            for atom2 in self.atoms:
                if atom1 in self._map and atom2 in self._map[atom1]:
                    matrix[i].append(1)
                else:
                    matrix[i].append(0)
                    
        return numpy.array(matrix)
    
    def draw(self, plot, color="black"):
        """
        Visualize this contact map.
        
        @param plot: L{csb.io.plots.Chart}'s plot to draw on
        @type plot: matplotlib.AxesSubplot
        @param color: pixel color (must be a matplotlib color constant)
        @type color: str
        """
        
        x, y = [], []
        
        for atom1 in self.atoms:
            for atom2 in self.atom_contacts(atom1):
                pos1 = self.position(atom1.residue.rank, atom1.name)
                pos2 = self.position(atom2.residue.rank, atom2.name)
            
                assert None not in (pos1, pos2), (atom1, atom2)
                x.append(pos1)
                y.append(pos2)
                    
        plot.plot(x, y, color=color, marker=",", linestyle='none')
                
        plot.set_xlim(0, self.chain.length)
        plot.set_ylim(0, self.chain.length)
        
        return plot      
    
    @staticmethod
    def compare(query, reference, min_distance=0):
        """
        Compare a query contact map against a reference.
        
        @type query: L{ContactMap}
        @type reference: L{ContactMap}
        
        @param min_distance: consider only contacts between atoms, separated by
                             the given minimum number of residues
        @type min_distance: int
        
        @return: precision and coverage
        @rtype: L{ContactMapComparsonInfo}  
        """
        if query.chain is not reference.chain:
            raise ValueError("Contact maps are not comparable")
        if not query._map and not reference._map:
            raise ValueError("Can't compare empty contact maps")
        
        true_pos = 0.0
        false_pos = 0.0
        false_neg = 0.0
        
        for a1, a2 in query.contacts:
            if abs(a1.residue.rank - a2.residue.rank) >= min_distance: 
                if reference.connected(a1, a2):
                    true_pos += 1.0
                else:
                    false_pos += 1.0

        for a1, a2 in reference.contacts:
            if abs(a1.residue.rank - a2.residue.rank) >= min_distance:            
                if not query.connected(a1, a2):
                    false_neg += 1.0
        
        try:
            precision = true_pos / (true_pos + false_pos)
            coverage = true_pos / (true_pos + false_neg)
            return ContactMapComparsonInfo(precision, coverage)
        
        except ZeroDivisionError:
            return ContactMapComparsonInfo(0, 0)
        
class ContactMapComparsonInfo(object):
    
    def __init__(self, precision, coverage):
        
        self.precision = precision
        self.coverage = coverage
        
        
class Label(object):
    """
    Utility class for working with chemical shift labels.
    
    @param residue: residue type
    @type residue: L{EnumItem}
    @param rank: residue position (1-based)
    @type rank: int
    @param atom_name: nucleus name
    @type atom_name: str 
    """
    
    @staticmethod
    def build(residue_type, position, atom_name):
        """
        Build a new string label by specifying its components.
        @rtype: str        
        """
        return '{0!s}#{1}:{2}'.format(residue_type, position, atom_name)

    @staticmethod    
    def from_shift(shift):
        """
        Build a new string label from a L{ChemShiftInfo}.
        @rtype: str        
        """
        return Label.build(shift.residue, shift.position, shift.name)

    @staticmethod    
    def from_atom(atom):
        """
        Build a new string label from an L{Atom}.
        @rtype: str        
        """        
        return Label.build(atom.residue.type, atom.residue.rank, atom.name)
    
    @staticmethod
    def match(shift, atom):
        """
        Return True if the labels of a L{ChemShiftInfo} and an L{Atom} match.
        @rtype: bool        
        """          
        
        l = Label.from_shift(shift)
        r = Label.from_atom(atom)
        
        return r == l
    
    @staticmethod
    def get_atom(chain, label):
        """
        Get the L{Atom} in a L{Chain}, designated by a given string label.
        @rtype: L{Atom}
        """
        dummy, rank, atom = Label.parse(label)
        return chain.residues[rank].atoms[atom]

    @staticmethod    
    def parse(label):
        """
        Parse the components of a string nucleus label.
        @return: (residue, rank, atom)
        @rtype: 3-tuple
        """        
        parts = label.split("#")
        residue = parts[0]
        
        subparts = parts[1].split(":")
        rank = int(subparts[0])
        atom = subparts[1]
        
        return (residue, rank, atom)
    
    @staticmethod
    def from_string(label):
        """
        Parse the a string nucleus label and create a new L{Label}.
        @rtype: L{Label}
        """           
        residue, rank, atom = Label.parse(label)
        return Label(residue, rank, atom)
    
    def __init__(self, residue, rank, atom_name):
        
        self._residue = residue
        self._rank = rank
        self._atom = atom_name
        
    @property
    def residue(self):
        """
        Residue type (a L{ProteinAlphabet} member)
        """
        return self._residue
    
    @property
    def rank(self):
        """
        Residue rank (1-based)
        """
        return self._rank
    
    @property
    def atom_name(self):
        """
        Nucleus name
        """        
        return self._atom
    
    def __str__(self):
        return Label.build(self._residue, self._rank, self._atom)
    

class ChemShiftInfo(object):
    """
    Chemical shift struct.
    
    @param position: residue rank (1-based)
    @type position: int
    @param residue: amino acid type (a member of L{ProteinAlphabet})
    @type residue: str or L{EnumItem}
    @param name: nucleus label
    @type name: str
    @param element: nucleus type (a member of L{ChemElements})
    @type element: str or L{EnumItem}
    @param shift: chemical shift value
    @type shift: float
    """
    
    def __init__(self, position, residue, name, element, shift):
        
        if not isinstance(residue, pu.EnumItem) or residue.enum is not ProteinAlphabet:
            residue = pu.Enum.parsename(ProteinAlphabet, str(residue))
            
        if not isinstance(element, pu.EnumItem) or element.enum is not ChemElements:
            element = pu.Enum.parsename(ChemElements, str(element))            
        
        self.position = int(position)
        self.residue = residue
        self.name = str(name)
        self.element = element
        self.shift = float(shift)
        
    def clone(self, name):
        """
        Clone the current shift and create a new one with the specified
        nucleus label.
        
        @rtype: L{ChemShiftInfo}
        """
        ni = self
        return ChemShiftInfo(ni.position, repr(ni.residue), name, repr(ni.element), ni.shift)
        
    def __str__(self):
        return "{0!s}#{1}:{2}".format(self.residue, self.position, self.name)
    
    @property
    def label(self):
        """
        String label representation
        @rtype: str
        """
        return str(self)

class ChemicalShiftNetwork(object):
    """
    Describes a network of covalently connected, chemical shift visible nuclei.
    
    @param shifts: chemical shift instances
    @type shifts: iterable of L{ChemShiftInfo}
    """
    
    def __init__(self, shifts):

        self._neighbors = {}
 
        labels = {}
        
        for cs in shifts:
            self._neighbors[cs] = set()
            id = Label.from_shift(cs)
            labels[id] = cs
        
        conn = AtomConnectivity.get()
        
        for cs in shifts:
            for atom_name in conn.connected_atoms(cs.residue, cs.name):
                target = Label.build(cs.residue, cs.position, atom_name)
                if target in labels:
                    self.connect(cs, labels[target])
    
    def connect(self, cs1, cs2):
        """
        Connect two nuclei.
        
        @param cs1: first chemical shift instance
        @type cs1: L{ChemShiftInfo}
        @param cs2: second chemical shift instance         
        @type cs2: L{ChemShiftInfo}
        """
        
        try:
            self._neighbors[cs1].add(cs2)
            self._neighbors[cs2].add(cs1)
        except KeyError:
            raise ValueError("Unknown chemical shift")
        
    def connected_shifts(self, source, element=None):
        """
        Return an iterator over all covalently connected neuclei to a given
        C{source}.
        
        @param source: source chemical shift
        @type source: L{ChemShiftInfo}
        
        @rtype: iterator of L{ChemShiftInfo}
        """
        
        
        if source not in self._neighbors:
            raise ValueError("No such chemical shift in this network")

        for cs in self._neighbors[source]:
            if element is None or cs.element == element:
                yield cs
                
    def __iter__(self):
        return iter(self._neighbors)
    
class ChemShiftScoringModel(object):
    """
    Chemical shift similarity scoring model. See C{ScoringModel.NUCLEI} for
    a list of supported chemical shift types. 
    """

    NUCLEI = ('CA', 'CB', 'C', 'N', 'HA')

    
    def __init__(self):
        
        self._pos = {}
        self._neg = {}
        
        self._pos['CA'] = GeneralizedNormal(0.02, 1.32, 1.1)
        self._neg['CA'] = GeneralizedNormal(-0.08, 4.23, 2.2)
                
        self._pos['CB'] = GeneralizedNormal(0.06, 1.32, 1.0)
        self._neg['CB'] = GeneralizedNormal(0.08, 2.41, 1.2)
                
        self._pos['C']  = GeneralizedNormal(0.12, 1.52, 1.4)
        self._neg['C']  = GeneralizedNormal(-0.13, 3.42, 2.1)
        
        self._pos['N']  = GeneralizedNormal(0.23, 4.39, 1.4)
        self._neg['N']  = GeneralizedNormal(0.17, 7.08, 1.9)
                
        self._pos['HA'] = GeneralizedNormal(0.00, 0.27, 1.0)
        self._neg['HA'] = GeneralizedNormal(-0.01, 0.66, 1.4)
        
        assert set(self._pos) == set(ChemShiftScoringModel.NUCLEI)
        assert set(self._neg) == set(ChemShiftScoringModel.NUCLEI) 

    def positive(self, nucleus, deltas):
        """
        Return the probability that a given chemical shift difference
        indicates structural similarity (true positive match).
        
        @param nucleus: chemical shift (a member of C{ScoringModel.NUCLEI})
        @type nucleus:  str
        @param deltas: chemical shift difference(s): q-s
        @type deltas:  float or list of floats
        
        @return: the raw value of the probability density function
        @rtype: float or array of floats
        """
        results = self._pos[nucleus].evaluate([deltas]) 
        return results[0]

    def negative(self, nucleus, deltas):
        """
        Return the probability that a given chemical shift difference
        indicates no structural similarity (true negative match).
        
        @param nucleus: chemical shift (a member of C{ScoringModel.NUCLEI})
        @type nucleus:  str
        @param deltas: chemical shift difference(s): q-s
        @type deltas:  float or list of floats
        
        @return: the raw value of the probability density function
        @rtype: float or array of floats
        """        
        results = self._neg[nucleus].evaluate([deltas]) 
        return results[0]
            
    def score(self, nucleus, deltas):
        """
        Return the bit score for a given chemical shift difference.
        
        @param nucleus: chemical shift (a member of C{ScoringModel.NUCLEI})
        @type nucleus:  str
        @param deltas: chemical shift difference(s): q-s
        @type deltas:  float or list of floats
        
        @return: bit score
        @rtype: float or array of floats
        """
        pos = self.positive(nucleus, deltas)
        neg = self.negative(nucleus, deltas)
        
        return numpy.log2(pos / neg)

            
class NOEPeak(object):
    """
    Describes a single NOE peak.
    
    @param intensity: peak intensity
    @type intensity: float
    @param dimensions: list of dimension values
    @type dimensions: iterable of float
    @param spectrum: owning NOE spectrum
    @type spectrum: L{NOESpectrum}
    """
    
    def __init__(self, intensity, dimensions, spectrum):
        
        self._dimensions = list(dimensions)
        self._intensity = float(intensity)
        self._spectrum = spectrum
        
    @property
    def intensity(self):
        """
        Peak intensity
        @rtype: float
        """
        return self._intensity
    
    @property
    def num_dimensions(self):
        """
        Number of dimensions
        @rtype: int
        """        
        return len(self._dimensions)
    
    def has_element(self, e):
        """
        Return True if the owning spectrum contains a dimension of the specified type
        
        @param e: element (dimension) type (see L{ChemElements})
        @type e: L{EnumItem}
        
        @rtype: bool
        """        
        return self._spectrum.has_element(e) 
    
    def __getitem__(self, column):
        return self.get(column)
    
    def __iter__(self):
        return iter(self._dimensions)
    
    def __str__(self):
        return '<NOEPeak: {0}, I={1}>'.format(self._dimensions, self._intensity)
    
    def element(self, i):
        """
        Return the dimension (nucleus) type at dimension index i
        
        @param i: dimension index (0-based)
        @type i: int
        
        @return: nucleus type
        @rtype: L{EnumItem}
        """
        return self._spectrum.element(i)
    
    def get(self, column):
        """
        Get the value of the specified dimension. 
        
        @param column: dimension index (0-based)
        @type column: int
        
        @return: dimension value        
        @rtype: float
        """
        if 0 <= column < len(self._dimensions):  
            return self._dimensions[column]
        else:
            raise IndexError("Dimension index out of range")
        
    def has_connected_dimensions(self, i):
        """
        Return True of dimension index C{i} has covalently connected dimensions.
        
        @param i: dimension index (0-based)
        @type i: int
        
        @rtype: bool
        """
        return self._spectrum.has_connected_dimensions(i)
        
    def connected_dimensions(self, i):
        """
        Return a list of all dimension indices, covalently connected to
        dimension C{i}.

        @param i: dimension index (0-based)
        @type i: int
        
        @rtype: iterable of L{EnumItem}         
        """
        return self._spectrum.connected_dimensions(i)


class NOESpectrum(object):
    """
    Describes an NOE spectrum.
    
    @param elements: list of dimension (nucleus) types for each dimension
    @type elements: iterable of L{EnumItem} (L{ChemElements}) or str 
    """
    def __init__(self, elements):
        
        self._elements = []
        self._elemset = set()        
        self._connected = {}
        self._protondim = set() 
        self._peaks = []
        self._min = float("inf")
        self._max = float("-inf")
        
        for i, n in enumerate(elements):
            
            if not isinstance(n, pu.EnumItem) or n.enum is not ChemElements:
                element = pu.Enum.parsename(ChemElements, n)
            else:
                element = n
            self._elements.append(element)
            
            if element == ChemElements.H:
                self._protondim.add(i)
            
        self._elemset = set(self._elements) 
        
    @staticmethod
    def join(spectrum, *spectra):
        """
        Merge multiple L{NOESpectrum} instances. All C{spectra} must have matching
        dimensions according to the master C{spectrum}.
        
        @return: merged spectrum
        @rtype: L{NOESpectrum}
        """
        elements = tuple(spectrum.dimensions)
        joint = NOESpectrum(map(repr, elements))
        
        for i, dummy in enumerate(elements):
            for j in spectrum.connected_dimensions(i):
                joint.connect(i, j)
        
        for s in [spectrum] + list(spectra):
            if tuple(s.dimensions) != elements:
                raise ValueError("Incompatible spectrum: {0}".format(s))
            for p in s: 
                joint.add(p.intensity, list(p))
                
        return joint 
    
        
    def __iter__(self):
        return iter(self._peaks)
    
    def __len__(self):
        return len(self._peaks)
    
    def __str__(self):
        return '<NOESpectrum: {0}>'.format(self._elements)
    
    def __getitem__(self, i):
        try:
            return self._peaks[i]
        except IndexError:
            raise IndexError("Peak index out of range")
    
    @property
    def min_intensity(self):
        """
        Minimum intensity
        @rtype: float
        """
        return self._min

    @property
    def max_intensity(self):
        """
        Maximum intensity
        @rtype: float
        """
        return self._max
            
    @property
    def dimensions(self):
        """
        Tuple of all dimensions (nucleus types)
        @rtype: tuple of L{EnumItem}
        """
        return tuple(self._elements)
    
    @property
    def proton_dimensions(self):
        """
        Tuple of all proton dimension indices
        @rtype: tuple of int
        """
        return tuple(self._protondim)    

    @property    
    def num_dimensions(self):
        """
        Number of dimensions
        @rtype: int
        """        
        return len(self._elements)
    
    @property    
    def num_proton_dimensions(self):
        """
        Number of proton dimensions
        @rtype: int
        """             
        return len(self._protondim)    
    
    def has_element(self, e):
        """
        Return True if the spectrum contains a dimension of the specified type
        
        @param e: element (dimension) type (see L{ChemElements})
        @type e: L{EnumItem}
        
        @rtype: bool
        """          
        return e in self._elemset
    
    def connect(self, i1, i2):
        """
        Mark dimensions with indices C{i1} and C{i2} as covalently connected.
        
        @param i1: dimension index 1 (0-based)
        @type i1: int
        @param i2: dimension index 2 (0-based)
        @type i2: int         
        """

        for i in [i1, i2]:
            if not 0 <= i < self.num_dimensions:
                raise IndexError("Dimension index out of range")
            
        if i1 == i2:
            raise ValueError("Can't connect a dimension to itself")
        if not self._can_connect(i1, i2):
            raise ValueError("Only proton-nonproton bonds are allowed")        
            
        self._connected.setdefault(i1, set()).add(i2)
        self._connected.setdefault(i2, set()).add(i1)
        
    def _can_connect(self, i1, i2):
        
        pair = set()

        for i in [i1, i2]:
            is_proton = self.element(i) == ChemElements.H
            pair.add(is_proton)
            
        if True in pair and False in pair:
            return True
        
        return False        
        
    def has_connected_dimensions(self, i):
        """
        Return True of dimension index C{i} has covalently connected dimensions.
        
        @param i: dimension index (0-based)
        @type i: int
        
        @rtype: bool
        """
        if i in self._connected:
            return len(self._connected[i]) > 0
        
        return False
        
    def connected_dimensions(self, i):
        """
        Return a list of all dimension indices, covalently connected to
        dimension C{i}.

        @param i: dimension index (0-based)
        @type i: int
        
        @rtype: iterable of int        
        """        
        if i in self._connected:
            return tuple(self._connected[i])
        
        return tuple()    
        
    def add(self, intensity, dimensions):
        """
        Add a new NOE peak.
        
        @param intensity: peak intensity
        @type intensity: float
        @param dimensions: list of dimension values
        @param dimensions: iterable of float
        """
        
        dimensions = list(dimensions)       
        if len(dimensions) != self.num_dimensions:
            raise ValueError("Invalid number of dimensions")
        
        peak = NOEPeak(intensity, dimensions, self)
        self._peaks.append(peak)
        
        if peak.intensity < self._min:
            self._min = peak.intensity            
        if peak.intensity > self._max:
            self._max = peak.intensity            
        
    def element(self, i):
        """
        Return the chemical element (nucleus) type at dimension index C{i}.
        @rtype: L{EnumItem}
        """
        return self._elements[i]


