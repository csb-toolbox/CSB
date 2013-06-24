"""
3D and secondary structure APIs.

This module defines some of the most fundamental abstractions in the library:
L{Structure}, L{Chain}, L{Residue} and L{Atom}. Instances of these objects may
exist independently and that is perfectly fine, but usually they are part of a
Composite aggregation. The root node in this Composite is a L{Structure} (or
L{Ensemble}). L{Structure}s are composed of L{Chain}s, and each L{Chain} is a
collection of L{Residue}s. The leaf node is L{Atom}. 

All of these objects implement the base L{AbstractEntity} interface. Therefore,
every node in the Composite can be transformed:
    
    >>> r, t = [rotation matrix], [translation vector]
    >>> entity.transform(r, t)
    
and it knows its immediate children:

    >>> entity.items
    <iterator>    # over all immediate child entities
    
If you want to traverse the complete Composite tree, starting at arbitrary level,
and down to the lowest level, use one of the L{CompositeEntityIterator}s. Or just
call L{AbstractEntity.components}:

    >>> entity.components()
    <iterator>   # over all descendants, of any type, at any level
    >>> entity.components(klass=Residue)
    <iterator>   # over all Residue descendants
    
Some of the inner objects in this hierarchy behave just like dictionaries
(but are not):

    >>> structure.chains['A']       # access chain A by ID
    <Chain A: Protein>
    >>> structure['A']              # the same
    <Chain A: Protein>
    >>> residue.atoms['CS']          
    <Atom: CA>                      # access an atom by its name
    >>> residue.atoms['CS']          
    <Atom: CA>                      # the same
        
Others behave like list collections:

    >>> chain.residues[10]               # 1-based access to the residues in the chain
    <ProteinResidue [10]: PRO 10>
    >>> chain[10]                        # 0-based, list-like access
    <ProteinResidue [11]: GLY 11>
    
Step-wise building of L{Ensemble}s, L{Chain}s and L{Residue}s is supported through
a number of C{append} methods, for example:

    >>> residue = ProteinResidue(401, ProteinAlphabet.ALA)
    >>> s.chains['A'].residues.append(residue)
    
See L{EnsembleModelsCollection}, L{StructureChainsTable}, L{ChainResiduesCollection}
and L{ResidueAtomsTable} for more details.

Some other objects in this module of potential interest are the self-explanatory
L{SecondaryStructure} and L{TorsionAngles}.     
"""

import os
import re
import copy
import math
import numpy

import csb.io
import csb.core
import csb.numeric
import csb.bio.utils

from abc import ABCMeta, abstractmethod, abstractproperty

from csb.bio.sequence import SequenceTypes, SequenceAlphabets, AlignmentTypes


class AngleUnits(csb.core.enum):
    """
    Torsion angle unit types
    """
    Degrees='deg'; Radians='rad'
    
class SecStructures(csb.core.enum):
    """
    Secondary structure types
    """
    Helix='H'; Strand='E'; Coil='C'; Turn='T'; Bend='S';
    Helix3='G'; PiHelix='I'; BetaBridge='B'; Gap='-'
    
class ChemElements(csb.core.enum):
    """
    Periodic table elements
    """
    H=1; He=2; Li=3; Be=4; B=5; C=6; N=7; O=8; F=9; Ne=10; Na=11; Mg=12; Al=13; Si=14; P=15; 
    S=16; Cl=17; Ar=18; K=19; Ca=20; Sc=21; Ti=22; V=23; Cr=24; Mn=25; Fe=26; Co=27; Ni=28; 
    Cu=29; Zn=30; Ga=31; Ge=32; As=33; Se=34; Br=35; Kr=36; Rb=37; Sr=38; Y=39; Zr=40; Nb=41; 
    Mo=42; Tc=43; Ru=44; Rh=45; Pd=46; Ag=47; Cd=48; In=49; Sn=50; Sb=51; Te=52; I=53; Xe=54;
    Cs=55; Ba=56; Hf=72; Ta=73; W=74; Re=75; Os=76; Ir=77; Pt=78; Au=79; Hg=80; Tl=81; Pb=82; 
    Bi=83; Po=84; At=85; Rn=86; Fr=87; Ra=88; Rf=104; Db=105; Sg=106; Bh=107; Hs=108; Mt=109; 
    Ds=110; Rg=111; La=57; Ce=58; Pr=59; Nd=60; Pm=61; Sm=62; Eu=63; Gd=64; Tb=65; Dy=66; 
    Ho=67; Er=68; Tm=69; Yb=70; Lu=71; Ac=89; Th=90; Pa=91; U=92; Np=93; Pu=94; Am=95; Cm=96; 
    Bk=97; Cf=98; Es=99; Fm=100; Md=101; No=102; Lr=103; x=-1 


class Broken3DStructureError(ValueError):
    pass

class Missing3DStructureError(Broken3DStructureError):
    pass   
    
class InvalidOperation(Exception):
    pass

class EntityNotFoundError(csb.core.ItemNotFoundError):
    pass

class ChainNotFoundError(EntityNotFoundError):
    pass

class AtomNotFoundError(EntityNotFoundError):
    pass

class EntityIndexError(csb.core.CollectionIndexError):
    pass

class DuplicateModelIDError(csb.core.DuplicateKeyError):
    pass

class DuplicateChainIDError(csb.core.DuplicateKeyError):
    pass

class DuplicateResidueIDError(csb.core.DuplicateKeyError):
    pass

class DuplicateAtomIDError(csb.core.DuplicateKeyError):
    pass

class AlignmentArgumentLengthError(ValueError):
    pass

class BrokenSecStructureError(ValueError):
    pass

class UnknownSecStructureError(BrokenSecStructureError):
    pass

class AbstractEntity(object):
    """
    Base class for all protein structure entities.
    
    This class defines uniform interface of all entities (e.g. L{Structure},
    L{Chain}, L{Residue}) according to the Composite pattern. 
    """
    
    __metaclass__ = ABCMeta

    @abstractproperty
    def items(self):
        """
        Iterator over all immediate children of the entity
        @rtype: iterator of L{AbstractEntity}
        """
        pass

    def components(self, klass=None):
        """
        Return an iterator over all descendants of the entity.
        
        @param klass: return entities of the specified L{AbstractEntity} subclass
                      only. If None, traverse the hierarchy down to the lowest level.
        @param klass: class
        """
        for entity in CompositeEntityIterator.create(self, klass):
            if klass is None or isinstance(entity, klass):
                yield entity
        
    def transform(self, rotation, translation):
        """
        Apply in place RotationMatrix and translation Vector to all atoms.
        
        @type rotation: numpy array
        @type translation: numpy array 
        """
        for node in self.items:
            node.transform(rotation, translation)
    
    def get_coordinates(self, what=None, skip=False):
        """
        Extract the coordinates of the specified kind(s) of atoms and return 
        them as a list.
        
        @param what: a list of atom kinds, e.g. ['N', 'CA', 'C']
        @type what: list or None
        
        @return: a list of lists, each internal list corresponding to the coordinates 
                 of a 3D vector
        @rtype: list
        
        @raise Broken3DStructureError: if a specific atom kind cannot be retrieved from a residue
        """
        coords = [ ]
        
        for residue in self.components(klass=Residue):
            for atom_kind in (what or residue.atoms):
                try:
                    coords.append(residue.atoms[atom_kind].vector)
                except csb.core.ItemNotFoundError:
                    if skip:
                        continue
                    raise Broken3DStructureError('Could not retrieve {0} atom from the structure'.format(atom_kind))
            
        return numpy.array(coords)
    
class CompositeEntityIterator(object):
    """
    Iterates over composite L{AbstractEntity} hierarchies.
    
    @param node: root entity to traverse
    @type node: L{AbstractEntity}
    """
    
    def __init__(self, node):
            
        if not isinstance(node, AbstractEntity):
            raise TypeError(node)
            
        self._node = node
        self._stack = csb.core.Stack()
        
        self._inspect(node)
                
    def __iter__(self):
        return self

    def __next__(self):
        return self.next()      
        
    def next(self):

        while True:
            if self._stack.empty():
                raise StopIteration()
            
            try:
                current = self._stack.peek()
                node = next(current)
                self._inspect(node)
                return node
            
            except StopIteration:
                self._stack.pop()
                
    def _inspect(self, node):
        """
        Push C{node}'s children to the stack.
        """
        self._stack.push(node.items)
        
    @staticmethod
    def create(node, leaf=None):
        """
        Create a new composite iterator.
        
        @param leaf: if not None, return a L{ConfinedEntityIterator}
        @type leaf: class
        @rtype: L{CompositeEntityIterator} 
        """
        if leaf is None:
            return CompositeEntityIterator(node)
        else:
            return ConfinedEntityIterator(node, leaf)
                
class ConfinedEntityIterator(CompositeEntityIterator):
    """
    Iterates over composite L{AbstractEntity} hierarchies, but terminates
    the traversal of a branch once a specific node type is encountered.
    
    @param node: root entity to traverse
    @type node: L{AbstractEntity}
    @param leaf: traverse the hierarchy down to the specified L{AbstractEntity}
    @type leaf: class
    """
    def __init__(self, node, leaf):
        
        if not issubclass(leaf, AbstractEntity):
            raise TypeError(leaf)
        
        self._leaf = leaf
        super(ConfinedEntityIterator, self).__init__(node)              
    
    def _inspect(self, node):
        
        if not isinstance(node, self._leaf):
            self._stack.push(node.items)
            
class Ensemble(csb.core.AbstractNIContainer, AbstractEntity):
    """
    Represents an ensemble of multiple L{Structure} models.
    Provides a list-like access to these models:
    
        >>> ensemble[0]
        <Structure Model 1: accn, x chains>
        >>> ensemble.models[1]
        <Structure Model 1: accn, x chains>
    """
    
    def __init__(self):
        self._models = EnsembleModelsCollection()
        
    def __repr__(self):
        return "<Ensemble: {0} models>".format(self.models.length)        
        
    @property
    def _children(self):
        return self._models
    
    @property
    def models(self):
        """
        Access Ensembles's models by model ID
        @rtype: L{EnsembleModelsCollection}
        """
        return self._models
    
    @property
    def items(self):
        return iter(self._models)
        
    @property
    def first_model(self):
        """
        The first L{Structure} in the ensemble (if available)
        @rtype: L{Structure} or None
        """
        if len(self._models) > 0:
            return self[0]
        return None
    
    def to_pdb(self, output_file=None):
        """
        Dump the ensemble in PDB format.
        
        @param output_file: output file name or open stream
        @type output_file: str or stream
        """
        from csb.bio.io.wwpdb import PDBEnsembleFileBuilder
        
        if self.models.length < 1:
            raise InvalidOperation("Can't dump an empty ensemble")
        
        temp = csb.io.MemoryStream()

        builder = PDBEnsembleFileBuilder(temp)        
        builder.add_header(self.first_model)

        for model in self.models:
            builder.add_structure(model)

        builder.finalize()
        
        data = temp.getvalue()        
        temp.close()
        
        if not output_file:
            return data
        else:
            with csb.io.EntryWriter(output_file, close=False) as out:
                out.write(data)  
        
class EnsembleModelsCollection(csb.core.CollectionContainer):
    
    def __init__(self):
        
        super(EnsembleModelsCollection, self).__init__(type=Structure, start_index=1)
        self._models = set()
        
    def append(self, structure):
        """
        Add a new model
        
        @param structure: model to append
        @type structure: L{Structure}
        """
        
        if not structure.model_id or not str(structure.model_id).strip():
            raise ValueError("Invalid model identifier: '{0.model_id}'".format(structure))
        if structure.model_id in self._models:
            raise DuplicateModelIDError(structure.model_id) 
        else:
            return super(EnsembleModelsCollection, self).append(structure)
        
    @property
    def _exception(self):
        return EntityIndexError
    

class Structure(csb.core.AbstractNIContainer, AbstractEntity):
    """
    Represents a single model of a PDB 3-Dimensional molecular structure.
    Provides access to the L{Chain} objects, contained in the model:
    
        >>> structure['A']
        <Chain A: Protein>
        >>> structure.chains['A']
        <Chain A: Protein>
        >>> structure.items
        <iterator of Chain-s>    
    
    @param accession: accession number of the structure
    @type accession: str
    """
    def __init__(self, accession):
                        
        self._accession = None
        self._chains = StructureChainsTable(self)
        self._model_id = None
        self._resolution = None
        
        self.accession = accession

    def __repr__(self):
        return "<Structure Model {0.model_id}: {0.accession}, {1} chains>".format(self, self.chains.length)

    @property
    def _children(self):
        return self._chains
    
    @property
    def chains(self):
        """
        Access chains by their chain identifiers
        @rtype: L{StructureChainsTable}
        """
        return self._chains
    
    @property
    def items(self):
        for chain in self._chains:
            yield self._chains[chain]
                
    @property
    def first_chain(self):
        """
        The first L{Chain} in the structure (if available)
        @rtype: L{Chain} or None
        """        
        if len(self._chains) > 0:
            return next(self.items)
        return None
        
    @property
    def accession(self):
        """
        Accession number
        @rtype: str
        """        
        return self._accession
    @accession.setter
    def accession(self, accession):
        if accession is None:
            raise ValueError(accession)
        self._accession = str(accession).strip().lower()
        for c in self.chains:
            self.chains[c]._accession = self._accession
            
    @property
    def model_id(self):
        """
        Model ID
        @rtype: int
        """        
        return self._model_id
    @model_id.setter
    def model_id(self, value):
        self._model_id = value
        
    @property
    def resolution(self):
        """
        Resolution in Angstroms
        """
        return self._resolution
    @resolution.setter
    def resolution(self, value):
        if value is not None:
            value = float(value)
        self._resolution = value
                    
    def to_fasta(self):
        """
        Dump the structure in FASTA format. 
        
        @return: FASTA-formatted string with all chains in the structure
        @rtype: str
        
        @deprecated: this method will be removed soon. Use
                     L{csb.bio.sequence.ChainSequence.create} instead
        """
        fasta = []
        
        for chain in self.items:

            if chain.length > 0:
                fasta.append('>{0}'.format(chain.header))
                fasta.append(chain.sequence)
        
        return os.linesep.join(fasta)

    def to_pdb(self, output_file=None):
        """
        Dump the whole structure in PDB format.
        
        @param output_file: output file name or open stream
        @type output_file: str or stream
        """
        from csb.bio.io.wwpdb import PDBFileBuilder
                
        temp = csb.io.MemoryStream()
        builder = PDBFileBuilder(temp)
        
        builder.add_header(self)
        builder.add_structure(self)
        builder.finalize()
        
        data = temp.getvalue()        
        temp.close()
        
        if not output_file:
            return data
        else:
            with csb.io.EntryWriter(output_file, close=False) as out:
                out.write(data)

    @staticmethod
    def from_chain(chain):
        """
        A Structure factory, which instantiates and returns a new Structure with 
        chain as deep cpoy of chain

        @param chain: the chain which will comprise the new structure
        @type chain: L{Chain}

        @rtype: L{Structure}
        """
        structure = Structure("NONE")
        structure.chains.append(chain.clone())

        return structure


class StructureChainsTable(csb.core.DictionaryContainer):
    
    def __init__(self, structure=None, chains=None):
        self.__container = structure
        super(StructureChainsTable, self).__init__()
        
        if chains is not None:
            for chain in chains:
                self.append(chain)
        
    def __repr__(self):
        if len(self) > 0:
            return "<StructureChains: {0}>".format(', '.join(self))
        else:
            return "<StructureChains: empty>"
        
    @property
    def _exception(self):
        return ChainNotFoundError        
    
    def append(self, chain):
        """
        Add a new Chain to the structure.
        
        @param chain: the new chain to be appended
        @type chain: L{Chain}
        
        @raise DuplicateChainIDError: if a chain with same ID is already defined
        @raise InvalidOperation: if the chain is already associated with a structure
        """
        
        if chain._structure and chain._structure is not self.__container:
            raise InvalidOperation('This chain is already part of another structure')
        if chain.id in self:
            raise DuplicateChainIDError('A chain with ID {0} is already defined'.format(chain.id))
            
        super(StructureChainsTable, self).append(chain.id, chain)
        
        if self.__container:
            chain._accession = self.__container.accession
            chain._structure = self.__container

    def remove(self, id):
        """
        Remove a chain from the structure.

        @param id: ID of the chain to be detached
        @type id: str
        @raise ChainNotFoundError: if C{id} is not a valid chain ID 
        """
        chain = self[id]
        self._remove(id)
        chain._structure = None   
    
    def _update_chain_id(self, chain, new_id):
        
        if chain.id not in self or self[chain.id] is not chain:
            raise InvalidOperation(chain)
        
        self._remove(chain.id)
        
        if new_id in self:
            raise DuplicateChainIDError('Chain ID {0} is already defined'.format(id))
        
        super(StructureChainsTable, self).append(new_id, chain)
        
class Chain(csb.core.AbstractNIContainer, AbstractEntity):
    """
    Represents a polymeric chain. Provides list-like and rank-based access to
    the residues in the chain:
    
        >>> chain[0]
        <ProteinResidue [1]: SER None>
        >>> chain.residues[1]
        <ProteinResidue [1]: SER None>
    
    You can also access residues by their PDB sequence number:
    
        >>> chain.find(sequence_number=5, insertion_code='A')
        <ProteinResidue [1]: SER 5A>
    
    @param chain_id: ID of the new chain
    @type chain_id: str
    @param type: sequence type (a member of the L{SequenceTypes} enum)
    @type type: L{csb.core.EnumItem}
    @param name: name of the chain
    @type name: str
    @param residues: initialization list of L{Residue}-s
    @type residues: list
    @param accession: accession number of the chain
    @type accession: str
    @param molecule_id: MOL ID of the chain, if part of a polymer
    
    """
    def __init__(self, chain_id, type=SequenceTypes.Protein, name='',           
                 residues=None, accession=None, molecule_id=None):       

        self._id = str(chain_id).strip()
        self._accession = None
        self._type = None
        self._residues = ChainResiduesCollection(self, residues)
        self._secondary_structure = None
        self._molecule_id = molecule_id
        self._torsion_computed = False
        self._name = str(name).strip()
        
        self._structure = None
        
        self.type = type
        if accession is not None:
            self.accession = accession
            
    @staticmethod
    def from_sequence(sequence, id="_"):
        """
        Create a new chain from an existing sequence.
        
        @param sequence: source sequence
        @type sequence: L{csb.bio.sequence.AbstractSequence}
        
        @rtype: L{Chain}
        """
        
        chain = Chain(id, type=sequence.type)
        
        for ri in sequence.residues:
            residue = Residue.create(sequence.type, ri.rank, ri.type, sequence_number=ri.rank)
            chain.residues.append(residue)
            
        return chain
            
    @property
    def _children(self):
        return self._residues

    def __repr__(self):
        return "<Chain {0.id}: {0.type!r}>".format(self)        

    def __len__(self):
        return self._residues.length

    @property
    def id(self):
        """
        Chain's ID
        @rtype: str
        """
        return self._id
    @id.setter
    def id(self, id):
        if not isinstance(id, csb.core.string):
            raise ValueError(id)
        id = id.strip()
        if self._structure:
            self._structure.chains._update_chain_id(self, id)
        self._id = id
    
    @property
    def accession(self):
        """
        Accession number
        @rtype: str
        """        
        return self._accession
    @accession.setter
    def accession(self, accession):
        if self._structure:
            raise InvalidOperation("Only the accession of the parent structure can be altered")
        if accession is None:
            raise ValueError(accession)
        self._accession = str(accession).strip()
        
    @property
    def type(self):
        """
        Chain type - any member of L{SequenceTypes}
        @rtype: enum item
        """
        return self._type
    @type.setter
    def type(self, type):
        if type.enum is not SequenceTypes:
            raise TypeError(type)
        self._type = type

    @property
    def residues(self):
        """
        Rank-based access to Chain's L{Residue}s
        @rtype: L{ChainResiduesCollection}
        """
        return self._residues
    
    @property
    def items(self):
        return iter(self._residues)
    
    @property
    def torsion(self):
        """
        Torsion angles
        @rtype: L{TorsionAnglesCollection}
        """
        if not self._torsion_computed:
            raise InvalidOperation('The correctness of the data is not guaranteed '
                                   'until chain.compute_torsion() is invoked.')
            
        torsion = TorsionAnglesCollection()
        
        for r in self.residues:
            if r.torsion is None:
                torsion.append(TorsionAngles(None, None, None))
            else:
                torsion.append(r.torsion)
                
        return torsion
    
    @property
    def has_torsion(self):
        """
        True if C{Chain.compute_torsion} had been invoked
        @rtype: bool
        """
        return self._torsion_computed

    @property
    def length(self):
        """
        Number of residues
        @rtype: int
        """
        return self._residues.length
    
    @property
    def entry_id(self):
        """
        Accession number + chain ID
        @rtype: str
        """
        if self._accession and self._id:
            return self._accession + self._id
        else:
            return None
    
    @property
    def name(self):
        """
        Chain name
        @rtype: str
        """
        return self._name
    @name.setter
    def name(self, value):
        if value is not None:
            value = str(value).strip()
        self._name = value

    @property
    def molecule_id(self):
        """
        PDB MOL ID of this chain
        @rtype: int
        """
        return self._molecule_id
    @molecule_id.setter
    def molecule_id(self, value):
        self._molecule_id = value
                
    @property
    def header(self):
        """
        FASTA header in PDB format
        @rtype: str
        """
        header = "{0._accession}_{0._id} mol:{1} length:{0.length} {0.name}"
        return header.format(self, str(self.type).lower())

    @property
    def sequence(self):
        """
        Chain sequence
        @rtype: str
        """    
        sequence = []
        gap = str(self.alphabet.GAP)
        
        for residue in self.residues:
            if residue and residue.type:
                sequence.append(str(residue.type))
            else:
                sequence.append(gap)
                
        return ''.join(sequence)
    
    @property
    def alphabet(self):
        """
        Sequence alphabet corresponding to the current chain type
        @rtype: L{csb.core.enum}
        """
        return SequenceAlphabets.get(self.type)
        
    @property
    def secondary_structure(self):
        """
        Secondary structure (if available)
        @rtype: L{SecondaryStructure}
        """
        return self._secondary_structure
    @secondary_structure.setter
    def secondary_structure(self, ss):
        if not isinstance(ss, SecondaryStructure):
            raise TypeError(ss)
        if len(ss) > 0:
            if (ss[ss.last_index].end > self._residues.last_index):
                raise ValueError('Secondary structure out of range')
        self._secondary_structure = ss        
        
    def clone(self):
        """
        Make a deep copy of the chain. If this chain is part of a structure, 
        detach from it.
        
        @return: a deep copy of self
        @rtype: L{Chain}
        """
        start, end = self.residues.start_index, self.residues.last_index
        return self.subregion(start, end, clone=True)
        
    def subregion(self, start, end, clone=False):
        """
        Extract a subchain defined by [start, end]. If clone is True, this
        is a deep copy of the chain. Otherwise same as:
        
            >>> chain.residues[start : end + 1]
        
        but coordinates are checked and a Chain instance is returned.
        
        @param start: start position of the sub-region
        @type start: int
        @param end: end position
        @type end: int
        @param clone: if True, a deep copy of the sub-region is returned, 
                      otherwise - a shallow one
        @type clone: bool
        
        
        @return: a new chain, made from the residues of the extracted region
        @rtype: L{Chain}
        
        @raise IndexError: if start/end positions are out of range
        """
        if start < self.residues.start_index or start > self.residues.last_index:
            raise IndexError('The start position is out of range {0.start_index} .. {0.last_index}'.format(self.residues))
        if end < self.residues.start_index or end > self.residues.last_index:
            raise IndexError('The end position is out of range {0.start_index} .. {0.last_index}'.format(self.residues))
                
        residues = self.residues[start : end + 1]
        
        if clone:
            residues = [r.clone() for r in residues]
        
        chain = Chain(self.id, accession=self.accession, name=self.name, 
                      type=self.type, residues=residues, molecule_id=self.molecule_id)
        if chain.secondary_structure:
            chain.secondary_structure = self.secondary_structure.subregion(start, end)
        chain._torsion_computed = self._torsion_computed
        
        return chain  
        
    def find(self, sequence_number, insertion_code=None):
        """
        Get a residue by its original Residue Sequence Number and Insertion Code.
        
        @param sequence_number: PDB sequence number of the residue
        @type sequence_number: str
        @param insertion_code: PDB insertion code of the residue (if any)
        @type insertion_code: str
        
        @return: the residue object with such an ID
        @rtype: L{Residue}
        
        @raise EntityNotFoundError: if no residue with that ID exists
        """
        res_id = str(sequence_number).strip()
        
        if insertion_code is not None:
            insertion_code = str(insertion_code).strip()
            res_id += insertion_code

        return self.residues._get_residue(res_id)
    
    def compute_torsion(self):
        """
        Iterate over all residues in the chain, compute and set their torsion property.
        
        @raise Missing3DStructureError: when a 3D structure is absent
        @raise Broken3DStructureError: when a given atom cannot be retrieved from any residue
        """
        if self.type != SequenceTypes.Protein:                          
            raise NotImplementedError()
               
        for i, residue in enumerate(self.residues, start=self.residues.start_index):
            
            prev_residue, next_residue = None, None            
            
            if i > self.residues.start_index:
                prev_residue = self.residues[i - 1]         
            if i < self.residues.last_index:
                next_residue = self.residues[i + 1] 
                
            residue.torsion = residue.compute_torsion(prev_residue, next_residue, strict=False)
            
        self._torsion_computed = True
    
    def superimpose(self, other, what=['CA'], how=AlignmentTypes.Global):                       
        """
        Find the optimal fit between C{self} and C{other}. Return L{SuperimposeInfo}
        (RotationMatrix, translation Vector and RMSD), such that:
        
            >>> other.transform(rotation_matrix, translation_vector)
            
        will result in C{other}'s coordinates superimposed over C{self}.
        
        @param other: the subject (movable) chain
        @type other: L{Chain}
        @param what: a list of atom kinds, e.g. ['CA']
        @type what: list
        @param how: fitting method (global or local) - a member of the L{AlignmentTypes} enum
        @type how: L{csb.core.EnumItem}
        
        @return: superimposition info object, containing rotation matrix, translation 
                 vector and computed RMSD
        @rtype: L{SuperimposeInfo}
        
        @raise AlignmentArgumentLengthError: when the lengths of the argument chains differ 
        """ 
        if self.length != other.length or self.length < 1:
            raise AlignmentArgumentLengthError('Both chains must be of the same and positive length')
        
        x = self.get_coordinates(what)
        y = other.get_coordinates(what) 
        assert len(x) == len(y)

        if how == AlignmentTypes.Global:                                    
            r, t = csb.bio.utils.fit(x, y)
        else:
            r, t = csb.bio.utils.fit_wellordered(x, y)
            
        rmsd = csb.bio.utils.rmsd(x, y) 
        
        return SuperimposeInfo(r, t, rmsd=rmsd)
                              
    def align(self, other, what=['CA'], how=AlignmentTypes.Global):         
        """
        Align C{other}'s alpha carbons over self in space and return L{SuperimposeInfo}. 
        Coordinates of C{other} are overwritten in place using the rotation matrix
        and translation vector in L{SuperimposeInfo}. Alias for::
        
            R, t = self.superimpose(other, what=['CA'])
            other.transform(R, t)
            
        @param other: the subject (movable) chain
        @type other: L{Chain}
        @param what: a list of atom kinds, e.g. ['CA']
        @type what: list
        @param how: fitting method (global or local) - a member of the L{AlignmentTypes} enum
        @type how: L{csb.core.EnumItem}
        
        @return: superimposition info object, containing rotation matrix, translation 
                 vector and computed RMSD
        @rtype: L{SuperimposeInfo}        
        """
        result = self.superimpose(other, what=what, how=how)
        other.transform(result.rotation, result.translation)
        
        return result
    
    def rmsd(self, other, what=['CA']):
        """
        Compute the C-alpha RMSD against another chain (assuming equal length).
        Chains are superimposed with Least Squares Fit / Singular Value Decomposition.
        
        @param other: the subject (movable) chain
        @type other: L{Chain}
        @param what: a list of atom kinds, e.g. ['CA']
        @type what: list
        
        @return: computed RMSD over the specified atom kinds
        @rtype: float
        """
        
        if self.length != other.length or self.length < 1:
            raise ValueError('Both chains must be of the same and positive length '
                             '(got {0} and {1})'.format(self.length, other.length))
        
        x = self.get_coordinates(what)
        y = other.get_coordinates(what)
        assert len(x) == len(y)

        return csb.bio.utils.rmsd(x, y) 
    
    def tm_superimpose(self, other, what=['CA'], how=AlignmentTypes.Global):                    
        """
        Find the optimal fit between C{self} and C{other}. Return L{SuperimposeInfo}
        (RotationMatrix, translation Vector and TM-score), such that:
        
            >>> other.transform(rotation_matrix, translation_vector)
            
        will result in C{other}'s coordinates superimposed over C{self}.
        
        @param other: the subject (movable) chain
        @type other: L{Chain}
        @param what: a list of atom kinds, e.g. ['CA']
        @type what: list
        @param how: fitting method (global or local) - a member of the L{AlignmentTypes} enum
        @type how: L{csb.core.EnumItem}
        
        @return: superimposition info object, containing rotation matrix, translation 
                 vector and computed TM-score
        @rtype: L{SuperimposeInfo}
        
        @raise AlignmentArgumentLengthError: when the lengths of the argument chains differ         
        """
        
        if self.length != other.length or self.length < 1:
            raise ValueError('Both chains must be of the same and positive length')
        
        x = self.get_coordinates(what)
        y = other.get_coordinates(what)
        assert len(x) == len(y)
        
        L_ini_min = 0
        if how == AlignmentTypes.Global:                                            
            fit = csb.bio.utils.fit
        elif how == AlignmentTypes.Local:
            fit = csb.bio.utils.fit_wellordered
        else:
            # TMscore.f like search (slow)
            fit = csb.bio.utils.fit
            L_ini_min = 4
                
        r, t, tm = csb.bio.utils.tm_superimpose(x, y, fit, None, None, L_ini_min)
        
        return SuperimposeInfo(r,t, tm_score=tm)         
    
    def tm_score(self, other, what=['CA']):
        """
        Compute the C-alpha TM-Score against another chain (assuming equal chain length
        and optimal configuration - no fitting is done).        
        
        @param other: the subject (movable) chain
        @type other: L{Chain}
        @param what: a list of atom kinds, e.g. ['CA']
        @type what: list
                
        @return: computed TM-Score over the specified atom kinds
        @rtype: float        
        """

        if self.length != other.length or self.length < 1:
            raise ValueError('Both chains must be of the same and positive length')
        
        x = self.get_coordinates(what)
        y = other.get_coordinates(what)
        assert len(x) == len(y)

        return csb.bio.utils.tm_score(x, y)             

class ChainResiduesCollection(csb.core.CollectionContainer):
    
    def __init__(self, chain, residues):
        super(ChainResiduesCollection, self).__init__(type=Residue, start_index=1)
        self.__container = chain
        self.__lookup = { }
        
        if residues is not None:
            for residue in residues:
                self.append(residue)        
        
    def __repr__(self):
        if len(self) > 0:
            return "<ChainResidues: {0} ... {1}>".format(self[self.start_index], self[self.last_index])
        else:
            return "<ChainResidues: empty>"
        
    @property
    def _exception(self):
        return EntityIndexError    
        
    def append(self, residue):
        """
        Append a new residue to the chain.
        
        @param residue: the new residue
        @type residue: L{Residue}
        
        @raise DuplicateResidueIDError: if a residue with the same ID already exists
        """
        if residue.id and residue.id in self.__lookup:
            raise DuplicateResidueIDError('A residue with ID {0} is already defined within the chain'.format(residue.id))
        index = super(ChainResiduesCollection, self).append(residue)
        residue._container = self
        self.__container._torsion_computed = False
        self._add(residue)        
        return index
        
    def _contains(self, id):
        return id in self.__lookup
    
    def _remove(self, id):
        if id in self.__lookup:
            del self.__lookup[id]

    def _add(self, residue):
        self.__lookup[residue.id] = residue
            
    def _get_residue(self, id):
        try:
            return self.__lookup[id]
        except KeyError:
            raise EntityNotFoundError(id)            
        
class Residue(csb.core.AbstractNIContainer, AbstractEntity):
    """
    Base class representing a single residue. Provides a dictionary-like
    access to the atoms contained in the residue:
    
        >>> residue['CA']
        <Atom [3048]: CA>
        >>> residue.atoms['CA']
        <Atom [3048]: CA>
        >>> residue.items
        <iterator of Atom-s>
    
    @param rank: rank of the residue with respect to the chain
    @type rank: int
    @param type: residue type - a member of any L{SequenceAlphabets}
    @type type: L{csb.core.EnumItem}
    @param sequence_number: PDB sequence number of the residue
    @type sequence_number: str
    @param insertion_code: PDB insertion code, if any
    @type insertion_code: str
    """            
    def __init__(self, rank, type, sequence_number=None, insertion_code=None):
        
        self._type = None    
        self._label = None
        self._rank = int(rank)
        self._atoms = ResidueAtomsTable(self) 
        self._secondary_structure = None
        self._torsion = None
        self._sequence_number = None
        self._insertion_code = None
        self._container = None
        
        self.type = type
        self.id = sequence_number, insertion_code
        self.label = repr(type)
        
    @property
    def _children(self):
        return self._atoms
        
    def __repr__(self):
        return '<{1} [{0.rank}]: {0.type!r} {0.id}>'.format(self, self.__class__.__name__)
    
    @property
    def label(self):
        """
        Original residue label (different from C{Residue.type} for modified
        residues)
        @rtype: str        
        """
        return self._label
    @label.setter
    def label(self, value):
        self._label = str(value)
        
    @property
    def is_modified(self):
        """
        Return True id this is a modified residue
        @rtype: bool        
        """        
        return self.label != repr(self.type)
        
    @property
    def type(self):
        """
        Residue type - a member of any sequence alphabet
        @rtype: enum item
        """
        return self._type
    @type.setter
    def type(self, type):
        if type.enum not in (SequenceAlphabets.Protein, SequenceAlphabets.Nucleic, SequenceAlphabets.Unknown):
            raise TypeError(type)
        self._type = type
        
    @property
    def rank(self):
        """
        Residue's position in the sequence (1-based)
        @rtype: int
        """
        return self._rank
    
    @property
    def secondary_structure(self):
        """
        Secondary structure element this residue is part of
        @rtype: L{SecondaryStructureElement}        
        """
        return self._secondary_structure
    @secondary_structure.setter
    def secondary_structure(self, structure):
        if not isinstance(structure, SecondaryStructureElement):
            raise TypeError(structure)
        self._secondary_structure = structure
        
    @property
    def torsion(self):
        """
        Torsion angles
        @rtype: L{TorsionAngles}
        """
        return self._torsion
    @torsion.setter
    def torsion(self, torsion):
        if not isinstance(torsion, TorsionAngles):
            raise TypeError(torsion)
        self._torsion = torsion
    
    @property
    def atoms(self):
        """
        Access residue's atoms by atom name
        @rtype: L{ResidueAtomsTable}
        """
        return self._atoms
    
    @property
    def items(self):
        for atom in self._atoms:
            yield self._atoms[atom]        

    @property
    def sequence_number(self):
        """
        PDB sequence number (if residue.has_structure is True)
        @rtype: int
        """
        return self._sequence_number

    @property
    def insertion_code(self):
        """
        PDB insertion code (if defined)
        @rtype: str
        """
        return self._insertion_code
    
    @property
    def id(self):
        """
        PDB sequence number [+ insertion code]
        @rtype: str
        """
        if self._sequence_number is None:
            return None
        elif self._insertion_code is not None:
            return str(self._sequence_number) + self._insertion_code
        else:
            return str(self._sequence_number)
    @id.setter
    def id(self, value):
        sequence_number, insertion_code = value
        old_id = self.id
        id = ''
        if sequence_number is not None:
            sequence_number = int(sequence_number)
            id = str(sequence_number)
        if insertion_code is not None:
            insertion_code = str(insertion_code).strip()
            id += insertion_code
            if sequence_number is None:
                raise InvalidOperation('sequence_number must be defined when an insertion_code is specified.')
        if old_id != id:
            if self._container:
                if self._container._contains(id):
                    raise DuplicateResidueIDError('A residue with ID {0} is already defined within the chain'.format(id))
                self._container._remove(old_id)
            self._sequence_number = sequence_number
            self._insertion_code = insertion_code
            if self._container:
                self._container._add(self)
    
    @property
    def has_structure(self):
        """
        True if this residue has any atoms
        @rtype: bool
        """
        return len(self.atoms) > 0
        
    def get_coordinates(self, what=None, skip=False):
        
        coords = []
        
        if not self.has_structure:
            if skip:
                return numpy.array([])
            raise Missing3DStructureError(self)
        
        for atom_kind in (what or self.atoms):
            if atom_kind in self.atoms:
                coords.append(self.atoms[atom_kind].vector)                 
            else:
                if skip:
                    continue
                raise Broken3DStructureError('Could not retrieve {0} atom'.format(atom_kind))

        return numpy.array(coords)
                    
    def clone(self):
        
        container = self._container
        self._container = None
        clone = copy.deepcopy(self)
        self._container = container
        
        return clone
        
    @staticmethod
    def create(sequence_type, *a, **k):
        """
        Residue factory method, which returns the proper L{Residue} instance based on 
        the specified C{sequence_type}. All additional arguments are used to initialize
        the subclass by passing them automatically to the underlying constructor. 
        
        @param sequence_type: create a Residue of that SequenceType 
        @type sequence_type: L{csb.core.EnumItem}
        
        @return: a new residue of the proper subclass
        @rtype: L{Residue} subclass
        
        @raise ValueError: if the sequence type is not known
        """        
        if sequence_type == SequenceTypes.Protein:                                      
            return ProteinResidue(*a, **k)
        elif sequence_type == SequenceTypes.NucleicAcid:                                
            return NucleicResidue(*a, **k)
        elif sequence_type == SequenceTypes.Unknown:
            return UnknownResidue(*a, **k)
        else:
            raise ValueError(sequence_type)        
           
class ProteinResidue(Residue):
    """
    Represents a single amino acid residue.
    
    @param rank: rank of the residue with respect to the chain
    @type rank: int
    @param type: residue type - a member of 
                 L{csb.bio.sequence.SequenceAlphabets.Protein}
    @type type: L{csb.core.EnumItem}
    @param sequence_number: PDB sequence number of the residue
    @type sequence_number: str
    @param insertion_code: PDB insertion code, if any
    @type insertion_code: str    
    """
    
    def __init__(self, rank, type, sequence_number=None, insertion_code=None):
          
        if isinstance(type, csb.core.string):
            try:
                if len(type) == 3:
                    type = csb.core.Enum.parsename(SequenceAlphabets.Protein, type)
                else:    
                    type = csb.core.Enum.parse(SequenceAlphabets.Protein, type)          
            except (csb.core.EnumMemberError, csb.core.EnumValueError):
                raise ValueError("'{0}' is not a valid amino acid".format(type))
        elif type.enum is not SequenceAlphabets.Protein:
            raise TypeError(type)
            
        super(ProteinResidue, self).__init__(rank, type, sequence_number, insertion_code)  
         
    def compute_torsion(self, prev_residue, next_residue, strict=True):
        """
        Compute the torsion angles of the current residue with neighboring residues
        C{prev_residue} and C{next_residue}. 
        
        @param prev_residue: the previous residue in the chain
        @type prev_residue: L{Residue}
        @param next_residue: the next residue in the chain
        @type next_residue: L{Residue}
        @param strict: if True, L{Broken3DStructureError} is raised if either C{prev_residue} 
                       or C{next_residue} has a broken structure, else the error is silently
                       ignored and an empty L{TorsionAngles} instance is created
        @type strict: bool
                        
        @return: a L{TorsionAngles} object, holding the phi, psi and omega values
        @rtype: L{TorsionAngles}
        
        @raise Broken3DStructureError: when a specific atom cannot be found 
        """       
        if prev_residue is None and next_residue is None:
            raise ValueError('At least one neighboring residue is required to compute the torsion.')
   
        angles = TorsionAngles(None, None, None, units=AngleUnits.Degrees)
        
        for residue in (self, prev_residue, next_residue):
            if residue is not None and not residue.has_structure:
                if strict:
                    raise Missing3DStructureError(repr(residue))
                elif residue is self:
                    return angles
        
        try:
            n = self._atoms['N'].vector
            ca = self._atoms['CA'].vector
            c = self._atoms['C'].vector
        except csb.core.ItemNotFoundError as missing_atom:
            if strict:
                raise Broken3DStructureError('Could not retrieve {0} atom from the current residue {1!r}.'.format(
                                                                                                missing_atom, self))
            else:
                return angles
        
        try:
            if prev_residue is not None and prev_residue.has_structure:
                prev_c = prev_residue._atoms['C'].vector
                angles.phi = csb.numeric.dihedral_angle(prev_c, n, ca, c)
        except csb.core.ItemNotFoundError as missing_prevatom:
            if strict:
                raise Broken3DStructureError('Could not retrieve {0} atom from the i-1 residue {1!r}.'.format(
                                                                                    missing_prevatom, prev_residue))    
        try:
            if next_residue is not None and next_residue.has_structure:    
                next_n = next_residue._atoms['N'].vector
                angles.psi = csb.numeric.dihedral_angle(n, ca, c, next_n)
                next_ca = next_residue._atoms['CA'].vector
                angles.omega = csb.numeric.dihedral_angle(ca, c, next_n, next_ca)
        except csb.core.ItemNotFoundError as missing_nextatom:
            if strict:
                raise Broken3DStructureError('Could not retrieve {0} atom from the i+1 residue {1!r}.'.format(
                                                                                    missing_nextatom, next_residue))              
                                
        return angles

class NucleicResidue(Residue):
    """
    Represents a single nucleotide residue.
    
    @param rank: rank of the residue with respect to the chain
    @type rank: int
    @param type: residue type - a member of 
                 L{csb.bio.sequence.SequenceAlphabets.Nucleic}
    @type type: L{csb.core.EnumItem}
    @param sequence_number: PDB sequence number of the residue
    @type sequence_number: str
    @param insertion_code: PDB insertion code, if any
    @type insertion_code: str        
    """
    
    def __init__(self, rank, type, sequence_number=None, insertion_code=None):
        
        if isinstance(type, csb.core.string):
            try:
                if len(type) > 1:
                    type = csb.core.Enum.parsename(SequenceAlphabets.Nucleic, type)
                else:    
                    type = csb.core.Enum.parse(SequenceAlphabets.Nucleic, type)
            except (csb.core.EnumMemberError, csb.core.EnumValueError):
                raise ValueError("'{0}' is not a valid nucleotide".format(type))
        elif type.enum is not SequenceAlphabets.Nucleic:
            raise TypeError(type)
            
        super(NucleicResidue, self).__init__(rank, type, sequence_number, insertion_code)  
        self.label = str(type)
        
    @property
    def is_modified(self):
        return self.label != str(self.type)        

class UnknownResidue(Residue):
    
    def __init__(self, rank, type, sequence_number=None, insertion_code=None):

        super(UnknownResidue, self).__init__(rank, SequenceAlphabets.Unknown.UNK,
                                             sequence_number, insertion_code)
            
class ResidueAtomsTable(csb.core.DictionaryContainer):
    """ 
    Represents a collection of atoms. Provides dictionary-like access,
    where PDB atom names are used for lookup.
    """
    def __init__(self, residue, atoms=None):
        
        self.__residue = residue
        super(ResidueAtomsTable, self).__init__()
        
        if atoms is not None:
            for atom in atoms:
                self.append(atom)
        
    def __repr__(self):
        if len(self) > 0:
            return "<ResidueAtoms: {0}>".format(', '.join(self.keys()))
        else:
            return "<ResidueAtoms: empty>"
        
    @property
    def _exception(self):
        return AtomNotFoundError
    
    def append(self, atom):
        """
        Append a new Atom to the catalog.
        
        If the atom has an alternate position, a disordered proxy will be created instead and the 
        atom will be appended to the L{DisorderedAtom}'s list of children. If a disordered atom 
        with that name already exists, the atom will be appended to its children only.
        If an atom with the same name exists, but it was erroneously not marked as disordered,
        that terrible condition will be fixed too.
        
        @param atom: the new atom to append
        @type atom: L{Atom}
        
        @raise DuplicateAtomIDError: if an atom with the same sequence number and 
                                     insertion code already exists in that residue
        """
        if atom.residue and atom.residue is not self.__residue:
            raise InvalidOperation('This atom is part of another residue')
        if atom.alternate or (atom.name in self and isinstance(self[atom.name], DisorderedAtom)):
            if atom.name not in self:
                atom._residue = self.__residue
                dis_atom = DisorderedAtom(atom)
                super(ResidueAtomsTable, self).append(dis_atom.name, dis_atom)
            else:
                if not isinstance(self[atom.name], DisorderedAtom):
                    buggy_atom = self[atom.name]
                    assert buggy_atom.alternate in (None, False)
                    buggy_atom.alternate = True
                    self.update(atom.name, DisorderedAtom(buggy_atom))
                if not atom.alternate:
                    atom.alternate = True 
                atom._residue = self.__residue
                self[atom.name].append(atom)          
        else:
            if atom.name in self:
                raise DuplicateAtomIDError('Atom {0} is already defined for {1}'.format(
                                                                        atom.name, self.__residue))
            else:                   
                super(ResidueAtomsTable, self).append(atom.name, atom)
                atom._residue = self.__residue
        
    def update(self, atom_name, atom):
        """ 
        Update the atom with the specified name.
        
        @param atom_name: update key
        @type atom_name: str
        @param atom: new value for this key
        @type atom: L{Atom}
        
        @raise ValueError: if C{atom} has a different name than C{atom_name}
        """
        if atom.name != atom_name:
            raise ValueError("Atom's name differs from the specified key.")
        if atom.residue is not self.__residue:
            atom._residue = self.__residue
        
        super(ResidueAtomsTable, self)._update({atom_name: atom})  
    
class Atom(AbstractEntity):
    """
    Represents a single atom in space.
    
    @param serial_number: atom's UID
    @type serial_number: int
    @param name: atom's name
    @type name: str
    @param element: corresponding L{ChemElements}
    @type element: L{csb.core.EnumItem}
    @param vector: atom's coordinates
    @type vector: numpy array
    @param alternate: if True, means that this is a wobbling atom with multiple alternative 
                      locations
    @type alternate: bool
    """        
    def __init__(self, serial_number, name, element, vector, alternate=False):
        
        self._serial_number = None
        self._name = None
        self._element = None
        self._residue = None
        self._vector = None
        self._alternate = False        
        self._bfactor = None
        self._occupancy = None
        self._charge = None

        if not isinstance(name, csb.core.string):
            raise TypeError(name)
        name_compact = name.strip()
        if len(name_compact) < 1:
            raise ValueError(name)
        self._name = name_compact
        self._full_name = name
            
        if isinstance(element, csb.core.string):
            element = csb.core.Enum.parsename(ChemElements, element)
        elif element is None:
            pass
        elif element.enum is not ChemElements:
            raise TypeError(element)
        self._element = element

        # pass type- and value-checking control to setters
        self.serial_number = serial_number
        self.vector = vector
        self.alternate = alternate
        
    def __repr__(self):
        return "<Atom [{0.serial_number}]: {0.name}>".format(self)
        
    def __lt__(self, other):
        return self.serial_number < other.serial_number
    
    def transform(self, rotation, translation):
        
        vector = numpy.dot(self.vector, numpy.transpose(rotation)) + translation
        self.vector = vector
    
    def get_coordinates(self, what=None, skip=False):
        
        if what is None:
            what = [self.name]
            
        if self.name in what:
            return numpy.array([self.vector.copy()])
        elif skip:
            return numpy.array([])
        else:
            raise Missing3DStructureError()
        
    def clone(self):
        
        residue = self._residue
        self._residue = None
        clone = copy.deepcopy(self)
        self._residue = residue
        
        return clone

    @property
    def serial_number(self):
        """
        PDB serial number
        @rtype: int
        """        
        return self._serial_number
    @serial_number.setter
    def serial_number(self, number):
        if not isinstance(number, int) or number < 1:
            raise TypeError(number)
        self._serial_number = number
    
    @property
    def name(self):
        """
        PDB atom name (trimmed)
        @rtype: str
        """
        return self._name
        
    @property
    def element(self):
        """
        Chemical element - a member of L{ChemElements}
        @rtype: enum item
        """
        return self._element
    
    @property
    def residue(self):
        """
        Residue instance that owns this atom (if available)
        @rtype: L{Residue}
        """
        return self._residue
    @residue.setter
    def residue(self, residue):
        if self._residue:
            raise InvalidOperation('This atom is already part of a residue.')
        if not isinstance(residue, Residue):
            raise TypeError(residue)
        self._residue = residue
    
    @property
    def vector(self):
        """
        Atom's 3D coordinates (x, y, z)
        @rtype: numpy.array
        """
        return self._vector
    @vector.setter
    def vector(self, vector):
        if numpy.shape(vector) != (3,):
            raise ValueError("Three dimensional vector expected")
        self._vector = numpy.array(vector)
        
    @property
    def alternate(self):
        """
        Alternative location flag
        @rtype: str
        """
        return self._alternate
    @alternate.setter
    def alternate(self, value):
        self._alternate = value
    
    @property
    def bfactor(self):
        """
        Temperature factor
        @rtype: float
        """
        return self._bfactor
    @bfactor.setter
    def bfactor(self, value):
        self._bfactor = value
    
    @property
    def occupancy(self):
        """
        Occupancy number
        @rtype: float
        """        
        return self._occupancy
    @occupancy.setter
    def occupancy(self, value):
        self._occupancy = value
    
    @property
    def charge(self):
        """
        Charge
        @rtype: int
        """         
        return self._charge
    @charge.setter
    def charge(self, value):
        self._charge = value
    
    @property
    def items(self):
        return iter([])
        
class DisorderedAtom(csb.core.CollectionContainer, Atom):
    """
    A wobbling atom, which has alternative locations. Each alternative is represented 
    as a 'normal' L{Atom}. The atom with a highest occupancy is selected as a representative,
    hence a DisorderedAtom behaves as a regular L{Atom} (proxy of the representative) as well
    as a collection of Atoms. 
    
    @param atom: the first atom to be appended to the collection of alternatives. It
                 is automatically defined as a representative, until a new atom with 
                 higher occupancy is appended to the collection
    @type atom: L{Atom}
    """  
        
    def __init__(self, atom):
        
        super(DisorderedAtom, self).__init__(type=Atom)
        
        self.__rep = None
        self.__alt = {}
        
        self.append(atom)

    def __getattr__(self, name):
        try:
            return object.__getattribute__(self, name)
        except AttributeError:
            subject = object.__getattribute__(self, '_DisorderedAtom__rep')
            return getattr(subject, name)             
            
    def append(self, atom):
        """
        Append a new atom to the collection of alternatives.
        
        @param atom: the new alternative
        @type atom: L{Atom}
        """
        self.__update_rep(atom)
        self.__alt[atom.alternate] = atom
        
        super(DisorderedAtom, self).append(atom)
        
    def find(self, altloc):
        """
        Retrieve a specific atom by its altloc identifier.
        
        @param altloc: alternative location identifier
        @type altloc: str
        
        @rtype: L{Atom}  
        """
        if altloc in self.__alt:
            return self.__alt[altloc]
        else:
            for atom in self:
                if atom.alternate == altloc:
                    return Atom
        
        raise EntityNotFoundError(altloc)
    
    def transform(self, rotation, translation):
        
        for atom in self:
            atom.transform(rotation, translation)
                
    def __update_rep(self, atom):
        
        if self.__rep is None or \
        ((self.__rep.occupancy != atom.occupancy) and (self.__rep.occupancy < atom.occupancy)):
        
            self.__rep = atom
            
    def __repr__(self):
        return "<DisorderedAtom: {0.length} alternative locations>".format(self)
            
class SuperimposeInfo(object):
    """
    Describes a structural alignment result.
    
    @type rotation: Numpy Array
    @type translation: L{Vector}
    @type rmsd: float
    """
    def __init__(self, rotation, translation, rmsd=None, tm_score=None):
        
        self.rotation = rotation
        self.translation = translation
        self.rmsd = rmsd
        self.tm_score = tm_score

class SecondaryStructureElement(object):
    """ 
    Describes a Secondary Structure Element.
    
    @param start: start position with reference to the chain
    @type start: float
    @param end: end position with reference to the chain
    @type end: float    
    @param type: element type - a member of the L{SecStructures} enum
    @type type: csb.core.EnumItem
    @param score: secondary structure prediction confidence, if available
    @type score: int
    
    @raise IndexError: if start/end positions are out of range
    """    
    def __init__(self, start, end, type, score=None):
        
        if not (0 < start <= end):
            raise IndexError('Element coordinates are out of range: 1 <= start <= end.')
                
        self._start = None
        self._end = None
        self._type = None
        self._score = None

        self.start = start
        self.end = end
        self.type = type
                
        if score is not None: 
            self.score = score
            
    def __lt__(self, other):
        return self.start < other.start
    
    def __eq__(self, other):
        return (self.type == other.type 
                and self.start == other.start 
                and self.end == other.end) 
    
    def __str__(self):
        return self.to_string()
    
    def __repr__(self):
        return "<{0.type!r}: {0.start}-{0.end}>".format(self)
    
    @property
    def start(self):
        """
        Start position (1-based)
        @rtype: int
        """
        return self._start
    @start.setter
    def start(self, value):
        if value is not None:
            value = int(value)
            if value < 1:
                raise ValueError(value)
        self._start = value
    
    @property
    def end(self):
        """
        End position (1-based)
        @rtype: int
        """        
        return self._end
    @end.setter
    def end(self, value):
        if value is not None:
            value = int(value)
            if value < 1:
                raise ValueError(value)            
        self._end = value 
    
    @property
    def type(self):
        """
        Secondary structure type - a member of L{SecStructures}
        @rtype: enum item
        """        
        return self._type
    @type.setter
    def type(self, value):
        if isinstance(value, csb.core.string):
            value = csb.core.Enum.parse(SecStructures, value)
        if not value.enum is SecStructures:
            raise TypeError(value)
        self._type = value
            
    @property
    def length(self):
        """
        Number of residues covered by this element
        @rtype: int
        """
        return self.end - self.start + 1
    
    @property
    def score(self):
        """
        Secondary structure confidence values for each residue
        @rtype: L{CollectionContainer}
        """
        return self._score
    @score.setter
    def score(self, scores):
        if not len(scores) == self.length:
            raise ValueError('There must be a score entry for each residue in the element.')        
        self._score = csb.core.CollectionContainer(
                                items=list(scores), type=int, start_index=self.start)
    
    def overlaps(self, other):
        """
        Return True if C{self} overlaps with C{other}.
        
        @type other: L{SecondaryStructureElement}
        @rtype: bool
        """
        this = set(range(self.start, self.end + 1))
        that = set(range(other.start, other.end + 1))
        return not this.isdisjoint(that)
    
    def merge(self, other):
        """
        Merge C{self} and C{other}.

        @type other: L{SecondaryStructureElement}
                
        @return: a new secondary structure element
        @rtype: L{SecondaryStructureElement}
        
        @bug: confidence scores are lost
        """
        if not self.overlaps(other):
            raise ValueError("Can't merge non-overlapping secondary structures")
        elif self.type != other.type:
            raise ValueError("Can't merge secondary structures of different type")            
        
        start = min(self.start, other.start)
        end = max(self.end, other.end)
        assert self.type == other.type
        
        return SecondaryStructureElement(start, end, self.type)    
    
    def to_string(self):
        """
        Dump the element as a string.
        
        @return: string representation of the element
        @rtype: str
        """
        return str(self.type) * self.length
    
    def simplify(self):
        """
        Convert to three-state secondary structure (Helix, Strand, Coil).
        """           
        if self.type in (SecStructures.Helix, SecStructures.Helix3, SecStructures.PiHelix):
            self.type = SecStructures.Helix
        elif self.type in (SecStructures.Strand, SecStructures.BetaBridge):
            self.type = SecStructures.Strand
        elif self.type in (SecStructures.Coil, SecStructures.Turn, SecStructures.Bend):
            self.type = SecStructures.Coil
        elif self.type == SecStructures.Gap or self.type is None:
            pass
        else:
            assert False, 'Unhandled SS type: ' + repr(self.type)    

class SecondaryStructure(csb.core.CollectionContainer):
    """
    Describes the secondary structure of a chain.
    Provides an index-based access to the secondary structure elements of the chain.
    
    @param string: a secondary structure string (e.g. a PSI-PRED output)
    @type string: str
    @param conf_string: secondary structure prediction confidence values, if available
    @type conf_string: str
    """
    def __init__(self, string=None, conf_string=None):

        super(SecondaryStructure, self).__init__(type=SecondaryStructureElement, start_index=1)
        
        self._minstart = None
        self._maxend = None
         
        if string is not None:
            for motif in SecondaryStructure.parse(string, conf_string):
                self.append(motif)
                
    def __str__(self):
        return self.to_string()
    
    def append(self, element):
        """
        Add a new SecondaryStructureElement. Then sort all elements by
        their start position.
        """
        super(SecondaryStructure, self).append(element)
        super(SecondaryStructure, self)._sort()
        
        if self._minstart is None or element.start < self._minstart:
            self._minstart = element.start
        if self._maxend is None or element.end > self._maxend:
            self._maxend = element.end            
                        
    @staticmethod  
    def parse(string, conf_string=None):
        """ 
        Parse secondary structure from DSSP/PSI-PRED output string.
        
        @param string: a secondary structure string (e.g. a PSI-PRED output)
        @type string: str
        @param conf_string: secondary structure prediction confidence values, if available
        @type conf_string: str
                
        @return: a list of L{SecondaryStructureElement}s.
        @rtype: list
        
        @raise ValueError: if the confidence string is not of the same length
        """
        if not isinstance(string, csb.core.string):
            raise TypeError(string)
                
        string = ''.join(re.split('\s+', string))
        if conf_string is not None:
            conf_string = ''.join(re.split('\s+', conf_string))
            if not len(string) == len(conf_string):
                raise ValueError('The confidence string has unexpected length.')
        motifs = [ ]

        if not len(string) > 0:
            raise ValueError('Empty Secondary Structure string')      
        
        currel = string[0]
        start = 0
                
        for i, char in enumerate(string + '.'):
            
            if currel != char:
                try:
                    type = csb.core.Enum.parse(SecStructures, currel)
                except csb.core.EnumValueError:
                    raise UnknownSecStructureError(currel)
                confidence = None
                if conf_string is not None:
                    confidence = list(conf_string[start : i])
                    confidence = list(map(int, confidence))
                motif = SecondaryStructureElement(start + 1, i, type, confidence)
                motifs.append(motif)
                
                currel = char
                start = i

        return motifs
    
    @property
    def start(self):
        """
        Start position of the leftmost element
        @rtype: int
        """
        return self._minstart
        
    @property
    def end(self):
        """
        End position of the rightmost element
        @rtype: int
        """        
        return self._maxend
    
    def clone(self):
        """
        @return: a deep copy of the object
        """
        return copy.deepcopy(self)
        
    def to_three_state(self):
        """
        Convert to three-state secondary structure (Helix, Strand, Coil).
        """           
        for e in self:
            e.simplify()
    
    def to_string(self, chain_length=None):
        """
        Get back the string representation of the secondary structure.
        
        @return: a string of secondary structure elements
        @rtype: str
        
        @bug: [CSB 0000003] If conflicting elements are found at a given rank,
              this position is represented as a coil.
        """  
        gap = str(SecStructures.Gap)
        coil = str(SecStructures.Coil)
        
        if chain_length is None:
            chain_length = max(e.end for e in self)

        ss = []
        
        for pos in range(1, chain_length + 1):
            elements = self.at(pos)
            if len(elements) > 0:
                if len(set(e.type for e in elements)) > 1:
                    ss.append(coil)                         # [CSB 0000003]                     
                else:    
                    ss.append(elements[0].to_string()) 
            else:
                ss.append(gap)        

        return ''.join(ss)
    
    def at(self, rank, type=None):
        """
        @return: all secondary structure elements covering the specifid position
        @rtype: tuple of L{SecondaryStructureElement}s 
        """
        return self.scan(start=rank, end=rank, filter=type, loose=True, cut=True)
    
    def scan(self, start, end, filter=None, loose=True, cut=True):
        """
        Get all secondary structure elements within the specified [start, end] region.
        
        @param start: the start position of the region, 1-based, inclusive
        @type start: int
        @param end: the end position of the region, 1-based, inclusive
        @type end: int     
        @param filter: return only elements of the specified L{SecStructures} kind
        @type filter: L{csb.core.EnumItem}
        @param loose: grab all fully or partially matching elements within the region.
                      if False, return only the elements which strictly reside within 
                      the region
        @type loose: bool
        @param cut: if an element is partially overlapping with the start..end region, 
                    cut its start and/or end to make it fit into the region. If False, 
                    return the elements with their real lengths
        @type cut: bool

        @return: a list of deep-copied L{SecondaryStructureElement}s, sorted by their 
                 start position
        @rtype: tuple of L{SecondaryStructureElement}s
        """        
        matches = [ ]
        
        for m in self:            
            if filter and m.type != filter:
                continue
            
            if loose:
                if start <= m.start <= end or start <= m.end <= end or (m.start <= start and m.end >= end):
                    partmatch = copy.deepcopy(m)
                    if cut:
                        if partmatch.start < start:
                            partmatch.start = start
                        if partmatch.end > end:
                            partmatch.end = end
                        if partmatch.score:  
                            partmatch.score = partmatch.score[start : end + 1]
                    matches.append(partmatch) 
            else:
                if m.start >= start and m.end <= end:
                    matches.append(copy.deepcopy(m))                                    

        matches.sort()
        return tuple(matches)
    
    def q3(self, reference, relaxed=True):
        """
        Compute Q3 score.
        
        @param reference: reference secondary structure
        @type reference: L{SecondaryStructure}
        @param relaxed: if True, treat gaps as coils
        @type relaxed: bool
        
        @return: the percentage of C{reference} residues with identical
                 3-state secondary structure.
        @rtype: float
        """
        
        this = self.clone()
        this.to_three_state()
        
        ref = reference.clone()
        ref.to_three_state()
        
        total = 0
        identical = 0
        
        def at(ss, rank):
            elements = ss.at(rank)
            if len(elements) == 0:
                return None
            elif len(elements) > 1:
                raise ValueError('Flat secondary structure expected')
            else:
                return elements[0] 
        
        for rank in range(ref.start, ref.end + 1):
            q = at(this, rank)
            s = at(ref, rank)            

            if s:
                if relaxed or s.type != SecStructures.Gap:
                    total += 1
                    if q:
                        if q.type == s.type:
                            identical += 1
                        elif relaxed:
                            pair = set([q.type, s.type])
                            match = set([SecStructures.Gap, SecStructures.Coil])
                            if pair.issubset(match):
                                identical += 1
                    
        if total == 0:
            return 0.0
        else:
            return identical * 100.0 / total
        
    def subregion(self, start, end):
        """
        Same as C{ss.scan(...cut=True)}, but also shift the start-end positions
        of all motifs and return a L{SecondaryStructure} instance instead of a list.
        
        @param start: start position of the subregion, with reference to the chain
        @type start: int
        @param end: start position of the subregion, with reference to the chain
        @type end: int
        
        @return: a deep-copy sub-fragment of the original L{SecondaryStructure}
        @rtype: L{SecondaryStructure}
        """
        sec_struct = SecondaryStructure()
        
        for motif in self.scan(start, end, loose=True, cut=True):
            
            motif.start = motif.start - start + 1
            motif.end = motif.end - start + 1
            if motif.score:
                motif.score = list(motif.score) # this will automatically fix the score indices in the setter
            sec_struct.append(motif) 
            
        return sec_struct
        
class TorsionAnglesCollection(csb.core.CollectionContainer):
    """
    Describes a collection of torsion angles. Provides 1-based list-like access.
    
    @param items: an initialization list of L{TorsionAngles}
    @type items: list
    """  
    def __init__(self, items=None, start=1):
        
        super(TorsionAnglesCollection, self).__init__(
                                items,type=TorsionAngles, start_index=start)

    def __repr__(self):
        if len(self) > 0:
            return "<TorsionAnglesList: {0} ... {1}>".format(self[self.start_index], self[self.last_index])
        else:
            return "<TorsionAnglesList: empty>"        
        
    @property
    def phi(self):
        """
        List of all phi angles
        @rtype: list
        """   
        return [a.phi for a in self]

    @property
    def psi(self):
        """
        List of all psi angles
        @rtype: list
        """           
        return [a.psi for a in self]            

    @property
    def omega(self):  
        """
        List of all omega angles
        @rtype: list
        """                
        return [a.omega for a in self] 
    
    def update(self, angles):
        self._update(angles)   
    
    def rmsd(self, other):
        """
        Calculate the Circular RSMD against another TorsionAnglesCollection.
        
        @param other: subject (right-hand-term)
        @type other: L{TorsionAnglesCollection}
        
        @return: RMSD based on torsion angles
        @rtype: float
        
        @raise Broken3DStructureError: on discontinuous torsion angle collections
        (phi and psi values are still allowed to be absent at the termini)
        @raise ValueError: on mismatching torsion angles collection lengths
        """                   
        if len(self) != len(other) or len(self) < 1:
            raise ValueError('Both collections must be of the same and positive length')
        
        length = len(self)
        query, subject = [], []
                
        for n, (q, s) in enumerate(zip(self, other), start=1):
            
            q = q.copy()
            q.to_radians()
            
            s = s.copy()
            s.to_radians()
            
            if q.phi is None or s.phi is None:
                if n == 1:
                    q.phi = s.phi = 0.0
                else:
                    raise Broken3DStructureError('Discontinuous torsion angles collection at {0}'.format(n))
                    
            if q.psi is None or s.psi is None:
                if n == length:
                    q.psi = s.psi = 0.0
                else:
                    raise Broken3DStructureError('Discontinuous torsion angles collection at {0}'.format(n))
                
            query.append([q.phi, q.psi])
            subject.append([s.phi, s.psi])
            
        return csb.bio.utils.torsion_rmsd(numpy.array(query), numpy.array(subject))
           
class TorsionAngles(object):
    """
    Describes a collection of phi, psi and omega backbone torsion angles.
    
    It is assumed that the supplied values are either None, or fitting into 
    the range of [-180, +180] for AngleUnites.Degrees and [0, 2pi] for Radians.  
    
    @param phi: phi angle value in C{units}
    @type phi: float
    @param psi: psi angle value in C{units}
    @type psi: float
    @param omega: omega angle value in C{units}
    @type omega: float    
    @param units: any of L{AngleUnits}'s enum members
    @type units: L{csb.core.EnumItem}
    
    @raise ValueError: on invalid/unknown units
    """
                     
    def __init__(self, phi, psi, omega, units=AngleUnits.Degrees):
        
        try:
            if isinstance(units, csb.core.string):
                units = csb.core.Enum.parse(AngleUnits, units, ignore_case=True)
            else:
                if units.enum is not AngleUnits:
                    raise TypeError(units)
                
        except ValueError:
            raise ValueError('Unknown angle unit type {0}'.format(units))                              

        self._units = units
        
        self._phi = None
        self._psi = None
        self._omega = None
                                        
        self.phi = phi
        self.psi = psi
        self.omega = omega        

    def __repr__(self):
        return "<TorsionAngles: phi={0.phi}, psi={0.psi}, omega={0.omega}>".format(self)
    
    def __nonzero__(self):
        return self.__bool__()

    def __bool__(self):
        return  self.phi is not None \
                or self.psi is not None \
                or self.omega is not None        

    @property
    def units(self):
        """
        Current torsion angle units - a member of L{AngleUnits}
        @rtype: enum item
        """
        return self._units
    
    @property
    def phi(self):
        return self._phi
    @phi.setter
    def phi(self, phi):
        TorsionAngles.check_angle(phi, self._units)
        self._phi = phi   

    @property
    def psi(self):
        return self._psi
    @psi.setter
    def psi(self, psi):
        TorsionAngles.check_angle(psi, self._units)
        self._psi = psi     
        
    @property
    def omega(self):
        return self._omega
    @omega.setter
    def omega(self, omega):
        TorsionAngles.check_angle(omega, self._units)
        self._omega = omega        
        
    def copy(self):
        """
        @return: a deep copy of C{self}
        """
        return TorsionAngles(self.phi, self.psi, self.omega, self.units)
        
    def to_degrees(self):
        """
        Set angle measurement units to degrees.
        Convert the angles in this TorsionAngles instance to degrees.
        """
        
        if self._units != AngleUnits.Degrees:
    
            phi = TorsionAngles.deg(self._phi)
            psi = TorsionAngles.deg(self._psi)
            omega = TorsionAngles.deg(self._omega)
            
            # if no ValueError is raised by TorsionAngles.check_angle in TorsionAngles.deg:
            # (we assign directly to the instance variables to avoid check_angle being invoked again in setters)
            self._phi, self._psi, self._omega = phi, psi, omega
            self._units = AngleUnits.Degrees

        
    def to_radians(self):
        """
        Set angle measurement units to radians.
        Convert the angles in this TorsionAngles instance to radians.
        """

        if self._units != AngleUnits.Radians:        

            phi = TorsionAngles.rad(self._phi)
            psi = TorsionAngles.rad(self._psi)
            omega = TorsionAngles.rad(self._omega)
            
            # if no ValueError is raised by TorsionAngles.check_angle in TorsionAngles.rad:
            # (we assign directly to the instance variables to avoid check_angle being invoked again in setters)
            self._phi, self._psi, self._omega = phi, psi, omega
            self._units = AngleUnits.Radians

    @staticmethod
    def check_angle(angle, units):
        """
        Check the value of a torsion angle expressed in the specified units.
        """
        if angle is None:
            return
        elif units == AngleUnits.Degrees:
            if not (-180 <= angle <= 180):
                raise ValueError('Torsion angle {0} is out of range -180..180'.format(angle))           
        elif units == AngleUnits.Radians:
            if not (0 <= angle <= (2 * math.pi)):
                raise ValueError('Torsion angle {0} is out of range 0..2pi'.format(angle))
        else:
            raise ValueError('Unknown angle unit type {0}'.format(units))  
                    
    @staticmethod
    def rad(angle):
        """ 
        Convert a torsion angle value, expressed in degrees, to radians.
        Negative angles are converted to their positive counterparts: rad(ang + 360deg). 
        
        Return the calculated value in the range of [0, 2pi] radians. 
        """
        TorsionAngles.check_angle(angle, AngleUnits.Degrees)
                       
        if angle is not None:
            if angle < 0:
                angle += 360.            
            angle = math.radians(angle)
        return angle               
    
    @staticmethod
    def deg(angle):    
        """ 
        Convert a torsion angle value, expressed in radians, to degrees.
        Negative angles are not accepted, it is assumed that negative torsion angles have been 
        converted to their ang+2pi counterparts beforehand.  
        
        Return the calculated value in the range of [-180, +180] degrees. 
        """    
        TorsionAngles.check_angle(angle, AngleUnits.Radians)
         
        if angle is not None:                   
            if angle > math.pi:
                angle = -((2. * math.pi) - angle)
            angle = math.degrees(angle)
            
        return angle
