"""
3D and secondary structure APIs. 
"""

import os
import re
import copy
import math
import numpy
import datetime

import csb.io
import csb.pyutils
import csb.numeric
import csb.bio.utils

from abc import ABCMeta, abstractmethod, abstractproperty

from csb.bio.sequence import SequenceTypes, SequenceAlphabets, AlignmentTypes


class AngleUnits(csb.pyutils.enum):
    """
    Torsion angle unit types
    """
    Degrees='deg'; Radians='rad'
    
class SecStructures(csb.pyutils.enum):
    """
    Secondary structure types
    """
    Helix='H'; Strand='E'; Coil='C'; Turn='T'; Bend='S';
    Helix3='G'; PiHelix='I'; BetaBridge='B'; Gap='-'
    
class ChemElements(csb.pyutils.enum):
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

class DuplicateModelIDError(csb.pyutils.DuplicateKeyError):
    pass

class DuplicateChainIDError(csb.pyutils.DuplicateKeyError):
    pass

class DuplicateResidueIDError(csb.pyutils.DuplicateKeyError):
    pass

class DuplicateAtomIDError(csb.pyutils.DuplicateKeyError):
    pass

class AlignmentArgumentLengthError(ValueError):
    pass

class BrokenSecStructureError(ValueError):
    pass

class UnknownSecStructureError(BrokenSecStructureError):
    pass

class Abstract3DEntity(object):
    """
    Base class for all protein structure entities.
    
    This class defines uniform interface of all entities (e.g. L{Structure},
    L{Chain}, L{Residue}) according to the Composite pattern. 
    """
    
    __metaclass__ = ABCMeta

    @abstractproperty
    def items(self):
        """
        Return an iterator over all immediate children of the entity.
        """
        pass

    def components(self, klass=None):
        """
        Return an iterator over all descendants of the entity.
        
        @param klass: return entities the specified L{Abstract3DEntity} subclass
                      only. If None, traverse the hierarchy down to the lowest level.
        @param klass: class
        """
        for entity in CompositeEntityIterator.create(self, klass):
            if klass is None or isinstance(entity, klass):
                yield entity
        
    def apply_transformation(self, rotation, translation):
        """
        Apply in place RotationMatrix and translation Vector to all atoms.
        
        @type rotation: numpy array
        @type translation: numpy array 
        """
        for node in self.items:
            node.apply_transformation(rotation, translation)
    
    def list_coordinates(self, what=None, skip=False):
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
                except csb.pyutils.ItemNotFoundError:
                    if skip:
                        continue
                    raise Broken3DStructureError('Could not retrieve {0} atom from the structure'.format(atom_kind))
            
        return numpy.array(coords)
    
class CompositeEntityIterator(object):
    """
    Iterates over composite L{Abstract3DEntity} hierarchies.
    
    @param node: root entity to traverse
    @type node: L{Abstract3DEntity}
    """
    
    def __init__(self, node):
            
        if not isinstance(node, Abstract3DEntity):
            raise TypeError(node)
            
        self._node = node
        self._stack = csb.pyutils.Stack()
        
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
    Iterates over composite L{Abstract3DEntity} hierarchies, but terminates
    the traversal of a branch once a specific node type is encountered.
    
    @param node: root entity to traverse
    @type node: L{Abstract3DEntity}
    @param leaf: traverse the hierarchy down to the specified L{Abstract3DEntity}
    @type leaf: class
    """
    def __init__(self, node, leaf):
        
        if not issubclass(leaf, Abstract3DEntity):
            raise TypeError(leaf)
        
        self._leaf = leaf
        super(ConfinedEntityIterator, self).__init__(node)              
    
    def _inspect(self, node):
        
        if not isinstance(node, self._leaf):
            self._stack.push(node.items)
            
class Ensemble(csb.pyutils.AbstractNIContainer, Abstract3DEntity):
    """
    Represents an ensemble of multiple L{Structure} models.
    Provides a list-like access to these models::
    
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
        return self._models
    
    @property
    def items(self):
        return iter(self._models)
        
    @property
    def first_model(self):
        if len(self._models) > 0:
            return self[0]
        return None
    
    def to_pdb(self, output_file=None):
        """
        Dump the ensemble in PDB format.
        
        @param output_file: output file name or open stream
        @type output_file: str or stream
        """
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
        
class EnsembleModelsCollection(csb.pyutils.CollectionContainer):
    
    def __init__(self):
        
        super(EnsembleModelsCollection, self).__init__(type=Structure, start_index=1)
        self._models = set()
        
    def append(self, structure):
        
        if not structure.model_id or not str(structure.model_id).strip():
            raise ValueError("Invalid model identifier: '{0.model_id}'".format(structure))
        if structure.model_id in self._models:
            raise DuplicateModelIDError(structure.model_id) 
        else:
            return super(EnsembleModelsCollection, self).append(structure)

class Structure(csb.pyutils.AbstractNIContainer, Abstract3DEntity):
    """
    Represents a single model of a PDB 3-Dimensional molecular structure.
    Provides access to the L{Chain} objects, contained in the model::
    
        >>> structure['A']
        <Chain A: Protein>
        >>> structure.chains['A']
        <Chain A: Protein>
        >>> structure.items[0]
        <Chain A: Protein> #(given a chain order A, B...)    
    
    @param accession: accession number of the structure
    @type accession: str
    """
    def __init__(self, accession):
                        
        self._accession = None
        self._chains = StructureChainsTable(self)
        self.model_id = None
        
        self.accession = accession

    def __repr__(self):
        return "<Structure Model {0.model_id}: {0.accession}, {1} chains>".format(self, self.chains.length)

    @property
    def _children(self):
        return self._chains
    
    @property
    def chains(self):
        return self._chains
    
    @property
    def items(self):
        for chain in self._chains:
            yield self._chains[chain]
                
    @property
    def first_chain(self):
        if len(self._chains) > 0:
            return next(self.items)
        return None
        
    @property
    def accession(self):
        return self._accession
    @accession.setter
    def accession(self, accession):
        if accession is None:
            raise ValueError(accession)
        self._accession = str(accession).strip().lower()
        for c in self.chains:
            self.chains[c]._accession = self._accession
            
    def define_molecule(self, molecule_id, chains):
        """
        Group C{chains}, defined by their IDs, to a polymer with id = C{molecule_id}.
        
        @param molecule_id: ID of the new polymer
        @param chains: a list of chain IDs
        @type chains: list
        
        @raise ValueError: if an invalid chain ID is provided
        """
        assert molecule_id is not None
        cs = [ ]
        
        for c in chains:
            if c not in self.chains:
                raise ValueError('No such chain {0} in structure'.format(c))
            cs.append(self.chains[c])
        for chain in cs:
            chain.molecule_id = molecule_id
                    
    def to_fasta(self):
        """
        Dump the structure in FASTA format. 
        
        @return: FASTA-formatted string with all chains in the structure
        @rtype: str
        
        @deprecated: this method will be removed soon. Use
                     L{csb.bio.sequence.ChainSequence.create} instead
        """
        fasta = []
        
        for chain_id in self.chains:

            chain = self.chains[chain_id]
            if chain.length > 0:
                fasta.append(chain.header)
                fasta.append(chain.sequence)
        
        return os.linesep.join(fasta)

    def to_pdb(self, output_file=None):
        """
        Dump the whole structure in PDB format.
        
        @param output_file: output file name or open stream
        @type output_file: str or stream
        """
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

class StructureChainsTable(csb.pyutils.DictionaryContainer):
    
    def __init__(self, structure=None, chains=None):
        self.__container = structure
        super(StructureChainsTable, self).__init__()
        
        if chains is not None:
            for chain in chains:
                self.append(chain)
        
    def __repr__(self):
        if len(self) > 0:
            return "<StructureChains: {0}>".format(', '.join(self.keys()))
        else:
            return "<StructureChains: empty>"     
    
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
        
class Chain(csb.pyutils.AbstractNIContainer, Abstract3DEntity):
    """
    Represents a polymeric chain. Provides list-like and rank-based access to
    the residues in the chain::
    
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
    @type type: L{csb.pyutils.EnumItem}
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
        self.molecule_id = molecule_id
        self._torsion_computed = False
        self.name = str(name).strip()
        
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
        return self._id
    @id.setter
    def id(self, id):
        if not isinstance(id, csb.pyutils.string):
            raise ValueError(id)
        id = id.strip()
        if self._structure:
            self._structure.chains._update_chain_id(self, id)
        self._id = id
    
    @property
    def accession(self):
        return self._accession
    @accession.setter
    def accession(self, accession):
        if self._structure:
            raise InvalidOperation("Only structure's accession can be altered in this case")
        if accession is None:
            raise ValueError(accession)
        self._accession = str(accession).strip()
        
    @property
    def type(self):
        return self._type
    @type.setter
    def type(self, type):
        if type.enum is not SequenceTypes:
            raise TypeError(type)
        self._type = type

    @property
    def residues(self):
        return self._residues
    
    @property
    def items(self):
        return iter(self._residues)
    
    @property
    def torsion(self):
        if not self._torsion_computed:
            raise InvalidOperation('The correctness of the data is not guaranteed until chain.compute_torsion() is invoked.')
            
        torsion = TorsionAnglesCollection()
        
        for r in self.residues:
            if r.torsion is None:
                torsion.append(TorsionAngles(None, None, None))
            else:
                torsion.append(r.torsion)
                
        return torsion
    
    @property
    def has_torsion(self):
        return self._torsion_computed

    @property
    def length(self):
        return self._residues.length
    
    @property
    def entry_id(self):
        if self._accession and self._id:
            return self._accession + self._id
        else:
            return None
    
    @property
    def header(self):
        return ">{0._accession}_{0._id} mol:{1} length:{0.length} {0.name}".format(self, str(self.type).lower())

    @property
    def sequence(self):        
        sequence = []
        for residue in self._residues:
            if residue and residue.type:
                sequence.append(str(residue.type))
            else:
                sequence.append('-')
        return ''.join(sequence)
    
    @property
    def alphabet(self):
        if self._type == SequenceTypes.Protein:                                                 
            return SequenceAlphabets.Protein
        elif self._type in (SequenceTypes.NucleicAcid, SequenceTypes.DNA, SequenceTypes.RNA):   
            return SequenceAlphabets.Nucleic            
        else:
            raise NotImplementedError()
        
    @property
    def secondary_structure(self):
        return self._secondary_structure
    @secondary_structure.setter
    def secondary_structure(self, ss):
        if not isinstance(ss, SecondaryStructure):
            raise TypeError(ss)
        if len(ss) > 0:
            if (ss[ss.last_index].end > self._residues.last_index):
                raise ValueError('The secondary structure is out of range')
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
    
    def format_residue(self, residue):
        """
        Return a PDB-friendly string representation of the residue. Safer than repr(residue).
        
        @param residue: the residue name to format (a member of any L{SequenceAlphabets})
        @type residue: L{csb.pyutils.EnumItem}
        
        @return: PDB-friendly string representation of the residue type
        @rtype: str
        """
        if self.type == SequenceTypes.NucleicAcid:                                  
            return str(residue.type)
        else:
            return repr(residue.type)        
        
    def find(self, sequence_number, insertion_code=None):
        """
        Get a residue by its original Residue Sequence Number and Insertion Code.
        
        @param sequence_number: PDB sequence number of the residue
        @type sequence_number: str
        @param insertion_code: PDB insertion code of the residue (if any)
        @type insertion_code: str
        
        @return: the residue object with such an ID
        @rtype: L{Residue}
        
        @raise csb.pyutils.ItemNotFoundError: if no residue with that ID exists
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
        
            >>> other.apply_transformation(rotation_matrix, translation_vector)
            
        will result in C{other}'s coordinates superimposed over C{self}.
        
        @param other: the subject (movable) chain
        @type other: L{Chain}
        @param what: a list of atom kinds, e.g. ['CA']
        @type what: list
        @param how: fitting method (global or local) - a member of the L{AlignmentTypes} enum
        @type how: L{csb.pyutils.EnumItem}
        
        @return: superimposition info object, containing rotation matrix, translation 
                 vector and computed RMSD
        @rtype: L{SuperimposeInfo}
        
        @raise AlignmentArgumentLengthError: when the lengths of the argument chains differ 
        """ 
        if self.length != other.length or self.length < 1:
            raise AlignmentArgumentLengthError('Both chains must be of the same and positive length')
        
        x = self.list_coordinates(what)
        y = other.list_coordinates(what) 
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
            other.apply_transformation(R, t)
            
        @param other: the subject (movable) chain
        @type other: L{Chain}
        @param what: a list of atom kinds, e.g. ['CA']
        @type what: list
        @param how: fitting method (global or local) - a member of the L{AlignmentTypes} enum
        @type how: L{csb.pyutils.EnumItem}
        
        @return: superimposition info object, containing rotation matrix, translation 
                 vector and computed RMSD
        @rtype: L{SuperimposeInfo}        
        """
        result = self.superimpose(other, what=what, how=how)
        other.apply_transformation(result.rotation, result.translation)
        
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
        
        x = self.list_coordinates(what)
        y = other.list_coordinates(what)
        assert len(x) == len(y)

        return csb.bio.utils.rmsd(x, y) 
    
    def tm_superimpose(self, other, what=['CA'], how=AlignmentTypes.Global):                    
        """
        Find the optimal fit between C{self} and C{other}. Return L{SuperimposeInfo}
        (RotationMatrix, translation Vector and TM-score), such that:
        
            >>> other.apply_transformation(rotation_matrix, translation_vector)
            
        will result in C{other}'s coordinates superimposed over C{self}.
        
        @param other: the subject (movable) chain
        @type other: L{Chain}
        @param what: a list of atom kinds, e.g. ['CA']
        @type what: list
        @param how: fitting method (global or local) - a member of the L{AlignmentTypes} enum
        @type how: L{csb.pyutils.EnumItem}
        
        @return: superimposition info object, containing rotation matrix, translation 
                 vector and computed TM-score
        @rtype: L{SuperimposeInfo}
        
        @raise AlignmentArgumentLengthError: when the lengths of the argument chains differ         
        """
        
        if self.length != other.length or self.length < 1:
            raise ValueError('Both chains must be of the same and positive length')
        
        x = self.list_coordinates(what)
        y = other.list_coordinates(what)
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
        
        x = self.list_coordinates(what)
        y = other.list_coordinates(what)
        assert len(x) == len(y)

        return csb.bio.utils.tm_score(x, y)             

class ChainResiduesCollection(csb.pyutils.CollectionContainer):
    
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
            raise csb.pyutils.ItemNotFoundError(id)
        
class Residue(csb.pyutils.AbstractNIContainer, Abstract3DEntity):
    """
    Base class representing a single residue. Provides a dictionary-like
    access to the atoms contained in the residue::
    
        >>> residue['CA']
        <Atom [3048]: CA>
        >>> residue.atoms['CA']
        <Atom [3048]: CA>
        >>> residue.items[1]
        <Atom [3047]: CA>
    
    @param rank: rank of the residue with respect to the chain
    @type rank: int
    @param type: residue type - a member of any L{SequenceAlphabets}
    @type type: L{csb.pyutils.EnumItem}
    @param sequence_number: PDB sequence number of the residue
    @type sequence_number: str
    @param insertion_code: PDB insertion code, if any
    @type insertion_code: str
    """            
    def __init__(self, rank, type, sequence_number=None, insertion_code=None):
        
        self._type = None    
        self._pdb_name = None
        self._rank = int(rank)
        self._atoms = ResidueAtomsTable(self) 
        self._secondary_structure = None
        self._torsion = None
        self._sequence_number = None
        self._insertion_code = None
        self._container = None
        
        self.type = type
        self.id = sequence_number, insertion_code
        self._pdb_name = repr(type)
        
    @property
    def _children(self):
        return self._atoms
        
    def __repr__(self):
        return '<{1} [{0.rank}]: {0.type!r} {0.id}>'.format(self, self.__class__.__name__)
        
    @property
    def type(self):
        return self._type
    @type.setter
    def type(self, type):
        if type.enum not in (SequenceAlphabets.Protein, SequenceAlphabets.Nucleic, SequenceAlphabets.Unknown):
            raise TypeError(type)
        self._type = type
        
    @property
    def rank(self):
        return self._rank
    
    @property
    def secondary_structure(self):
        return self._secondary_structure
    @secondary_structure.setter
    def secondary_structure(self, structure):
        if not isinstance(structure, SecondaryStructureElement):
            raise TypeError(structure)
        self._secondary_structure = structure
        
    @property
    def torsion(self):
        return self._torsion
    @torsion.setter
    def torsion(self, torsion):
        if not isinstance(torsion, TorsionAngles):
            raise TypeError(torsion)
        self._torsion = torsion
    
    @property
    def atoms(self):
        return self._atoms
    
    @property
    def items(self):
        for atom in self._atoms:
            yield self._atoms[atom]        

    @property
    def sequence_number(self):
        return self._sequence_number

    @property
    def insertion_code(self):
        return self._insertion_code
    
    @property
    def id(self):
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
        return len(self.atoms) > 0
        
    def list_coordinates(self, what=None, skip=False):
        
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
    def create(sequence_type, *arguments, **keyword_arguments):
        """
        Residue factory method, which returns the proper L{Residue} instance based on 
        the specified C{sequence_type}. All additional arguments are used to initialize
        the subclass by passing them automatically to the underlying constructor. 
        
        @param sequence_type: create a Residue of that SequenceType 
        @type sequence_type: L{csb.pyutils.EnumItem}
        
        @return: a new residue of the proper subclass
        @rtype: L{Residue} subclass
        
        @raise ValueError: if the sequence type is not known
        """        
        if sequence_type == SequenceTypes.Protein:                                      
            return ProteinResidue(*arguments, **keyword_arguments)
        elif sequence_type == SequenceTypes.NucleicAcid:                                
            return NucleicResidue(*arguments, **keyword_arguments)
        elif sequence_type == SequenceTypes.Unknown:
            return UnknownResidue(*arguments, **keyword_arguments)
        else:
            raise ValueError(sequence_type)        
           
class ProteinResidue(Residue):
    """
    Represents a single amino acid residue.
    
    @param rank: rank of the residue with respect to the chain
    @type rank: int
    @param type: residue type - a member of 
                 L{csb.bio.sequence.SequenceAlphabets.Protein}
    @type type: L{csb.pyutils.EnumItem}
    @param sequence_number: PDB sequence number of the residue
    @type sequence_number: str
    @param insertion_code: PDB insertion code, if any
    @type insertion_code: str    
    """
    
    def __init__(self, rank, type, sequence_number=None, insertion_code=None):
          
        if isinstance(type, csb.pyutils.string):
            try:
                if len(type) == 3:
                    type = csb.pyutils.Enum.parsename(SequenceAlphabets.Protein, type)
                else:    
                    type = csb.pyutils.Enum.parse(SequenceAlphabets.Protein, type)          
            except (csb.pyutils.EnumMemberError, csb.pyutils.EnumValueError):
                raise ValueError("What is '{0}'?".format(type))
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
        except csb.pyutils.ItemNotFoundError as missing_atom:
            if strict:
                raise Broken3DStructureError('Could not retrieve {0} atom from the current residue {1!r}.'.format(missing_atom, self))
            else:
                return angles
        
        try:
            if prev_residue is not None and prev_residue.has_structure:
                prev_c = prev_residue._atoms['C'].vector
                angles.phi = csb.numeric.dihedral_angle(prev_c, n, ca, c)
        except csb.pyutils.ItemNotFoundError as missing_prevatom:
            if strict:
                raise Broken3DStructureError('Could not retrieve {0} atom from the i-1 residue {1!r}.'.format(missing_prevatom, prev_residue))    
        try:
            if next_residue is not None and next_residue.has_structure:    
                next_n = next_residue._atoms['N'].vector
                angles.psi = csb.numeric.dihedral_angle(n, ca, c, next_n)
                next_ca = next_residue._atoms['CA'].vector
                angles.omega = csb.numeric.dihedral_angle(ca, c, next_n, next_ca)
        except csb.pyutils.ItemNotFoundError as missing_nextatom:
            if strict:
                raise Broken3DStructureError('Could not retrieve {0} atom from the i+1 residue {1!r}.'.format(missing_nextatom, next_residue))              
                                
        return angles

class NucleicResidue(Residue):
    """
    Represents a single nucleotide residue.
    
    @param rank: rank of the residue with respect to the chain
    @type rank: int
    @param type: residue type - a member of 
                 L{csb.bio.sequence.SequenceAlphabets.Nucleic}
    @type type: L{csb.pyutils.EnumItem}
    @param sequence_number: PDB sequence number of the residue
    @type sequence_number: str
    @param insertion_code: PDB insertion code, if any
    @type insertion_code: str        
    """
    
    def __init__(self, rank, type, sequence_number=None, insertion_code=None):
        
        if isinstance(type, csb.pyutils.string):
            try:
                if len(type) > 1:
                    type = csb.pyutils.Enum.parsename(SequenceAlphabets.Nucleic, type)
                else:    
                    type = csb.pyutils.Enum.parse(SequenceAlphabets.Nucleic, type)
            except (csb.pyutils.EnumMemberError, csb.pyutils.EnumValueError):
                raise ValueError("What is '{0}'?".format(type))
        elif type.enum is not SequenceAlphabets.Nucleic:
            raise TypeError(type)
            
        super(NucleicResidue, self).__init__(rank, type, sequence_number, insertion_code)  
        self._pdb_name = str(type)

class UnknownResidue(Residue):
    
    def __init__(self, rank, type, sequence_number=None, insertion_code=None):
        super(UnknownResidue, self).__init__(rank, SequenceAlphabets.Unknown.UNK, sequence_number, insertion_code)
            
class ResidueAtomsTable(csb.pyutils.DictionaryContainer):
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
    
    def append(self, atom):
        """
        Append a new Atom to the catalog.
        
        If the atom has an alternate position, a disordered proxy will be created instead and the 
        atom will be appended to the L{DisorderedAtom}'s list of children. If a disordered atom 
        with that name already exists, the atom will be appended to its children only.
        If an atom with the same name exists, but it was erroneously not marked as disordered,
        that terrible condition will be fixed too :-(
        
        @param atom: the new atom to append
        @type atom: L{Atom}
        
        @raise DuplicateAtomIDError: if an atom with the same sequence number and 
                                     insertion code already exists in that residue
        """
        if atom._residue and atom._residue is not self.__residue:
            raise InvalidOperation('This atom is part of another residue')
        if atom.alternate or (atom.name in self and isinstance(self[atom.name], DisorderedAtom)):
            if atom.name not in self:
                atom.residue = self.__residue
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
                self[atom.name].append(atom)          
        else:
            if atom.name in self:
                raise DuplicateAtomIDError('Atom {0} is already defined for {1}'.format(atom.name, self.__residue))
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
            raise ValueError('Atom\'s name differs from the specified update key.')
        if atom.residue is not self.__residue:
            atom._residue = self.__residue
        
        super(ResidueAtomsTable, self)._update({atom_name: atom})  
    
class Atom(Abstract3DEntity):
    """
    Represents a single atom in space.
    
    @param serial_number: atom's UID
    @type serial_number: int
    @param name: atom's name
    @type name: str
    @param element: corresponding L{ChemElements}
    @type element: L{csb.pyutils.EnumItem}
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
        self.alternate = False        
        self.temperature = None
        self.occupancy = None
        self.charge = None

        if not isinstance(name, csb.pyutils.string):
            raise TypeError(name)
        name_compact = name.strip()
        if len(name_compact) < 1:
            raise ValueError(name)
        self._name = name_compact
        self._full_name = name
            
        if isinstance(element, csb.pyutils.string):
            element = csb.pyutils.Enum.parsename(ChemElements, element)
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
    
    def apply_transformation(self, rotation, translation):
        
        vector = numpy.dot(self.vector, numpy.transpose(rotation)) + translation
        self.vector = vector
    
    def list_coordinates(self, what=None, skip=False):
        
        if what is None:
            what = [self.name]
            
        if self.name in what:
            return numpy.array([self.vector.copy()])
        elif skip:
            return numpy.array([])
        else:
            raise Missing3DStructureError()        

    @property
    def serial_number(self):
        return self._serial_number
    @serial_number.setter
    def serial_number(self, number):
        if not isinstance(number, int) or number < 1:
            raise TypeError(number)
        self._serial_number = number
    
    @property
    def name(self):
        return self._name
        
    @property
    def element(self):
        return self._element
    
    @property
    def residue(self):
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
        return self._vector
    @vector.setter
    def vector(self, vector):
        if numpy.shape(vector) != (3,):
            raise ValueError("Three dimensional vector expected")
        self._vector = numpy.array(vector)
    
    @property
    def items(self):
        return iter([])
        
class DisorderedAtom(csb.pyutils.CollectionContainer, Atom):
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
        self._rep = None
        self.append(atom)
            
    def append(self, atom):
        """
        Append a new atom to the collection of alternatives.
        
        @param atom: the new alternative
        @type atom: L{Atom}
        """
        self.__update_rep(atom)
        super(DisorderedAtom, self).append(atom)
    
    def apply_transformation(self, rotation, translation):
        
        for atom in self:
            atom.apply_transformation(rotation, translation)
        self._vector = self._rep._vector
        
    def __update_rep(self, atom):
        
        if self._rep is None or \
        ((self._rep.occupancy != atom.occupancy) and (self._rep.occupancy < atom.occupancy)):
        
            self._rep = atom      

            self._serial_number = self._rep.serial_number
            self._name = self._rep.name
            self._full_name = self._rep._full_name            
            self._element = self._rep.element
            self.alternate = self._rep.alternate
            self._residue = self._rep.residue
            self._vector = self._rep.vector
            self.temperature = self._rep.temperature
            self.occupancy = self._rep.occupancy
            self.charge = self._rep.charge
            
    def __repr__(self):
        return "<DisorderedAtom: {0.length} alternative locations>".format(self)
        
class FileBuilder(object):
    """
    Base abstract files for all structure file formatters.
    Defines a common step-wise interface according to the Builder pattern.
    
    @param output: output stream (this is where the product is constructed)
    @type param: stream
    """
    
    __metaclass__ = ABCMeta

    def __init__(self, output):

        if not hasattr(output, 'write'):
            raise TypeError(output)
        
        def isnull(this, that, null=None):
            if this is null:
                return that
            else:
                return this        

        self._out = output
        self._isnull = isnull
        
    @property
    def output(self):
        return self._out
    
    @property
    def isnull(self):
        return self._isnull
    
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
    def add_header(self, master_structure):
        pass
    
    @abstractmethod
    def add_structure(self, structure):
        pass
    
    def finalize(self):
        pass

class PDBFileBuilder(FileBuilder):
    """
    PDB file format builder.
    """
        
    def writeline(self, text):
        self.write('{0:80}\n'.format(text))
        
    def add_header(self, master):
        """
        Write the HEADER of the file using C{master}
        
        @type master: L{Structure}
        """

        isnull = self.isnull
        
        header = 'HEADER    {0:40}{1:%d-%b-%y}   {2:4}'
        self.writeline(header.format('.', datetime.datetime.now(), master.accession.upper()))
        
        molecules = { }
        
        for chain_id in master.chains:
            chain = master.chains[chain_id]
            if chain.molecule_id not in molecules:
                molecules[chain.molecule_id] = [ ]
            molecules[chain.molecule_id].append(chain_id)
        
        k = 0
        for mol_id in sorted(molecules):
            
            chains = molecules[mol_id]
            first_chain = master.chains[ chains[0] ]            
            
            self.writeline('COMPND {0:3} MOL_ID: {1};'.format(k + 1, isnull(mol_id, '0')))
            self.writeline('COMPND {0:3} MOLECULE: {1};'.format(k + 2, isnull(first_chain.name, '')))
            self.writeline('COMPND {0:3} CHAIN: {1};'.format(k + 3, ', '.join(chains)))
            k += 3
            
        for chain_id in master.chains:
            
            chain = master.chains[chain_id]
            res = [ r._pdb_name for r in chain.residues ]

            rn = 0
            for j in range(0, chain.length, 13):
                rn += 1
                residues = [ '{0:>3}'.format(r) for r in res[j : j + 13] ]
                self.writeline('SEQRES {0:>3} {1} {2:>4}  {3}'.format(
                                            rn, chain.id, chain.length, ' '.join(residues) ))
                
    def add_structure(self, structure):
        """
        Append a new model to the file
        
        @type structure: L{Structure}
        """

        isnull = self.isnull
         
        for chain_id in structure.chains:
        
            chain = structure.chains[chain_id]
            for residue in chain.residues:
        
                atoms = [ ]
                for an in residue.atoms:
                    atom = residue.atoms[an]
                    if type(atom) is DisorderedAtom:
                        for dis_atom in atom: atoms.append(dis_atom)
                    else:
                        atoms.append(atom)
                atoms.sort()
                
                for atom in atoms:

                    alt = atom.alternate
                    if alt is True:
                        alt = 'A'
                    elif alt is False:
                        alt = ' '
                    
                    if atom.element:
                        element = repr(atom.element)
                    else:
                        element = ' '
                    self.writeline('ATOM  {0:>5} {1:>4}{2}{3:>3} {4}{5:>4}{6}   {7:>8.3f}{8:>8.3f}{9:>8.3f}{10:>6.2f}{11:>6.2f}{12:>12}{13:2}'.format(
                                        atom.serial_number, atom._full_name, isnull(alt, ' '), 
                                        residue._pdb_name, chain.id, 
                                        isnull(residue.sequence_number, residue.rank), isnull(residue.insertion_code, ' '), 
                                        atom.vector[0], atom.vector[1], atom.vector[2], isnull(atom.occupancy, 0.0), isnull(atom.temperature, 0.0), 
                                        element, isnull(atom.charge, ' ') ))        

            self.writeline('TER')
        
    def finalize(self):
        """
        Add the END marker
        """
        self.writeline('END')
        self._out.flush()     

class PDBEnsembleFileBuilder(PDBFileBuilder):
    """
    Supports serialization of NMR ensembles.
    
    Functions as a simple decorator, which wraps C{add_structure} with
    MODEL/ENDMDL records.
    """

    def add_structure(self, structure):
        
        model_id = self.isnull(structure.model_id, 1)
        
        self.writeline('MODEL     {0:>4}'.format(model_id))       
        super(PDBEnsembleFileBuilder, self).add_structure(structure)
        self.writeline('ENDMDL')
            
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
    @type type: csb.pyutils.EnumItem
    @param score: secondary structure prediction confidence, if available
    @type score: int
    
    @raise IndexError: if start/end positions are out of range
    """    
    def __init__(self, start, end, type, score=None):
        
        if not (0 < start <= end):
            raise IndexError('Element coordinates are out of range: 0 < start <= end.')
        if isinstance(type, csb.pyutils.string):
            type = csb.pyutils.Enum.parse(SecStructures, type)
        if type.enum is not SecStructures:
            raise TypeError(type)
                
        self.start = start
        self.end = end
        self.type = type
        self._score = None
        
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
    def length(self):
        return self.end - self.start + 1
    
    @property
    def score(self):
        return self._score
    @score.setter
    def score(self, scores):
        if not len(scores) == self.length:
            raise ValueError('There must be a score entry for each residue in the element.')        
        self._score = csb.pyutils.CollectionContainer(items=list(scores), type=int, start_index=self.start)
    
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

class SecondaryStructure(csb.pyutils.CollectionContainer):
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
        Add a new SecondaryStructureElement. Then sort all added elements by
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
        if not isinstance(string, csb.pyutils.string):
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
                    type = csb.pyutils.Enum.parse(SecStructures, currel)
                except csb.pyutils.EnumValueError:
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
        return self._minstart
        
    @property
    def end(self):
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
        return self.scan(start=rank, end=rank, filter=type, loose=True, cut=True)
    
    def scan(self, start, end, filter=None, loose=True, cut=True):
        """
        Get all secondary structure elements within the specified start..end region.
        
        @param start: the start position of the region, 1-based, inclusive
        @type start: int
        @param end: the end position of the region, 1-based, inclusive
        @type end: int     
        @param filter: return only elements of the specified L{SecStructures} kind
        @type filter: L{csb.pyutils.EnumItem}
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
        @rtype: list
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
        return matches
    
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
        
class TorsionAnglesCollection(csb.pyutils.CollectionContainer):
    """
    Describes a collection of torsion angles. Provides 1-based list-like access.
    
    @param items: an initialization list of L{TorsionAngles}
    @type items: list
    """  
    def __init__(self, items=None, start=1):
        super(TorsionAnglesCollection, self).__init__(items, type=TorsionAngles, start_index=start)

    def __repr__(self):
        if len(self) > 0:
            return "<TorsionAnglesList: {0} ... {1}>".format(self[self.start_index], self[self.last_index])
        else:
            return "<TorsionAnglesList: empty>"        
        
    @property
    def phi(self):       
        return [a.phi for a in self]

    @property
    def psi(self):       
        return [a.psi for a in self]            

    @property
    def omega(self):       
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
    @type units: L{csb.pyutils.EnumItem}
    
    @raise ValueError: on invalid/unknown units
    """
                     
    def __init__(self, phi, psi, omega, units=AngleUnits.Degrees):
        
        try:
            if isinstance(units, csb.pyutils.string):
                units = csb.pyutils.Enum.parse(AngleUnits, units, ignore_case=True)
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
        return  self.phi is not None \
                or self.psi is not None \
                or self.omega is not None

    @property
    def units(self):
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
