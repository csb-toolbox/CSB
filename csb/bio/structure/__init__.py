import os
import re
import copy
import math
import numpy
import csb.pyutils
import csb.bio.utils

from csb.bio.sequence import SequenceTypes, SequenceAlphabets, AlignmentTypes
from itertools import izip

AngleUnits = csb.pyutils.enum( Degrees='deg', Radians='rad' )
"""
@var AngleUnits: Torsion angle unit types
"""

SecStructures = csb.pyutils.enum( Helix='H', Strand='E', Coil='C', Turn='T',
                                  Bend='S', Helix3='G', PiHelix='I', BetaBridge='B', Gap='-' )
"""
@var SecStructures: Secondary structure types
"""

ChemElements = csb.pyutils.enum( H=1, He=2, Li=3, Be=4, B=5, C=6, N=7, O=8, F=9, Ne=10, Na=11, Mg=12, Al=13, Si=14, P=15, 
                                 S=16, Cl=17, Ar=18, K=19, Ca=20, Sc=21, Ti=22, V=23, Cr=24, Mn=25, Fe=26, Co=27, Ni=28, 
                                 Cu=29, Zn=30, Ga=31, Ge=32, As=33, Se=34, Br=35, Kr=36, Rb=37, Sr=38, Y=39, Zr=40, Nb=41, 
                                 Mo=42, Tc=43, Ru=44, Rh=45, Pd=46, Ag=47, Cd=48, In=49, Sn=50, Sb=51, Te=52, I=53, Xe=54,
                                 Cs=55, Ba=56, Hf=72, Ta=73, W=74, Re=75, Os=76, Ir=77, Pt=78, Au=79, Hg=80, Tl=81, Pb=82, 
                                 Bi=83, Po=84, At=85, Rn=86, Fr=87, Ra=88, Rf=104, Db=105, Sg=106, Bh=107, Hs=108, Mt=109, 
                                 Ds=110, Rg=111, La=57, Ce=58, Pr=59, Nd=60, Pm=61, Sm=62, Eu=63, Gd=64, Tb=65, Dy=66, 
                                 Ho=67, Er=68, Tm=69, Yb=70, Lu=71, Ac=89, Th=90, Pa=91, U=92, Np=93, Pu=94, Am=95, Cm=96, 
                                 Bk=97, Cf=98, Es=99, Fm=100, Md=101, No=102, Lr=103,
                                 x=-1 ) 
"""
@var ChemElements: Periodic table elements
"""

class Broken3DStructureError(ValueError):
    pass
class Missing3DStructureError(Broken3DStructureError):
    pass       
class InvalidOperation(Exception):
    pass
class DuplicateChainIDError(csb.pyutils.DuplicateKeyError):
    pass
class DuplicateResidueIDError(csb.pyutils.DuplicateKeyError):
    pass
class DuplicateAtomIDError(csb.pyutils.DuplicateKeyError):
    pass
class AlignmentArgumentLengthError(ValueError):
    pass

class Structure(object):
    """
    Represents a PDB 3-Dimensional molecular structure.
    
    @param accession: accession number of the structure
    @type accession: str
    """
    def __init__(self, accession):
                        
        self._accession = None
        self._chains = StructureChainsTable(self)
        self.model_id = None
        
        self.accession = accession

    def __repr__(self):
        return "<Structure: {0.accession}, {1} chains>".format(self, self.chains.length)

    def __getitem__(self,key):
        return self._chains[key]

    def __iter__(self):
        return iter(self._chains)
    
    @property
    def chains(self):
        return self._chains
    
    @property
    def first_chain(self):
        chains = self.chains.keys()
        if len(chains) > 0:
            return self.chains[chains[0]]
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

    def apply_transformation(self, rotation, translation):
        """
        Apply in place RotationMatrix and translation Vector to all atoms in the structure.
        
        @type rotation: L{RotationMatrix}
        @type translation: L{Vector} 
        """
        for chain_id in self.chains:
            self.chains[chain_id].apply_transformation(rotation, translation)
                    
    def to_fasta(self):
        """
        Dump the structure in FASTA format. 
        
        @return: FASTA-formatted string with all chains in the structure
        @rtype: str
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
        
        @param output_file: output file name
        @type output_file: str 
        """
        from datetime import datetime
        import StringIO
        
        def isnull(this, that, null=None):
            if this is null: return that
            return this
            
        class MyStringIO(StringIO.StringIO):
            def writeline(self, data):
                self.write('{0:80}\n'.format(data))
        
        stream = MyStringIO()
        
        stream.writeline('HEADER    {0:40}{1:%d-%b-%y}   {2:4}'.format('.', datetime.now(), self.accession.upper()))
        
        molecules = { }
        for chain_id in self.chains:
            chain = self.chains[chain_id]
            if chain.molecule_id not in molecules:
                molecules[chain.molecule_id] = [ ]
            molecules[chain.molecule_id].append(chain_id)
        
        k = 0
        for mol_id in sorted(molecules):
            
            chains = molecules[mol_id]
            first_chain = self.chains[ chains[0] ]            
            
            stream.writeline('COMPND {0:3} MOL_ID: {1};'.format(k + 1, isnull(mol_id, '0')))
            stream.writeline('COMPND {0:3} MOLECULE: {1};'.format(k + 2, isnull(first_chain.name, '')))
            stream.writeline('COMPND {0:3} CHAIN: {1};'.format(k + 3, ', '.join(chains)))
            k += 3
            
        for chain_id in self.chains:
            
            chain = self.chains[chain_id]
            res = [ chain.format_residue(r) for r in chain.residues ]

            rn = 0
            for j in range(0, chain.length, 13):
                rn += 1
                residues = [ '{0:>3}'.format(r) for r in res[j : j + 13] ]
                stream.writeline('SEQRES {0:>3} {1} {2:>4}  {3}'.format(
                                            rn, chain.id, chain.length, ' '.join(residues) ))
        
        for chain_id in self.chains:
        
            chain = self.chains[chain_id]
            for residue in chain.residues:
        
                atoms = [ ]
                for an in residue.structure:
                    atom = residue.structure[an]
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
                    stream.writeline('ATOM  {0:>5} {1:>4}{2}{3:3} {4}{5:>4}{6}   {7:>8.3f}{8:>8.3f}{9:>8.3f}{10:>6.2f}{11:>6.2f}{12:>12}{13:2}'.format(
                                        atom.serial_number, atom._full_name, isnull(alt, ' '), 
                                        chain.format_residue(residue), chain.id, 
                                        isnull(residue.sequence_number, residue.rank), isnull(residue.insertion_code, ' '), 
                                        atom.vector.x, atom.vector.y, atom.vector.z, isnull(atom.occupancy, 0.0), isnull(atom.temperature, 0.0), 
                                        element, isnull(atom.charge, ' ') ))
            stream.writeline('TER')
        stream.writeline('END')
        
        data = stream.getvalue()        
        stream.close()
        
        if not output_file:
            return data
        else:
            with open(output_file, 'w') as out:
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
                
    def update(self):
        raise NotImplementedError()
    def set(self):
        raise NotImplementedError()
        
class Chain(object):
    """
    Represents a polymeric chain. Provides indexed access to the residues in the chain.
    
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
    def __init__(self, chain_id, type=SequenceTypes.Protein, name='', residues=None, accession=None, molecule_id=None):

        self._id = str(chain_id).strip()
        self._accession = None
        self._type = None
        self._residues = ChainResiduesCollection(self, residues)
        self.molecule_id = molecule_id
        self._torsion_computed = False
        self.name = str(name).strip()
        
        self._structure = None
        
        self.type = type
        if accession is not None:
            self.accession = accession

    def __repr__(self):
        return "<Chain {0.id}: {0.type!r}>".format(self)        

    def __len__(self):
        return self._residues.length

    def __getitem__(self,key):
        return self._residues[key]

    def __iter__(self):
        return iter(self._residues)

    @property
    def id(self):
        return self._id
    @id.setter
    def id(self, id):
        if not isinstance(id, basestring):
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
        if not isinstance(type, csb.pyutils.EnumItem):
            raise TypeError(type)
        self._type = type

    @property
    def residues(self):
        return self._residues
    
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
        elif self._type == SequenceTypes.NucleicAcid:
            return SequenceAlphabets.Nucleic            
        else:
            raise NotImplementedError()
        
    def clone(self):
        """
        Make a deep copy of the chain. If this chain is part of a structure, 
        detach from it.
        
        @return: a deep copy of self
        @rtype: L{Chain}
        """
        
        new_chain = csb.pyutils.deepcopy(self) # self.subregion(self.residues.start_index, self.residues.last_index, clone=True)
        new_chain.type = self.type
        new_chain._structure = None
        
        return new_chain       
        
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
            residues = csb.pyutils.deepcopy(residues)
        
        chain = Chain(self.id, accession=self.accession, name=self.name, 
                      type=self.type, residues=residues, molecule_id=self.molecule_id)
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
            
    def apply_transformation(self, rotation, translation):
        """
        Apply in place RotationMatrix and translation Vector to all atoms in the chain.
        
        @type rotation: L{RotationMatrix}
        @type translation: L{Vector}         
        """
        for residue in self.residues:
            for atom_name in residue.structure:
                residue.structure[atom_name]._transform_vector(rotation, translation)
                
    def list_coordinates(self, what):
        """
        Extract the coordinates of the specified kind(s) of atoms and return 
        them as a list.
        
        @param what: a list of atom kinds, e.g. ['N', 'CA', 'C']
        @type what: list
        
        @return: a list of lists, each internal list corresponding to the coordinates 
                 of a 3D vector
        @rtype: list
        
        @raise Broken3DStructureError: if a specific atom kind cannot be retrieved from a residue
        """
        coords = [ ]
        
        for residue in self.residues:
            if not residue.has_structure:
                raise Missing3DStructureError('Discontinuous structure at residue {0}'.format(residue))
            try:
                for atom_kind in what:
                    coords.append(residue.structure[atom_kind].vector.row())
            except csb.pyutils.ItemNotFoundError:
                raise Broken3DStructureError('Could not retrieve {0} atom from the structure'.format(atom_kind))
            
        return numpy.array(coords)
    
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
        
        r = [ Vector(*row) for row in r ]
        
        return SuperimposeInfo(RotationMatrix(*r), Vector(*t), rmsd=rmsd)
                              
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
        
        if how == AlignmentTypes.Global:
            fit = csb.bio.utils.fit
        else:
            fit = csb.bio.utils.fit_wellordered
            
        r, t, tm = csb.bio.utils.tm_superimpose(x, y, fit)
        r = [ Vector(*row) for row in r ]
        
        return SuperimposeInfo(RotationMatrix(*r), Vector(*t), tm_score=tm)         
    
    def tm_score(self, other, what=['CA']):
        """
        Compute the C-alpha TM-Score against another chain (assuming equal chain length
        and optimal configuration - no fitting is done).        
        
        @param other: the subject (movable) chain
        @type other: L{Chain}
        @param what: a list of atom kinds, e.g. ['CA']
        @type what: list
        @param how: fitting method (global or local) - a member of the L{AlignmentTypes} enum
        @type how: L{csb.pyutils.EnumItem}
                
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
    
    def update(self):
        raise NotImplementedError()
        
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
        
class Residue(object):
    """
    Base class representing a single residue.
    
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
        self._rank = int(rank)
        self._structure = ResidueAtomsTable(self)    
        self._secondary_structure = None
        self._torsion = None
        self._sequence_number = None
        self._insertion_code = None
        self._container = None
        
        self.type = type
        self.id = sequence_number, insertion_code
        
    def __repr__(self):
        return '<{1} [{0.rank}]: {0.type!r} {0.id}>'.format(self, self.__class__.__name__)

    def __getitem__(self, key):
        return self._structure[key]
    
    def __iter__(self):
        return iter(self._structure)

        
    @property
    def type(self):
        return self._type
    @type.setter
    def type(self, type):
        if not isinstance(type, csb.pyutils.EnumItem):
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
    def structure(self):
        return self._structure

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
    def id(self, (sequence_number, insertion_code)):
        old_id = self.id
        id = ''
        if sequence_number is not None:
            sequence_number = int(sequence_number)
            id = str(sequence_number)
        if insertion_code is not None:
            insertion_code = str(insertion_code).strip()
            id += insertion_code
            if sequence_number is None:
                raise ValueError('sequence_number must be defined when an insertion_code is specified.')
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
        if hasattr(self, 'structure') and self.structure is not None:
            return len(self.structure) > 0
        else:
            return False
        
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
    @param type: residue type - a member of L{SequenceAlphabets.Protein}
    @type type: L{csb.pyutils.EnumItem}
    @param sequence_number: PDB sequence number of the residue
    @type sequence_number: str
    @param insertion_code: PDB insertion code, if any
    @type insertion_code: str    
    """
    
    def __init__(self, rank, type, sequence_number=None, insertion_code=None):
          
        if isinstance(type, basestring):
            try:
                if len(type) == 3:
                    type = csb.pyutils.Enum.parsename(SequenceAlphabets.Protein, type)
                else:    
                    type = csb.pyutils.Enum.parse(SequenceAlphabets.Protein, type)          
            except (csb.pyutils.EnumMemberError, csb.pyutils.EnumValueError):
                raise ValueError("What is '{0}'?".format(type))
        elif not csb.pyutils.Enum.ismember(type, SequenceAlphabets.Protein):
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
            n = self._structure['N'].vector
            ca = self._structure['CA'].vector
            c = self._structure['C'].vector
        except csb.pyutils.ItemNotFoundError as missing_atom:
            if strict:
                raise Broken3DStructureError('Could not retrieve {0} atom from the current residue {1!r}.'.format(missing_atom, self))
            else:
                return angles
        
        try:
            if prev_residue is not None and prev_residue.has_structure:
                prev_c = prev_residue._structure['C'].vector
                angles.phi = Vector.dihedral(prev_c, n, ca, c)
        except csb.pyutils.ItemNotFoundError as missing_prevatom:
            if strict:
                raise Broken3DStructureError('Could not retrieve {0} atom from the i-1 residue {1!r}.'.format(missing_prevatom, prev_residue))    
        try:
            if next_residue is not None and next_residue.has_structure:    
                next_n = next_residue._structure['N'].vector
                angles.psi = Vector.dihedral(n, ca, c, next_n)
                next_ca = next_residue._structure['CA'].vector
                angles.omega = Vector.dihedral(ca, c, next_n, next_ca)
        except csb.pyutils.ItemNotFoundError as missing_nextatom:
            if strict:
                raise Broken3DStructureError('Could not retrieve {0} atom from the i+1 residue {1!r}.'.format(missing_nextatom, next_residue))              
                                
        return angles

class NucleicResidue(Residue):
    """
    Represents a single nucleotide residue.
    
    @param rank: rank of the residue with respect to the chain
    @type rank: int
    @param type: residue type - a member of L{SequenceAlphabets.Nucleic}
    @type type: L{csb.pyutils.EnumItem}
    @param sequence_number: PDB sequence number of the residue
    @type sequence_number: str
    @param insertion_code: PDB insertion code, if any
    @type insertion_code: str        
    """
    
    def __init__(self, rank, type, sequence_number=None, insertion_code=None):
        
        if isinstance(type, basestring):
            try:
                if len(type) > 1:
                    type = csb.pyutils.Enum.parsename(SequenceAlphabets.Nucleic, type)
                else:    
                    type = csb.pyutils.Enum.parse(SequenceAlphabets.Nucleic, type)
            except (csb.pyutils.EnumMemberError, csb.pyutils.EnumValueError):
                raise ValueError("What is '{0}'?".format(type))
        elif not csb.pyutils.Enum.ismember(type, SequenceAlphabets.Nucleic):
            raise TypeError(type)
            
        super(NucleicResidue, self).__init__(rank, type, sequence_number, insertion_code)  

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
                raise DuplicateAtomIDError('Atom {0} is already defined for the current residue.'.format(atom.name))
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
        
        super(ResidueAtomsTable, self).update({atom_name: atom})  
    
class Atom(object):
    """
    Represents a single atom in space.
    
    @param serial_number: atom's UID
    @type serial_number: int
    @param name: atom's name
    @type name: str
    @param element: corresponding L{ChemElements}
    @type element: L{csb.pyutils.EnumItem}
    @param vector: atom's coordinates
    @type vector: L{Vector}
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

        if not isinstance(name, basestring):
            raise TypeError(name)
        name_compact = name.strip()
        if len(name_compact) < 1:
            raise ValueError(name)
        self._name = name_compact
        self._full_name = name
            
        if isinstance(element, basestring):
            element = csb.pyutils.Enum.parsename(ChemElements, element)
        elif element is None:
            pass
        elif not csb.pyutils.Enum.ismember(element, ChemElements):
            raise TypeError(element)
        self._element = element

        # pass type- and value-checking control to setters
        self.serial_number = serial_number
        self.vector = vector
        
    def __repr__(self):
        return "<Atom [{0.serial_number}]: {0.name}>".format(self)
        
    def __cmp__(self, other):
        return cmp(self.serial_number, other.serial_number)
    
    def _transform_vector(self, rotation, translation):
        self.vector = self.vector.transform(rotation, translation)
        
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
        if not isinstance(vector, Vector):
            raise TypeError(vector)
        self._vector = vector
        
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
        atom._proxy = self
    
    def update(self, new_items):
        raise NotImplementedError()
    
    def _transform_vector(self, rotation, translation):
        for atom in self:
            atom._transform_vector(rotation, translation)
        self._vector = self._rep._vector
        
    def __update_rep(self, atom):
        
        if self._rep is None or (self._rep.occupancy < atom.occupancy):
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
        
class Vector(object):
    """
    Simple 3D vector. Provides basic vector calculations: cross, dot, +, -, angle.
    
    @param x: x coordinate
    @type x: float
    @param y: y coordinate
    @type y: float
    @param z: z coordinate
    @type z: float        
    """  
    def __init__(self, x, y, z):
        
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)  
        
    def __repr__(self):
        return '<Vector: x={0.x}, y={0.y}, z={0.z}>'.format(self)
        
    def __getitem__(self, i):
        return (self.x, self.y, self.z)[i]
        
    def __neg__(self):
        return Vector(-self.x, -self.y, -self.z)
                
    def __add__(self, other):   
        return Vector(self.x + other.x, 
                      self.y + other.y, 
                      self.z + other.z)

    def __sub__(self, other):
        return Vector(self.x - other.x, 
                      self.y - other.y, 
                      self.z - other.z)       
        
    def row(self):
        """
        @return: a list representation of the vector: [x, y, z]
        @rtype: list
        """
        return [self.x, self.y, self.z] 
    
    def clone(self):
        """
        @return: deep copy of the vector
        @rtype: L{Vector}
        """
        return Vector(self.x, self.y, self.z)
    
    def transform(self, rotation, translation=None):
        """
        Apply L{RotationMatrix} and translation L{Vector}.
        
        @type rotation: L{RotationMatrix}
        @type translation: L{Vector}
        
        @return: transformed copy of the vector
        @rtype: L{Vector}
        """
        if translation is None:
            translation = Vector(0, 0, 0)
            
        v = Vector(0, 0, 0)    
        v.x = self.dot(rotation.x) + translation.x
        v.y = self.dot(rotation.y) + translation.y
        v.z = self.dot(rotation.z) + translation.z
        
        return v
    
    def normalized(self):
        """
        @return: normalized copy of the vector
        @rtype: L{Vector}
        """
        n = self.norm
        return Vector(self.x/n, self.y/n, self.z/n)
        
    def dot(self, other):
        """
        @param other: right-hand-side term
        @type other: L{Vector}
                
        @return: dot product of C{self} and C{other}
        @rtype: float
        """
        return self.x * other.x + self.y * other.y + self.z * other.z 
    
    def cross(self, other):
        """
        @param other: right-hand-side term
        @type other: L{Vector}
        
        @return: cross product of C{self} and C{other}
        @rtype: L{Vector}
        """
        return Vector(self.y * other.z - self.z * other.y,
                      self.z * other.x - self.x * other.z, 
                      self.x * other.y - self.y * other.x)
    
    def angle(self, other, units=AngleUnits.Degrees):
        """
        @param other: right-hand-side term
        @type other: L{Vector}
        @param units: target L{AngleUnits}
        @type units: L{csb.pyutils.EnumItem}
                
        @return: angle between C{self} and C{other} in C{units}
        @rtype: float
        
        @raise ValueError: if any of the vectors is a null vector
        @raise ValueError: when the target units are not valid 
        """
        
        if self.is_null or other.is_null:
            raise ValueError('Zero (Null) vector is not a valid argument.')
        
        cos_angle = self.dot(other) / float(self.norm * other.norm)
        
        if cos_angle < -1.0:
            cos_angle = -1.0
        if cos_angle > 1.0:
            cos_angle = 1.0
            
        angle =  math.acos(cos_angle)
                
        if units == AngleUnits.Degrees:
            return math.degrees(angle)
        elif units == AngleUnits.Radians:
            return angle 
        else:
            raise ValueError('Unknown angle units {0}'.format(units))
        
    @property
    def norm(self):
        return self.dot(self) ** 0.5
    
    @property
    def length(self):
        return self.norm
    
    @property
    def is_null(self):
        return self.x == self.y == self.z == 0
    
    @staticmethod
    def dihedral(a, b, c, d):
        """
        Compute the dihedral angle from a set of 4 connected atoms (a-b-c-d)
        around the b-c bond. 
        
        @param a: first vector
        @type a: L{Vector}
        @param b: second vector
        @type b: L{Vector}
        @param c: third vector
        @type c: L{Vector}
        @param d: fourth vector
        @type d: L{Vector}
                                
        @return: the computed angle in degrees: [-180, 180]
        @rtype: float
        """                  
        v = b - c
        m = (a - b).cross(v).normalized()
        n = (d - c).cross(v).normalized()
        
        cos = m.dot(n)
        sin = n.cross(m).dot(v) / v.norm
        angle = math.degrees(math.atan2(sin, cos))        
            
        while angle > 180:
            angle -= 360
        while angle < -180:
            angle += 360 
        
        return angle                    

    @staticmethod
    def dihedral2(a, b, c, d):
        """
        Alternative implementation of Vector.dihedral().
        
        @deprecated: use Vector.dihedral(), although the results should be identical
        """           
        abc_normal = (a - b).cross(c - b)
        bcd_normal = (b - c).cross(d - c)
         
        angle = abc_normal.angle(bcd_normal, units=AngleUnits.Degrees)
         
        if abc_normal.dot(d - c) >= 0:
            angle = -angle      
            
        while angle > 180:
            angle -= 360
        while angle < -180:
            angle += 360
        
        return angle 

class RotationMatrix(object):
    """
    Represents a 3D vector rotation matrix. x, y, z  are the rows of the matrix, 
    each row holds one L{Vector}, such that for some vector::
    
        i' = vector.dot(matrix.i)        (i=x,y,z)
        
    gives the rotated components of the vector.
    
    @type x: L{Vector}
    @type y: L{Vector}
    @type z: L{Vector}        
    """    
    
    def __init__(self, x, y, z):
        
        for i in [x, y, z]:
            if type(i) is not Vector:
                raise TypeError(i)
        
        self.x = x
        self.y = y
        self.z = z 
        
    def matrix(self):
        """
        @return: a list (3x3 'matrix') representation of the matrix
        @rtype: list
        """
        return [ self.x.row(), 
                 self.y.row(), 
                 self.z.row() ]
        
class SuperimposeInfo(object):
    """
    Describes a structural alignment result.
    
    @type rotation: L{RotationMatrix}
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
        if type not in csb.pyutils.Enum.members(SecStructures):
            raise TypeError(type)
        if isinstance(type, basestring):
            type = csb.pyutils.Enum.parse(SecStructures, type)
        
        self.start = start
        self.end = end
        self.type = type
        self._score = None
        
        if score is not None: 
            self.score = score
            
    def __cmp__(self, other):
        return cmp(self.start, other.start)
    
    def __str__(self):
        return self.to_string()
            
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
        
    def to_string(self):
        """
        Dump the element as a string.
        
        @return: string representation of the element
        @rtype: str
        """
        return str(self.type) * self.length

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
                       
        if string is not None:
            for motif in SecondaryStructure.parse(string, conf_string):
                self.append(motif)
                
    def __str__(self):
        return self.to_string()
            
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
        if not isinstance(string, basestring):
            raise TypeError(string)
                
        string = ''.join(re.split('\s+', string))
        if conf_string is not None:
            conf_string = ''.join(re.split('\s+', conf_string))
            if not len(string) == len(conf_string):
                raise ValueError('The confidence string has unexpected length.')
        motifs = [ ]

        if not len(string) > 0:
            raise ValueError(string)      
        
        currel = string[0]
        start = 0
                
        for i, char in enumerate(string + '.'):
            
            if currel != char:
                type = csb.pyutils.Enum.parse(SecStructures, currel)
                confidence = None
                if conf_string is not None:
                    confidence = list(conf_string[start : i])
                    confidence = map(int, confidence)
                motif = SecondaryStructureElement(start + 1, i, type, confidence)
                motifs.append(motif)
                
                currel = char
                start = i

        return motifs   
    
    def to_string(self):
        """
        Get back the string representation of the secondary structure.
        
        @return: a string of secondary structure elements
        @rtype: str
        """
        return ''.join([ str(m) for m in self ])   
        
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
    def __init__(self, items=None):
        super(TorsionAnglesCollection, self).__init__(items, type=TorsionAngles, start_index=1)

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
    
    def rmsd(self, other):
        """
        Calculate the Circular RSMD against another TorsionAnglesCollection.
        
        @param other: subject (right-hand-term)
        @type other: L{TorsionAnglesCollection}
        
        @return: RMSD based on torsion angles
        @rtype: float
        
        @raise Broken3DStructureError: on discontinuous torsion angle collections 
        @raise ValueError: on mismatching torsion angles collection lengths
        """           
        if len(self) != len(other) or len(self) < 1:
            raise ValueError('Both collections must be of the same and positive length (got {0} and {1})'.format(len(self), len(other)))
                
        from numpy import array, sin, cos, sqrt
        
        query = TorsionAnglesCollection()
        subject = TorsionAnglesCollection()
        
        for q, s in izip(iter(self), iter(other)):
            if q.phi is None or q.psi is None or s.phi is None or s.psi is None:
                raise Broken3DStructureError('Should we compute RMSD if the structure is discontinuous?')
            else:
                q = copy.copy(q)
                q.to_radians()
                
                s = copy.copy(s)
                s.to_radians()
                
                query.append(q)
                subject.append(s)
                
        phi_diff = array(query.phi) - array(subject.phi)
        psi_diff = array(query.psi) - array(subject.psi)

        window = min([len(self), len(other)])
        if window == 0:
            raise Missing3DStructureError("Cannot compute RMSD with/against an empty TorsionAnglesCollection.")    
        elif len(phi_diff) < int(0.95 * window):                                                                             # useless now since exception is thrown near izip?
            raise Broken3DStructureError("Got {0} out of {1} phi-psi pairs - less than 95%.".format(len(phi_diff), window))

        assert len(phi_diff) == len(psi_diff)
        
        r = sin(phi_diff).sum()**2 + cos(phi_diff).sum()**2 + sin(psi_diff).sum()**2 + cos(psi_diff).sum()**2
        
        rmsd = 1 - (1.0/len(phi_diff)) * sqrt(r/2.0)
        
        return rmsd    
           
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
            if isinstance(units, csb.pyutils.EnumItem):
                if not csb.pyutils.Enum.ismember(units, AngleUnits):
                    raise ValueError()
            else:
                units = csb.pyutils.Enum.parse(AngleUnits, units, ignore_case=True)
                
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
