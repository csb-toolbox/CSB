"""
NMR related objects.
"""

import os
import csb.io.tsv
import csb.pyutils as pu

from csb.bio.sequence import SequenceAlphabets


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
            residue = pu.Enum.parsename(SequenceAlphabets.Protein, row[0])
            nucleus, value = row[1:]
            
            if residue not in self._reference:
                self._reference[residue] = {}
            
            self._reference[residue][nucleus] = value
        
        header = 'Residue:str Nucleus:str CS1:float CS2:float CS3:float CS4:float'
        
        for row in csb.io.tsv.Table.from_tsv(cor, header):   
            residue = pu.Enum.parsename(SequenceAlphabets.Protein, row[0])
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
                    residue = pu.Enum.parse(SequenceAlphabets.Protein, residue)
                else:
                    residue = pu.Enum.parsename(SequenceAlphabets.Protein, residue)
            else:
                if residue.enum is not SequenceAlphabets.Protein:
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
            