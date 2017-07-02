## Basic Tutorial

Perhaps one of the most frequently used parts of the library is 
the _csb.bio.structure_ module, which provides the _Structure_, 
_Chain_, _Residue_ and _Atom_ objects. You could easily build a 
_Structure_ from scratch, but a far more common scenario is parsing 
a structure from a PDB file using one of the _AbstractStructureParsers_. 
All bio IO objects, including the _StructureParser_ factory, are 
defined in _csb.bio.io_ and sub-packages:

```python
>>> from csb.bio.io.wwpdb import StructureParser
>>> p = StructureParser("/some/file/pdb1x80.ent")
>>> s = p.parse_structure()
>>> print(s)
<Structure: 1x80, 2 chains>
```
     
The last statement will return a _csb.bio.structure.Structure_ instance, 
which is a composite hierarchical object:

```python
>>> for chain_id in s.chains:
    chain = s.chains[chain_id](chain_id)
    for residue in chain.residues:
        for atom_id in residue.atoms:
            atom = residue.atoms[atom_id](atom_id)
            print(atom.vector)
```

Some of the inner objects in this hierarchy behave just like dictionaries (but are not):

```python
>>> s.chains['A']('A')        # access chain A by ID
<Chain A: Protein>
>>> s['A']('A')               # the same
<Chain A: Protein>
```

Others behave like collections:

```python
>>> chain.residues[10](10)               # 1-based access to the residues in the chain
<ProteinResidue [10](10): PRO 10>
>>> chain[10](10)                        # 0-based, list-like access
<ProteinResidue [11](11): GLY 11>
```

But all entities are iterable because they inherit the items iterator from 
_AbstractEntity_. The above loop can be shortened:

```python
>>> for chain in s.items:
    for residue in chain.items:
        for atom in residue.items:
            print(atom.vector)
```

or even more:

```python
>>> from csb.bio.structure import Atom
>>> for atom in s.components(klass=Atom):
    print(atom.vector)
```
 
You may also be interested in extracting a sub-chain from this structure:

```python
>>> s.chains['B']('B').subregion(3, 20)    # from positions 3 to 20, inclusive
<Chain B: Protein>
```

or modifying it in some way, for example, in order to append a new residue, try:

```python
>>> from csb.bio.structure import ProteinResidue
>>> from csb.bio.sequence import ProteinAlphabet
>>> residue = ProteinResidue(401, ProteinAlphabet.ALA)
>>> s.chains['A']('A').residues.append(residue)
```

Finally, you would probably want to save your structure back to a PDB file:

```python
>>> s.to_pdb('/some/file/name.pdb')
```

## Where to go from here

If you want to dive into statistics, you could peek inside _csb.statistics_ 
and its sub-packages. For example, _csb.statistics.pdf_ contains a collection 
of probability density objects, like _Gaussian_ or _Gamma_.

But chances are you would first like to try reading some files, so you could start 
exploring _csb.bio.io_ right now. As we have already seen, _csb.bio.io.wwpdb_ 
provides PDB Structure parsers, for example _csb.bio.io.wwpdb.RegularStructureParser_ 
and _csb.bio.io.wwpdb.LegacyStructureParser_.

_csb.bio.io.fasta_ is all about reading FASTA Sequences and _SequenceAlignment_-s. 
Be sure to check out _csb.bio.io.fasta.SequenceParser_, 
_csb.bio.io.fasta.SequenceAlignmentReader_ and _csb.bio.io.fasta.StructureAlignmentFactory_.

If you are working with HHpred (_ProfileHMM_-s, _HHpredHit_-s), then 
_csb.bio.io.hhpred_ is for you. This package provides _csb.bio.io.hhpred.HHProfileParser_ 
and _csb.bio.io.hhpred.HHOutputParser_, which are used to read _.hhm_ and _.hhr_ files.

Finally, if you want to make some nice plots with matplotlib, you may like the clean 
object-oriented interface of our _Chart_. See _csb.io.plots_ and maybe also 
_csb.io.tsv_ to get started.

