## Basic Tutorial

Perhaps one of the most frequently used parts of the library is 
the ``csb.bio.structure`` module, which provides the ``Structure``, 
``Chain``, ``Residue`` and ``Atom`` objects. You could easily build a 
``Structure`` from scratch, but a far more common scenario is parsing 
a structure from a PDB file using one of the ``AbstractStructureParsers``. 
All bio IO objects, including the ``StructureParser`` factory, are 
defined in ``csb.bio.io`` and sub-packages:

```python
>>> from csb.bio.io.wwpdb import StructureParser
>>> p = StructureParser("/some/file/pdb1x80.ent")
>>> s = p.parse_structure()
>>> print(s)
<Structure: 1x80, 2 chains>
```
     
The last statement will return a ``csb.bio.structure.Structure`` instance, 
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
``AbstractEntity``. The above loop can be shortened:

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

If you want to dive into statistics, you could peek inside ``csb.statistics`` 
and its sub-packages. For example, ``csb.statistics.pdf`` contains a collection 
of probability density objects, like ``Gaussian`` or ``Gamma``.

But chances are you would first like to try reading some files, so you could start 
exploring ``csb.bio.io`` right now. As we have already seen, ``csb.bio.io.wwpdb`` 
provides PDB Structure parsers, for example ``csb.bio.io.wwpdb.RegularStructureParser`` 
and ``csb.bio.io.wwpdb.LegacyStructureParser``.

``csb.bio.io.fasta`` is all about reading FASTA Sequences and ``SequenceAlignment``-s. 
Be sure to check out ``csb.bio.io.fasta.SequenceParser``, 
``csb.bio.io.fasta.SequenceAlignmentReader`` and ``csb.bio.io.fasta.StructureAlignmentFactory``.

If you are working with HHpred (``ProfileHMM``-s, ``HHpredHit``-s), then 
``csb.bio.io.hhpred`` is for you. This package provides ``csb.bio.io.hhpred.HHProfileParser`` 
and ``csb.bio.io.hhpred.HHOutputParser``, which are used to read ``.hhm`` and ``.hhr`` files.

Finally, if you want to make some nice plots with matplotlib, you may like the clean 
object-oriented interface of our ``Chart``. See ``csb.io.plots`` and maybe also 
``csb.io.tsv`` to get started.

