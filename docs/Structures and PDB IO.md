## Protein Structures

The ``csb.bio.structure`` module defines some of the most fundamental 
abstractions in the library: ``Structure``, ``Chain``, ``Residue`` and 
``Atom``. Instances of these objects may exist independently and that is 
perfectly fine, but usually they are part of a ``Composite`` aggregation. 
The root node in this Composite is a ``Structure`` (or ``Ensemble``). 
``Structure``-s are composed of ``Chain``-s, and each ``Chain`` is a collection 
of ``Residue``-s. The leaf nodes are ``Atom``-s.

All of these objects implement the base ``AbstractEntity`` interface. 
Therefore, every node in the Composite can be transformed:

```python
>>> r, t = [rotation matrix](rotation-matrix), [translation vector](translation-vector)
>>> entity.transform(r, t)
```

and it knows its immediate children:

```python
>>> entity.items
<iterator>    # over all immediate child entities
```

If you want to traverse the complete Composite tree, starting at arbitrary 
level, and down to the lowest level, use one of the ``CompositeEntityIterators``. 
Or just call ``AbstractEntity.components()``:

```python
>>> entity.components()
<iterator>   # over all descendants, of any type, at any level
>>> entity.components(klass=Residue)
<iterator>   # over all Residue descendants
```
Some of the inner objects in this hierarchy behave just like dictionaries (but are not):

```python
>>> structure.chains['A']('A')       # access chain A by ID
<Chain A: Protein>
>>> structure['A']('A')              # the same
<Chain A: Protein>
>>> residue.atoms['CS']('CS')          
<Atom: CA>                      # access an atom by its name
>>> residue.atoms['CS']('CS')          
<Atom: CA>                      # the same
```

Others behave like list collections:

```python
>>> chain.residues[10](10)               # 1-based access to the residues in the chain
<ProteinResidue [10](10): PRO 10>
>>> chain[10](10)                        # 0-based, list-like access
<ProteinResidue [11](11): GLY 11>
```

Step-wise building of ``Ensemble``-s, ``Chain``-s and ``Residue``-s is supported 
through a number of append methods, for example:

```python
>>> residue = ProteinResidue(401, ProteinAlphabet.ALA)
>>> s.chains['A']('A').residues.append(residue)
```

See ``EnsembleModelsCollection``, ``StructureChainsTable``, 
``ChainResiduesCollection`` and ``ResidueAtomsTable`` in our API docs for more 
details.

Some other objects in this module of potential interest are the 
self-explanatory ``SecondaryStructure`` and ``TorsionAngles``.


## PDB I/O

CSB comes with a number of PDB structure parsers, format builders and 
database providers, all defined in the ``csb.bio.io.wwpdb`` package. 
The most basic usage is:

```python
>>> parser = StructureParser('structure.pdb')
>>> parser.parse_structure()
<Structure>     # a Structure object (model)
```

or if this is an NMR ensemble:

```python
>>> parser.parse_models()
<Ensemble>      # an Ensemble object (collection of alternative Structure-s)
```

This module introduces a family of PDB file parsers. The common interface 
of all parsers is defined in ``AbstractStructureParser``. This class has 
several implementations:

* ``RegularStructureParser`` - handles normal PDB files with SEQRES fields
* ``LegacyStructureParser`` - reads structures from legacy or malformed PDB 
  files, which are lacking SEQRES records (initializes all residues from 
  the ATOMs instead)
* ``PDBHeaderParser`` - reads only the headers of the PDB files and produces 
  structures without coordinates. Useful for reading metadata (e.g. 
  ccession numbers or just plain SEQRES sequences) with minimum overhead

Unless you have a special reason, you should use the ``StructureParser`` 
factory, which returns a proper ``AbstractStructureParser`` implementation, 
depending on the input PDB file. If the input file looks like a regular 
PDB file, the factory returns a ``RegularStructureParser``, otherwise it 
instantiates ``LegacyStructureParser``. ``StructureParser`` is in fact an 
alias for ``AbstractStructureParser.create_parser``.

Writing your own, customized PDB parser is easy. Suppose that you are 
trying to parse a PDB-like file which misuses the charge column to store 
custom info. This will certainly crash ``AbstractStructureParser`` 
(for good), but you can create your own parser as a workaround. All you 
need to do is override the virtual ``_read_charge_field`` hook method:

```python
class CustomParser(RegularStructureParser):

    def _read_charge(self, line):
        try:
            return super(CustomParser, self)._read_charge(line)
        except StructureFormatError:
            return None
```

Another important abstraction in this module is ``StructureProvider``. 
It has several implementations which can be used to retrieve PDB 
Structures from various sources: file system directories, remote 
URLs, etc. You can easily create your own provider as well. See 
``StructureProvider`` for details.

Finally, this module gives you some ``FileBuilder``-s, used for text 
serialization of ``Structure``-s and ``Ensemble``-s:

```python
>>> builder = PDBFileBuilder(stream)
>>> builder.add_header(structure)
>>> builder.add_structure(structure)
```

where stream is any Python stream, e.g. an open file or sys.stdout.


 
