## Protein Structures

The _csb.bio.structure_ module defines some of the most fundamental 
abstractions in the library: _Structure_, _Chain_, _Residue_ and 
_Atom_. Instances of these objects may exist independently and that is 
perfectly fine, but usually they are part of a _Composite_ aggregation. 
The root node in this Composite is a _Structure_ (or _Ensemble_). 
_Structure_s are composed of _Chain_s, and each _Chain_ is a collection 
of _Residue_s. The leaf node is _Atom_.

All of these objects implement the base _AbstractEntity_ interface. 
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
level, and down to the lowest level, use one of the _CompositeEntityIterators_. 
Or just call _AbstractEntity.components()_:

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

Step-wise building of _Ensemble_-s, _Chain_-s and _Residue_-s is supported 
through a number of append methods, for example:

```python
>>> residue = ProteinResidue(401, ProteinAlphabet.ALA)
>>> s.chains['A']('A').residues.append(residue)
```

See _EnsembleModelsCollection_, _StructureChainsTable_, 
_ChainResiduesCollection_ and _ResidueAtomsTable_ in our API docs for more 
details.

Some other objects in this module of potential interest are the 
self-explanatory _SecondaryStructure_ and _TorsionAngles_.


## PDB I/O

CSB comes with a number of PDB structure parsers, format builders and 
database providers, all defined in the _csb.bio.io.wwpdb_ package. 
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
of all parsers is defined in _AbstractStructureParser_. This class has 
several implementations:

* _RegularStructureParser_ - handles normal PDB files with SEQRES fields
* _LegacyStructureParser_ - reads structures from legacy or malformed PDB 
  files, which are lacking SEQRES records (initializes all residues from 
  the ATOMs instead)
* _PDBHeaderParser_ - reads only the headers of the PDB files and produces 
  structures without coordinates. Useful for reading metadata (e.g. 
  ccession numbers or just plain SEQRES sequences) with minimum overhead

Unless you have a special reason, you should use the _StructureParser_ 
factory, which returns a proper _AbstractStructureParser_ implementation, 
depending on the input PDB file. If the input file looks like a regular 
PDB file, the factory returns a _RegularStructureParser_, otherwise it 
instantiates _LegacyStructureParser_. _StructureParser_ is in fact an 
alias for ``AbstractStructureParser.create_parser``.

Writing your own, customized PDB parser is easy. Suppose that you are 
trying to parse a PDB-like file which misuses the charge column to store 
custom info. This will certainly crash ``AbstractStructureParser`` 
(for good), but you can create your own parser as a workaround. All you 
need to to is to override the virtual ``_read_charge_field`` hook method:

```python
class CustomParser(RegularStructureParser):

    def _read_charge(self, line):
        try:
            return super(CustomParser, self)._read_charge(line)
        except StructureFormatError:
            return None
```

Another important abstraction in this module is _StructureProvider_. 
It has several implementations which can be used to retrieve PDB 
Structures from various sources: file system directories, remote 
URLs, etc. You can easily create your own provider as well. See 
_StructureProvider_ for details.

Finally, this module gives you some _FileBuilder_s, used for text 
serialization of _Structure_s and _Ensemble_-s:

```python
>>> builder = PDBFileBuilder(stream)
>>> builder.add_header(structure)
>>> builder.add_structure(structure)
```

where stream is any Python stream, e.g. an open file or sys.stdout.


 