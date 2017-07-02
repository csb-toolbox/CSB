## Overview

The _csb.bio.sequence_ module defines the base interfaces of our sequence 
and sequence alignment objects: _AbstractSequence_ and _AbstractAlignment_. 
This module provides also a number of useful enumerations, like 
_SequenceTypes_ and _SequenceAlphabets_.

## Sequences

_AbstractSequence_ has a number of implementations. These are of course 
interchangeable, but have different intents and may differ significantly 
in performance. The standard _Sequence_ implementation is what you are 
after if all you need is high performance and efficient storage (e.g. 
when you are parsing big files). _Sequence_ objects store their underlying 
sequences as strings. _RichSequence_-s on the other hand will store their 
residues as _ResidueInfo_ objects, which have the same basic interface as 
the _csb.bio.structure.Residue_ objects. This of course comes at the 
expense of degraded performance. A _ChainSequence_ is a special case of a 
rich sequence, whose residue objects are actually real 
_csb.bio.structure.Residue_-s.

Basic usage:

```python
>>> seq = RichSequence('id', 'desc', 'sequence', SequenceTypes.Protein)
>>> seq.residues[1](1)
<ResidueInfo [1](1): SER>
>>> seq.dump(sys.stdout)
>desc
SEQUENCE
``` 

See _AbstractSequence_ in the API docs for details.

## Alignments

_AbstractAlignment_ defines a table-like interface to access the data 
in an alignment:

```python
>>> ali = SequenceAlignment.parse(">a\nABC\n>b\nA-C")
>>> ali[0, 0](0,-0)
<SequenceAlignment>   # a new alignment, constructed from row #1, column #1
>>> ali[0, 1:3](0,-1_3)
<SequenceAlignment>   # a new alignment, constructed from row #1, columns #2..#3
```

which is just a shorthand for using the standard 1-based interface:

```python
>>> ali.rows[1](1)
<AlignedSequenceAdapter: a, 3>                        # row #1 (first sequence)
>>> ali.columns[1](1)
(<ColumnInfo a [1](1)(1): ALA>, <ColumnInfo b [1](1)(1): ALA>)    # residues at column #1
```

See _AbstractAlignment_ in our API docs for all details and more examples.

There are a number of _AbstractAlignment_ implementations defined here. 
_SequenceAlignment_ is the default one, nothing surprising. _A3MAlignment_ 
is a more special one: the first sequence in the alignment is a master 
sequence. This alignment is usually used in the context of HHpred. More 
important is the StructureAlignment, which is an alignment of 
_csb.bio.structure.Chain_ objects. The residues in every aligned sequence 
are really the _csb.bio.structure.Residue_ objects taken from those chains.

## Sequence and Alignment I/O

CSB provides parsers and writers for sequences and alignments in FASTA 
format, defined in _csb.bio.io.fasta_. The most basic usage is:

```python
>>> parser = SequenceParser()
>>> parser.parse_file('sequences.fa')
<SequenceCollection>   # collection of L{AbstractSequence}s
```

This will load all sequences in memory. If you are parsing a huge file, 
then you could efficiently read the file sequence by sequence:

```python
>>> for seq in parser.read('sequences.fa'):
        ...            # seq is an L{AbstractSequence}
```

_BaseSequenceParser_ is the central class in this module, which defines 
a common infrastructure for all sequence readers. _SequenceParser_ is a 
standard implementation, and _PDBSequenceParser_ is specialized to read 
FASTA sequences with PDB headers.

For parsing alignments, have a look at _SequenceAlignmentReader_ and 
_StructureAlignmentFactory_.

Finally, this module provides a number of _OutputBuilder_-s, which know 
how to write _AbstractSequence_ and _AbstractAlignment_ objects to FASTA 
files:

```python
>>> with open('file.fa', 'w') as out:
        builder = OutputBuilder.create(AlignmentFormats.FASTA, out)
        builder.add_alignment(alignment)
        builder.add_sequence(sequence)
        ...
```

or you could instantiate any of the _OutputBuilder_-s directly.
