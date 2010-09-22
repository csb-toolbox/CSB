"""
FASTA format parsers.
"""

import csb.io

from abc import abstractmethod, ABCMeta
from csb.bio.sequence import Sequence, PDBSequence


class BaseSequenceParser(object):
    """
    FASTA parser interface. Subclasses are obliged to implement their proper
    L{Sequence} C{_factory} method.
    """
    __metaclass__ = ABCMeta

    @abstractmethod
    def _factory(self):
        """
        @return: an object which implements L{Sequence} construction from a
                 single FASTA string, such that:

                    >>> obj.from_string(string)
                    L{Sequence}

        @rtype: object

        @note: This method must be implemented by each individual subclass
        """
        pass

    def parse_string(self, fasta_string):
        """
        Read FASTA sequences from a (m)FASTA-formatted string

        @param fasta_string: FASTA string to parse
        @type fasta_string: str

        @return: a list of L{Sequence}s
        @rtype: list
        """
        import cStringIO

        stream = cStringIO.StringIO()
        stream.write(fasta_string)

        return self.parse_file(stream)

    def parse_file(self, fasta_file):
        """
        Read FASTA sequences from a (m)FASTA file

        @param fasta_file: input FASTA file name or opened stream
        @type fasta_file: str, file

        @return: a list of L{Sequence}s
        @rtype: list
        """
        if isinstance(fasta_file, basestring):
            stream = open(fasta_file)
        else:
            stream = fasta_file

        seqs = []

        reader = csb.io.EntryReader(stream, '>', None)
        factory = self._factory()

        for entry in reader.entries():
            seqs.append(factory.from_string(entry))

        return seqs


class SequenceParser(BaseSequenceParser):
    """
    Generic FASTA format parser
    """

    def _factory(self):
        return Sequence


class PDBSequenceParser(BaseSequenceParser):
    """
    PDB SEQRES FASTA format parser
    """

    def _factory(self):
        return PDBSequence
