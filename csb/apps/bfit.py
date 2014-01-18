"""
Python application for robust structure superposition of two structures.
bfit models non-rigid displacements in protein ensembles with outlier-tolerant
probability distributions.
"""
import numpy

import csb.apps
import csb.bio.utils

from csb.bio.io.wwpdb import LegacyStructureParser
from csb.bio.sequence import SequenceAlignment


class ExitCodes(csb.apps.ExitCodes):
    IO_ERROR = 2
    INPUT_ERROR = 3

class AppRunner(csb.apps.AppRunner):

    @property
    def target(self):
        return BFitApp

    def command_line(self):
        
        cmd = csb.apps.ArgHandler(self.program, __doc__)

        # Input structures
        cmd.add_positional_argument('pdb1', str,
                                    'full path to the first structure')

        cmd.add_positional_argument('pdb2', str,
                                    'full path to the second structure')

        # Optional arguments 
        cmd.add_scalar_option('chain1', 'c', str,
                              'Chain of the first structure',
                              default='A')

        cmd.add_scalar_option('chain2', 'd', str,
                              'Chain of the second structure',
                              default='A')
        
        cmd.add_scalar_option('scalemixture', 's', str,
                              'Scale mixture distribution',
                              default='student',
                              choices=['student', 'k'])

        
        cmd.add_scalar_option('alignment', 'a', str,
                              'Alignment in fasta format defining equivalent positions\n'
                              + 'Assumes that chain1 is the first sequence of '
                              + 'the alignment and chain2 the second sequence')

        cmd.add_scalar_option('outfile', 'o', str,
                              'file to which the rotated second ' + 
                              'structure will be written',
                              default='bfit.pdb')

        cmd.add_scalar_option('niter', 'n', int,
                              'Number of optimization steps',
                              default=200)

        cmd.add_boolean_option('em', None,
                               'Use the EM algorithm for optimsation',
                               default = False)

        return cmd



class BFitApp(csb.apps.Application):
    """
    Python application for robust structure superposition of two protein structures
    """

    def main(self):
        try:
            parser = LegacyStructureParser(self.args.pdb1)
            r = parser.parse()

            parser = LegacyStructureParser(self.args.pdb2)
            m = parser.parse()
        except IOError as e:
            self.exit('PDB file parsing failed\n' + str(e.value), ExitCodes.IO_ERROR)

        X = numpy.array(r[self.args.chain1].get_coordinates(['CA'], True))
        Y = numpy.array(m[self.args.chain2].get_coordinates(['CA'], True))

        if self.args.alignment is not None:
            align = SequenceAlignment.parse(file(self.args.alignment).read())
            align = align[:2, :]
            
            matches = []
            for i in range(1, align.length + 1):
                if not align.gap_at(i):
                    matches.append([align.columns[i][0].rank - 1,
                                    align.columns[i][1].rank - 1])
            matches = numpy.array(matches)
            X = X[matches[:, 0], :]
            Y = Y[matches[:, 1], :]

        
        if len(X) != len(Y):
            self.exit('Structures are of different lengths,' + 
                      ' please specify an alignment',
                      ExitCodes.INPUT_ERROR)

        R, t = csb.bio.utils.bfit(X, Y, self.args.niter,
                self.args.scalemixture, self.args.em)

        m.transform(R, t)
        m.to_pdb(self.args.outfile)
        

def main():
    AppRunner().run()
    
    
if __name__ == '__main__':
    main()
    
