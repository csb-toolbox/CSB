"""
Python application for robust structure superposition of two structures
"""
import numpy

import csb.apps
import csb.bio.utils

from csb.bio.io.wwpdb import LegacyStructureParser
from csb.bio.utils import probabilistic_fit
from csb.statistics.scalemixture import ScaleMixture, GammaPrior, InvGammaPrior
from csb.bio.sequence import SequenceAlignment


class ExitCodes(csb.apps.ExitCodes):
    IO_ERROR = 2
    INPUT_ERROR = 3

class AppRunner(csb.apps.AppRunner):

    @property
    def target(self):
        return BFitApp

    def command_line(self):
        __doc__ = "bFit is a probabilistic method for robust superposition" \
                  + "and comparison of protein structures. To do so, "\
                  + "non-rigid displacements in protein structures are "\
                  + "modelled with outlier-tolerant probability "\
                  + "distributions."
        
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
                              'Number of Gibbs sampling steps',
                              default=200)

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

        X = numpy.array(r[self.args.chain1].list_coordinates(['CA'], True))
        Y = numpy.array(m[self.args.chain2].list_coordinates(['CA'], True))

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
        
        if self.args.scalemixture == 'student':
            prior = GammaPrior()
        elif self.args.scalemixture == 'k':
            prior = InvGammaPrior()
        
        mixture = ScaleMixture(scales=X.shape[0],
                               prior=prior, d=3)

        R, t = csb.bio.utils.fit(X, Y)

        # gibbs sampling cycle
        for i in range(self.args.niter):
            # apply rotation
            data = numpy.sum((X - numpy.dot(Y, numpy.transpose(R)) - t) ** 2,
                             axis= -1) ** (1. / 2)
            # sample scales
            mixture.estimate(data)
            # sample rotations
            R, t = probabilistic_fit(X, Y, mixture.scales)
            

        m.apply_transformation(R, t)
        m.to_pdb(self.args.outfile)
        


if __name__ == '__main__':
    AppRunner().run()
