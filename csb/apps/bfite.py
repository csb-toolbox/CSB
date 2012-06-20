"""
Python application for robust structure superposition of an ensemble of structures.
bfite models non-rigid displacements in protein ensembles with outlier-tolerant
probability distributions.
"""
import numpy

import csb.apps
import csb.bio.structure

from csb.bio.io.wwpdb import LegacyStructureParser
from csb.bio.utils import average_structure, fit, wfit
from csb.statistics.scalemixture import ScaleMixture, GammaPrior


class ExitCodes(csb.apps.ExitCodes):
    IO_ERROR = 2

class AppRunner(csb.apps.AppRunner):

    @property
    def target(self):
        return BFitApp

    def command_line(self):
        
        cmd = csb.apps.ArgHandler(self.program, __doc__)

        # Input structures
        cmd.add_positional_argument('pdb', str,
                                    'full path to the ensemble')

        # Optional arguments 
        cmd.add_scalar_option('chain', 'c', str,
                              'Chain',
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
        
        return cmd



class BFitApp(csb.apps.Application):
    """
    Python application for robust structure superposition of two protein structures
    """

    def main(self):
        try:
            parser = LegacyStructureParser(self.args.pdb)
            models = parser.models()

        except IOError as e:
            self.exit('PDB file parsing failed\n' + str(e.value), ExitCodes.IO_ERROR)

        if len(models) < 2:
            self.exit('PDB file contains only one model', ExitCodes.USAGE_ERROR)

        ensemble = parser.parse_models(models)
        X = numpy.array([model[self.args.chain].list_coordinates(['CA'], True) for model in ensemble])
        x_mu = average_structure(X)
        #n = X.shape[1]
        m = X.shape[0]
        R = numpy.zeros((m, 3, 3))
        t = numpy.ones((m, 3))


        prior = GammaPrior()
        mixture = ScaleMixture(scales=X.shape[1],
                               prior=prior, d=3)

        for i in range(m):
            R[i, :, :], t[i, :] = fit(x_mu, X[i])
        
        # gibbs sampling cycle
        for j in range(self.args.niter):
            # apply rotation
            data = numpy.array([numpy.sum((x_mu - numpy.dot(X[i], numpy.transpose(R[i])) - t[i]) ** 2, -1) ** 0.5
                                for i in range(m)]).T
            # sample scales
            mixture.estimate(data)
            # sample rotations
            for i in range(m):
                R[i, :, :], t[i, :] = wfit(x_mu, X[i], mixture.scales)


        out_ensemble = csb.bio.structure.Ensemble()

        for i, model in enumerate(ensemble):
            model.transform(R[i], t[i])
            out_ensemble.models.append(model)

        out_ensemble.to_pdb(self.args.outfile)
        


if __name__ == '__main__':
    AppRunner().run()
