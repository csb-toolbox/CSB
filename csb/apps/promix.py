"""
ProMix: Take a protein structure ensemble and find a mixture of rigid
segments or a mixture of conformers. Writes K copies of the ensemble
(for segments) or K subsets of the ensemble (for conformers) as PDB
files, each superposed on different components.

Reference: Hirsch M, Habeck M. - Bioinformatics. 2008 Oct 1;24(19):2184-92
"""

import numpy

import csb.apps
import csb.bio.structure

from csb.bio.io.wwpdb import LegacyStructureParser
from csb.statistics import mixtures


class ExitCodes(csb.apps.ExitCodes):
    IO_ERROR = 2

class AppRunner(csb.apps.AppRunner):

    @property
    def target(self):
        return ProMixApp

    def command_line(self):
        cmd = csb.apps.ArgHandler(self.program, __doc__)

        cmd.add_scalar_option('components', 'K', int, 'Number of components', -1)
        cmd.add_scalar_option('type', 't', str, 'Type of mixture', 'segments', ('segments', 'conformers'))
        cmd.add_positional_argument('infile', str, 'input PDB file')

        return cmd

    def initapp(self, args):
        app = self.target
        return app(args)

class ProMixApp(csb.apps.Application):

    def main(self):
        try:
            parser = LegacyStructureParser(self.args.infile)
            models = parser.models()
        except:
            self.exit('PDB file parsing failed', ExitCodes.IO_ERROR)

        if len(models) < 2:
            self.exit('PDB file contains only one model', ExitCodes.USAGE_ERROR)

        ensemble = parser.parse_models(models)
        X = numpy.array([model.get_coordinates(['CA'], True) for model in ensemble])

        if self.args.type == 'segments':
            self.main_segments(ensemble, X)
        elif self.args.type == 'conformers':
            self.main_conformers(ensemble, X)
        else:
            raise ValueError('type must be "segments" or "conformers"')

    def main_segments(self, ensemble, X):

        mixture = mixtures.SegmentMixture.new(X, self.args.components)
        self.log('Number of segments: {0}'.format(mixture.K))

        for k,(sigma,w) in enumerate(zip(mixture.sigma, mixture.w)):
            outfile = 'promix_segment_{0}.pdb'.format(k+1)
            self.log('  {0}: sigma = {1:6.3f}, w = {2:.3f}, file = {3}'.format(k+1, sigma, w, outfile))

            for model, R, t in zip(ensemble, mixture.R, mixture.t):
                if k > 0:
                    model.transform(R[k-1], t[k-1])
                R = R[k].T
                t = -numpy.dot(R, t[k])
                model.transform(R, t)

            ensemble.to_pdb(outfile)

    def main_conformers(self, ensemble, X):

        mixture = mixtures.ConformerMixture.new(X, self.args.components)
        self.log('Number of conformers: {0}'.format(mixture.K))

        membership = mixture.membership

        for k,(sigma,w) in enumerate(zip(mixture.sigma, mixture.w)):
            outfile = 'promix_conformer_{0}.pdb'.format(k+1)
            self.log('  {0}: sigma = {1:6.3f}, w = {2:.3f}, file = {3}'.format(k+1, sigma, w, outfile))

            ek = csb.bio.structure.Ensemble()

            for model, R, t, mk in zip(ensemble, mixture.R, mixture.t, membership):
                if mk != k:
                    continue
                R = R[k].T
                t = -numpy.dot(R, t[k])
                model.transform(R, t)
                ek.models.append(model)

            ek.to_pdb(outfile)


def main():
    AppRunner().run()
    
    
if __name__ == '__main__':
    main()