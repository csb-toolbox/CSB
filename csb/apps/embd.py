"""
Sharpening of EM maps by non-negative blind deconvolution.
For details see:

Hirsch M, Schoelkopf B and Habeck M (2010)
A New Algorithm for Improving the Resolution of Cryo-EM Density Maps.
"""

import os
import numpy
import csb.apps

from numpy import sum, sqrt

from csb.numeric import convolve, correlate, trim
from csb.bio.io.mrc import DensityMapReader, DensityMapWriter, DensityInfo, DensityMapFormatError


class ExitCodes(csb.apps.ExitCodes):
    
    IO_ERROR = 2
    INVALID_DATA = 3
    ARGUMENT_ERROR = 4


class AppRunner(csb.apps.AppRunner):
    
    @property
    def target(self):
        return DeconvolutionApp
     
    def command_line(self):
        
        
        cmd = csb.apps.ArgHandler(self.program, __doc__)
        
        cmd.add_scalar_option('psf-size', 's', int, 'size of the point spread function', default=15)
        cmd.add_scalar_option('output', 'o', str, 'output directory of the sharpened maps', default='.')
        cmd.add_scalar_option('iterations', 'i', int, 'number of iterations', default=1000)
        cmd.add_scalar_option('output-frequency', 'f', int, 'create a map file each f iterations', default=50)
        cmd.add_boolean_option('verbose', 'v', 'verbose mode')
                        
        cmd.add_positional_argument('mapfile', str, 'Input Cryo EM file in CCP4 MRC format')
                        
        return cmd
    

class DeconvolutionApp(csb.apps.Application):
    
    def main(self):
        
        if not os.path.isfile(self.args.mapfile):
            DeconvolutionApp.exit('Input file not found.', code=ExitCodes.IO_ERROR)
                    
        if not os.path.isdir(self.args.output):
            DeconvolutionApp.exit('Output directory does not exist.', code=ExitCodes.IO_ERROR)
    
        if self.args.psf_size < 1:
            DeconvolutionApp.exit('PSF size must be a positive number.', code=ExitCodes.ARGUMENT_ERROR)
                         
        if self.args.iterations < 1:
            DeconvolutionApp.exit('Invalid number of iterations.', code=ExitCodes.ARGUMENT_ERROR)
            
        if self.args.output_frequency < 1:
            DeconvolutionApp.exit('Output frequency must be a positive number.', code=ExitCodes.ARGUMENT_ERROR)
            
        if self.args.iterations < self.args.output_frequency:
            DeconvolutionApp.exit('Output frequency is too low.', code=ExitCodes.ARGUMENT_ERROR)            
            
        self.args.output = os.path.abspath(self.args.output)
                                
        self.run()
    
    def run(self):

        writer = DensityMapWriter()

        self.log('Reading input density map...')
        try:
            input = DensityMapReader(self.args.mapfile).read()
            embd = Deconvolution(input.data, self.args.psf_size)
            
        except DensityMapFormatError as e:
            msg = 'Error reading input MRC file: {0}'.format(e)  
            DeconvolutionApp.exit(msg, code=ExitCodes.INVALID_DATA)

        self.log('Running {0} iterations...'.format(self.args.iterations))
        self.log(' Iteration             Loss Correlation  Output')
                
        for i in range(1, self.args.iterations + 1):
            embd.run_once()

            if i % self.args.output_frequency == 0:
                output = OutputPathBuilder(self.args, i)
                
                density = DensityInfo(embd.data, None, None, header=input.header)
                writer.write_file(output.fullpath, density)
                
                self.log('{0:>9}.  {1:15.2f}  {2:10.4f}  {3}'.format(
                                    i, embd.loss, embd.correlation, output.filename))

        self.log('Done: {0}.'.format(output.fullpath))
                            
    def log(self, *a, **k):
        
        if self.args.verbose:
            super(DeconvolutionApp, self).log(*a, **k)
      

class OutputPathBuilder(object):
    
    def __init__(self, args, i):
        
        basename = os.path.basename(args.mapfile)
        file, extension = os.path.splitext(basename)        
        
        self._newfile = '{0}.{1}{2}'.format(file, i, extension)
        self._path = os.path.join(args.output, self._newfile)
        
    @property
    def fullpath(self):
        return self._path
    
    @property
    def filename(self):
        return os.path.basename(self._newfile)
      
class Util(object):
    
    @staticmethod                    
    def corr(x, y, center=False):
    
        if center:
            x = x - x.mean()
            y = y - y.mean()
    
        return sum(x * y) / sqrt(sum(x * x)) / sqrt(sum(x * x))

class Deconvolution(object):
    """
    Blind deconvolution for n-dimensional images.
    
    @param data: EM density map data (data field of L{csb.bio.io.mrc.DensityInfo})
    @type data: array
    @param psf_size: point spread function size
    @type psf_size: ints
    @param beta_x: hyperparameters of sparseness constraints
    @type beta_x: float
    @param beta_f: hyperparameters of sparseness constraints
    @type beta_f: float
    """
    
    def __init__(self, data, psf_size, beta_x=1e-10, beta_f=1e-10, cache=True):

        self._f = []
        self._x = []
        self._y = numpy.array(data)
        self._loss = []
        self._corr = []
        
        self._ycache = None
        self._cache = bool(cache)

        self._beta_x = float(beta_x)
        self._beta_f = float(beta_f)
                
        shape_psf = (psf_size, psf_size, psf_size)
        self._initialize(shape_psf)
        
    @property
    def beta_x(self):
        return self._beta_x

    @property
    def beta_f(self):
        return self._beta_f
    
    @property
    def loss(self):
        """
        Current loss value.
        """        
        if len(self._loss) > 0:
            return float(self._loss[-1])
        else:
            return None
        
    @property
    def correlation(self):
        """
        Current correlation value.
        """
        if len(self._corr) > 0:
            return float(self._corr[-1])
        else:
            return None
        
    @property
    def data(self):
        return trim(self._x, self._f.shape)
            
    def _initialize(self, shape_psf):
        """
        Initialize with flat image and psf.
        """
        self._f = numpy.ones(shape_psf)
        self._x = numpy.ones(numpy.array(self._y.shape) + numpy.array(shape_psf) - 1)

        self._normalize_psf()
                
    def _normalize_psf(self):
        self._f /= self._f.sum()
        
    def _calculate_image(self):
        return convolve(self._f, self._x)

    def calculate_image(self, cache=False):

        if cache and self._ycache is not None:
            return self._ycache
        else:
            y = self._calculate_image()
            if self._cache:
                self._ycache = y
            return y

    def _update_map(self):

        y = self.calculate_image()

        N = correlate(self._f, self._y) - self.beta_x
        D = correlate(self._f, y)

        self._x *= numpy.clip(N, 1e-300, 1e300) / numpy.clip(D, 1e-300, 1e300)

    def _update_psf(self):

        y = self.calculate_image()

        N = correlate(self._x, self._y) - self.beta_f
        D = correlate(self._x, y)

        self._f *= numpy.clip(N, 1e-300, 1e300) / numpy.clip(D, 1e-300, 1e300)
        self._normalize_psf()

    def eval_loss(self, cache=False):

        y = self.calculate_image(cache=cache)
        
        return 0.5 * ((self._y - y) ** 2).sum() + \
                + self.beta_f * self._f.sum() + self.beta_x * self._x.sum()

    def eval_corr(self, cache=False):
        
        y = self.calculate_image(cache=cache)
        return Util.corr(self._y, y)

    def run_once(self):
        """
        Run a single iteration.
        """

        self._loss.append(self.eval_loss(cache=True))
        self._corr.append(self.eval_corr(cache=True))

        self._update_map()                
        self._update_psf()
            
    def run(self, iterations):
        """
        Run multiple iterations.
        
        @param iterations: number of iterations to run
        @type iterations: int 
        """
        for i in range(iterations):
            self.run_once()


def main():
    AppRunner().run()
    
    
if __name__ == '__main__':
    main()