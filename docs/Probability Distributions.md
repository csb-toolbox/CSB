## Probability Density Functions

The _csb.statistics.pdf_ module defines _AbstractDensity_: a common interface 
for all PDFs. Each _AbstractDensity_ describes a specific type of probability
distribution, for example _Normal_ is an implementation of the Gaussian 
distribution:

```python
>>> pdf = Normal(mu=10, sigma=1.1)
>>> pdf.mu, pdf['sigma']('sigma')
10.0, 1.1
``` 

Every PDF provides an implementation of the _AbstractDensity.evaluate()_ 
method, which evaluates the PDF for a list of input data points:

```python
>>> pdf.evaluate([10, 9, 11, 12](10,-9,-11,-12))
array([ 0.3626748 ,  0.2399147 ,  0.2399147 ,  0.06945048](-0.3626748-,--0.2399147-,--0.2399147-,--0.06945048))
``` 

PDF instances also behave like functions:

```python
>>> pdf(data)    # the same as pdf.evaluate(data)
``` 

Some AbstractDensity implementations may support drawing random numbers 
from the distribution (or raise an exception otherwise):

```python
>>> pdf.random(2)
array([ 9.86257083,  9.73760515](-9.86257083,--9.73760515))
``` 

Each implementation of AbstractDensity may support infinite number of 
estimators, used to estimate and re-initialize the PDF parameters from a 
set of observed data points:

```python
>>> pdf.estimate([5, 5, 10, 10](5,-5,-10,-10))
>>> pdf.mu, pdf.sigma
(7.5, 2.5)
>>> pdf.estimator
<csb.statistics.pdf.GaussianMLEstimator>
``` 

Estimators implement the _AbstractEstimator_ interface. They are treated 
as pluggable tools, which can be exchanged through the 
_AbstractDensity.estimator_ property (you could create, initialize and plug in 
your own estimator as well). This is a classic Strategy.

