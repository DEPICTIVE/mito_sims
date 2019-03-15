
import numpy as np
from scipy.optimize import bisect

# =============================================
# =============================================
# =============================================
# =============================================
class find:
    def __init__(self, f, xlims, args, nsamples=100, max_iter=500):
        '''
        Find the fixed points of a dynamic system

        Input :
        f :  function, model to compute the fixed points
        xlims : python list, the boundary in which to search for values
        args :  tuple, parameters for f
        nsamples : int, the number of partitions of the range specified in xlims.  In each partition the code will use the bisect method to search for roots of f by the bisect method.
        max_iter : max number of iterations
        '''
        self.xlims = xlims
        self.args = args
        self.nsamples = nsamples
        self.max_iter = max_iter
        self._find(f)

    # =============================================
    # =============================================
    def _find(self, f, xo=None):
        '''
        Numerically solve for the fixed points and estimate stability.  To ensure that all fixed points are found, solve for roots with nsamples randomly choosen initial values between range xlims.

        Input:
        f = model dydt function
        xlims = list or numpy array indicating bounds of x
        args : arguments passed to model dydt function
        '''
        # set partition values
        xo = np.linspace(self.xlims[0], self.xlims[1], self.nsamples)
        # instantiate fixed point list
        fp = []
        for w in range(1, self.nsamples):
            # compute the flux at points (xo[w], xo[w-1])
            dy = np.array([f(xo[w], self.args[0], self.args[1],
                           self.args[2], self.args[3]),
                           f(xo[w-1], self.args[0], self.args[1],
                           self.args[2], self.args[3])])
            # if there is a sign change then their must be a root between xo[w] and xo[w-1]
            if np.sum(np.abs(dy)) != np.abs(np.sum(dy)):
                # solve for the root using the bisect function in scipy
                fp += [bisect(f, xo[w-1], xo[w],
                        args=self.args)]
        # store fixed points
        self.fp = np.array(fp)
        # compute the stability of the fixed point
        self._get_stability(f)

    # =============================================
    # =============================================
    def _get_stability(self, f):
        '''
        Compute the stability of a fixed point by measuring the response of random small perturbations.
        '''
        # instantiate stability array
        self.stability = np.zeros(self.fp.size)
        # loop over fixed points
        for w in range(self.fp.size):
            # sample small perturbation of the independent variable
            x = self.fp[w] + 0.01*np.random.rand(10)
            # compute the flux for each perturbation
            y = np.zeros(x.size)
            for wy in range(x.size):
                y[wy] = f(np.abs(x[wy]), self.args[0], self.args[1],
                            self.args[2], self.args[3])
            # find the slope of the flux about the perturbations
            slope = compute_slope(x, y)
            # if the slope is less than 0 then the fp is stable, the fixed point is unstable
            if slope < 0:
                self.stability[w] = 1
            else:
                self.stability[w] = 0

# =============================================
# =============================================
# =============================================
# =============================================
def compute_slope(x, y):
    '''
    Compute slope using linear regression

    Input :
    x : numpy array of floats representing independent variable
    y : numpy array of floats representing dependent variable

    Return :
    float, representing the slope of the line that minimizes sum-squared error.
    '''
    c = np.cov(x,y)
    return c[0, 1] / c[0, 0]
