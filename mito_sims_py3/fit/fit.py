

import numpy as np
from scipy.optimize import fmin
from ..simulate import dynamics



# ====================================================
# SUM SQUARED ERRORS
# ====================================================

# Error
def chi(k, time, pars, s, fit_pars, y, sem):
    """
    k : ndarray
        parameters to infer, in order of fit_pars
    time : ndarray
        (N, ) array of time points
    pars : dictionary
        Model parameter values, see model.dynamics for required values
    fit_pars : list
        [parameter 1 to update, parameter 2 to update, ...]
    s : float
        stimulus strength
    """
    # update non pl parameters
    error = 0
    j = 0
    for w in fit_pars:
        pars[w] = k[j]
        # values must be greater than 0
        if k[j] < 0:
            error = 1e7
        j += 1
    # compute chi-square errors
    max_stim_chi = np.sum((y - dynamics(time, pars, s))**2 / sem) + error
    no_stim_chi = np.sum((pars['yo'] - dynamics(time, pars, 0))**2 / pars['yo'])
    return max_stim_chi + no_stim_chi + error


# ====================================================
# FITDATA
# ====================================================

def fit(time, y, sem, s, pars, fit_pars):
    """
    Compute the profile likilhood from fitting the dynamic model to dynamic data

    Input
    -----
    time : ndarray
        (N,) array of time points
    y : ndarray
        (N, ) array of y values to fit, represent mean over replicate set
    sem : ndarray
        (N, ) standard error of y for each time point
    s : float
        stimulus strength
    pars : dictionary
        model parameters, see model.dynamics for explanation
    fit_pars : list
        names of parameters to update in the fit

    Return
    ------
        pars : dictionary
            inferred parameters as specified by fit_pars and fopt
        warnflag : int
            1 if max fun evals, 2 if max num iterations made, in both failure to converge
    """
    # run nelder mead optimization
    out = fmin(chi,
            [pars[w] for w in fit_pars],
            args=(time, pars, s, fit_pars, y, sem),
            full_output=True,
            maxfun=2500,
            maxiter=2500,
            disp=False)
    # return solution
    outPars = {}
    for j in range(len(fit_pars)):
        outPars[fit_pars[j]] = [out[0][j]]
    outPars['chi'] = [out[1]]
    return [outPars,
            out[-1]]
