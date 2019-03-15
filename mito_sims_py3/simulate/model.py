from multiprocessing import Pool

import numpy as np
from scipy.optimize import bisect


# ========================================
# ========================================
# deterministic dynamics
# ========================================
# ========================================

def dynamics(time, par, s, b_i=0., m=1.):
    '''
    Fourth order Runge-Kutta integrator of the dynamic equation.

    Input:
    time : numpy array of timepoints to evaluate
    par : dictionary of parameters
        y0 : float, initial abundance of y
        yt : float, effective max concenctration of caspase
        n : int, average number of monomers per pore
        r : float, effective spontaneous production rate of IC
        mu : float, effective feedback strenth from EC
        gamma : float, effective autocatalysis and TRAIL induction amplitude
    s : float, stimulus strength, s \in [0, 1]
    b_i : float, abundance of bcl-2 inhibitor normalized to half-max constant, i.e. b_i = 1. is assumed to be the half-max constant.
    m : float, fold change over mean mitochondria abundance, <m> = 1.

    Returns:
    numpy array of length time.size, for which each element is the abundance of y at each time point
    '''
    # instantiate the numpy array of IC abundance y
    y = np.zeros(time.size)
    # initial abundance of y is from input parameters
    y[0] = par['yo']
    # compute y for each time point
    for w in range(1, time.size):
        # compute delta t
        dt = time[w] - time[w-1]
        # first order approximation
        k1 = dt * dydt(y[w-1], par, s, b_i, m)
        # second order approximation
        k2 = dt * dydt(y[w-1] + 0.5*k1, par, s, b_i, m)
        # third order approximation
        k3 = dt * dydt(y[w-1] + 0.5*k2, par, s, b_i, m)
        # fourth order approximation
        k4 = dt * dydt(y[w-1] + k3, par, s, b_i, m)
        # combine the approximation of each order to the abundance of y
        y[w] = y[w-1] + k1/6. + k2/3. + k3/3. + k4/6.
    return y



# ========================================
# ========================================

def dydt(y, par, s, b_i, m):
    '''
    Compute net flux given the abundance of y, kinetic parameters, stimulus strength, the abundance of Bcl-2 inhibitor, and the fold change of mitochondria abundance.

    Input:
    y : numpy array, current evaluation of y
    par : python dictionary of parameters, see dynamics functions
    s : float, stimulus stength, defined as in dynamics function
    b_i : float, defined as in dynamics function
    m : flaot, defined as in dynamics function

    Returns:
    float, of the net flux of y
    '''
    return _fp(y, par, s, b_i, m) - _fn(y)

# ========================================
# ========================================

def _fp(y, par, s, b_i, m):
    '''
    Positive flux

    Input:
    y : numpy array, current evaluation of y
    par : python dictionary of parameters, see dynamics functions
    s : float, stimulus stength, defined as in dynamics function
    b_i : float, defined as in dynamics function
    m : flaot, defined as in dynamics function

    Returns:
    float of positive flux of y
    '''
    # compute number of bax monomers
    x = _get_x(y, b_i, m, par)
    # compute the positive flux
    return (par['yt'] - y) * (par['r'] + s*par['gamma']*y + par['mu']*m*(x/m)**par['n'])

# ========================================
# ========================================

def _fn(y):
    '''
    Negative flux

    Input:
    y : float, representing the abundance of effective IC

    Return:
    float, represents the negative flux of y
    '''
    return y

# ========================================
# ========================================

def _get_x(y, b_i, m, par):
    '''
    Use conservation of Bcl-2 and mitochondria associated Bax to numerically solve for the abundance of Bax monomers.

    Input:
    y : scalar, the effective abundance of IC
    b_i : float, defined as in dynamics function
    m : float, defined as in dynamics function
    par : python dictionary of parameters, see dynamics functions

    Return:
    float, the abundance of Bax monomers
    '''
    x_mito = y / par['alpha']
    try:
        out = bisect(_mass_conserve, 0, x_mito,
                    args=(x_mito, b_i, m, par))
    except ValueError:
        print('b_i : {:0.3}'.format(b_i))
        print('rho : {:0.3}'.format(m))
        print('y : {:0.3}'.format(x_mito))
        print('Using absolute value')
        out = bisect(_mass_conserve, 0, np.abs(x_mito),
                    args=(np.abs(x_mito), b_i, m, par))
        print(out)
    return out

# ========================================
# ========================================

def _mass_conserve(x, x_mito, b_i, m, par):
    '''
    The conservation of mass equation for the numerical solver.

    Input :
    x : float, the amount of bax monomers
    x_mito : float, the total amount of Bax on the mitochondria surface
    b_i : float, defined as in dynamics function
    m : float, defined as in dynamics function
    par : python dictionary of parameters, see dynamics functions
    '''
    # bcl2 sequestration of bax
    xt = par['bt'] * x / (par['kbb']*m*(1 + b_i) + x)
    # d.b. bax polymerization
    for w in range(1, par['n']+1):
        xt += m * w * (x/m)**w
    return x_mito - xt


# ========================================
# ========================================
# Stochastic simulations
# ========================================
# ========================================

def noise(y, par, s, b_i, m, dt):
    '''
    Compute the magnitude of the noise from the Langevin equation.

    Input:
    y : float, the effective abundance of IC
    par : python dictionary of parameters, see dynamics functions
    s : stimulation stength
    b_i : float, defined as in dynamics function
    m : float, defined as in dynamics function
    dt : float, the width of time interval of simulation

    Returns:
    float, representing the standard deviation of the fluctuating force term
    '''
    return par['Omega']*np.sqrt((_fp(y, par, s, b_i, m) + _fn(y))*dt)

# ========================================
# ========================================

def _update(y, par, s, b_i, m, dt):
    '''
    Time step updater using modiefied Runge-Kutta method.  The drift term updates as 4 order Runge-Kutta.

    Input :
    y : float, the effective abundance of IC
    par : python dictionary of parameters, see dynamics functions
    s : stimulation stength
    b_i : float, defined as in dynamics function
    m : float, defined as in dynamics function
    dt : float, the width of time interval of simulation

    Returns:
    float, representing the next time point stochastic prediction for y
    '''
    # random noise term
    eta = np.random.randn()
    k1 = dt * dydt(y, par, s, b_i, m)
    k2 = dt * dydt(y+0.5*k1, par, s, b_i, m)
    k3 = dt * dydt(y+0.5*k2, par, s, b_i, m)
    k4 = dt * dydt(y+k3, par, s, b_i, m)
    e1 = noise(y, par, s, b_i, m, dt) * eta
    e2 = noise(y+0.5*k1, par, s, b_i, m, dt) * eta
    e3 = noise(y+0.5*k2, par, s, b_i, m, dt) * eta
    e4 = noise(y+k3, par, s, b_i, m, dt) * eta
    return y + (k1+e1)/6 + (k2+e2)/3 + (k3+e3)/3 + (k4+e4)/6

# ========================================
# stochastic simulator
# ========================================


def stochastic(time, par, s, b_i=0, m=1.):
    '''
    Simulate apoptosis by chemical Langevin equation.
    Input:
    time : numpy array
    par : dictionary of model parameters
    s : stimulus strength [0, 1]
    b_i : abundance of BCL2 inhibitor [0, \infty)
    m : mitochondria abundance (0, \infty)

    Return:
    y : numpy array of stochastic simulation results
    '''
    # instantiate the numpy array of the effective IC, y
    y = np.zeros(time.size)
    # the initial abundance of y is an input parameter
    y[0] = par['yo']
    # iterate over each time point
    for w in range(1, time.size):
        # time interval
        dt = time[w] - time[w-1]
        # compute the y at next time point from the stochastic simulation
        tmp = _update(y[w-1], par, s, b_i, m, dt)
        # apply the boundary condition that y >= 0.
        if (tmp < 0) | (np.isnan(tmp) == True) | (tmp > par['yt']):
            y[w] = y[w-1]
        else:
            y[w] = tmp
    return y

# ========================================
# parallelize simulations
# ========================================

def _call_sims(args):
    np.random.seed(args['counter'])
    return stochastic(args['time'], args['pars'], args['s'],
            b_i=args['bi'], m=args['m'])

def parallel_sims(time, pars, s, N, b_i=0., m=1.):
    # prepare arguments dictionary
    args = []
    for wn in range(N):
        args += [{'time':time, 'pars':pars, 's':s,
                    'bi':b_i, 'm':m, 'counter':wn}]
    # instantiate pool object
    pool = Pool()
    # run simulations in parallel
    xsim = pool.map(_call_sims, args)
    pool.close()
    return np.array(xsim)
