
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable as csm

import fixed_points as fpoints
import model

# ================================================================
# ================================================================
# bifurcation diagrams
# ================================================================
# ================================================================

def bifurcation_diagram(pars, s, bbox, bi=0.,
            nPar=100, ylims=[0, 1.], sname=None):
    '''
    Bifurcation diagram for changing mitochondria levels

    Input :
    stable_points
    '''
    # compute stable and unstable fixed points
    mito = np.logspace(np.log10(bbox[0]), np.log10(bbox[1]), nPar)
    stable, unstable = _get_organized_fixed_points(mito, pars, s, ylims, bi=bi)

    # prepare arrorws
    # set arrows posotions
    Xm, Xsig = np.meshgrid(np.logspace(np.log10(bbox[0]),
                                       np.log10(bbox[1]), 7),
                           np.linspace(ylims[0],ylims[1], 7))
    # magnitude of arrow in the x-axis direction, which is of course 0.
    U = np.zeros(shape=Xm.shape)

    # magnitude of the arrow in the y-axis direction
    V = []
    for wm in Xm[0, :]:
        tmp = []
        for wsig in Xsig[:, 0]:
            tmp += [model.dydt(wsig, pars, s, bi, wm)]
        V += [tmp]
    V = np.array(V).transpose()

    # plot
    fig = plt.figure(figsize=(4.5, 3.5))
    plt.quiver(Xm, Xsig, U, V / np.abs(V), headwidth=10, headlength=6,
               alpha=0.5, width=0.005, pivot='tip',
               scale=15)
    for w in range(stable.shape[1]):
        plt.plot(mito, stable[:, w], '.', color='k')
    for w in range(unstable.shape[1]):
        plt.plot(mito, unstable[:, w], 'o', mec='k', mfc='none', ms=5)
    plt.ylim(ylims)
    plt.xlabel(r'$\rho_{Mitochondria}$', fontsize=15)
    plt.ylabel(r'y', fontsize=15)
    ax = plt.gca()
    ax.set_position([0.175, 0.175, 0.8, 0.725])
    plt.xscale('log')
    save(fig, sname)
    plt.show(block=False)


# ================================================================
# ================================================================
# stable, unstable = _get_organized_fixed_points(mito, pars, s, ylims, bi=bi)
def _get_organized_fixed_points(mito, pars, s, ylims, bi=0):
    # compute the fixed points for each mito value
    # in addition, store the maximum number of fixed points
    n_stable = 0
    n_unstable = 0
    fp = []
    for wm in mito:
        # store the fixed point class
        fp += [fpoints.find(model.dydt, ylims,
                    (pars, s, bi, wm), nsamples=20)]
        nfp = len(fp[-1].fp)
        sum_stable = int(np.sum(fp[-1].stability))
        # find the maximum number of stable fixed points
        if sum_stable > n_stable:
            n_stable = sum_stable
        # compute the maximum number of unstable fixed points
        if nfp - sum_stable > n_unstable:
            n_unstable = nfp - sum_stable
    # instantiate the stable fixed point array
    stable = np.zeros(shape=(mito.size, n_stable)) + np.nan
    # instantiate the unstble fixed point array
    unstable = np.zeros(shape=(mito.size, n_unstable)) + np.nan

    for wfp in range(len(fp)):
        i_stable = 0
        i_unstable = 0
        for w in range(len(fp[wfp].fp)):
            if fp[wfp].stability[w] == 1:
                stable[wfp, i_stable] = fp[wfp].fp[w]
                i_stable += 1
            else:
                unstable[wfp, i_unstable] = fp[wfp].fp[w]
                i_unstable += 1
    return [stable, unstable]

# ================================================================
# ================================================================
# ================================================================
# ================================================================

def save(fig, sname):
    if sname is not None:
        fig.savefig(sname, fmt=sname.split('.')[-1])


def get_coluer(N, cmap='viridis'):
    cm = csm(cmap=cmap)
    return cm.to_rgba(range(N))

# ================================================================
# ================================================================
# deterministic dynamics
# ================================================================
# ================================================================

def dynamics(pars, s=1., bi=0., m=1.,
        sname=None, Nt = 500, fsize=(4.5, 4),
        cmap='viridis', iMax = 100):
    time = np.linspace(pars['to'], pars['tf'], Nt)
    int_count = 0
    if (type(s) == int) | (type(s) == float):
        s = [s]
        int_count +=1
    if (type(bi) == int) | (type(bi) == float):
        bi = [bi]
        int_count += 1
    if (type(m) == int) | (type(m) == float):
        m = [m]
        int_count += 1

    if (int_count >= 2) & (sum((len(m), len(bi), len(s))) < iMax+2):
        if int_count == 3:
            coluer = get_coluer(10, cmap='tab10')
        else:
            coluer = get_coluer(sum((len(m), len(bi), len(s)))-2,
                cmap=cmap)
        fig = plt.figure(figsize=fsize)
        count = 0
        for ws in s:
            for wb in bi:
                for wm in m:
                    plt.plot(time,
                        model.dynamics(time, pars, ws, b_i=wb, m=wm),
                        '-', linewidth=2.5,
                        color=coluer[count, :])
                    count += 1
        plt.xlabel(r'Time', fontsize=15)
        plt.ylabel('Active IC', fontsize=15)
        plt.ylim(0, 1)
        plt.tight_layout()
        plt.show(block=False)
        save(fig, sname)
    else:
        print '================================'
        print '1)  You are trying to plots sequences of {} variables'.format(3 - int_count)
        print '2)  You are trying to plot a sequence with more than {} iterations'.format(iMax)
        print '================================'
        print 'For an interpretable plot you can only plot a sequence of 1 variable whose length is less than or equal to iMax.  If you would like to plot a sequence of more than iMax entries, please include in input.'

# ================================================================
# ================================================================
#  Stochastic dynamics
# ================================================================
# ================================================================
# need to parallelize



def stochastic(pars, N=10, s=1., bi=0., m=1.,
        Nt = 500, cmap='tab10', n_max = 100,
        fsize=(4.5, 4), sname=None):
        '''
        Simulate the stochastic dynamics of IC in N cells and plot
        Input:
        pars : python dictionary of parameter values
        N : int, number of virtual cells to simulate

        '''
        if N <= n_max:
            time = np.linspace(pars['to'], pars['tf'], Nt)
            coluer = get_coluer(N,
                    cmap=cmap)

            x = model.parallel_sims(time, pars, s, N, b_i=bi, m=m)

            fig = plt.figure(figsize=fsize)
            for w in range(N):
                plt.plot(time, x[w, :],
                    '-', linewidth=2.5,
                    color=coluer[w, :])
            plt.xlabel(r'Time', fontsize=15)
            plt.ylabel('Active IC', fontsize=15)
            plt.ylim(0, 1)
            plt.tight_layout()
            save(fig, sname)
            plt.show(block=False)
        else:
            print 'Simulating N cells greather than {} will be much slower to process.  If you are sure you want to simulate this many cells please explicitly change n_max in the function arguments.'.format(n_max)


# ================================================================
# ================================================================
# Phase plane
# ================================================================
# ================================================================

def phase_plane(pars, s, bi=0, m=1., cmap='viridis',
        fsize=(4.5, 4), yopt=None, sname=None):
    # open circles for unstable fixed points, closed circles for stable
    stability_markers = [{'marker':'o', 'ms':15, 'mfc':'w', 'mec':'k'},
                    {'marker':'o', 'ms':12.5, 'mfc':'k', 'mec':'none'}]
    # set ymin to slightly negative so we can find the zeroth fixed point
    ymin = -0.01
    # find the fixed points
    fp = fpoints.find(model.dydt, [0, 1.], (pars, s, bi, m))
    # find optimal plotting range
    # compute y values and dydt
    if yopt is None:
        if len(fp.fp) > 1:
            y = np.linspace(ymin, pars['yt'], 250)
            dydt = np.array([model.dydt(wy, pars, s, bi, m) for wy in y])
            # expandd high fixed point and compute slope
            xtilde = np.max(fp.fp) + 0.025*np.random.randn(10)
            ytilde = [model.dydt(wy, pars, s, bi, m) for wy in xtilde]
            # slope is computed from the covariance matrix
            c = np.cov(xtilde, ytilde)
            # yopt is such that |dydt| <= |dydt[y <= fp_max]|
            slope = c[0, 1] / c[0, 0]
            yint = np.mean(ytilde) - np.mean(xtilde) * slope
            yopt = -(np.max(np.abs(dydt[y < np.max(fp.fp)])) + yint) / slope
        else:
            yopt = pars['yt']
    # using yopt generate plotting values
    y = np.linspace(ymin, yopt, 250)
    dydt = [model.dydt(wy, pars, s, bi, m) for wy in y]

    fig = plt.figure(figsize=fsize)
    plt.plot([0, 1.], [0, 0], '-', color='k', linewidth=1)
    plt.plot(y, dydt, '-', color='k', linewidth=2.5)
    for w in range(fp.fp.size):
        idx = int(fp.stability[w])
        plt.plot(fp.fp[w], 0, stability_markers[idx]['marker'],
                ms=stability_markers[idx]['ms'],
                mec=stability_markers[idx]['mec'],
                mfc=stability_markers[idx]['mfc'])
    plt.xlabel('y', fontsize=15)
    plt.ylabel('dy/dt', fontsize=15)
    plt.xlim([ymin, yopt])
    plt.tight_layout()
    save(fig, sname)
    plt.show(block=False)
