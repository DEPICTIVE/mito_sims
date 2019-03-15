MITO_SIM: A coarse grained model for simulating initiator Caspase dynamics during TRAIL induced apoptosis
===========================================================

This package provides simple tools to simulate the dynamics of apoptosis and to fit model parameters to data.  Among these tools includes:

1. deterministic dynamic simulator,
2. stability analysis: finding fixed points and assessing their stability,
3. stochastic dynamic simulator,
4. parallelization for many experiments, and
5. plotting tools.

To see how these tools are used please refer to the our wiki.

Installation
------------

### Dependencies

To run simulations you will need:

- Python (>=3.5)
- Numpy (>= 1.12.0)
- SciPy (>= 0.18.1)
- matplotlib (>= 2.0.0)

We recommend using the MITO_SIM package in a virtual environment.  Which can simply be done by:

```bash
$ python -m venv env
```

where `env` can be replaced with any name that you desire.  Then to activate,

```bash
$ source env/bin/activate
(env) $
```

Installing the MITO_SIM package dependencies in the virtual environemnt is easy with the `requirements.txt` file:

```bash
(env) $ pip install -r requirements.txt
```

Then lastly migrate to the parent directory in which you have saved the MITO_SIM package:

```bash
(env) $ pip install mito_sims_py3/
```

and you are done!

Using the package
-----------------

Check out the examples in our wiki!


Citation
--------

Please don't forget to cite us at:

```
  @article {2017mito_sims,
    author = {Santos, L.C. and Vogel, R. and Chipuk, J.E. and
      Birtwistle, M.R. and Stolovitzky, G. and Meyer, P.},
    title = {Origins of fractional control in regulated cell death},
    year = {2017},
    doi = {10.1101/201160},
    publisher = {Cold Spring Harbor Laboratory},
    URL = {https://www.biorxiv.org/content/early/2017/10/10/201160},
    eprint = {https://www.biorxiv.org/content/early/2017/10/10/201160.full.pdf},
    journal = {bioRxiv}
  }
```
