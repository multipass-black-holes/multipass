# BHMF inference code

This is the code to the paper

> Population analysis with split mass functions
>
> Yannick Ulrich, Djuna Croon, Jeremy Sakstein, ...
>
> arXiv:23xx.xxxxx

## Usage

You can clone the repository by using
```shell
# read only
$ git clone --recursive https://github.com/samueldmcdermott/bhmf.git
# read & write
$ git clone --recursive git@github.com:samueldmcdermott/bhmf.git
```

### Install requirements

You will need a reasonably modern fortran compiler (tested with gcc 10.2) and at least python 3.6.
You may want to use a virtual environment for the required python packags
```shell
$ python -m venv bhmf
$ source bhmf/bin/activate
(bhmf) $ pip install -r requirements.txt
```
Next, you need to obtain the data from LVC.
We use the [GW Open Science](https://www.gw-openscience.org/) API for this:
```shell
(bhmf) $ python get_data.py
```

### Compile the code
Go to the directory `multinest` and run
```shell
$ cd multinest/
$ make
```
This should produce the multinest executable `main`.

### Preparing the data
Before you can run our code, you need convert the data from the LVC format into our own.
To do this, run
```shell
(bhmf) $ python prepare.py
```
This will load and merge all event files as well as the injection and create `data.rec` and `inj.rec`.
These files are sufficient for inference and the full dataset from LVC is no longer required

### Running the code
To set up a job, you need to prepare a run card describing the model, the priors, and the MCMC parameters.
For example
```
m ppisn+planck+trivial
p  2.,  0.,  20.,  0.0, -4., -10., -7.0, -7.0,  20.,-4
P 10., 10., 120.,  0.5,  0.,   0., -0.3, -0.3, 150.,12.
n 100
t 0.5
e 0.3
o ppisn-h0/short
s
```
See below for a list of models and parameters.
You can now run the MCMC
```shell
$ mpirun -np <number of jobs> ./main /path/to/run.card
$ python analyse.py /path/to/output
```

## Models
 * `plp+flat+trivial+trivial`:
    powerlaw+peak for the primary, flat for the secondary, no redshift or spin.
    Parameters:
          $m_{\rm min}$,
          $\delta m$,
          $m_{\rm max}$,
          $\mu$,
          $\sigma$,
          $\alpha$,
          $\lambda_p$
 * `plp+plp+trivial+trivial`:
    powerlaw+peak for the primary and secondary, no redshift or spin.
    Primary and secondary are coupled through $q^\beta_q$.
    Parameters:
          $m_{\rm min}$,
          $\delta m$,
          $m_{\rm max}$,
          $\mu$,
          $\sigma$,
          $\alpha$,
          $\lambda_p$,
          $\beta_q$
 * `plp+pow+trivial+trivial`:
    powerlaw+peak for the primary, power law for the secondary, no redshift or spin.
    Parameters:
          $m_{\rm min}$,
          $\delta m$,
          $m_{\rm max}$,
          $\mu$,
          $\sigma$,
          $\alpha$,
          $\lambda_p$,
          $k$
 * `plp+flat+planck+trivial`:
    powerlaw+peak for the primary, flat for the secondary, fitting for $H_0$ but no spin.
    Parameters:
          $m_{\rm min}$,
          $\delta m$,
          $m_{\rm max}$,
          $\mu$,
          $\sigma$,
          $\alpha$,
          $\lambda_p$,
          $H_0$
 * `plp+pow+trivial+trivial`:
    powerlaw+peak for the primary, power law for the secondary, fitting for $H_0$ but no spin.
    Parameters:
          $m_{\rm min}$,
          $\delta m$,
          $m_{\rm max}$,
          $\mu$,
          $\sigma$,
          $\alpha$,
          $\lambda_p$,
          $k$,
          $H_0$
 * `plp+flat+trivial+beta`:
    powerlaw+peak for the primary, flat for the secondary, no redshift and $\beta$ distribution for spin.
    Parameters:
          $m_{\rm min}$,
          $\delta m$,
          $m_{\rm max}$,
          $\mu$,
          $\sigma$,
          $\alpha$,
          $\lambda_p$,
          $\alpha$,
          $\beta$
 * `plp+pow+trivial+beta`:
    powerlaw+peak for the primary, power law for the secondary, no redshift and $\beta$ distribution for spin.
    Parameters:
          $m_{\rm min}$,
          $\delta m$,
          $m_{\rm max}$,
          $\mu$,
          $\sigma$,
          $\alpha$,
          $\lambda_p$,
          $k$,
          $\alpha$,
          $\beta$

 * `ppisn+flat+trivial+trivial`:
    PPISN for the primary, flat for the secondary, no spin or redshift
    Parameters:
          $m_{\rm min}$,
          $\delta m$,
          $m_{\rm BHMG}$,
          $a$,
          $b$,
          $d$,
          $\log_{10}\lambda_{21}$,
          $\log_{10}\lambda_{12}$
 * `ppisn+trivial+trivial`:
    PPISN for the primary, physical for the secondary, no spin or redshift.
    Primary and secondary are coupled using $q^\beta_q$,
    Parameters:
          $m_{\rm min}$,
          $\delta m$,
          $m_{\rm BHMG}$,
          $a$,
          $b$,
          $d$,
          $\log_{10}\lambda_{21}$,
          $\log_{10}\lambda_{12}$,
          $\beta_q$
 * `ppisn2P+trivial+trivial`:
    PPISN for the primary, physical for the secondary, no spin or redshift.
    Primary and secondary are coupled using $q^\beta_q$ with two different $\beta$,
    Parameters:
          $m_{\rm min}$,
          $\delta m$,
          $m_{\rm BHMG}$,
          $a$,
          $b$,
          $d$,
          $\log_{10}\lambda_{21}$,
          $\log_{10}\lambda_{12}$,
          $\beta_q^{(0)}$,
          $\beta_q^{(1)}$
 * `ppisn+planck+trivial`:
    PPISN for the primary, physical for the secondary, fitting for $H_0$ but not for spin.
    Primary and secondary are coupled using $q^\beta_q$ with two different $\beta$,
    Parameters:
          $m_{\rm min}$,
          $\delta m$,
          $m_{\rm BHMG}$,
          $a$,
          $b$,
          $d$,
          $\log_{10}\lambda_{21}$,
          $\log_{10}\lambda_{12}$,
          $\beta_q$,
          $H_0$
 * `ppisn+trivial+beta`:
    PPISN for the primary, physical for the secondary, no redshift and $\beta$ distributions for spin.
    Primary and secondary are coupled using $q^\beta_q$ with two different $\beta$,
    Parameters:
          $m_{\rm min}$,
          $\delta m$,
          $m_{\rm BHMG}$,
          $a$,
          $b$,
          $d$,
          $\log_{10}\lambda_{21}$,
          $\log_{10}\lambda_{12}$,
          $\beta_q$,
          $\alpha_1$,
          $\beta_1$,
          $\alpha_2$,
          $\beta_2$
