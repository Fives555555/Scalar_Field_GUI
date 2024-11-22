# Scientific Programming in Python – submission 3

### Project title: $\Lambda$ CDM cosmology + a single axion scalar field

### Student name: Lauren Gaughan

This repository contains code to investigate simple changes to parameter values in  $\Lambda$ CDM cosmology with an added single axion scalar field.

## Environment

The provided (docs/environment.yml) file can be used to create the conda environment required to run this code:

```
conda create -n new_env -f environment.yml
```

## Code

### LCDM.py

The LCDM.py file contains simple code to produce a plot of scale factor against conformal time for a flat  $\Lambda$ CDM universe. 

![LCDM_figure](https://user-images.githubusercontent.com/108680435/206687677-7b474d8e-b749-4b18-a231-968d81bc9800.png)

A simple GUI allows the user to change the value of $\Omega_{m}$ (and so $\Omega_{\Lambda}$) to investigate the effects on the scale factor over time.

Additionally the user can change whether they view the analytic solution, numerical solution, or both. 


### Axion.py

The Axion.py file contains code with a flat $\Lambda$ CDM universe with an additional scalar field responsible for producing the axion.

![Axion_figure](https://user-images.githubusercontent.com/108680435/206687744-3456fea8-05a1-4a40-9ad4-8d434e3bdf9c.png)

Again, a GUI allows the user to investigate the effects of changing parameter values on various quantities over time.

The user can see the effects of altering cosmological parameters:
* $\Omega_{m}$ and $\Omega_{\Lambda}$ = 1 - $\Omega_{m}$
* $\Omega_{r}$
* $h_{0}$
* G
Or the axion scalar field parameters:
* $\mu$
* $f$ - axion decay constant
* $\phi_{0}$ - the initial field value, with $\psi$ being the field velocity

The user can again alter which solutions shown for the scale factor as well as energy densities.

The plots are (left to right):
* Scale factor against conformal time
* Scalar field potential against scale factor
* $\Lambda$ CDM universe energy densities against redshift z
* Hubble ratio against z
* Scalar field velocity against conformal time
* $\Lambda$ CDM + single axion scalar field universe energy densities against redshift z

### module.py

The module.py file contains all relevant functions ultilised. 

These functions have been created for ease of use when calling functions such as solve_ivp or solving analytic equations.


