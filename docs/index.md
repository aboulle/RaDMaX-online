_RaDMaX online_ is a web-based program that allows to retrieve strain and disorder depth-profiles in ion-irradiated materials from the simulation of X-ray diffraction data recorded in symmetric $\theta/2\theta$ geometry.

_RaDMaX online_ is written in [Python](https://www.python.org/), using the [NumPy](https://numpy.org/) and [SciPy](https://scipy.org/) libraries. The graphical user interface is written within a [Jupyter](https://jupyter.org/) notebook using [ipywidgets](https://github.com/jupyter-widgets/ipywidgets) for interactive widgets and [bqplot](https://github.com/bloomberg/bqplot) for interactive plots. The html/css/javascript rendering is achieved with [voilà](https://github.com/voila-dashboards/voila). Some crystallographic calculations are performed using [xrayutilities](https://xrayutilities.sourceforge.io/).

**Data privacy**: **no data is stored on the server**. All data files uploaded to _RaDMaX online_ and all calculations performed with _RaDMaX online_ are stored in RAM and are definitely lost when uploading a new data set or when closing the program. Session saving capabilities are available if you run the Jupyter notebook in offline mode (see [below](#offline-mode)).

**Help and support**: new crystal structures can be added upon request. Bug reports and improvement suggestions are welcome. Contact info: [alexandre.boulle@unilim.fr](mailto:alexandre.boulle@unilim.fr)

## Program usage

_RaDMaX online_ is a web application hosted at Universtiy of Limoges (the page may take a few seconds to load):

[<img src="https://www.unilim.fr/wp-content/uploads/sites/8/2015/09/logo-ul@2x.png" alt="unilim" width="110"/>](https://radmax.unilim.fr/)

If for any reason the previous link doesn't work properly, a Binder instance of _RaDMaX online_ can be launched here:

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/aboulle/RaDMaX-online/master?urlpath=voila%2Frender%2FRaDMaX.ipynb)

At first launch Binder converts the github repository into a Docker image which might take some time (up to a few minutes). See the [Binder website](https://mybinder.org/) for further details.

## Video tutorial

<iframe width="560" height="315" src="https://www.youtube.com/embed/53y2MJijKps" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

## Offline mode

The Jupyter notebook can also be executed locally. Clone or download this repository and install all required dependencies. Using pip :

```bash
pip install numpy scipy jupyter ipywidgets bqplot voila xrayutilities
```

If you are using [anaconda](https://www.anaconda.com/distribution/), the scipy stack and jupyter are already installed:

```bash
conda install -c conda-forge ipywidgets bqplot voila
```

then

```bash
pip install xrayutilities
```

When working in local mode change the first line of RaDMaX.ipynb to ```local = True``` to benefit from automatic session saving capabilities.

## Scientific background

In *RaDMaX online* the diffracted X-ray intensity is computed within the framework of the dynamical theory of diffraction from distorted crystals as first theorized by [S. Takagi](https://doi.org/10.1143/JPSJ.26.1239) and [D. Taupin](https://doi.org/10.3406/bulmi.1964.5769) which is the state of the art approach for the analysis of irradiated crystals using XRD. More specifically, *RaDMaX online* implements the 1D solution to the Takagi - Taupin equations as derived by [Bartels, Hornstra and Lobeek](https://doi.org/10.1107/S0108767386098768). Within this approach, the irradiated region is divided in a number of sub-layers each having distinct, but fixed, compostion/strain/disorder. The scattering from the ensemble is computed iteratively, starting from the film-substrate interface, up to the surface. The X-ray amplitude ratio $X_{n+1}$ at the top of a layer is related to the amplitude ratio at the bottom, $X_n$, by :

$$
X_{n+1} = \eta + \sqrt { \eta^2 - 1} {S_1 + S_2 \over S_1 - S_2}
$$

with

$$
S_1 = \left(X_n - \eta + \sqrt {\eta ^2 -1} \right) \exp \left( - i T \sqrt {\eta ^2 -1} \right)
$$

$$
S_2 = \left(X_n - \eta - \sqrt {\eta ^2 -1}  \right) \exp \left( i T \sqrt {\eta ^2 -1} \right)
$$

The amplitude ratio X is given by (in order not to overly complexify the equations, the index n is ommited below):

$$
X = \sqrt {\left( F_{\bar H} \right / F_{H}) \left| \gamma_H  / \gamma_0\right|} D_H / D_0

$$

$\eta$ is dynamical theory's deviation parameter:

$$
\eta = \left[ -b (\theta - \theta_B) \sin 2\theta_B - \frac{1}{2} \Gamma F_0 (1-b) \right] / \left(C \Gamma \sqrt {|b| F_H F_{\bar H}} \right)
$$

and $T$ is related to the sub-layer thickness: 

$$
T = \pi C \Gamma \sqrt{F_H F_{\bar H}} t / \left( \lambda \sqrt{\gamma_0 \gamma_H} \right)
$$

where $\Gamma = r_e \lambda / \pi V$, $r_e = e^2 / 4 \pi \epsilon_0 m c^2$ and $b = \gamma_0 / \gamma_H$. $C$ is the polarization of the incident beam, $D_0$ and $D_H$ are the incident and diffracted beam amplitudes, $\gamma_0$ and $\gamma_H$ are the direction cosines of the incident and diffracted beam with respect to the surface normal. In the above equations, $t$ is the sub-layer thickness, $F_H$ and $F_{\bar H}$ are the structure factor of the $hkl$ and $\overline{hkl}$ reflections and $\theta_B$ is the Bragg angle within the given layer. This angle is related to the strain component $e_{zz}$ via:

$$
\theta_B = \arcsin \left[ \lambda / 2 d_0 (1+e_{zz})\right]
$$

where $d_0$ is the lattice spacing of the bulk (unstrained) material. The structure factor $F_H$ is sensitive to the atomic disorder via:

$$
F_H = DW \times F_{0,H}
$$

where $F_{0,H}$ is he structure factor of the bulk material and $DW$ is the Debye-Waller factor which ranges between 0 and 1. The $DW$ is related to random atomic displacements $\delta \bold u$ via:

$$
DW = \left\langle  \exp\left[  i \bold{H} \delta \bold{u}  \right]  \right\rangle
$$

The equations above are used to generated the diffracted intensity for all $\theta$ values provided by the user. **The fitting procedure consist in finding the best $e_{zz}$ and $DW$ values so that the calculated curve matches the experimental data provided by the user. This can be done either manually (via the interactive strain/DW plots), or automatically using a least-squares fitting procedure**.

In order to limit the number of adjustable parameters and avoid numerical instabilities, the values of $e_{zz}$ and $DW$ are constrained to exhibit a *cubic spline* shape. This is performed using *cubic B-spline* functions:

$$
f (z) = \sum_{i = 1}^{N_w^{S,D}} w_i^{S,D} B_{i,3}(z)
$$

where $N_w^{S,D}$ is the number of B-splines used to describe the strain (‘S’) and disorder (‘D’) profiles (typical values are in the 5-15 range), $w_i^{S,D}$ are the weights to be determined in the fitting procedure and $B_{i,3}(z)$
is the third-degree basis function, and $z$ is the depth coordinate.

Furter details regarding B-spline functions and the theoretical background can be found in the following references: [Boulle & Debelle, 2010](https://www.unilim.fr/pages_perso/alexandre.boulle/files/Boulle%20JAC%202010.pdf), [Souilah, Boulle & Debelle, 2016](https://www.unilim.fr/pages_perso/alexandre.boulle/files/Souilah%20JAC%2015.pdf) and [Boulle & Mergnac, 2020](https://arxiv.org/pdf/1911.06521.pdf)

The python code implementing these equations for the case of a single crystal is given below:

```python
def f_Refl_Default(th, param, cst):
    offset = cst["offset"]*np.pi/360
    G = cst["G"]
    thB_S = cst["thB_S"]
    wl = cst["wl"]
    t = cst["t"]
    N = cst["N"]
    resol = cst["resol"]
    b_S = cst["b_S"]
    phi = cst["phi"]
    t_l = cst["t_l"]
    z = cst["z"]
    FH = cst["FH"]
    FmH = cst["FmH"]
    F0 = cst["F0"]
    Nspline = cst["sdw_basis"]
    model = cst["sdw_model"]
    bkg = cst["bkg"]

    th = th + offset
    strain = f_strain(z, param[:Nspline:], t, model)
    DW = f_DW(z, param[Nspline:2*Nspline:], t, model)
    thB = thB_S - strain * np.tan(thB_S)

    eta = (-b_S*(th-thB_S)*np.sin(2*thB_S) - 0.5*G*F0*(1-b_S))\
    /(b_S)**0.5)*G*(FH*FmH)**0.5 )
    res = (eta - np.sign(eta.real)*((eta*eta - 1)**0.5))

    n = 1
    while (n<=N):
        g0 = np.sin(thB[n] - phi)
        gH = -np.sin(thB[n] + phi)
        b = g0 / gH
        T = np.pi * G * ((FH*FmH)**0.5) * t_l * DW[n]/ (wl\
        *(abs(g0*gH)**0.5) )
        eta = (-b*(th-thB[n])*np.sin(2*thB_S) - 0.5*G*F0*(1-b))\
        /((abs(b)**0.5)*G*DW[n]*(FH*FmH)**0.5)
        sqrt_eta2 = (eta*eta-1)**0.5

        S1 = (res - eta + sqrt_eta2)*np.exp(-1j*T*sqrt_eta2)
        S2 = (res - eta - sqrt_eta2)*np.exp(1j*T*sqrt_eta2)

        res = (eta + sqrt_eta2*((S1+S2)/(S1-S2)))
        n += 1

    ical = np.convolve(abs(res)**2, resol, mode='same')
    return (ical/ical.max())+bkg
```

## License

This program is licensed under the  CeCILL license. See [LICENSE](https://github.com/aboulle/RaDMaX-online/blob/master/LICENSE.txt) file.
