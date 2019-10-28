_RaDMaX online_ is a web-based program that allows to retrieve strain and disorder profiles in ion-irradiated materials from the simulation of X-ray diffraction data recorded in symmetric thêta-2thêta geometry.

_RaDMaX online_ is written in [Python](https://www.python.org/), using the [NumPy](https://numpy.org/) and [SciPy](https://scipy.org/) libraries. The graphical user interface is written within a [Jupyter](https://jupyter.org/) notebook using [ipywidgets](https://github.com/jupyter-widgets/ipywidgets) for interactive widgets and [bqplot](https://github.com/bloomberg/bqplot) for interactive plots. The html/css/javascript rendering is achieved with [voilà](https://github.com/voila-dashboards/voila).

Some crystallographic calculations are performed using [xrayutilities](https://xrayutilities.sourceforge.io/).

## Program usage

_RaDMaX online_ is a web application hosted at Universtiy of Limoges: 
[<img src="https://www.unilim.fr/wp-content/uploads/sites/8/2015/09/logo-ul@2x.png" alt="unilim" width="110"/>](https://radmax.unilim.fr/)

If the previous link doesn't work or if the server is saturated (no or slow access), a Binder instance of _RaDMaX online_ can be launched here:
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/aboulle/RaDMaX-webapp/master?urlpath=voila%2Frender%2FRaDMaX.ipynb)

At first launch Binder converts the github repository into a Docker image which might take some time (up to a few minutes). See the [Binder website](https://mybinder.org/) for further details.

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
When working in local mode change the first line of RaDMaX.ipynb to ```local = 1``` to benefit from automatic session saving capabilities.

## License
This program is licensed under the  CeCILL license. See [LICENSE](https://github.com/aboulle/RaDMaX-webapp/blob/master/LICENSE.txt) file.
