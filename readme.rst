General Information
===================
Implementation of an algorithm for determining differential flat ouputs
for nonlinear dynamic systems [1,2]. For further information on how
to use this package, see doc/user_guide.pdf.

An example has been implemented in the following IPython Notebook:
`Gantry crane <http://nbviewer.jupyter.org/github/klim-/uc_algorithm/blob/master/examples/gantry_crane_notebook/brueckenkran.ipynb>`_



[1] M. Franke and K. Röbenack: "On the Computation of Flat Outputs for Nonlinear Control Systems". In European Proc. of  the European Control Conference (ECC), Zurich, July 2013, p. 167-172 (https://ieeexplore.ieee.org/document/6669771/)

[2] K. Fritzsche, M. Franke, C. Knoll and K. Röbenack: "Zur systematischen Bestimmung flacher Ausgänge nichtlinearer Mehrgrößensysteme". In at - Automatisierungstechnik, 64(12), 2016, pp. 948-960, (http://dx.doi.org/10.1515/auto-2016-0079)


Installation
============
Make sure you have the following dependencies installed:

- sympy
- numpy
- symbtools
- pycartan

It is recommended to install these using the Python Package Index PyPI::

    $ pip install sympy numpy symbtools pycartan

Install uc_algorithm by cloning the repo::

    $ git clone https://github.com/klim-/uc_algorithm.git
