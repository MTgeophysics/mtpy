MTpy: A Python Toolbox for Magnetotelluric (MT) Data Processing, Analysis, Modelling and Visualization
==================================

|Build Status|  |Documentation Status|


Overview
========

A Python Toolbox for Magnetotelluric (MT) Data Processing, Analysis, Modelling and Visualization

- Documentation: http://mtpy2.readthedocs.io/en/develop/

- Home Page: https://github.com/MTgeophysics/mtpy

- Issue tracking: https://github.com/MTgeophysics/mtpy/issues

- Wiki Pages: https://github.com/MTgeophysics/mtpy/wiki


Note that this repository has superseded the `geophysics/mtpy <https://github.com/geophysics/mtpy/tree/beta>`_
and `GeoscienceAustralia/mtpy2 <https://github.com/GeoscienceAustralia/mtpy2/tree/develop>`_


Contacts
==========

| **Alison Kirkby**
| Alison.Kirkby@ga.gov.au

| **Fei Zhang**
| fei.zhang@ga.gov.au

| **Jared Peacock**
| peacock.jared@gmail.com

| **Rakib Hassan**
| Rakib.Hassan@ga.gov.au

| **Jinming Duan**
| Jingming.Duan@ga.gov.au



System Requirements
==========================

-  Python 2.7
-  Python 3.6+

Setup Guide for Developers
==========================

1. Install Python environment and dependency packages. See Wiki Pages: https://github.com/MTgeophysics/mtpy/wiki

2. Obtain the source code from https://github.com/MTgeophysics/mtpy:

-  ``git clone https://github.com/MTgeophysics/mtpy.git``
- ``cd mtpy``

   - ``pip install -v --user -e .`` (into user's own home ~/.local/lib/python2.7/site-packages/mtpy.egg-link)
   OR 
   
   - ``python setup.py develop --user``   
   OR 
   
   - ``pip install -v -e .``  (into python lib's dir site-packages, write-permission required)   
   OR 
   
   - `` export  PYTHONPATH=/Path2/mtpy:$PYTHONPATH `` (Only valid for each session)
   
   
To verify the install: 

- ``pip list | grep mtpy``

- ``pip show mtpy``

To uninstall the mtpy package: 

- ``pip uninstall -v mtpy``



License
===============

MTpy is licensed under the GPL version 3

The license agreement is contained in the repository and should be kept together with the code.


Conventions
===============

1. MTpy uses E- and B-fields (although the sensors may be confusingly named as H-sensors in EDI files)
2. [E] = microvolts/meter (muV/m)
3. [B] = nanotesla (nT)
4. [Z] = [E]/[B] = km/s
5. Apparent resistivty rho = 0.2 * T * |Z|^2  (in Ohm m)
6. Angles are given in degrees (mod 360)
7. EDI files can contain data in Z- or rho/phi-form
8. EDI files contain data from one station only
9. Coordinates are handled in decimal degrees (converted when reading)
10. Time stamps refer to UTC
11. Internal coordinates: X = North-South, Y = East-West
12. Rotations are interpreted clockwise (mathematically negative)
13. 0 degrees azimuth = North



.. |Build Status| image:: https://travis-ci.org/MTgeophysics/mtpy.svg?branch=develop
   :target: https://travis-ci.org/MTgeophysics/mtpy

.. |Coverage Status| image:: https://coveralls.io/repos/github/MTgeophysics/mtpy/badge.svg?branch=develop
   :target: https://coveralls.io/github/MTgeophysics/mtpy?branch=develop

.. |Documentation Status| image:: https://readthedocs.org/projects/mtpy2/badge/?version=develop
   :target: http://mtpy2.readthedocs.io/en/develop/


