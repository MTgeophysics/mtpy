MTpy: A Python Toolbox for Magnetotelluric (MT) Data Processing, Analysis, Modelling and Visualization
==================================

|Build Status|  |Documentation Status|

How to Cite
===========

If you use this software in a scientific publication, we'd very much appreciate if you could cite the following papers:

- Kirkby, A.L., Zhang, F., Peacock, J., Hassan, R., Duan, J., 2019. The MTPy software package for magnetotelluric data analysis and visualisation. Journal of Open Source Software, 4(37), 1358. https://doi.org/10.21105/joss.01358
   
- Krieger, L., and Peacock, J., 2014. MTpy: A Python toolbox for magnetotellurics. Computers and Geosciences, 72, p167-175. https://doi.org/10.1016/j.cageo.2014.07.013

Overview
========

A Python Toolbox for Magnetotelluric (MT) Data Processing, Analysis, Modelling and Visualization

- Home Page: https://github.com/MTgeophysics/mtpy

- API Documentation: http://mtpy2.readthedocs.io/en/develop/

- Issue tracking: https://github.com/MTgeophysics/mtpy/issues

- Installation Guide (Wiki Pages): https://github.com/MTgeophysics/mtpy/wiki

- User Guide: https://github.com/MTgeophysics/mtpy/blob/develop/docs/MTPy%20User%20Guide.pdf


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

| **Bren Moushall**
| brenainn.moushall@ga.gov.au

| **Rakib Hassan**
| Rakib.Hassan@ga.gov.au

| **Jingming Duan**
| Jingming.Duan@ga.gov.au



System Requirements
==========================

-  Python 2.7
-  Python 3.6+


License
===============

MTpy is licensed under the GPL version 3

The license agreement is contained in the repository and should be kept together with the code.


Conventions Used in the MTPy Software
=====================================

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


