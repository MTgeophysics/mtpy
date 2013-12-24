MTpy
====

A Python toolbox for standard MT data processing

Installation
------------

    $ python setup.py install

If you have [setuptools](http://pypi.python.org/pypi/setuptools) installed, we
suggest you use:

    $ python setup.py develop

Documentation
-------------

See http://geophysics.github.com/mtpy



License
-------

MTpy is licensed under the GPL version 3

The license agreement is contained in the repository and should be kept together with the code.


Copyright: Lars Krieger, Jared Peacock, Kent Inverarity; Adelaide, 2013



Conventions
-----------

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
11. Internal coordinates: X = North, Y = East
12. Rotations are interpreted clockwise (mathematically negative)
13. 0 degrees azimuth = North





