MTpy Version 2 Proposed Updates
===================================

General Changes
-----------------
	* use Path instead of os and glob wherever possible
	* Remove compatiblity with Python 2
	* Use f"{variable}" for print statements and writing a string
	* Use Git Actions for testing, updating version, updating PyPi and Conda
  * Refactor some of the modules to new repositories
  * Suggest putting mt_metadata under MTGeophysics
  * Create a GIS module for all the functions that we have for plotting, making shapefiles, etc, and maybe leverage cartopy?

Core
-------
z
~~
	* Create a TransferFunction class that can hold covariance matrices. The covariance matrices are needed to make sure the errors are correct when the data are rotated.  There could also be different types of transfer functions besides MT and Tipper, like interstation transfer functions.  We need a robust way of dealing with them.   
		- TransferFunction would be under Z and Tipper
		- Could be inherited by them, though that might create redunant storage
		- Z and Tipper could be attributes of TransferFunction
		- Alternatively, could have Z and Tipper be able to handle covariance matrices.
MT
~~~
	- Be able to read/write more file types for the transfer functions into a standardized container with standard metadata.
		- File types:
			- **.edi**: most commonly used exchangeable format, does not store covariances. Decent amount of metadata
			- **.zmm**, **.zrr**: Output of Egbert's EMTF code, has covariances. Minimal metadata
			- **.j**: Output of BIRRP, does not natively store covariances, those are in different files, minimal metadata.
			- **.avg**: output of Zonge processing codes, no covariances, minimal metadata.
			- **.xml**: A. Kelbert's EMTF XML format, has covariances, has metadata based on her standard.
		- Have a plug-in type setup for each file type that will translate into an MT object and vice-versa.  
		- Metadata standards are based on mt_metadata.transfer_function
    - Add to_xarray/from_xarray functions for easier manipulation in a collection of transfer functions.
      - Metadata is stored in xarray.DataArray.attrs
      - Dimensions are frequencies
      - Coordinates are z, z_err, tipper, tipper_err, resisitivity, phase, phase_tensor, resistivity_tensor, covariances, etc
		
	
MTCollection
~~~~~~~~~~~~~
	- Updated EDICollection to be all inclusive of any type of transfer function file supported
  - Use xarray.DataSet as the container
    - Keeps metadata of each xarray.DataArray
    - Can easily interpolate or pick frequencies
    - Can easily query for location, period range, date, data type, etc.
  - This should be the base for nearl all functionality for imaging, modeling, analysis.  This will likely reduce the code base as well because it will remove redundancies.
  - Could add to/from_modem, to/from_occam, etc here or in the modeling modules.
  - Add to/from NetCDF or HDF5 so that a user would only have to read in the various transfer function files once.
  - Add functions to add/remove stations
  - Add functions to use geopandas to make GIS files of station locations with metadata, phase tensors, induction vectors, etc.
    
		
Analysis
---------
  - Update modules to use MTCollection
  - Update staticshift to use MTCollection to estimate a local static shift from nearby stations.

Imaging
---------
  - Updated to use MTCollection
  - Many of these files are long because of all the 

GUI
----
  - Suggest moving this to a new repository that calls MTpy. 

Processing
-----------
  - Develop tools to work with IRIS' new Aurora processing code (releasing in Fall, 2021)
  - Develop tools to work with `Razorback <https://github.com/BRGM/razorback/>`_?
  - Develop tools to work with `Resistics <https://github.com/resistics/resistics>`_?
  - Develop tools to read covariances from BIRRP output files (New in BIRRP 5.3).

Modeling
-----------
  - Update Data to use MTCollection
  - Tools for other modeling programs that might become available in the future?
  - Develop tools to store outputs in NetCDF for other programs to read.
  - Develop tools to read other geophyiscal models from NetCDF and input into models
  - Think about integrating with `Discretize <https://github.com/simpeg/discretize>`_ for making meshes.  
  - Develop tools to use `SimPEG <https://github.com/simpeg/simpeg>`_ for 1-D, 2-D, and 3-D modeling (in development).


UofA
-----
  - Suggest moving this to a new repository

USGS
------
  - Suggest moving this to a new repository


Utils
---------
  - Suggest moving files to more logical places.  There are some plotting tools in there   

Documentation
---------------
  - Suggest updating to the "sphinx_rtd_theme"
  - Suggest adding examples, usage, history, introduction in the docs.
  - Update doc strings, never ending!

Logging
---------
  - Update how the loggers are initiated and where the logs go.
  - Update functions and classes to have comprehensive logging

Tests
-------
  - Test for backwards compatibility to version 1.*
  - Update test, , also never ending!

Examples
----------
  - Add as many examples as possible
  - Suggest adding an mtpy_examples respository
    - Put most example data here, which would lighten the size of the mtpy distribution
    - Have the same folder structure as mtpy with an example for each.
    - Have examples of where MTpy was used in published studies.
    
  
  
