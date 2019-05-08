---
title: "The MTPy software package for magnetotelluric data analysis and visualisation"
tags:
  - Python
  - MTPy
  - magnetotellurics
authors:
  - name: Alison Louise Kirkby
    orcid: 0000-0003-1361-440X
    affiliation: 1
  - name: Fei Zhang
    affiliation: 1
  - name: Jared Peacock
    affiliation: 2
  - name: Rakib Hassan
    affiliation: 1
  - name: Jingming Duan
    affiliation: 1
affiliations:
  - name: Geoscience Australia
    index: 1
  - name: United States Geological Survey
    index: 2
date: 1 April 2019
bibliography: paper.bib
---


# Introduction

The magnetotelluric (MT) method is increasingly being applied to a wide variety of geoscience problems. However, the software available for MT data analysis and interpretation is still very limited in comparison to many of the more mature geophysical methods such as the gravity, magnetic or seismic reflection methods.

MTPy is an open source Python package to assist with MT data processing, analysis, modelling, visualization and interpretation. It was initiated at the University of Adelaide in 2013 as a means to store and share Python code amongst the MT community [@kriegerpeacock2014]. Here we provide an overview of the software and describe recent developments to MTPy. These include new functionality and a clean up and standardisation of the source code, as well as the addition of an integrated testing suite, documentation, and examples in order to facilitate the use of MT in the wider geophysics community.

# Key Functionality

MTPy contains modules to carry out much of the processing and analysis that can be applied to MT data. MTPy works primarily with the industry standard Electronic Data Interchange (EDI) file format however it is compatible with other formats including the Jfile format exported by the Bounded Influence, Remote Reference Processing (BIRRP) code [@chaveetal1987; @chavethomson2004], and xml. The MTPy structure is based on key steps in working with MT data, namely processing, analysis, modelling, and imaging. These sub-packages are based on modules within the core sub-package [Figure 1, @kriegerpeacock2014].

![MTPy package structure showing the key modules that are under active development. The colours in the workflow diagram represent which parts of MTPy are used at each step. For example, modules in the analysis sub-package are used on the impedance tensor data, while visualisation is carried out on time series data, processed data, and the results of analysis and modelling or inversion. This diagram is modified after @kriegerpeacock2014.](paper_figures/mtpy_diagram.png){width="5in"}

The core sub-package contains functionality to read and write MT data from industry standard formats such as EDI. When data are read in, an MT object is created. The MT object contains site metadata from the header of the EDI file, including the location. It also contains the impedance tensor (Z) and the vertical transfer functions (Tipper).

The processing sub-package is designed to facilitate working with time series data and generating inputs for existing third-party processing codes, for example the BIRRP code [@chaveetal1987; @chavethomson2004].

The analysis sub-package includes functions that are commonly used to analyse MT data prior to modelling and inversion. These include: dimensionality analysis, strike angle calculation, phase tensor analysis [@caldwelletal2004], and calculation and removal of distortion [@bibbyetal2005]. In addition, there are functions to calculate and visualize depth of investigation using the Niblett-Bostick transform [@niblett1960; @bostick1977; @jones1983].

MTPy also contains the modeling sub-package, which allows users to create inputs to, and visualize outputs from, several of the commonly-used MT inversion codes. These include the Occam 1D and 2D inversion codes [@constableetal1987; @degroothedlin1990], and the ModEM 3D inversion code [@egbertkelbert2012; @kelbertetal2014].

![Example of a single station plot generated using the plotresponse module within the imaging sub-package in MTPy. The plot shows Z~XY~ and Z~YX~ resistivity and phase, induction vectors [Tipper; @parkinson1962], and phase tensors [@caldwelletal2004] as a function of period.](paper_figures/Synth00.png){width="4.4015748031496065in"}

The imaging sub-package contains functionality to visualise MT data, as well as many of the outputs from the analysis and modelling modules. These include simple plots of resistivity and phase at a single station (e.g. Figure 2), resistivity and phase pseudosections and maps, phase tensor sections and maps (e.g. Figure 3), penetration depth and geoelectric strike plots, and images of resistivity models and the model responses.

![Example of a phase tensor and induction vector map created by the phase_tensor_maps module in the imaging sub-package in MTPy. Phase tensor ellipses are calculated using the method detailed by @caldwelletal2004 and coloured by skew angle in degree.](paper_figures/phase_tensor_map100s.png){width="5in"}

The functions in utils underpin these modules, by handling basic underlying functionality such as coordinate transformations, common filehandling operations and calculations that are commonly applied to MT data, such as conversion from impedance tensor (Z) to apparent resistivity and phase.

# Recent development

In order to make MTPy a comprehensive toolkit that can be easily installed and used by the wider geophysics community, we are applying software engineering best practices and techniques to improve software quality and usability. New developments include incorporation of a unit-testing suite, continuous build, automatic documentation generation (<http://mtpy2.readthedocs.io/en/develop/>) and addition of installation guide wiki pages (<https://github.com/MTgeophysics/mtpy/wiki>).

Recent functionality development to MTPy includes significant expansion of the modules that generate inputs to, and allow visualiation of outputs from, the ModEM 3D inversion code [@egbertkelbert2012; @kelbertetal2014]. This includes new functionality to extract vertical slices from 3D models for display on profiles. Both of these support overlay of other datasets in both profile and plan view. These can be used with the seismic module in imaging for plotting 2D seismic reflection data (in profile) and the geology module for plotting shapefiles. We have now also implemented functions to export ModEM 3D models to GOCAD^TM^ sgrid format, and to export depth slices to a column based text format containing x, y, z resistivity.

The analysis and imaging sub-packages have been expanded through the addition of modules that calculate and visualise depth of investigation of MT measurements, for a single site, a profile, or a grid of sites.

We have also added an edi_collection module that works specifically with collections of edi files and enables phase tensors from a collection of stations to be exported to csv or ArcView shapefiles.

Projections and transformations are now handled using the open-source Pyproj library, or if users prefer, the more comprehensive GDAL library.

Finally, we have compiled a set of jupyter notebooks, example scripts, and workshop material within the examples folder in the repository, which demonstrate the key functionality contained in MTPy.


# Availability

The software is distributed under a GNU General Public License v3.0 and is available from <https://github.com/MTgeophysics/mtpy>.

# Acknowledgements

MTPy was initiated at the University of Adelaide and is now developed by Geoscience Australia and the United States Geological Survey. This paper is published with the permission of the CEO, Geoscience Australia.

# References


