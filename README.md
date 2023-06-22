scatterBrains: an open database of human head models and companion optode locations for realistic Monte Carlo photon simulations <!-- omit in toc -->
===================================================

Author: Melissa M. Wu (wu dot melissa dot m at gmail.com)

## Table of contents <!-- omit in toc -->

- [Introduction](#introduction)
- [Dependencies](#dependencies)
- [Subject models](#subject-models)
- [Example scripts](#example-scripts)
- [Running MMC](#running-mmc)
- [Acknowledgements](#acknowledgements)
- [References](#references)
  
## Introduction
===================================================

This is an open data set of MRI-based human head models, along with some associated functions, compiled to generate input data for mesh-based Monte Carlo light transport simulations. The data set consists of volume and mesh models from sixteen subject MRI scans, as well as approximately 800 systematically placed head surface locations for each model. The models themselves are divided into five tissue layers, representing the scalp, skull, cerebrospinal fluid, grey matter, and white matter. On any given model, users are able to select one of the pre-defined head surface locations (or another location of their own choosing) and input any desired source-detector distances into a probe-placing function - the function will then return 3D source and detector locations, along with a source launch vector, appropriate for running a mesh-based Monte Carlo simulation.

## Dependencies
===================================================

Required dependencies are the mesh-based Monte Carlo simulation software (version 1.9, v2020, Moon Cake - beta, https://sourceforge.net/projects/mcx/files/mmc/mmc%20binary/1.9%20%28v2020%20Moon%20Cake%20-%20beta%29/) and iso2mesh package (version 1.9.0, Century Egg, https://iso2mesh.sourceforge.net/cgi-bin/index.cgi) written by Qianqian Fang. Please download them into the same folder as this package.

## Subject models
===================================================

For details on subject demographics, please see the associated paper in the References section below.

Each subject folder contains several types of files:

> .dat files 
 
contain the subject mesh information required to run the Monte Carlo simulations. They can be loaded into Matlab using the "readmmc" functions provided in the MMC software (see example.m lines 19-22 for an example). The models have been de-identified and also have had their ears removed for more accurate probe placement in certain areas of the head. If models with the ears are needed, please send Melissa an email and she can provide them to you.

>"subject locs" .mat files 

contain the ~800 points that have been systematically placed on each subject's head surface. The example.m script has a visualization of these for subject 3. These files also contain rz and lz, two other points which are needed as inputs to the probe placing function.

> "subject mesh" .mat files 

contain the node, elem, and face variables for the models. They also contain layer-by-layer surface mesh variables for user visualization. An example can be found in "example_visualization.m". Since MMC requires the .dat files for model specification, these subject mesh files are primarily for user visualization.

> "subject volume" .mat files 

contain the volume files which were used to generate the mesh files. 

## Example scripts
===================================================

The main folder contains several example scripts:

> Example.m:

runs through the entire pipeline with subject 3. It reads the subject mesh, visualizes the pre-defined points, and selects one of them. It then calls the probe-placement function at that location for several source-detector distances. This function automatically places sources and detectors at specified source-detector distances on a given location on the head. After ensuring proper source placement, the example script writes this information to a .json file used as an input to MMC. The user can then run the MMC simulation in their preferred manner (an example command line call is given).

For users working with diffuse correlation spectroscopy, an additional post-processing example is given at the end of the script to calculate and fit autocorrelations.

> Example_visualization.m:

allows the user to load the subject mesh files and visualize the model structure.

> Example_mesh_generation.m:
 
shows how the volume and surface meshes were generated from the original subject volume files.

> Example_distance2brain.m:

calculates the distance from a specified surface location to the nearest brain tissue voxel.

## Running MMC
===================================================

Please remember to add an optical property file (e.g. "prop_subject03.dat") to the relevant subject folder before running MMC. More information about this can be found on the MMC help page: http://mcx.space/wiki/index.cgi?MMC

Sometimes the source is placed a little too far away from the surface of the head, which prevents MMC from running. To address this, the source can be moved incrementally along the direction of the source launch vector until the source is detected to be touching an element. This is shown in lines 64-72 of the example.m script.

## Acknowledgements
===================================================

This dataset and associated functions are dependent on and/or primarily written to be used in conjunction with MMC and iso2mesh, both written by Qianqian Fang, who heads the Computational Optics and Translational Imaging Lab at Northeastern University.

> MMC:

- Summary: a mesh-based Monte Carlo photon transport simulation software
- License: GNU General Public License version 3.0
- URL: http://mcx.space/wiki/index.cgi?MMC
  
> iso2mesh:

- Summary: a 3D mesh generator for Matlab/Octave
- License: GNU General Public License version 2.0 or later
- URL: https://iso2mesh.sourceforge.net/cgi-bin/index.cgi

## References
===================================================