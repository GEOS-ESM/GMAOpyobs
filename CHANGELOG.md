# ChangeLog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

# [Unreleased]

### Added 
- Added cubed-sphere binning capability
- pm class to aop.py - with some additional comments
- G2GAOP can now take a string as the config file variable input
- add a function to sampler that can append a vertical coordinate to a sampled dataset
- example Jupyter notebooks on using pyobs utilities to sample GEOS and compare to CALIPSO and DC-8 obs

### Changed 

- Updated README to document lite_install
- Update `components.yaml` to match that of AeroApps
  - Mainly for newer ESMA_env that allows building on RHEL8 GMAO
    machines (e.g., calculon)
  - Use postfix-@ for subrepos to match AeroApps
- Allow ability to not build f2py code for CI purposes

### Fixed 
- import of IGBP_
- missing conversion from the sulfate ion to ammonium sulfate
- aop.py *getAOPrt* phase function now being correctly normalized 
  by total scattering
- aop.py - protect against divide by zero in getAOPrt and getAOPext when doing calculation for an individual species
- aop.py - remove dependency on having 'DU' as a species in your yaml optics table definition
- calipso_l1p5.py - took out extinction from list of variables.  L1.5 files don't have this
- calipso_l2.py - use variable attributes to mask missing data, read the altitude coordinate from metadata
- constants.py - fix typo in units of gravity and calipso_l2 scripts 

### Removed 

# [v1.1.0] 2024-03-20

### Added 

- Module *waypoint* to handle waypoint files for flight modules 
- Module *sampler*, in pure python, including both station and 
- Module *mietable* for handling GEOSmie tables 
- Module *aop*, first draft of a replacement to the old aod_calculator 


### Changed 

- pyobs *__init__* method no longer loads submodules by default 
- stn_sample command line utility rewritten in terms of the new 
  *sampler* module 
- trj_sample command line utility rewritten in terms of the new 
  *sampler* module 
  - module *icartt* extebded with to_xarray() method.
  - mietables.py wavelength bug fix (nm to m unit conversion)

### Fixed 

- Modernized mcd43.py by making use of xarray and cartopy map 
  transforms. 

### Removed 

# [v1.0.8]

### Added

### Fixed
- fixed granules search in vx04.py to use updates to VIIRS path logic
- added a DB_DEEP retrieval.  splits up land retrievals into the 2 pathways - one that uses the 412 surface channel and one that does not

### Removed

# [v1.0.7]

### Fixed

### Changed 

- converted active_aeronet.py to py3



# [v1.0.6] 2023-09-08

### Added

- mxd04 and vx04 will write out gridded angstrom exponent if available



# [v1.0.5] 2023-06-15

### Added

- mxd04 Level3 gridded product writes out the number of observations that went into the grid cell
- man.py function to create a concatenated of all the seperate cruises

### Fixed

- reading date/time variables from MAN concatenated file needs to be unicode not byte

### Removed

# [v.1.0.4] 2023-06-07


### Fixed

- removed all instances of using the 'types' module in pyobs



# [v1.0.3] 2023-05-25

### Added

- add interpolated 490 AOD to aeronet.py. lines up with VIIRS wavelengths
- subroutine binObsCnt3D to binObs_py to counts obs in NNR L3 files
- vx04.py VIIRS reader


# [v1.0.2] 2023-05-17


### Fixed

- mxd04.writeods now writes a 'post_anal' file. this saves the original retrieved AOD in the ods files


## [1.0.1] - 2023-05-16

### Fixed

- bug in config.py. integer divide needs // in python3
- bug in icarrtt.py. Config call not updated when going to python3

### Added

- Added changelog enforcer and yaml linter

## [1.0.0] - 2023-05-11

### Added

- Initial release
- Constructing GMAOpyobs from bits and pieces of GMAO_Shared/GMAO_pyobs

