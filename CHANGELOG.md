# ChangeLog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

# [Unreleased] - yyyy-mm-dd

### Fixed

### Added

### Changed

### Removed

### Deprecated

# [v1.4.0] - 2025-09-25

### Fixed

- fixed error when duplicate variables occur in an icartt file. ICARTT now reads the first instance of the duplicate variable.
- fixed bug in aop.py that broke if you used v1.X.X optics tables (that don't include pmatrix)
- fixed bug in aop.py that broke if you had lowercase DELP in your model files

### Added

- add KX for lunar AERONET obs
- add KX for NOAA-20 and NOAA-21 AOD obs


# [v1.3.0] - 2025-07-28

### Fixed
- Fixed some misspelled `ValueError` exceptions
- added a min dim to aeronet reads.  fixes error when there is only 1 observation in the file
  
### Added
- now getAOPrt returns the phase matrix as detault instead of pmom (the pmatrix expansion)

### Changed
- Update to ESMA_cmake v3.62.1


# [v1.2.6] - 2025-07-28

### Fixed
- hotfix in aeronet.py for double decoding


# [v1.2.5] - 2025-06-05

### Fixed
- fixed mietable runtime warning for pback calculation. don't do this recursively anymore.
- commented out "hack" in xrctl. it wasn't doing anything but added memory overhead and runtime
### Added
- can pass additional keywork arguments to TRAJECTORY and SAMPLER
### Changed

### Removed

### Deprecated


# [v1.2.4] - 2025-05-29

### Fixed
- Added list parsing for variables in trajectory sampler
- fixed byte string bug in aeronet.py
- - use a local copy of RH in aop calculator.  otherwise it overwrites when fixRH is used
### Added
- MPL reader and plot curtain 
- calculation of total backscatter coefficient in aop.py
- xrctl supports providing a list of control files
- parse time in MPL reader to return datetimes
- sampler notebook that uses station sampler at an MPL
- add option for vacuum aerodynamic size cutoff
### Changed
- add auto chunking to TRAJECTORY and STATION. This enables dask
- preload some key variables in aop.py so you don't hit dask repeatedly in for loop
- rename p11 and p22 to pback11 & 22
- updates aop aback calculation to use name indexing. generalizes to any dimensionality
- updates station sampler to sample at an input time frequency. defaults to model time frequency
- updates xrctl to read multiple GEOS collections at once
- updates the compile and build instructions in the readme
### Removed

### Deprecated

# [v1.2.3] - 2025-03-20

### Fixed
- bug/typo in aop.py for vector AOP calculation

### Added
- ability to limit number of moments requested from getAOPrt

### Changed

### Removed

### Deprecated

# [v1.2.2] - 2025-03-03 

### Fixed 
- bad reshape import in kde.py
### Added 
- Added TROPOMI Level 2 reader	
### Changed 

### Removed 

### Deprecated 


# [v1.2.1] - 2025-02-04 

### Fixed 
- Bug fix in aop.py for returning PMOM

### Added 
	
### Changed 

### Removed 

### Deprecated 


# [v1.2.0] - 2025-01-16

### Added 
- Added module mcbef for handling plume rise/fire related functionality.	
- Added cubed-sphere binning capability
- pm class to aop.py - with some additional comments
- G2GAOP can now take a string as the config file variable input
- add a function to sampler that can append a vertical coordinate to a sampled dataset
- example Jupyter notebooks on using pyobs utilities to sample GEOS
   and compare to CALIPSO and DC-8 obs
- Added combustion module.

### Changed 
- Updated README to document lite_install
- Update `components.yaml` to match that of AeroApps
  - Mainly for newer ESMA_env that allows building on RHEL8 GMAO
    machines (e.g., calculon)
  - Use postfix-@ for subrepos to match AeroApps
- Allow ability to not build f2py code for CI purposes
	
### Fixed 
- import of IGBP_ in igbp.py
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

