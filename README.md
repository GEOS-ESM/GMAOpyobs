# GMAOpyobs

## What is it?

**GMAOpyobs** is a set of python modules (some built with f2py) for
implementing interfaces to several of satellite and ground-based
observing systems. The main design philosophy is to make it as *pure
Python* as possible, with *f2py* only when absolutely necessary. In
particular, it should not contain any dependence on external libraries
such as ODS or GFIO.

## Python 3 only

GMAOpyobs no longer supports Python 2.7.


## How to build GMAOpyobs on a GMAO supported system 

### Preliminary steps

#### Load build modules

In your `.bashrc` or `.tcshrc` or other rc file add a line:

##### NCCS (SLES12)

```
module use -a /discover/swdev/gmao_SIteam/modulefiles-SLES12
```

##### NAS
```
module use -a /nobackup/gmao_SIteam/modulefiles
```

##### GMAO Desktops
On the GMAO desktops, the SI Team modulefiles should automatically be
part of running `module avail` but if not, they are in:

```
module use -a /ford1/share/gmao_SIteam/modulefiles
```

Also do this in any interactive window you have. This allows you to get module files needed to correctly checkout and build the model.

Now load the `GEOSenv` module:
```
module load GEOSenv
```
which obtains the latest `git`, `CMake`, etc. modules needed to build.

### Use mepo to clone the repository

[Mepo](https://github.com/GEOS-ESM/mepo) is a multiple repository tool available on github.

```
mepo clone git@github.com:GEOS-ESM/GMAOpyobs.git
```


### Building and installing GMAOpyobs

#### Load compiler, MPI stack, and baselibs
On tcsh:
```
source @env/g5_modules
```
or on bash:
```
source @env/g5_modules.sh
```
#### One-step default build
Run the cmake_it script
```
./cmake_it
```

#### Multi-step build

Follow these instructions for managing multi-version builds, or custom builds.

##### Create build directory

We currently do not allow in-source builds of GEOSgcm. So we must make a directory:
```
mkdir build
```
The advantages of this is that you can build both a Debug and Release version with the same clone if desired.

##### Run cmake

CMake generates the Makefiles needed to build the model.
```
cmake -B build -DBASEDIR=$BASEDIR/Linux -DCMAKE_Fortran_COMPILER=ifort --install-prefix=$(pwd)/install
```
This will install to a directory `install` parallel to your `build` directory. If you prefer to install elsewhere change
the `--install-prefix` option  and CMake will install there.

##### Building with debugging flags

To build with debugging flags add:
```
-DCMAKE_BUILD_TYPE=Debug
```
to the cmake line.

##### Disabling f2py

To disable building of f2py modules add:
```
-DUSE_F2PY=OFF
```
NOTE: This is really only used for systems that do not (yet) support f2py like those used in CI. As the f2py
code in GMAOpyobs is essential for the code to work, this should not be used in general.

### Compile and install with make

```
cmake --build build --target install -j 6
```

## How to build GMAOpyobs on other systems

Building of f2py codes is not implemented for systems not supported by GMAO.

#### Use git to clone the repository

```
git clone git@github.com:GEOS-ESM/GMAOpyobs.git
```

## Contributing

Please check out our [contributing guidelines](CONTRIBUTING.md).

## License

All files are currently licensed under the Apache-2.0 license, see [`LICENSE`](LICENSE).

Previously, the code was licensed under the [NASA Open Source Agreement, Version 1.3](LICENSE-NOSA).
