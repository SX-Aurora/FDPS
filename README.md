# About
FDPS (https://github.com/FDPS/FDPS) ported and partially optimized version for SX-Aurora TSUBASA Vector Engine. This is based on FDPS version 7.1.

This is based on "VEOffload" function.
On FDPS applications, main program is run on VH and user defined interaction functions are run on VE.

See [HOWTO_SAMPLES](HOWTO_SAMPLES.md) about how to make/run samples and [HOWTO_IMPLEMENTATION](HOWTO_IMPLEMENTATION.md) about how to implement FDPS applications with VEOffload.

# Directories
Directory tee is as a following.

```
FDPS-SX-Aurora/
  |
  +--  src/: Modified source code.
  |      |
  |      +-- fortran_interface/
  |             |
  |             +-- blueprints_new/: Fortran I/F template with time profile modifications.
  |             +-- blueprints_veo/: Fortran I/F template for VEOffload.
  |
  |
  +--  src_base/: Original source code.
  |
  +--  scripts/
  |      |
  |      +--  gen_ftn_if.py: Original Fortran I/F generation script.
  |      +--  gen_ftn_if_new.py: Fortran I/F generation script with time profile modifications.
  |      +--  gen_ftn_if_veo.py: Fortran I/F generation script for VEOffload.
  |
  +--  sample/
         |
         +--  c++/
         |       |
         |       +--  sph/: Original SPH sample.
         |       +--  sph_veoffload/: SPH sample for VEOffload.
         |       +--  nbody/: Original N-body sample.
         |       +--  nbody_veoffload/: N-body sample for VEOffload.
         |
         +--  fortran/
                 |
                 +--  sph/: Original SPH sample.
                 +--  sph_veoffload/: SPH sample for VEOffload.
                 +--  sph_velib/: Shared library of user defined interaction function run on VE.
                 +--  nbody/: Original N-body sample.
                 +--  nbody_veoffload/: N-body sample for VEOffload.
                 +--  nbody_velib/: Shared library of use defined interaction function run on VE.
```

The following source files are modified.

1\. ps_defs.hpp

2\. tree_for_force_impl_force.hpp

3\. veo_defs.hpp
