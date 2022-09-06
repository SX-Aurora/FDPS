# How to implement FDPS applications with VEOffload
The followings are how to implement FDPS applications with VEOffload using C++ SPH sample as an example.

## Shared library runs on VE
1\. Preparing EPJ particles

Function to get particles utilized on user defined interaction function from all "EssentialParticleJ" particles is necessary. Which function is utilized depends on the type of "EPJ" particles.

For example, the following is the function for density "EPJ" particles.

```
extern "C" int calc_density_prepare_epj_for_force(const EP* const sorted_epj,
                                                  const PS::S32* const adr_epj_for_force,
                                                  const PS::S32 Njp,
                                                  EP* const epj_for_force) {
      for(PS::S32 i = 0; i < Njp ; i++ ) {
          epj_for_force[i] = sorted_epj[adr_epj_for_force[i]];
      }

      return (int)Njp;
}
```

Here,

```
Argument 1, sorted_epj: Array of sorted whole EPJ particles.
Argument 2, adr_epj_for_force: Array of gotten EPJ particle index.
Argument 3, Njp: Number of gotten EPJ particles.
Argument 4, epj_for_force: Array to store gotten EPJ particles.
```

2\. User defined interaction function

User defined interaction function is defined as follows.

```
extern "C" int calc_density(const EP* const ep_i, const PS::S32 Nip,
                 const EP* const ep_j, const PS::S32 Njp, Dens* const dens);

extern "C" int calc_hydro(const EP* const ep_i, const PS::S32 Nip,
                 const EP* const ep_j, const PS::S32 Njp, Hydro* const hydro);
```

Here,

```
Argument 1, ep_i: Array of EPI particles.
Argument 2, Nip: Number of EPI particles.
Argument 3, ep_j: Array of EPJ particles.
Argument 4, Njp: Number of EPJ particles.
Argument 5, dens/hydro: Array of force as result.
```

Force array has to be cleared in advance.

3\. Compiling

To make shared library runs on VE, compile the sources with NEC compiler (ncc). See "Makefile_g++" and etc. about how to make.

## Main program
1\. Caller function for user defined interaction function

Define the class inherits "VEO::FunctorEpEp".

```
class CalcDensity : public VEO::FunctorEpEp<EP,EP> {
    public:
    CalcDensity(const char* pFuncPrepareEpjSymbol, const char* pFuncEpEpSymbol)
             : VEO::FunctorEpEp<EP,EP>(pFuncPrepareEpjSymbol, pFuncEpEpSymbol) {}
};

class CalcHydroForce : public VEO::FunctorEpEp<EP,EP> {
    public:
    CalcHydroForce(const char* pFuncPrepareEpjSymbol, const char* pFuncEpEpSymbol)
             : VEO::FunctorEpEp<EP,EP>(pFuncPrepareEpjSymbol, pFuncEpEpSymbol) {}
};
```

Definition for each tree is necessary. The first argument of "VEO::FunctorEpEp" is the type of "EPI" particles and the second argument is the type of "EPJ" particles.

2\. Initialization of VEOffload interface

The following function call is necessary just after "PS::Initialize" is called.

```
#if defined(_SX_VEOFFLOAD)
    VEO::VeoInterface::initialize((const char *)argv[0],"libvesph.so");
#endif
```

The name of shared library runs on VE is specified on the second argument.

3\. Calling offloaded functions

Instantiate functions to call offloaded functions. 

```
#if defined(_SX_VEOFFLOAD)
    CalcDensity funcCalcDensity("calc_density_prepare_epj_for_force","calc_density");
    CalcHydroForce funcCalcHydroForce("calc_hydro_prepare_epj_for_force","calc_hydro");
#endif
```

The first argument is the name of the function to get "EPJ" particles and the second is the name of the function of user defined interaction function.

These instances are specified on the first argument of "calcForceAllAndWriteBack".

```
dens_tree.calcForceAllAndWriteBack(funcCalcDensity, sph_system, dinfo);
hydr_tree.calcForceAllAndWriteBack(funcCalcHydroForce, sph_system, dinfo);
```

4\. Cleanup user defined interaction functions

After the main loop, "cleanup" function has to be called.

```
#if defined(_SX_VEOFFLOAD)
    funcCalcDensity.cleanup();
    funcCalcHydroForce.cleanup();
#endif
```

5\. Disposing VEOffload interface

After the main loop, "disposeInterface" function has to be called.

```
#if defined(_SX_VEOFFLOAD)
    VEO::VeoInterface::disposeInterface();
#endif
```

6\. Compiling

To compile, specify "-D_SX_VEOFFLOAD" and include "FDPS-SX-Aurora/src" and "/opt/nec/ve/veos/include" into include path. The following option is necessary to link.

```
-L/opt/nec/ve/veos/lib64 -Wl,-rpath=/opt/nec/ve/veos/lib64 -lveo
```

See "Makefile_g++" and etc. about how to make.

## Functions for long-range interaction
SPH sample is of short-range interaction without super-particle. But N-body sample is of long-range interaction  with super-particle.

The followings are about the differences from SPH (short-range interaction) case.

1\. User defined interaction functions

User defined interaction functions are defined as follows in N-Body sample.

```
extern "C" int calc_gravity(const FPGrav * ep_i,
             const PS::S32 n_ip,
             const FPGrav * ep_j,
             const PS::S32 n_jp,
             FPGrav * force) ;

extern "C" int calc_gravity_spj(const FPGrav * ep_i,
             const PS::S32 n_ip,
             const PS::SPJMonopole * ep_j,
             const PS::S32 n_jp,
             FPGrav * force) ;
```

The first function is of short-range interaction and the second is of long-range interaction. Here,

```
Argument 1, ep_i: Array of EPI particles.
Argument 2, Nip: Number of EPI particles.
Argument 3, ep_j: Array of EPJ particles on the first function. Array of "SuperParticleJ" particles on the second function.
Argument 4, Njp: Number of EPJ particles.
Argument 5, dens/hydro: Array of force as result.
```

2\. User defined interaction function

Define the class inherits "VEO::FunctorEpEp" for short-range interaction.

```
class CalcGravity : public VEO::FunctorEpEp<FPGrav,FPGrav> {
    public:
    CalcGravity(const char* pFuncPrepareEpjSymbol, const char* pFuncEpEpSymbol)
            : VEO::FunctorEpEp<FPGrav,FPGrav>(pFuncPrepareEpjSymbol, pFuncEpEpSymbol) {}
};
```

The first argument of "VEO::FunctorEpEp" is the type of "EPI" particles and the second argument is the type of "EPJ" particles.

On the other hand, define the class inherits "VEO::FunctorEpSp" for long-range interaction.

```
class CalcGravitySpj : public VEO::FunctorEpSp<FPGrav,FPGrav,PS::SPJMonopole> {
    public:
    CalcGravitySpj(const char* pFuncPrepareSpjSymbol, const char* pFuncEpEpSymbol, CalcGravity& calcGravity)
            : VEO::FunctorEpSp<FPGrav,FPGrav,PS::SPJMonopole>(pFuncPrepareSpjSymbol, pFuncEpEpSymbol, calcGravity) {}
};
```

The first argument of "VEO::FunctorEpSp" is the type of "EPI" particles and the second argument is the type of "SuperParticleJ" particles.

Specify the reference to the function for short-range interaction force calculation into the third argument.

3\. Calling offloaded functions

Instantiate functions to call offloaded functions. 

```
#if defined(_SX_VEOFFLOAD)
    CalcGravity          funcCalcGravity("prepare_epj_for_force","calc_gravity");
    CalcGravitySpj       funcCalcGravitySpj("prepare_spj_for_force","calc_gravity_spj",funcCalcGravity);
#endif
```

These instances are specified on the first and the second arguments of "calcForceAllAndWriteBack".

```
tree_grav.calcForceAllAndWriteBack(funcCalcGravity,
                                       funcCalcGravitySpj,
                                       system_grav,
                                       dinfo);
```
