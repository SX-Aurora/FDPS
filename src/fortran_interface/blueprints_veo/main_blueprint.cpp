extern "C" void f_main_();

/* Standard headers */
#include <iostream>
#include <fstream>
/* FDPS headers */
#define _SX_FORTRAN_VEO
#include <particle_simulator.hpp> 
/* User-defined headers */
#include "FDPS_Manipulators.h"

// VeoInterfaceクラスのシングルトンインスタンス
VEO::VeoInterface* VEO::VeoInterface::pInstance_ = NULL;

int main(int argc, char *argv[])
{
   
   //* Initialize fdps_manip
   FDPS_Manipulators::Initialize(argc,argv);
#if defined(_SX_VEOFFLOAD)
   FDPS_Manipulators::PS_Initialize();
/*
   // Initialize VEOffloading Interface
   VEO::VeoInterface::initialize((const char *)argv[0],"xxx.so");
*/
#endif
   //* Call Fortran main subroutine
   f_main_();

#if defined(_SX_VEOFFLOAD)
   VEO::VeoInterface::disposeInterface();
#endif

   return 0;

}
