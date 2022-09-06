// Include FDPS header
#include <particle_simulator.hpp>
// Include the standard C++ headers
#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>
#include <sys/stat.h>
//#include <ftrace.h>
#include <veo_hmem.h>

#include "user-defined.hpp"


//#define F_EPS (1.0/32.0)
PS::F64 FPGrav::eps = 1.0/32.0;

extern "C" int prepare_epj_for_force(const FPGrav* const sorted_epj, const PS::S32* const adr_epj_for_force, const PS::S32 Njp,
                                FPGrav* const epj_for_force) {
      for(PS::S32 i = 0; i < Njp ; i++ ) {
          epj_for_force[i] = sorted_epj[adr_epj_for_force[i]];
      }

      return (int)Njp;
}

extern "C" int prepare_spj_for_force(const PS::SPJMonopole* const sorted_spj, const PS::S32* const adr_spj_for_force, const PS::S32 Njp,
                                PS::SPJMonopole* const spj_for_force) {
      for(PS::S32 i = 0; i < Njp ; i++ ) {
          spj_for_force[i] = sorted_spj[adr_spj_for_force[i]];
      }

      return (int)Njp;
}

extern "C" int calc_gravity(const FPGrav * ep_i,
                 const PS::S32 n_ip,
                 const FPGrav * ep_j,
                 const PS::S32 n_jp,
                 FPGrav * force) {

	const PS::F64 eps2 = FPGrav::eps * FPGrav::eps;
	for(PS::S32 i = 0; i < n_ip; i++){
		//PS::F64vec xi = ep_i[i].getPos();
		PS::F64vec ai = 0.0;
		PS::F64 poti = 0.0;
		//force[i].acc = 0.0;
		//force[i].pot = 0.0;
		for(PS::S32 j = 0; j < n_jp; j++){
			PS::F64vec rij    = ep_i[i].getPos() - ep_j[j].getPos();
			//PS::F64vec rij    = xi - ep_j[j].getPos();
			//PS::F64    r3_inv = rij * rij + FPGrav::eps * FPGrav::eps;
			PS::F64    r3_inv = rij * rij + eps2;
			PS::F64    r_inv  = 1.0/sqrt(r3_inv);
			r3_inv  = r_inv * r_inv;
			r_inv  *= ep_j[j].getCharge();
			r3_inv *= r_inv;
			ai -= r3_inv * rij;
            		poti -= r_inv;
        	}
		force[i].acc = ai;
		force[i].pot = poti;
	}

	return (int)n_ip;
}


extern "C" int calc_gravity_spj(const FPGrav * ep_i,
                 const PS::S32 n_ip,
                 const PS::SPJMonopole * ep_j,
                 const PS::S32 n_jp,
                 FPGrav * force) {

	const PS::F64 eps2 = FPGrav::eps * FPGrav::eps;
	for(PS::S32 i = 0; i < n_ip; i++){
		//PS::F64vec xi = ep_i[i].getPos();
		PS::F64vec ai = 0.0;
		PS::F64 poti = 0.0;
		for(PS::S32 j = 0; j < n_jp; j++){
			PS::F64vec rij    = ep_i[i].getPos() - ep_j[j].getPos();
			//PS::F64vec rij    = xi - ep_j[j].getPos();
			PS::F64    r3_inv = rij * rij + eps2;
			PS::F64    r_inv  = 1.0/sqrt(r3_inv);
			r3_inv  = r_inv * r_inv;
			r_inv  *= ep_j[j].getCharge();
			r3_inv *= r_inv;
			ai += r3_inv * rij;
            		poti += r_inv;
        	}
		force[i].acc -= ai;
		force[i].pot -= poti;
	}

	return (int)n_ip;
}
