// Include FDPS header
#include <particle_simulator.hpp>
// Include the standard C++ headers
#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>
#include <sys/stat.h>

/* Parameters */
const short int Dim = 3;
const PS::F64 SMTH = 1.2;
const PS::U32 OUTPUT_INTERVAL = 10;
const PS::F64 C_CFL = 0.3;

/* Kernel Function */
const PS::F64 pi = atan(1.0) * 4.0;
const PS::F64 kernelSupportRadius = 2.5;

PS::F64 W(const PS::F64vec dr, const PS::F64 h){
   const PS::F64 H = kernelSupportRadius * h;
   const PS::F64 s = sqrt(dr * dr) / H;
   const PS::F64 s1 = (1.0 - s < 0) ? 0 : 1.0 - s;
   const PS::F64 s2 = (0.5 - s < 0) ? 0 : 0.5 - s;
   PS::F64 r_value = pow(s1, 3) - 4.0 * pow(s2, 3);
   //if # of dimension == 3
   r_value *= 16.0 / pi / (H * H * H);
   return r_value;
}

PS::F64vec gradW(const PS::F64vec dr, const PS::F64 h){
   const PS::F64 H = kernelSupportRadius * h;
   const PS::F64 s = sqrt(dr * dr) / H;
   const PS::F64 s1 = (1.0 - s < 0) ? 0 : 1.0 - s;
   const PS::F64 s2 = (0.5 - s < 0) ? 0 : 0.5 - s;
   PS::F64 r_value = - 3.0 * pow(s1, 2) + 12.0 * pow(s2, 2);
   //if # of dimension == 3
   r_value *= 16.0 / pi / (H * H * H);
   //return dr * r_value / (sqrt(dr * dr) * H + 1.0e-6 * h);
   return dr * (r_value / (sqrt(dr * dr) * H + 1.0e-6 * h));
}

/* Class Definitions */
//** Force Class (Result Class)
class Dens{
   public:
   PS::F64 dens;
   PS::F64 smth;
   void clear(){
      dens = 0;
   }
};
class Hydro{
   public:
   PS::F64vec acc;
   PS::F64 eng_dot;
   PS::F64 dt;
   void clear(){
      acc = 0;
      eng_dot = 0;
   }
};

//** Full Particle Class
struct FP{
   PS::F64 mass;
   PS::F64vec pos;
   PS::F64vec vel;
   PS::F64vec acc;
   PS::F64 dens;
   PS::F64 eng;
   PS::F64 pres;
   PS::F64 smth;
   PS::F64 snds;
   PS::F64 eng_dot;
   PS::F64 dt;
   PS::S64 id;
   PS::F64vec vel_half;
   PS::F64 eng_half;
   void copyFromForce(const Dens& dens){
      this->dens = dens.dens;
   }
   void copyFromForce(const Hydro& force){
      this->acc     = force.acc;
      this->eng_dot = force.eng_dot;
      this->dt      = force.dt;
   }
   PS::F64 getCharge() const{
      return this->mass;
   }
   PS::F64vec getPos() const{
      return this->pos;
   }
   PS::F64 getRSearch() const{
      return kernelSupportRadius * this->smth;
   }
   void setPos(const PS::F64vec& pos){
      this->pos = pos;
   }
   void writeAscii(FILE* fp) const{
      fprintf(fp,
              "%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t"
              "%lf\t%lf\t%lf\t%lf\t%lf\n",
              this->id, this->mass,
              this->pos.x, this->pos.y, this->pos.z,
              this->vel.x, this->vel.y, this->vel.z,
              this->dens, this->eng, this->pres);
   }
   void readAscii(FILE* fp){
      fscanf(fp,
             "%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t"
             "%lf\t%lf\t%lf\t%lf\t%lf\n",
             &this->id, &this->mass,
             &this->pos.x, &this->pos.y, &this->pos.z,
             &this->vel.x, &this->vel.y, &this->vel.z,
             &this->dens, &this->eng, &this->pres);
   }
   void setPressure(){
      const PS::F64 hcr = 1.4;
      pres = (hcr - 1.0) * dens * eng;
      snds = sqrt(hcr * pres / dens);
   }
};

//** Essential Particle Class
struct EP{
   PS::F64vec pos;
   PS::F64vec vel;
   PS::F64    mass;
   PS::F64    smth;
   PS::F64    dens;
   PS::F64    pres;
   PS::F64    snds;
   void copyFromFP(const FP& rp){
      this->pos  = rp.pos;
      this->vel  = rp.vel;
      this->mass = rp.mass;
      this->smth = rp.smth;
      this->dens = rp.dens;
      this->pres = rp.pres;
      this->snds = rp.snds;
   }
   PS::F64vec getPos() const{
      return this->pos;
   }
   PS::F64 getRSearch() const{
      return kernelSupportRadius * this->smth;
   }
   void setPos(const PS::F64vec& pos){
      this->pos = pos;
   }
};

extern "C" int calc_density_prepare_epj_for_force(const EP* const sorted_epj, const PS::S32* const adr_epj_for_force, const PS::S32 Njp,
                                EP* const epj_for_force) {
      for(PS::S32 i = 0; i < Njp ; i++ ) {
          epj_for_force[i] = sorted_epj[adr_epj_for_force[i]];
      }

      return (int)Njp;
}

extern "C" int calc_hydro_prepare_epj_for_force(const EP* const sorted_epj, const PS::S32* const adr_epj_for_force, const PS::S32 Njp,
                                EP* const epj_for_force) {
      for(PS::S32 i = 0; i < Njp ; i++ ) {
          epj_for_force[i] = sorted_epj[adr_epj_for_force[i]];
      }

      return (int)Njp;
}

extern "C" int calc_density(const EP* const ep_i, const PS::S32 Nip,
                     const EP* const ep_j, const PS::S32 Njp, Dens* const dens){
      for(PS::S32 i = 0 ; i < Nip ; ++i){
         dens[i].clear();
         for(PS::S32 j = 0 ; j < Njp ; ++j){
            const PS::F64vec dr = ep_j[j].pos - ep_i[i].pos;
            dens[i].dens += ep_j[j].mass * W(dr, ep_i[i].smth);
         }
      }

      return (int)Nip;
}

extern "C" int calc_hydro(const EP* const ep_i, const PS::S32 Nip,
                     const EP* const ep_j, const PS::S32 Njp, Hydro* const hydro){
      for(PS::S32 i = 0; i < Nip ; ++ i){
         hydro[i].clear();
         PS::F64 v_sig_max = 0.0;
         for(PS::S32 j = 0; j < Njp ; ++j){
            const PS::F64vec dr = ep_i[i].pos - ep_j[j].pos;
            const PS::F64vec dv = ep_i[i].vel - ep_j[j].vel;
            const PS::F64 w_ij = (dv * dr < 0) ? dv * dr / sqrt(dr * dr) : 0;
            const PS::F64 v_sig = ep_i[i].snds + ep_j[j].snds - 3.0 * w_ij;
            v_sig_max = std::max(v_sig_max, v_sig);
            const PS::F64 AV = - 0.5 * v_sig * w_ij / (0.5 * (ep_i[i].dens + ep_j[j].dens));
            const PS::F64vec gradW_ij = 0.5 * (gradW(dr, ep_i[i].smth) + gradW(dr, ep_j[j].smth));
            hydro[i].acc     -= ep_j[j].mass * (ep_i[i].pres / (ep_i[i].dens * ep_i[i].dens) + ep_j[j].pres / (ep_j[j].dens * ep_j[j].dens) + AV) * gradW_ij;
            hydro[i].eng_dot += ep_j[j].mass * (ep_i[i].pres / (ep_i[i].dens * ep_i[i].dens) + 0.5 * AV) * dv * gradW_ij;
         }
         hydro[i].dt = C_CFL * 2.0 * ep_i[i].smth / v_sig_max;
      }

      return (int)Nip;
}

