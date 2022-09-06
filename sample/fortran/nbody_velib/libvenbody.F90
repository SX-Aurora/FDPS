!===============================
!   MODULE: User defined types
!===============================
module user_defined_types
   use, intrinsic :: iso_c_binding
   !use fdps_vector
   !use fdps_super_particle
   implicit none

   !* Public variables
   real(kind=c_double), public :: eps_grav ! gravitational softening

   !**** PS::F32vec
   type, public, bind(c) :: fdps_f32vec
      real(kind=c_float) :: x,y,z
   end type fdps_f32vec

   !**** PS::F64vec
   type, public, bind(c) :: fdps_f64vec
      real(kind=c_double) :: x,y,z
   end type fdps_f64vec

   !**** PS::SPJMonopole
   type, public, bind(c) :: fdps_spj_monopole
      !real(kind=c_float) :: mass
      !type(fdps_f32vec) :: pos
      real(kind=c_double) :: mass
      type(fdps_f64vec) :: pos
   end type fdps_spj_monopole

   !**** Full particle type
   type, public, bind(c) :: full_particle !$fdps FP,EPI,EPJ,Force
      !$fdps copyFromForce full_particle (pot,pot) (acc,acc)
      !$fdps copyFromFP full_particle (id,id) (mass,mass) (pos,pos) 
      !$fdps clear id=keep, mass=keep, pos=keep, vel=keep
      integer(kind=c_long_long) :: id
      real(kind=c_double)  mass !$fdps charge
      type(fdps_f64vec) :: pos !$fdps position
      type(fdps_f64vec) :: vel !$fdps velocity
      real(kind=c_double) :: pot
      type(fdps_f64vec) :: acc
   end type full_particle

   !* The following types are used in PIKG-generated kenrels
   type, public, bind(c) :: epi_grav
      type(fdps_f32vec) :: pos
   end type epi_grav

   type, public, bind(c) :: epj_grav
      type(fdps_f32vec) :: pos
      real(kind=c_float) :: mass
   end type epj_grav

   type, public, bind(c) :: force_grav
      type(fdps_f32vec) :: acc
      real(kind=c_float) :: pot
   end type force_grav
  
   contains

   !**** Interaction function (particle-particle)
   subroutine prepare_epj_for_force(sorted_epj, adr_epj_for_force, n_jp, epj_for_force) bind(c)
      integer(kind=c_int), intent(in), value :: n_jp
      type(full_particle), dimension(n_jp), intent(in) :: sorted_epj
      integer(kind=c_int), dimension(n_jp), intent(in) :: adr_epj_for_force
      type(full_particle), dimension(n_jp), intent(inout) :: epj_for_force
      integer(kind=c_int) i
      !integer(kind=c_int) j

      !print '("prepare_epj_for_force n_jp=",i0)', n_jp
      do i=1,n_jp
         !print '("adr_epj_for_force(", i0, ")=",i0)', i , adr_epj_for_force(i)
         epj_for_force(i) = sorted_epj( adr_epj_for_force(i) + 1 )
      end do

      !print '("epj_for_force(", i0, ")%pos%x=",f10.8)',  n_jp, epj_for_force(n_jp)%pos%x
   end subroutine prepare_epj_for_force

   subroutine prepare_spj_for_force(sorted_spj, adr_spj_for_force, n_jp, sp_j) bind(c)
      integer(kind=c_int), intent(in), value :: n_jp
      type(fdps_spj_monopole), dimension(n_jp), intent(in) :: sorted_spj
      integer(kind=c_int), dimension(n_jp), intent(in) :: adr_spj_for_force
      type(fdps_spj_monopole), dimension(n_jp), intent(inout) :: sp_j
      integer(kind=c_int) i
      !integer(kind=c_int) j

      !print '("prepare_spj_for_force n_jp=",i0)', n_jp
      do i=1,n_jp
         !print '("adr_spj_for_force(", i0, ")=",i0)', i , adr_spj_for_force(i)
         sp_j(i) = sorted_spj( adr_spj_for_force(i) + 1 )
         !spj_for_force(i)%pos%x = sorted_spj( adr_spj_for_force(i) + 1 )%pos%x
         !spj_for_force(i)%pos%y = sorted_spj( adr_spj_for_force(i) + 1 )%pos%y
         !spj_for_force(i)%pos%z = sorted_spj( adr_spj_for_force(i) + 1 )%pos%z
      	!print '("spj_for_force(", i0, ")%pos%x=",f10.8," y=",f10.8)',  i, spj_for_force(i)%pos%x, spj_for_force(i)%pos%y
         sp_j(i)%mass = sorted_spj( adr_spj_for_force(i) + 1 )%mass
      end do

   end subroutine prepare_spj_for_force

   subroutine calc_gravity_ep_ep(ep_i,n_ip,ep_j,n_jp,f) bind(c)
      implicit none
      integer(c_int), intent(in), value :: n_ip,n_jp
      type(full_particle), dimension(n_ip), intent(in) :: ep_i
      type(full_particle), dimension(n_jp), intent(in) :: ep_j
      type(full_particle), dimension(n_ip), intent(inout) :: f
      !* Local variables
      integer(c_int) :: i,j
      real(c_double) :: eps2,poti,r3_inv,r_inv
      type(fdps_f64vec) :: xi,ai,rij

      eps_grav = 1.0d0/32.0d0
      eps2 = eps_grav * eps_grav
      do i=1,n_ip
         f(i)%pot = 0.0d0
         f(i)%acc%x = 0.0d0
         f(i)%acc%y = 0.0d0
         f(i)%acc%z = 0.0d0
         xi = ep_i(i)%pos
         !ai = 0.0d0
         ai%x = 0.0d0
         ai%y = 0.0d0
         ai%z = 0.0d0
         poti = 0.0d0
         do j=1,n_jp
            rij%x  = xi%x - ep_j(j)%pos%x
            rij%y  = xi%y - ep_j(j)%pos%y
            rij%z  = xi%z - ep_j(j)%pos%z
            r3_inv = rij%x*rij%x &
                   + rij%y*rij%y &
                   + rij%z*rij%z &
                   + eps2
            r_inv  = 1.0d0/sqrt(r3_inv)
            r3_inv = r_inv * r_inv
            r_inv  = r_inv * ep_j(j)%mass
            r3_inv = r3_inv * r_inv
            ai%x   = ai%x - r3_inv * rij%x
            ai%y   = ai%y - r3_inv * rij%y
            ai%z   = ai%z - r3_inv * rij%z
            poti   = poti - r_inv
            ! [IMPORTANT NOTE]
            !   In the innermost loop, we use the components of vectors
            !   directly for vector operations because of the following
            !   reasion. Except for intel compilers with `-ipo` option,
            !   most of Fortran compilers use function calls to perform
            !   vector operations like rij = x - ep_j(j)%pos.
            !   This significantly slow downs the speed of the code.
            !   By using the components of vector directly, we can avoid 
            !   these function calls.
         end do
         f(i)%pot = f(i)%pot + poti
         !f(i)%acc = f(i)%acc + ai
         f(i)%acc%x = f(i)%acc%x + ai%x
         f(i)%acc%y = f(i)%acc%y + ai%y
         f(i)%acc%z = f(i)%acc%z + ai%z
      end do

   end subroutine calc_gravity_ep_ep

   !**** Interaction function (particle-super particle)
   subroutine calc_gravity_ep_sp(ep_i,n_ip,sp_j,n_jp,f) bind(c)
      implicit none
      integer(c_int), intent(in), value :: n_ip,n_jp
      type(full_particle), dimension(n_ip), intent(in) :: ep_i
      type(fdps_spj_monopole), dimension(n_jp), intent(in) :: sp_j
      type(full_particle), dimension(n_ip), intent(inout) :: f
      !* Local variables
      integer(c_int) :: i,j
      real(c_double) :: eps2,poti,r3_inv,r_inv
      type(fdps_f64vec) :: xi,ai,rij

      eps_grav = 1.0d0/32.0d0
      eps2 = eps_grav * eps_grav

      do i=1,n_ip
         xi = ep_i(i)%pos
         ai%x = 0.0d0
         ai%y = 0.0d0
         ai%z = 0.0d0
         poti = 0.0d0
         do j=1,n_jp
            rij%x  = xi%x - sp_j(j)%pos%x
            rij%y  = xi%y - sp_j(j)%pos%y
            rij%z  = xi%z - sp_j(j)%pos%z
            r3_inv = rij%x*rij%x &
                   + rij%y*rij%y &
                   + rij%z*rij%z &
                   + eps2
            r_inv  = 1.0d0/sqrt(r3_inv)
            r3_inv = r_inv * r_inv
            r_inv  = r_inv * sp_j(j)%mass
            r3_inv = r3_inv * r_inv
            ai%x   = ai%x - r3_inv * rij%x
            ai%y   = ai%y - r3_inv * rij%y
            ai%z   = ai%z - r3_inv * rij%z
            poti   = poti - r_inv
         end do
         f(i)%pot = f(i)%pot + poti
         f(i)%acc%x = f(i)%acc%x + ai%x
         f(i)%acc%y = f(i)%acc%y + ai%y
         f(i)%acc%z = f(i)%acc%z + ai%z
      end do

   end subroutine calc_gravity_ep_sp

end module user_defined_types
