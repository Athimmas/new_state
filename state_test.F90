program state_dummy

      use kinds_mod
      use constants
      use omp_lib
      use state_mod

implicit none

      integer (int_kind), parameter :: &
      nx_block = 84, &
      ny_block = 196, &
      km = 60

      integer (int_kind) :: k,this_block,kk  
      
      real (r8), dimension(nx_block,ny_block,km,3) :: TMIX, TRCR  

      real (r8), dimension(nx_block,ny_block,km) :: RHOK1,RHOK2,RHOK3,RHOK4

      !integer (int_kind) :: &
      !k=1 ,this_block=1
      kk=1
      call random_number(TMIX)
      call random_number(TRCR) 

      !dir$ offload begin target(mic:0) 
      print *,"here"
      call state(kk,kk,TRCR(:,:,:,1), TRCR(:,:,:,2), this_block,RHOOUT=RHOK1,RHOFULL=RHOK2, DRHODT=RHOK3, DRHODS=RHOK4) 
      print *,"ended state"
      !dir$ end offload

      end program state_dummy

