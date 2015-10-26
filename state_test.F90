program state_dummy

      use kinds_mod
      use constants
      use omp_lib

implicit none
 
      integer (int_kind) k,this_block  
      
      real (r8), dimension(nx_block,ny_block,km,3) :: TMIX, TRCR  

      real (r8), dimension(nx_block,ny_block,km) :: RHOK1 

      !integer (int_kind) :: &
      !k=1 ,this_block=1

      call state(k,1,TRCR(:,:,:,1), TRCR(:,:,:,2), this_block, RHOFULL=RHOK1) 

contains

      subroutine state(k, kk, TEMPK, SALTK, this_block, RHOOUT, RHOFULL, DRHODT, DRHODS)

      implicit none 

      integer (int_kind), parameter :: &
      nx_block = 84, & 
      ny_block = 196, &
      km = 60

      integer (int_kind) i,j,this_block  

      real (r8) start_time, end_time

      real (r8), parameter ::                  &
      mwjfnp0s0t0 =   9.99843699e+2_r8 * p001, &
      mwjfnp0s0t1 =   7.35212840e+0_r8 * p001, &
      mwjfnp0s0t2 =  -5.45928211e-2_r8 * p001, &
      mwjfnp0s0t3 =   3.98476704e-4_r8 * p001, &
      mwjfnp0s1t0 =   2.96938239e+0_r8 * p001, &
      mwjfnp0s1t1 =  -7.23268813e-3_r8 * p001, &
      mwjfnp0s2t0 =   2.12382341e-3_r8 * p001, &
      mwjfnp1s0t0 =   1.04004591e-2_r8 * p001, &
      mwjfnp1s0t2 =   1.03970529e-7_r8 * p001, &
      mwjfnp1s1t0 =   5.18761880e-6_r8 * p001, &
      mwjfnp2s0t0 =  -3.24041825e-8_r8 * p001, &
      mwjfnp2s0t2 =  -1.23869360e-11_r8* p001

   !*** these constants will be used to construct the denominator

      real (kind=r8), parameter ::       &
      mwjfdp0s0t0 =   1.0e+0_r8,         &
      mwjfdp0s0t1 =   7.28606739e-3_r8,  &
      mwjfdp0s0t2 =  -4.60835542e-5_r8,  &
      mwjfdp0s0t3 =   3.68390573e-7_r8,  &
      mwjfdp0s0t4 =   1.80809186e-10_r8, &
      mwjfdp0s1t0 =   2.14691708e-3_r8,  &
      mwjfdp0s1t1 =  -9.27062484e-6_r8,  &
      mwjfdp0s1t3 =  -1.78343643e-10_r8, &
      mwjfdp0sqt0 =   4.76534122e-6_r8,  &
      mwjfdp0sqt2 =   1.63410736e-9_r8,  &
      mwjfdp1s0t0 =   5.30848875e-6_r8,  &
      mwjfdp2s0t3 =  -3.03175128e-16_r8, &
      mwjfdp3s0t1 =  -1.27934137e-17_r8
 

      integer (int_kind) ,intent(in) :: &
      k,                    &! depth level index
      kk                     ! level to which water is adiabatically displaced

      real (r8), dimension(nx_block,ny_block,km), optional, intent(out) :: & 
      RHOOUT,  &! perturbation density of water
      RHOFULL, &! full density of water
      DRHODT,  &! derivative of density with respect to temperature
      DRHODS    ! derivative of density with respect to salinity


      real (r8), dimension(nx_block,ny_block,km), intent(in) :: & 
      TEMPK,             &! temperature at level k
      SALTK,             &! salinity    at level k
      TQ,SQ,             &! adjusted T,S
      SQR,DENOMK,        &
      WORK1, WORK2, WORK3, WORK4  

      real (r8), dimension(km) :: pressz 
        

      real (r8) ::      & 
      tmin, tmax,       &! valid temperature range for level k
      smin, smax        ! valid salinity    range for level k
                        ! ref pressure (bars) at each level

     real (r8) :: p, p2, &! temporary pressure scalars 
      mwjfnums0t0, mwjfnums0t1, mwjfnums0t2, mwjfnums0t3,              &
      mwjfnums1t0, mwjfnums1t1, mwjfnums2t0,                           &
      mwjfdens0t0, mwjfdens0t1, mwjfdens0t2, mwjfdens0t3, mwjfdens0t4, &
      mwjfdens1t0, mwjfdens1t1, mwjfdens1t3,                           &
      mwjfdensqt0, mwjfdensqt2
 
 
      tmin =  -2.0_r8  ! limited   on the low  end
      tmax = 999.0_r8  ! unlimited on the high end
      smin =   0.0_r8  ! limited   on the low  end
      smax = 0.999_r8  ! unlimited on the high end


      call random_number(TQ)
      call random_number(SQ)
      call random_number(TEMPK)
      call random_number(SALTK)
      call random_number(pressz)
 
      !***
      !*** first calculate numerator of MWJF density [P_1(S,T,p)]
      !***

      mwjfnums0t1 = mwjfnp0s0t1 
      mwjfnums0t3 = mwjfnp0s0t3
      mwjfnums1t1 = mwjfnp0s1t1
      mwjfnums2t0 = mwjfnp0s2t0

      mwjfdens0t2 = mwjfdp0s0t2
      mwjfdens0t4 = mwjfdp0s0t4
      mwjfdens1t0 = mwjfdp0s1t0
      mwjfdens1t1 = mwjfdp0s1t1
      mwjfdens1t3 = mwjfdp0s1t3
      mwjfdensqt0 = mwjfdp0sqt0
      mwjfdensqt2 = mwjfdp0sqt2

      mwjfnums0t0 = mwjfnp0s0t0 + p*(mwjfnp1s0t0 + p*mwjfnp2s0t0)
      mwjfnums0t2 = mwjfnp0s0t2 + p*(mwjfnp1s0t2 + p*mwjfnp2s0t2)
      mwjfnums1t0 = mwjfnp0s1t0 + p*mwjfnp1s1t0
     
      mwjfdens0t0 = mwjfdp0s0t0 + p*mwjfdp1s0t0
      mwjfdens0t1 = mwjfdp0s0t1 + p**3 * mwjfdp3s0t1
      mwjfdens0t3 = mwjfdp0s0t3 + p**2 * mwjfdp2s0t3 

      start_time = omp_get_wtime()

      do k=1,km
      do j=1,ny_block
      do i=1,nx_block
      TQ(i,j,k) = min(TEMPK(i,j,k),tmax)
      TQ(i,j,k) = max(TQ(i,j,k),tmin)
      SQ(i,j,k) = min(SALTK(i,j,k),smax)
      SQ(i,j,k) = max(SQ(i,j,k),smin)

      p   = c10*pressz(kk)

      SQ(i,j,k)  = c1000 * SQ(i,j,k)
      SQR(i,j,k) = sqrt( SQ(i,j,k) )

      WORK1(i,j,k) = mwjfnums0t0 + TQ(i,j,k) * (mwjfnums0t1 + TQ(i,j,k) * (mwjfnums0t2 + &
              mwjfnums0t3 * TQ(i,j,k) )) + SQ(i,j,k) * (mwjfnums1t0 +              &
              mwjfnums1t1 * TQ(i,j,k) + mwjfnums2t0 * SQ(i,j,k) )

      WORK2(i,j,k) = mwjfdens0t0 + TQ(i,j,k) * (mwjfdens0t1 + TQ(i,j,k) * (mwjfdens0t2 +    &
           TQ(i,j,k) * (mwjfdens0t3 + mwjfdens0t4 * TQ(i,j,k) ))) +                   &
           SQ(i,j,k) * (mwjfdens1t0 + TQ(i,j,k) * (mwjfdens1t1 + TQ(i,j,k) * TQ(i,j,k) *mwjfdens1t3)+ &
           SQR(i,j,k) * (mwjfdensqt0 + TQ(i,j,k) * TQ(i,j,k) *mwjfdensqt2))


      DENOMK(i,j,k) = c1/WORK2(i,j,k)

      if (present(RHOOUT)) then
      RHOOUT(i,j,k)  = WORK1(i,j,k) * DENOMK(i,j,k)
      endif 

      if (present(RHOFULL)) then
      RHOFULL(i,j,k) = WORK1(i,j,k) * DENOMK(i,j,k)
      endif

      if(present(DRHODT))then 
      WORK3(i,j,k) = &! dP_1/dT
                 mwjfnums0t1 + TQ(i,j,k) * (c2*mwjfnums0t2 +    &
                 c3*mwjfnums0t3 * TQ(i,j,k)) + mwjfnums1t1 * SQ(i,j,k)

      WORK4(i,j,k) = &! dP_2/dT
                 mwjfdens0t1 + SQ(i,j,k) * mwjfdens1t1 +               &
                 TQ(i,j,k) * (c2*(mwjfdens0t2 + SQ(i,j,k) * SQR(i,j,k) * mwjfdensqt2) +  &
                 TQ(i,j,k) * (c3*(mwjfdens0t3 + SQ(i,j,k) * mwjfdens1t3) +    &
                 TQ(i,j,k) *  c4*mwjfdens0t4))

      DRHODT(i,j,k) = (WORK3(i,j,k) - WORK1(i,j,k) * DENOMK(i,j,k) * WORK4(i,j,k)) * DENOMK(i,j,k)
      endif 

      if(present(DRHODS))then
      WORK3(i,j,k) = &! dP_1/dS
                 mwjfnums1t0 + mwjfnums1t1 * TQ(i,j,k) + c2*mwjfnums2t0 * SQ(i,j,k)

      WORK4(i,j,k) = mwjfdens1t0 +   &! dP_2/dS
                 TQ(i,j,k) * (mwjfdens1t1 + TQ(i,j,k) * TQ(i,j,k) * mwjfdens1t3) +   &
                 c1p5* SQR(i,j,k) * (mwjfdensqt0 + TQ(i,j,k) * TQ(i,j,k) *mwjfdensqt2)

      DRHODS(i,j,k) = (WORK3(i,j,k) - WORK1(i,j,k) * DENOMK(i,j,k) * WORK4(i,j,k)) * DENOMK(i,j,k) * c1000
      endif 

         enddo
         enddo
         enddo
 
         end_time = omp_get_wtime() 

         print *, end_time - start_time 

         print *,sum(WORK1)
         print *,sum(WORK2)
         print *,sum(WORK3)
         print *,sum(WORK4)
         print *,sum(DRHODT)
         print *,sum(DRHODS)
         print *,sum(RHOFULL)
         print *,sum(RHOOUT)

         end subroutine state  

end program state_dummy

