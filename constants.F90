!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module constants

!BOP
! !MODULE: constants
!
! !DESCRIPTION:
!  This module defines a variety of physical and numerical constants
!  used throughout the Parallel Ocean Program.
!
! !REVISION HISTORY:
!  SVN:$Id: constants.F90 24379 2010-08-13 19:54:51Z njn01 $

! !USES:

   use kinds_mod

   implicit none
   private
   save

! !DEFINED PARAMETERS:

   ! numbers

   real (r8), parameter, public :: &
      c0     =    0.0_r8   ,&
      c1     =    1.0_r8   ,&
      c2     =    2.0_r8   ,&
      c3     =    3.0_r8   ,&
      c4     =    4.0_r8   ,&
      c5     =    5.0_r8   ,&
      c8     =    8.0_r8   ,&
      c10    =   10.0_r8   ,&
      c16    =   16.0_r8   ,&
      c1000  = 1000.0_r8   ,&
      c10000 =10000.0_r8   ,&
      c1p5   =    1.5_r8   ,&
      p33    = c1/c3       ,&
      p5     = 0.500_r8    ,&
      p25    = 0.250_r8    ,&
      p125   = 0.125_r8    ,&
      p001   = 0.001_r8    ,&
      eps    = 1.0e-10_r8  ,&
      eps2   = 1.0e-20_r8  ,&
      bignum = 1.0e+30_r8


 end module constants

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
