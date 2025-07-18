! --------------------------------------------------------------------
! Copyright (C) LAPACK
! Copyright (C) 2007 - 2025 - EDF R&D - www.code-aster.org
! This file is part of code_aster.
!
! code_aster is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! code_aster is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
! --------------------------------------------------------------------
!
! ===============================================================
! THIS LAPACK 2.0 ROUTINE IS DEPRECATED
! DO NOT USE IT : YOU SHOULD PREFER UP-TO-DATE LAPACK ROUTINE
!
! BUT DO NOT REMOVE IT :
! THE PRESENT ROUTINE IS MANDATORY FOR ARPACK LIBRARY
! WHICH STICKS TO LAPACK 2.0 VERSION
! ==============================================================
subroutine ar_dlartg(f, g, cs, sn, r)
!
!     SUBROUTINE LAPACK GENERANT UNE ROTATION PLANE QUI EST UNE
!     VERSION PLUS PRECISE QUE LA ROUTINE BLAS1 DROTG.
!---------------------------------------------------------------------
!  -- LAPACK AUXILIARY ROUTINE (VERSION 2.0) --
!     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
!     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
!     SEPTEMBER 30, 1994
!
!  PURPOSE
!  =======
!
!  DLARTG GENERATE A PLANE ROTATION SO THAT
!
!     (  CS  SN  )  .  ( F )  =  ( R )   WHERE CS**2 + SN**2 = 1.
!     ( -SN  CS  )     ( G )     ( 0 )
!
!  THIS IS A SLOWER, MORE ACCURATE VERSION OF THE BLAS1 ROUTINE DROTG,
!  WITH THE FOLLOWING OTHER DIFFERENCES:
!     F AND G ARE UNCHANGED ON RETURN.
!     IF G=0, THEN CS=1 AND SN=0.
!     IF F=0 AND (G .NE. 0), THEN CS=0 AND SN=1 WITHOUT DOING ANY
!        FLOATING POINT OPERATIONS (SAVES WORK IN DBDSQR WHEN
!        THERE ARE ZEROS ON THE DIAGONAL).
!
!  IF F EXCEEDS G IN MAGNITUDE, CS WILL BE POSITIVE.
!
!  ARGUMENTS
!  =========
!
!  F       (INPUT) REAL*8
!          THE FIRST COMPONENT OF VECTOR TO BE ROTATED.
!
!  G       (INPUT) REAL*8
!          THE SECOND COMPONENT OF VECTOR TO BE ROTATED.
!
!  CS      (OUTPUT) REAL*8
!          THE COSINE OF THE ROTATION.
!
!  SN      (OUTPUT) REAL*8
!          THE SINE OF THE ROTATION.
!
!  R       (OUTPUT) REAL*8
!          THE NONZERO COMPONENT OF THE ROTATED VECTOR.
!
! ASTER INFORMATION
! 11/01/2000 TOILETTAGE DU FORTRAN SUIVANT LES REGLES ASTER,
!            REMPLACEMENT DE DLAMCH PAR R8PREM(), R8MIEM() ET
!            ISBAEM().
! 28/01/2000 RAJOUT DE LA VARIABLE BASE.
! INTRINSIC FUNCTIONS
!            ABS, INT, LOG, MAX, SQRT, DBLE.
!----------------------------------------------------------------------
! CORPS DU PROGRAMME
    implicit none
!
!     .. SCALAR ARGUMENTS ..
#include "asterf_types.h"
#include "asterc/isbaem.h"
#include "asterc/r8miem.h"
#include "asterc/r8prem.h"
    real(kind=8) :: cs, f, g, r, sn
!     ..
!     .. PARAMETERS ..
    real(kind=8) :: zero
    parameter(zero=0.0d0)
    real(kind=8) :: one
    parameter(one=1.0d0)
    real(kind=8) :: two
    parameter(two=2.0d0)
!     ..
!     .. LOCAL SCALARS ..
    aster_logical :: first
    integer(kind=8) :: count, i
    real(kind=8) :: eps, f1, g1, safmin, safmn2, safmx2, scale, base
!     ..
!     .. EXTERNAL FUNCTIONS ..
!     ..
!     .. SAVE STATEMENT ..
    save first, safmx2, safmin, safmn2
!     ..
!     .. DATA STATEMENTS ..
    data first/.true./
!     ..
!     .. EXECUTABLE STATEMENTS ..
!
    if (first) then
        first = .false.
        safmin = r8miem()
        eps = r8prem()*0.5d0
        base = dble(isbaem())
        safmn2 = base**int(log(safmin/eps)/log(base)/two)
        safmx2 = one/safmn2
    end if
    if (g .eq. zero) then
        cs = one
        sn = zero
        r = f
    else if (f .eq. zero) then
        cs = zero
        sn = one
        r = g
    else
        f1 = f
        g1 = g
        scale = max(abs(f1), abs(g1))
        if (scale .ge. safmx2) then
            count = 0
10          continue
            count = count+1
            f1 = f1*safmn2
            g1 = g1*safmn2
            scale = max(abs(f1), abs(g1))
            if (scale .ge. safmx2) goto 10
            r = sqrt(f1**2+g1**2)
            cs = f1/r
            sn = g1/r
            do i = 1, count
                r = r*safmx2
            end do
        else if (scale .le. safmn2) then
            count = 0
30          continue
            count = count+1
            f1 = f1*safmx2
            g1 = g1*safmx2
            scale = max(abs(f1), abs(g1))
            if (scale .le. safmn2) goto 30
            r = sqrt(f1**2+g1**2)
            cs = f1/r
            sn = g1/r
            do i = 1, count
                r = r*safmn2
            end do
        else
            r = sqrt(f1**2+g1**2)
            cs = f1/r
            sn = g1/r
        end if
        if (abs(f) .gt. abs(g) .and. cs .lt. zero) then
            cs = -cs
            sn = -sn
            r = -r
        end if
    end if
!
!     END OF DLARTG
!
end subroutine
