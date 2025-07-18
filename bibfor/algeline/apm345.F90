! --------------------------------------------------------------------
! Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
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
subroutine apm345(nbtetc, typcon, rayonc, centrc, nk, &
                  k24rc, pivot2, ltest, typcha, lraide, &
                  lmasse, ldynam, solveu, lamor, lc, &
                  impr, ifapm)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8pi.h"
#include "asterc/r8rddg.h"
#include "asterfort/apchar.h"
#include "asterfort/assert.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: nbtetc, nk, pivot2, lraide, lmasse, ldynam, lamor, ifapm
    real(kind=8) :: rayonc
    complex(kind=8) :: centrc
    aster_logical :: ltest, lc
    character(len=3) :: impr
    character(len=8) :: typcon, typcha
    character(len=19) :: solveu
    character(len=24) :: k24rc
!     STEPS 3/4/5 OF THE ARGUMENT PRINCIPAL METHOD THAT COUNT THE
!     EIGENVALUES WITHIN A GIVEN SHAPE OF THE COMPLEX PLANE
!     ------------------------------------------------------------------
! IN NBTETC : IN : NUMBER OF CHECKED POINTS THAT DISCRETIZE THE SHAPE
! IN TYPCON : K8 : TYPE OF SHAPE
! IN RAYONC : R8 : RADIUS OF THE DISC OF THE GIVEN SHAPE
! IN CENTRC : C16: CENTRE OF THE DISC
! IN NK     : IN : SIZE OF THE EIGENVALUE PROBLEM
! IN K24RC  : K24: JEVEUX NAME OF THE COEFFICIENTS OF THE CHARACTERIS
!                  TIC POLYNOMIAL (ONLY NEEDED IF TYPCHA='ROMBOUT')
! OUT PIVOT2: IN : NUMBER OF EIGENVALUES IN THE SHAPE
! IN LTEST  : LOG: RUN OF INTERNAL TEST FOR DEBUGGING IF LTEST=.TRUE.
! IN TYPCHA : K8 : TYPE OF METHOD TO EVALUATE THE CHARAC. POLYNOME
! IN LRAIDE : IN : JEVEUX ATTRIBUT OF THE STIFFNESS MATRIX
! IN LMASSE : IN : JEVEUX ATTRIBUT OF THE MASS MATRIX
! IN LDYNAM : IN : JEVEUX ATTRIBUT OF THE DYNAMIC MATRIX
! IN SOLVEU : K19: JEVEUX SD OF THE LINEAR SOLVER
! IN LAMOR  : IN : JEVEUX ATTRIBUT OF THE DAMPING MATRIX
! IN LC     : LOG: FLAG THAT INDICATES IF THE PB IS QUADRATIC OR NOT
! IN IMPR/IFAPM: IN/K3 : PRINT PARAMETERS FOR DEBBUGING
!     ------------------------------------------------------------------
! person_in_charge: olivier.boiteau at edf.fr
!
!
    integer(kind=8) :: i, jcont, jtheta, ifm, niv
    real(kind=8) :: raux1, pi, raux2, rauxx, rauxy, theta, thetao, thetap
    real(kind=8) :: thetam, raddeg, prec2, piprec
!
!   --- MISCELLANEOUS ---
    call jemarq()
    call infniv(ifm, niv)
    pi = r8pi()
    raddeg = r8rddg()
!
!   --- VERBOSE MODE FOR APM STEPS ---
!      NIV=2
!    --- PERTURBATION TO BE ADDED TO THE EVALUATION OF THE ARGUMENT ---
!    --- FOR NOT MISSING A LOOP WHEN ARG IS JUST UNDER 2*PI RADIENT ---
    prec2 = 1.d-3
    piprec = 2.d0*pi-prec2
!
!   --- STEP 3: DISCRETISATION OF THE SHAPE ---
    call wkvect('&&APM345.CONTOUR.DIS', 'V V C', nbtetc, jcont)
    if (typcon(1:6) .eq. 'CERCLE') then
        ASSERT(nbtetc .gt. 1)
        raux1 = 2.d0*pi/(nbtetc-1)
        do i = 1, nbtetc
            raux2 = (i-1)*raux1
            rauxx = rayonc*cos(raux2)
            rauxy = rayonc*sin(raux2)
            zc(jcont+i-1) = centrc+dcmplx(rauxx, rauxy)
        end do
    end if
!
!   --- STEP 4: EVALUATING THE ARGUMENT VALUE OF P(CONTOUR) ---
    call wkvect('&&APM345.CONTOUR.THETA', 'V V R', nbtetc, jtheta)
    do i = 1, nbtetc
        call apchar(typcha, k24rc, nk, zc(jcont+i-1), theta, &
                    lraide, lmasse, ldynam, solveu, lamor, &
                    lc, impr, ifapm, i)
        zr(jtheta+i-1) = theta
    end do
    call jedetr('&&APM345.CONTOUR.DIS')
!
!   --- STEP 5: COUNTING THE DIFFERENCE BETWEEN ANGLE         ---
    if (ltest .or. (niv .ge. 2)) then
        write (ifm, *) 'DISPLAY I/ARG(P(LAMBDA_I))/NB_EIGENVALUE_I'
        write (ifm, *) '------------------------------------------'
    end if
    pivot2 = 0
    thetao = zr(jtheta)
    if (ltest .or. (niv .ge. 2)) write (ifm, *) 'STEP 5: ', 1, thetao*raddeg, pivot2
    do i = 2, nbtetc
        theta = zr(jtheta+i-1)
!   --- CORRECTION IF THETA GOES NEAR 2*PI TO NOT MISS A LOOP ---
        if (theta .gt. piprec) theta = prec2
!
!   --- HEURISTIC BASED ON THE WORK OF O.BERTRAND (PHD INRIA) ---
!   --- WE ADD ONLY THE CONDITION I>2                         ---
        thetap = thetao+pi
        thetam = thetao-pi
        if ((thetao .lt. pi) .and. ((theta .gt. thetap) .or. abs(theta) .lt. prec2) .and. &
            (i .gt. 2)) then
            pivot2 = pivot2-1
        else if ((thetao .gt. pi) .and. ((theta .lt. thetam) .or. abs(theta) &
                                         .lt. prec2)) then
            pivot2 = pivot2+1
        end if
        thetao = theta
!
        if (ltest .or. (niv .ge. 2)) write (ifm, *) 'STEP 5: ', i, theta*raddeg, pivot2
    end do
    call jedetr('&&APM345.CONTOUR.THETA')
!
    call jedema()
!
end subroutine
