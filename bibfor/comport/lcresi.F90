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

subroutine lcresi(fami, kpg, ksp, rela_comp, typmod, &
                  imat, nmat, materd, materf, &
                  nbcomm, cpmono, pgl, nfs, nsg, &
                  toutms, hsr, nr, nvi, vind, &
                  vinf, itmax, toler, timed, timef, &
                  yd, yf, deps, epsd, dy, &
                  r, iret, crit)
! aslint: disable=W1504
    implicit none
!       CALCUL DES TERMES DU SYSTEME NL A RESOUDRE = R(DY)
!       IN  FAMI   :  FAMILLE DU POINT DE GAUSS
!           KPG    :  POINT DE GAUSS
!           KSP    :  SOUS-POINT DE GAUSS
!           LOI    :  MODELE DE COMPORTEMENT
!           TYPMOD    :  TYPE DE MODELISATION
!           IMAT   :  NOM DU MATERIAU
!           NMAT   :  DIMENSION MATER
!           MATERD :  COEFFICIENTS MATERIAU A T
!           MATERF :  COEFFICIENTS MATERIAU A T+DT
!           TIMED  :  INSTANT  T
!           TIMEF  :  INSTANT  T+DT
!           DEPS   :  INCREMENT DE DEFORMATION
!           EPSD   :  DEFORMATION A T
!           VIND   :  VARIABLES INTERNES A T
!           VINF   :  VARIABLES INTERNES A T+DT
!           YD     :  VARIABLES A T      =    ( SIGD  VIND  (EPSD3)  )
!           YF     :  VARIABLES A T + DT =    ( SIGF  VINF  (EPS3F)  )
!           DY     :  SOLUTION           =    ( DSIG  DVIN  (DEPS3)  )
!       OUT R      :  SYSTEME NL A T + DT
!       ----------------------------------------------------------------
!
#include "asterfort/cvmres.h"
#include "asterfort/hayres.h"
#include "asterfort/irrres.h"
#include "asterfort/lcmmre.h"
#include "asterfort/lcresa.h"
#include "asterfort/lkresi.h"
#include "asterfort/srresi.h"
    integer(kind=8) :: imat, nmat, nr, nvi, kpg, ksp, itmax, iret
    integer(kind=8) :: nfs, nsg
    real(kind=8) :: deps(6), epsd(6), vind(*), toler
    real(kind=8) :: r(*), yd(*), yf(*), dy(*), vinf(*)
    real(kind=8) :: materd(nmat, 2), materf(nmat, 2)
    real(kind=8) :: timed, timef, crit(*)
    character(len=8) :: typmod
    character(len=16) :: rela_comp
    real(kind=8) :: toutms(nfs, nsg, 6), hsr(nsg, nsg)
    character(len=*) :: fami
!
    integer(kind=8) :: nbcomm(nmat, 3)
    real(kind=8) :: pgl(3, 3)
    character(len=24) :: cpmono(5*nmat+1)
!
!       ----------------------------------------------------------------
!
    iret = 0
    if (rela_comp .eq. 'VISCOCHAB') then
        call cvmres(typmod, nmat, materd, materf, timed, &
                    timef, yd, yf, epsd, deps, &
                    dy, r)
!
    else if (rela_comp .eq. 'MONOCRISTAL') then
        call lcmmre(typmod, nmat, materd, materf, &
                    nbcomm, cpmono, pgl, nfs, nsg, &
                    toutms, hsr, nr, nvi, vind, &
                    itmax, toler, timed, timef, yd, &
                    yf, deps, dy, r, iret)
    else if (rela_comp .eq. 'IRRAD3M') then
        call irrres(fami, kpg, ksp, typmod, nmat, &
                    materd, materf, yd, yf, deps, &
                    dy, r)
    else if (rela_comp .eq. 'LETK') then
        call lkresi(typmod, nmat, materf, timed, timef, &
                    nvi, vind, vinf, yd, yf, &
                    deps, nr, r)
    else if (rela_comp .eq. 'LKR') then
        call srresi(nmat, materf, timed, timef, &
                    nvi, vind, vinf, yd, yf, deps, nr, r)
    else if (rela_comp .eq. 'HAYHURST') then
        call hayres(typmod, nmat, materd, materf, timed, &
                    timef, yd, yf, deps, dy, &
                    r, crit, iret)
    else
        call lcresa(fami, kpg, ksp, typmod, imat, &
                    nmat, materd, materf, rela_comp, nr, &
                    nvi, timed, timef, deps, epsd, &
                    yf, dy, r, iret, yd, &
                    crit)
!
    end if
!
end subroutine
