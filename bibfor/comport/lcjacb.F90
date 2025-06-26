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
! aslint: disable=W1504
subroutine lcjacb(fami, kpg, ksp, rela_comp, mod, &
                  nmat, materf, timed, timef, &
                  yf, deps, itmax, toler, nbcomm, &
                  cpmono, pgl, nfs, nsg, toutms, &
                  hsr, nr, nvi, vind, &
                  vinf, epsd, yd, dy, &
                  crit, &
                  drdy, iret)

    implicit none
!       CALCUL DU JACOBIEN DU SYSTEME NL A RESOUDRE = DRDY(DY)
!       IN  FAMI   :  FAMILLE DES POINTS DE GAUSS
!           KPG    :  NUMERO DU POINT DE GAUSS
!           KSP    :  NUMERO DU SOUS POINT DE GAUSS
!           LOI    :  MODELE DE COMPORTEMENT
!           MOD    :  TYPE DE MODELISATION
!           IMAT   :  ADRESSE DU MATERIAU CODE
!           NMAT   :  DIMENSION MATER
!           MATERF :  COEFFICIENTS MATERIAU A T+DT
!           YF     :  VARIABLES A T + DT =    ( SIGF  VINF  (EPS3F)  )
!           TOLER  :  TOLERANCE DE CONVERGENCE LOCALE
!           DEPS   :  INCREMENT DE DEFORMATION
!           ITMAX  :  NOMBRE MAXI D'ITERATIONS LOCALES
!           TIMED  :  INSTANT  T
!           TIMEF  :  INSTANT  T+DT
!           NBCOMM :  INCIDES DES COEF MATERIAU monocristal
!           CPMONO :  NOM DES COMPORTEMENTS monocristal
!           PGL    :  MATRICE DE PASSAGE
!           TOUTMS :  TENSEURS D'ORIENTATION monocristal
!           HSR    :  MATRICE D'INTERACTION monocristal
!           NR     :  DIMENSION DECLAREE DRDY
!           COMP   :  COMPOR - LOI ET TYPE DE DEFORMATION
!           NVI    :  NOMBRE DE VARIABLES INTERNES
!           VIND   :  VARIABLE INTERNES A T
!           VINF   :  VARIABLE INTERNES A T+DT
!           EPSD   :  DEFORMATION A T
!           YD     :  VARIABLES A T   = ( SIGD  VARD  ) A T
!           DY     :  SOLUTION           =    ( DSIG  DVIN  (DEPS3)  )
!           CRIT   :  CRITERES LOCAUX
!       OUT DRDY   :  JACOBIEN DU SYSTEME NON LINEAIRE
!           IRET   :  CODE RETOUR
!       ----------------------------------------------------------------
!
#include "asterf_types.h"
#include "asterfort/cvmjac.h"
#include "asterfort/hayjac.h"
#include "asterfort/irrjac.h"
#include "asterfort/lcmmja.h"
#include "asterfort/lkijac.h"
#include "asterfort/srijac.h"
    integer(kind=8) :: nr, nmat, kpg, ksp, itmax, iret, nvi, nfs, nsg
    real(kind=8) :: deps(*), epsd(*), toler, crit(*)
    real(kind=8) :: drdy(nr, nr), yf(nr), dy(nr), yd(nr)
!
    real(kind=8) :: materf(nmat, 2)
    real(kind=8) :: timed, timef, vind(*), vinf(*)
    real(kind=8) :: toutms(nfs, nsg, 6), hsr(nsg, nsg)
!
    character(len=*) :: fami
    character(len=8) :: mod
    character(len=16) :: rela_comp
!
    integer(kind=8) :: nbcomm(nmat, 3)
    real(kind=8) :: pgl(3, 3)
    character(len=24) :: cpmono(5*nmat+1)
!
!       ----------------------------------------------------------------
!
    iret = 0
    if (rela_comp .eq. 'VISCOCHAB') then
        call cvmjac(mod, nmat, materf, timed, timef, &
                    yf, dy, nr, epsd, deps, &
                    drdy)
!
    else if (rela_comp .eq. 'MONOCRISTAL') then
        call lcmmja(mod, nmat, materf, timed, &
                    timef, itmax, toler, nbcomm, cpmono, &
                    pgl, nfs, nsg, toutms, hsr, &
                    nr, nvi, vind, deps, yf, &
                    yd, dy, drdy, iret)
!
    else if (rela_comp .eq. 'IRRAD3M') then
        call irrjac(fami, kpg, ksp, mod, nmat, &
                    materf, yf, dy, nr, drdy)
!
    else if (rela_comp .eq. 'LETK') then
        call lkijac(mod, nmat, materf, timed, timef, &
                    yf, deps, nr, nvi, vind, &
                    vinf, yd, dy, drdy, iret)
!
    else if (rela_comp .eq. 'LKR') then
        call srijac(nmat, materf, timed, timef, &
                    yf, deps, nr, nvi, vind, &
                    vinf, yd, drdy)
!
    else if (rela_comp .eq. 'HAYHURST') then
        call hayjac(mod, nmat, materf(1, 1), materf(1, 2), timed, &
                    timef, yf, deps, nr, nvi, &
                    vind, vinf, yd, dy, crit, &
                    drdy, iret)
!
    end if
!
end subroutine
