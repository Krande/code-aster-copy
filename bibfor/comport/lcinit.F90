! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
subroutine lcinit(fami, kpg, ksp, rela_comp, typess, &
                  essai, mod, nmat, materf, &
                  timed, timef, intg, nr, nvi, &
                  yd, epsd, deps, dy, compor, &
                  nbcomm, cpmono, pgl, nfs, nsg, &
                  toutms, vind, sigd, sigf, epstr, &
                  bnews, mtrac, indi, iret)

    implicit none
!       ROUTINE AIGUILLAGE
!       ----------------------------------------------------------------
!       CALCUL DE LA SOLUTION INITIALE ESSAI DY = ( DSIG DVIN (DEPS3) )
!       IN  FAMI   :  FAMILLE DE POINT DE GAUSS
!           KPG    :  NUMERO DU POINT DE GAUSS
!           KSP    :  NUMERO DU SOUS-POINT DE GAUSS
!           LOI    :  MODELE DE COMPORTEMENT
!           TYPESS :  TYPE DE SOLUTION D ESSAI POUR DY(DEPEND DU MODELE)
!                      > VOIR XXXCVG ET XXXINI
!           ESSAI  :  SOLUTION D ESSAI
!           MOD    :  TYPE DE MODELISATION
!           NMAT   :  DIMENSION MATER
!           MATERD :  COEFFICIENTS MATERIAU A T
!           MATERF :  COEFFICIENTS MATERIAU A T+DT
!           TIMED  :  INSTANT  T
!           TIMEF  :  INSTANT  T+DT
!           INTG   :  NOMBRE DE TENTATIVES D'INTEGRATION
!           EPSD   :  DEFORMATION A T
!           YD     :  VARIABLES A T   = ( SIG  VIN  (EPS3)  )
!           NVI    :  NOMBRE VARIABLES INTERNES
!           VIND   :  VECTEUR VARIABLES INTERNES A T
!           NR     :  DIMENSION VECTEUR INCONNUES
!           SIGF   :  PREDICTION ELASTIQUE DES CONTRAINTES (LCELAS)
!           BNEWS  :  GESTION TRACTION AVEC HUJEUX (IN/OUT)
!           MTRAC  :  GESTION TRACTION AVEC HUJEUX (BIS) (IN/OUT)
!           IRET   :  IRET = 2 - RELANCE DU PROCESSUS DE RESOLUTION
!       VAR DEPS   :  INCREMENT DE DEFORMATION
!       OUT DY     :  SOLUTION ESSAI  = ( DSIG DVIN (DEPS3) )
!           INDI   :  INDICATEURS DES MECANISMES POT. ACTIFS (HUJEUX)
!       ----------------------------------------------------------------
#include "asterf_types.h"
#include "asterfort/cvmini.h"
#include "asterfort/hujini.h"
#include "asterfort/irrini.h"
#include "asterfort/lcmmin.h"
#include "asterfort/lklini.h"
#include "asterfort/srlini.h"
#include "asterfort/Behaviour_type.h"
    integer :: typess, nmat, nr, nvi, kpg, ksp, nfs, nsg
    integer :: nbcomm(nmat, 3), iret, indi(7), intg
    real(kind=8) :: deps(6), epsd(6), essai
    real(kind=8) :: yd(*), dy(*)
    real(kind=8) :: materf(nmat, 2)
    real(kind=8) :: timed, timef
    real(kind=8) :: pgl(3, 3)
    real(kind=8) :: vind(*), sigd(6), epstr(6)
    real(kind=8) :: toutms(nfs, nsg, 6), sigf(6)
    character(len=*) :: fami
    character(len=8) :: mod
    character(len=16), intent(in) :: rela_comp
    character(len=16), intent(in) :: compor(COMPOR_SIZE)
    character(len=24) :: cpmono(5*nmat+1)
    aster_logical :: bnews(3), mtrac
!       ----------------------------------------------------------------
!
    iret = 0
!
    if (rela_comp .eq. 'VISCOCHAB') then
        call cvmini(typess, essai, mod, nmat, materf, &
                    timed, timef, yd, epsd, deps, &
                    dy)
!
    else if (rela_comp .eq. 'MONOCRISTAL') then
        call lcmmin(typess, essai, mod, nmat, materf, &
                    nr, nvi, yd, deps, dy, &
                    compor, nbcomm, cpmono, pgl, nfs, &
                    nsg, toutms, timed, timef, vind, &
                    sigd, epstr)
!
    else if (rela_comp .eq. 'IRRAD3M') then
        call irrini(fami, kpg, ksp, typess, essai, &
                    mod, nmat, materf, yd, deps, &
                    dy)
!
    else if (rela_comp .eq. 'LETK') then
        call lklini(sigf, nr, yd, dy)
!
    else if (rela_comp .eq. 'LKR') then
        call srlini(sigf, nr, yd, dy)
!
    else if (rela_comp .eq. 'HUJEUX') then
        call hujini(mod, nmat, materf, intg, deps, &
                    nr, yd, nvi, vind, sigd, &
                    sigf, bnews, mtrac, dy, indi, &
                    iret)
!
    else
!        SOLUTION INITIALE = ZERO
        dy(1:nr) = 0.d0
        if (mod(1:6) .eq. 'C_PLAN') then
            deps(3) = 0.d0
        end if
    end if
!
end subroutine
