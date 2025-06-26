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

subroutine lccnvx(fami, kpg, ksp, rela_comp, &
                  imat, nmat, materf, sigd, &
                  sigf, deps, vind, vinf, nbcomm, &
                  cpmono, pgl, nvi, vp, vecp, &
                  hsr, nfs, nsg, toutms, timed, &
                  timef, &
                  seuil)
! aslint: disable=W1504
    implicit none
! --- BUT : CONVEXE ELASTO PLASTIQUE A T+DT POUR (SIGF , VIND) DONNES --
! ----------------------------------------------------------------------
! IN  : FAMI   :  FAMILLE DES POINTS DE GAUSS  -------------------------
! --- : KPG    :  NUMERO DU POINT DE GAUSS  ----------------------------
! --- : KSP    :  NUMERO DU SOUS POINT DE GAUSS ------------------------
! --- : LOI    :  NOM DU MODELE DE COMPORTEMENT ------------------------
! --- : SIGD   :  CONTRAINTE A T ---------------------------------------
! --- : SIGF   :  CONTRAINTE A T+DT ------------------------------------
! --- : DEPS   :  INCRMEENT DE DEFORMATION -----------------------------
! --- : VIND   :  VARIABLES INTERNES A T -------------------------------
! --- : VINF   :  VARIABLES INTERNES A T +DT----------------------------
! --- : IMAT   :  ADRESSE DU MATERIAU CODE -----------------------------
! --- : NMAT   :  DIMENSION MATER --------------------------------------
! --- : TOLER  :  TOLERANCE DE CONVERGENCE LOCALE-----------------------
! --- : MATERD :  COEFFICIENTS MATERIAU A T ----------------------------
! --- : MATERF :  COEFFICIENTS MATERIAU A T+DT -------------------------
! --- : NBCOMM :  INDICES DES COEF MATERIAU ----------------------------
! --- : TIMED  :  INSTANT T --------------------------------------------
! --- : TIMEF  :  INSTANT T+DT -----------------------------------------
! OUT : VP     :  VALEURS PROPRES DU DEVIATEUR ELASTIQUE (HOEK-BROWN) --
! --- : VECP   :  VECTEURS PROPRES DU DEVIATEUR ELASTIQUE (HOEK-BROWN) -
! --- : SEUIL  :  SEUIL  ELASTICITE  A T+DT ----------------------------
! --- : YD     :  VECTEUR INCONNUES A T --------------------------------
! --- : YF     :  VECTEUR INCONNUES A T+DT -----------------------------
! ----------------------------------------------------------------------
! ======================================================================
#include "asterfort/cvmcvx.h"
#include "asterfort/hbrcvx.h"
#include "asterfort/irrcvx.h"
#include "asterfort/lcmmvx.h"
#include "asterfort/lglcvx.h"
#include "asterfort/lkcnvx.h"
#include "asterfort/srcnvx.h"
#include "asterfort/rslcvx.h"
    integer(kind=8) :: nmat, imat, nvi, kpg, ksp, nfs, nsg
    character(len=*) :: fami
    real(kind=8) :: materf(nmat, 2), seuil
    real(kind=8) :: timed, timef, deps(6), vinf(*)
    real(kind=8) :: sigd(6), sigf(6), vind(*), hsr(nsg, nsg)
    character(len=16), intent(in) :: rela_comp
    integer(kind=8) :: nbcomm(nmat, 3)
    real(kind=8) :: pgl(3, 3), vp(3), vecp(3, 3), toutms(nfs, nsg, 6)
    character(len=24) :: cpmono(5*nmat+1)
! ======================================================================
    if (rela_comp .eq. 'ROUSS_PR') then
        call rslcvx(fami, kpg, ksp, imat, nmat, &
                    materf, sigf, vind, seuil)
! ======================================================================
    else if (rela_comp .eq. 'ROUSS_VISC') then
        call rslcvx(fami, kpg, ksp, imat, nmat, &
                    materf, sigf, vind, seuil)
! ======================================================================
    else if (rela_comp .eq. 'VISCOCHAB') then
        call cvmcvx(nmat, materf, sigf, vind, seuil)
! ======================================================================
    else if (rela_comp .eq. 'LAIGLE') then
        call lglcvx(sigf, vind, nmat, materf, seuil)
! ======================================================================
    elseif ((rela_comp .eq. 'HOEK_BROWN') .or. (rela_comp .eq. 'HOEK_BROWN_EFF')) then
        call hbrcvx(sigf, vind, nmat, materf, seuil, &
                    vp, vecp)
! ======================================================================
    else if (rela_comp .eq. 'MONOCRISTAL') then
        call lcmmvx(sigf, vind, nmat, materf, nbcomm, &
                    cpmono, pgl, nvi, hsr, nfs, &
                    nsg, toutms, timed, timef, deps, &
                    seuil)
! ======================================================================
    else if (rela_comp .eq. 'IRRAD3M') then
        call irrcvx(fami, kpg, ksp, nmat, materf, &
                    sigf, vind, seuil)
! ======================================================================
    else if (rela_comp .eq. 'LETK') then
        call lkcnvx(sigd, sigf, nvi, vind, nmat, &
                    materf, seuil, vinf)
! ======================================================================
    else if (rela_comp .eq. 'LKR') then
        call srcnvx(sigd, sigf, nvi, vind, nmat, materf, seuil, vinf)
! ======================================================================
    end if
! ======================================================================
end subroutine
