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
subroutine lcmate(fami, kpg, ksp, comp, mod, &
                  imat, nmat, tempd, tempf, tref, impexp, &
                  typma, hsr, materd, materf, matcst, &
                  nbcomm, cpmono, angmas, pgl, itmax, &
                  toler, ndt, ndi, nr, crit, &
                  nvi, vind, nfs, nsg, toutms, &
                  nhsr, numhsr, sigd, mult_comp_)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/cvmmat.h"
#include "asterfort/haymat.h"
#include "asterfort/hbrmat.h"
#include "asterfort/irrmat.h"
#include "asterfort/lcmatt.h"
#include "asterfort/lcmmap.h"
#include "asterfort/lcmmat.h"
#include "asterfort/lglmat.h"
#include "asterfort/lkimat.h"
#include "asterfort/srimat.h"
#include "asterfort/matect.h"
#include "asterfort/rslmat.h"
#include "asterfort/rsvmat.h"
#include "asterfort/vecmat.h"
!
!       RECUPERATION DU MATERIAU A TEMPF ET TEMPD
!       IN  FAMI   :  FAMILLE DE POINT DE GAUSS (RIGI,MASS,...)
!           KPG,KSP:  NUMERO DU (SOUS)POINT DE GAUSS
!           COMP   :  COMPORTEMENT
!           MOD    :  TYPE DE MODELISATION
!           IMAT   :  ADRESSE DU MATERIAU CODE
!           NMAT   :  DIMENSION 1 DE MATER
!           TEMPD  :  TEMPERATURE A T
!           TEMPF  :  TEMPERATURE A T + DT
!           IMPEXP : 0 IMPLICITE, 1 EXPLICITE
!          ANGMAS  :  LES TROIS ANGLES DU MOT_CLEF MASSIF
!           SIGD   :  ETAT DE CONTRAINTES A T
!       OUT MATERD :  COEFFICIENTS MATERIAU A T    (TEMPD )
!           MATERF :  COEFFICIENTS MATERIAU A T+DT (TEMPF )
!                     MATER(*,I) = CARACTERISTIQUES MATERIAU
!                                    I = 1  CARACTERISTIQUES ELASTIQUES
!                                    I = 2  CARACTERISTIQUES PLASTIQUES
!           MATCST :  'OUI' SI  MATERIAU A T = MATERIAU A T+DT
!                     'NON' SINON OU 'NAP' SI NAPPE DANS 'VECMAT.F'
!           NBCOMM : POSITION DES COEF POUR CHAQUE LOI DE CHAQUE SYSTEME
!           CPMONO : NOMS DES LOIS POUR CHAQUE FAMILLE DE SYSTEME
!           PGL    : MATRICE DE PASSAGE
!           NDT    :  NB TOTAL DE COMPOSANTES TENSEURS
!           NDI    :  NB DE COMPOSANTES DIRECTES  TENSEURS
!           NR     :  NB DE COMPOSANTES SYSTEME NL
!           NVI    :  NB DE VARIABLES INTERNES
!           TOUTMS :  TOUS LES TENSEURS MS
!           HSR    : MATRICE D'INTERACTION POUR L'ECROUISSAGE ISOTROPE
!                    UTILISEE SEULEMENT POUR LE MONOCRISTAL IMPLICITE
!       ----------------------------------------------------------------

    integer(kind=8), intent(in):: nvi
    integer(kind=8) :: imat, nmat, ndt, ndi, nr, i, itmax, kpg, ksp, impexp
    real(kind=8) :: materd(nmat, 2), materf(nmat, 2), tempd, tempf, tref
    real(kind=8) :: vind(*), pgl(3, 3), angmas(3), toler, crit(*), sigd(6)
    character(len=16) :: rela_comp, comp(*), mult_comp
    character(len=8) :: mod, typma
    character(len=3) :: matcst
    character(len=*) :: fami
!     SPECIFIQUE MONOCRISTAL
    integer(kind=8) :: numhsr(*), nbcomm(*), nfs, nsg, nhsr
    real(kind=8) :: hsr(*), toutms(*)
    character(len=24) :: cpmono(*)
    character(len=16), optional, intent(in) :: mult_comp_
!       ----------------------------------------------------------------
!
! -     INITIALISATION DE MATERD ET MATERF A 0.
!
    do i = 1, nmat
        materd(i, 1) = 0.d0
        materd(i, 2) = 0.d0
        materf(i, 1) = 0.d0
        materf(i, 2) = 0.d0
    end do
! - For number of phases when is not a crystal behaviour (issue30310)
    nbcomm(1) = 1
!
    mult_comp = ' '
    if (present(mult_comp_)) then
        mult_comp = mult_comp_
    end if
    rela_comp = comp(1)
    if (rela_comp .eq. 'ROUSS_PR') then
        call rslmat(fami, kpg, ksp, mod, imat, &
                    nmat, materd, materf, matcst, ndt, &
                    ndi, nr, nvi, vind)
!
    else if (rela_comp .eq. 'ROUSS_VISC') then
        call rsvmat(fami, kpg, ksp, mod, imat, &
                    nmat, materd, materf, matcst, ndt, &
                    ndi, nr, nvi, vind)
!
    else if (rela_comp .eq. 'VISCOCHAB') then
        call cvmmat(fami, kpg, ksp, mod, imat, &
                    nmat, materd, materf, matcst, typma, &
                    ndt, ndi, nr, crit, vind, &
                    nvi, sigd)
!
    else if (rela_comp .eq. 'VENDOCHAB' .or. rela_comp .eq. 'VISC_ENDO_LEMA') then
        call vecmat(fami, kpg, ksp, mod, rela_comp, &
                    imat, nmat, materd, materf, matcst, &
                    typma, ndt, ndi, nr, nvi)
!
    else if (rela_comp(1:6) .eq. 'LAIGLE') then
        call lglmat(mod, imat, nmat, tempd, materd, &
                    materf, matcst, ndt, ndi, nr, &
                    nvi)
!
    elseif ((rela_comp .eq. 'HOEK_BROWN') .or. (rela_comp .eq. 'HOEK_BROWN_EFF')) then
        call hbrmat(mod, imat, nmat, tempd, materd, &
                    materf, matcst, ndt, ndi, nr, &
                    nvi)
!
    else if (rela_comp .eq. 'MONOCRISTAL') then
        ASSERT(mult_comp .ne. ' ')
        call lcmmat(fami, kpg, ksp, mult_comp, mod, &
                    imat, nmat, angmas, pgl, materd, &
                    materf, matcst, nbcomm, cpmono, ndt, &
                    ndi, nr, nvi, hsr, nfs, &
                    nsg, toutms, vind, impexp)
        typma = 'COHERENT'
        if (mod .ne. '3D') then
            sigd(5) = 0.d0
            sigd(6) = 0.d0
        end if
!
    else if (rela_comp .eq. 'POLYCRISTAL') then
        ASSERT(mult_comp .ne. ' ')
        call lcmmap(fami, kpg, ksp, mult_comp, mod, &
                    imat, nmat, angmas, pgl, materd, &
                    materf, matcst, nbcomm, cpmono, ndt, &
                    ndi, nr, nvi, nfs, nsg, &
                    nhsr, numhsr, hsr)
        typma = 'COHERENT'
!
    else if (rela_comp .eq. 'IRRAD3M') then
        call irrmat(fami, kpg, ksp, mod, imat, &
                    nmat, itmax, toler, materd, materf, &
                    matcst, ndt, ndi, nr, nvi)
!
    else if (rela_comp .eq. 'LETK') then
        call lkimat(mod, imat, nmat, materd, materf, &
                    matcst, ndt, ndi, nvi, nr)
        typma = 'COHERENT'
!
    else if (rela_comp .eq. 'LKR') then
        call srimat(mod, imat, nmat, tempd, tempf, tref, materd, materf, &
                    matcst, ndt, ndi, nvi, nr)
        typma = 'COHERENT'
!
    else if (rela_comp .eq. 'HAYHURST') then
        call haymat(fami, kpg, ksp, mod, imat, &
                    nmat, '-', materd(1, 1), materd(1, 2), nvi, &
                    nr)
        call haymat(fami, kpg, ksp, mod, imat, &
                    nmat, '+', materf(1, 1), materf(1, 2), nvi, &
                    nr)
        call matect(materd, materf, nmat, matcst)
        typma = 'COHERENT'
!
    else
!
! CAS GENERAL
!
        call lcmatt(fami, kpg, ksp, mod, imat, &
                    nmat, '-', rela_comp, materd(1, 1), materd(1, 2), &
                    typma, ndt, ndi, nr, nvi)
        call lcmatt(fami, kpg, ksp, mod, imat, &
                    nmat, '+', rela_comp, materf(1, 1), materf(1, 2), &
                    typma, ndt, ndi, nr, nvi)
!
        call matect(materd, materf, nmat, matcst)
!
    end if
!
!     - DANS LCPLNL ON DIMENSIONNE DES TABLES AVEC (NDT+NVI) QUI SONT
!       ENSUITE UTILISEES PAR NEWTON
!     - LA DIMENSION DU SYSTEME DIFFERENTIEL EST NR
!     ==> IL FAUT DONC NDT+NVI >= NR
    ASSERT((ndt+nvi) .ge. nr)
end subroutine
