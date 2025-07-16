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

subroutine nmorth(fami, kpg, ksp, ndim, phenom, &
                  imate, poum, deps, sigm, option, &
                  angl_naut, sigp, dsidep)
    implicit none
#include "asterf_types.h"
#include "asterc/r8vide.h"
#include "asterfort/d1ma3d.h"
#include "asterfort/d1mamc.h"
#include "asterfort/dmat3d.h"
#include "asterfort/dmatmc.h"
#include "asterfort/lteatt.h"
#include "asterfort/matrot.h"
#include "asterfort/utmess.h"
#include "asterfort/verift.h"
#include "asterfort/verifh.h"
#include "asterfort/verifs.h"
#include "asterfort/verifepsa.h"
#include "asterfort/utbtab.h"
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in)          :: kpg
    integer(kind=8), intent(in)          :: ksp
    integer(kind=8), intent(in)          :: ndim
    character(len=16), intent(in):: phenom
    integer(kind=8), intent(in)          :: imate
    character(len=*), intent(in) :: poum
    real(kind=8), intent(in)      :: deps(2*ndim)
    real(kind=8), intent(in)      :: sigm(2*ndim)
    character(len=16), intent(in):: option
    real(kind=8), intent(in)      :: angl_naut(3)
    real(kind=8), intent(out)     :: sigp(2*ndim)
    real(kind=8), intent(out)     :: dsidep(2*ndim, 2*ndim)
!
! --------------------------------------------------------------------
!  IN    FAMI   : FAMILLE DE POINT DE GAUSS
!  IN    KPG    : NUMERO DU POINT DE GAUSS
!  IN    KSP    : NUMERO DU SOUS POINT DE GAUSS
!  IN    NDIM   : DIMENSION DU PROBLEME
!  IN    PHENOM : PHENOMENE (ELAS_ORTH OU ELAS_ISTR)
!  IN    TYPMOD : TYPE DE MODELISATION
!  IN    IMATE  : ADRESSE DU MATERIAU
!  IN    POUM   : '+' INSTANT SUIVANT OU '-' INSTANT COURANT
!                 OU 'T' = '+' - '-' INCREMENT
!  IN    EPSM   : DEFORMATION A L INSTANT T-
!  IN    DESPS  : INCREMENT DE DEFORMATION
!  IN    SIGM   : CONTRAINTE A L INSTANT T-
!  IN    OPTION : OPTION A CALCULER
!  IN    ANGL_NAUT : ANGLE DU REPERE LOCAL D ORTHOTROPIE
!  OUT   SIGP   : CONTRAINTE A L INSTANT T+
!                 CAR IL EN EXISTE FORCEMENT UNE)
!  OUT   DSIDEP : MATRICE DE RIGIDITE TANGENTE
!
! --------------------------------------------------------------------
    real(kind=8) :: p(3, 3)
    real(kind=8) :: rbid, hookf(36), mkooh(36)
    real(kind=8) :: depstr(6)
    real(kind=8) :: epsth_anis(3), deplth(6), depgth(6)
    real(kind=8) :: depghy, depgse, depgepsa(6)
    real(kind=8) :: depsme(6), rac2, epsm2(6)
    real(kind=8) :: work(3, 3), deplth_mat(3, 3), depgth_mat(3, 3)
    real(kind=8) :: sigm2(2*ndim)
    integer(kind=8) :: ndimsi, i, j
    aster_logical :: vrai
! --------------------------------------------------------------------
!
    rbid = r8vide()
!
    if (phenom .eq. 'ELAS_ISTR' .and. ndim .eq. 2) then
        call utmess('F', 'ELEMENTS3_2')
    end if

    rac2 = sqrt(2.d0)
    ndimsi = ndim*2
    dsidep = 0
    depgth = 0
!
    if (option .eq. 'FULL_MECA' .or. option .eq. 'RAPH_MECA') then
        do i = 1, ndimsi
            if (i .le. 3) then
                depstr(i) = deps(i)
            else
                depstr(i) = deps(i)*rac2
            end if
        end do
    end if
!
    if (angl_naut(1) .eq. r8vide()) then
        call utmess('F', 'ALGORITH8_20')
    end if
!
! - VERIFICATION DE L'ELEMENT
!
    vrai = .false.
    if (fami .eq. 'PMAT') then
!        ON VIENT DE OP0033
        vrai = .true.
    else
        if (lteatt('DIM_TOPO_MAILLE', '3')) then
            vrai = .true.
        else if (lteatt('C_PLAN', 'OUI')) then
            vrai = .true.
        else if (lteatt('D_PLAN', 'OUI')) then
            vrai = .true.
        else if (lteatt('AXIS', 'OUI')) then
            vrai = .true.
        end if
    end if
!
    if (.not. vrai) then
        call utmess('F', 'ALGORITH8_22')
    end if
!
! - MATRICES TANGENTES
!
    if (fami .eq. 'PMAT') then
!        ON VIENT DE OP0033
        if (option .eq. 'RIGI_MECA_TANG') then
            call dmat3d(fami, imate, rbid, '-', kpg, &
                        ksp, angl_naut, hookf)
        else
            call d1ma3d(fami, imate, rbid, '-', kpg, &
                        ksp, angl_naut, mkooh)
            call dmat3d(fami, imate, rbid, '+', kpg, &
                        ksp, angl_naut, hookf)
        end if
!
    else
        if (option .eq. 'RIGI_MECA_TANG') then
            call dmatmc(fami, imate, rbid, '-', kpg, &
                        ksp, angl_naut, ndimsi, hookf)
        else
            call d1mamc(fami, imate, rbid, '-', kpg, &
                        ksp, angl_naut, ndimsi, mkooh)
            call dmatmc(fami, imate, rbid, '+', kpg, &
                        ksp, angl_naut, ndimsi, hookf)
        end if
    end if
!
    if (option(1:10) .eq. 'RIGI_MECA_' .or. option(1:9) .eq. 'FULL_MECA') then
        do i = 1, ndimsi
            do j = 1, ndimsi
                dsidep(i, j) = hookf(ndimsi*(j-1)+i)
            end do
        end do
    end if
!
! - INTEGRATION
!
    if (option .eq. 'FULL_MECA' .or. option .eq. 'RAPH_MECA') then
!
!   DEFORMATION THERMIQUES
!       DANS LE REPERE LOCAL
        if (phenom .eq. 'ELAS_ORTH') then
!
            call verift(fami, kpg, ksp, poum, imate, &
                        epsth_anis_=epsth_anis)
            deplth(1) = epsth_anis(1)
            deplth(2) = epsth_anis(2)
            deplth(3) = epsth_anis(3)
!
        else if (phenom .eq. 'ELAS_ISTR') then
!
            call verift(fami, kpg, ksp, poum, imate, &
                        epsth_anis_=epsth_anis)
            deplth(1) = epsth_anis(1)
            deplth(2) = epsth_anis(1)
            deplth(3) = epsth_anis(2)
!
        end if
!
        deplth(4) = 0.d0
        deplth(5) = 0.d0
        deplth(6) = 0.d0
!
!       RECUPERATION DE LA MATRICE DE PASSAGE
!
        call matrot(angl_naut, p)
!
!       PASSAGE DU TENSEUR DES DEFORMATIONS THERMIQUES DANS LE REPERE GLOBAL
!
        do i = 1, 3
            deplth_mat(i, i) = deplth(i)
        end do
        deplth_mat(1, 2) = deplth(4)
        deplth_mat(1, 3) = deplth(5)
        deplth_mat(2, 3) = deplth(6)
        deplth_mat(2, 1) = deplth_mat(1, 2)
        deplth_mat(3, 1) = deplth_mat(1, 3)
        deplth_mat(3, 2) = deplth_mat(2, 3)
!
        call utbtab('ZERO', 3, 3, deplth_mat, p, &
                    work, depgth_mat)
!
        depgth(1) = depgth_mat(1, 1)
        depgth(2) = depgth_mat(2, 2)
        depgth(3) = depgth_mat(3, 3)
        depgth(4) = depgth_mat(1, 2)
        depgth(5) = depgth_mat(1, 3)
        depgth(6) = depgth_mat(2, 3)
!
!   RETRAIT ENDOGENE ET RETRAIT DE DESSICCATION (SCALAIRE)
!       IDENTIQUES DANS LES 2 REPERES L ET G CAR ISOTROPES
        call verifh(fami, kpg, ksp, poum, imate, &
                    depghy)
        call verifs(fami, kpg, ksp, poum, imate, &
                    depgse)
!
!   DEFORMATIONS ANELASTIQUES EPSAXX, EPSAYY, EPSAZZ, EPSAXY, EPSAXZ, EPSAYZ
!       DEFINIES DANS LE REPERE GLOBAL
!
        call verifepsa(fami, kpg, ksp, poum, depgepsa)

!
! CALCUL DES DEFORMATIONS MECANIQUES
!
        do i = 1, ndimsi
            if (i .le. 3) then
                depsme(i) = depstr(i)-depgth(i)-depghy-depgse-depgepsa(i)
            else
                depsme(i) = depstr(i)-2.0*depgth(i)-2.0*depgepsa(i)
            end if
        end do
!
! CONTRAINTE A L ETAT +
        sigm2 = sigm
        do i = 4, ndimsi
            sigm2(i) = sigm2(i)/rac2
        end do

! MODIFICATION DE SIGM POUR PRENDRE EN COMPTE LA VARIATION DE
! COEF ELASTIQUES AVEC LA TEMPERATURE
!
        do i = 1, ndimsi
            epsm2(i) = 0.d0
            do j = 1, ndimsi
                epsm2(i) = epsm2(i)+mkooh(ndimsi*(j-1)+i)*sigm2(j)
            end do
        end do
!
        do i = 1, ndimsi
            sigp(i) = 0.d0
            do j = 1, ndimsi
                sigp(i) = sigp(i)+hookf(ndimsi*(j-1)+i)*(depsme(j)+ &
                                                         epsm2(j))
            end do
        end do
!
! REMISE AU FORMAT ASTER DES VALEURS EXTRA DIAGONALES
        do i = 4, ndimsi
            sigp(i) = sigp(i)*rac2
        end do
    end if
!
    if (option(1:10) .eq. 'RIGI_MECA_' .or. option(1:9) .eq. 'FULL_MECA') then
        do i = 1, ndimsi
            do j = 4, ndimsi
                dsidep(i, j) = dsidep(i, j)*sqrt(2.d0)
            end do
        end do
        do i = 4, ndimsi
            do j = 1, ndimsi
                dsidep(i, j) = dsidep(i, j)*sqrt(2.d0)
            end do
        end do
    end if
!
end subroutine
