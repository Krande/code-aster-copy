! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
subroutine dis_choc_frot_syme(for_discret, icodma, ulp, xg, klv, &
                              dvl, dpe, dve, Predic, force, &
                              varmo, varpl)
!
    use te0047_type
    implicit none
!
#include "asterf_types.h"
#include "asterc/r8prem.h"
#include "asterfort/diraidklv.h"
#include "asterfort/diklvraid.h"
#include "asterfort/rcvala.h"
#include "asterfort/ut2vgl.h"
#include "asterfort/utpvgl.h"
#include "blas/dcopy.h"
!
    type(te0047_dscr), intent(in) :: for_discret
    integer :: icodma
    real(kind=8) :: ulp(*), dvl(*)
    real(kind=8) :: dpe(*), dve(*)
    real(kind=8) :: klv(*), xg(*)
    real(kind=8) :: varmo(*), varpl(*), force(*)
    aster_logical :: Predic
!
! --------------------------------------------------------------------------------------------------
!
!     RELATION DE COMPORTEMENT "DIS_CHOC"
!
! --------------------------------------------------------------------------------------------------
! in :
!       icodma  : adresse du materiau code
!       ulp     : deplacement
!       xg      : coordonnees des noeuds repere global
!       dvl     : vitesse
!       dpe     : déplacement d'entrainement
!       dve     : vitesse d'entrainement
!       varmo   : variables internes (temps moins)
! in/out :
!       klv     : matrice de raideur symétrique initiale     (triangulaire supérieure)
!       klv     : matrice de raideur symétrique              (triangulaire supérieure)
! out :
!       force   : efforts
!       varpl   : variables internes (temps plus)
! --------------------------------------------------------------------------------------------------
! person_in_charge: jean-luc.flejou at edf.fr
!
    integer, parameter :: nbre1 = 8
    real(kind=8) :: valre1(nbre1)
    integer :: codre1(nbre1)
    character(len=8) :: nomre1(nbre1)
!   Index des variables internes
    integer, parameter :: idepx = 1, idepy = 2, idepz = 3, iidic = 4, idepyp = 5, idepzp = 6
    integer, parameter :: ifx = 7, ify = 8, ifz = 9, icalc = 10
!   État du discret : adhérent, glissant, décollé
    integer, parameter :: EtatAdher = 0, EtatGliss = 1, EtatDecol = 2
    integer, parameter :: EnVitesse = 1, EnPlasticite = 2
!
    integer :: ii
    real(kind=8) :: xl(6), xd(3), dirl(6), raide(6), rignor, rigtan
    real(kind=8) :: coulom, dist12, utot, vit2, vit3, depx, depy, depz, psca
    real(kind=8) :: vitt, vity, vitz, fort
    blas_int :: b_incx, b_incy, b_n
!
    data nomre1/'RIGI_NOR', 'RIGI_TAN', 'AMOR_NOR', 'AMOR_TAN', &
        'COULOMB', 'DIST_1', 'DIST_2', 'JEU'/
! ----------------------------------------------------------------------
!
!   Définition des parametres
    xl(:) = 0.d0; dirl(:) = 0.d0; xd(:) = 0.d0
!   Coordonnees dans le repere local
    if (for_discret%ndim .eq. 3) then
        call utpvgl(for_discret%nno, 3, for_discret%pgl, xg, xl)
    else
        call ut2vgl(for_discret%nno, 2, for_discret%pgl, xg, xl)
    end if
!   Raideurs du discret
!       ==> Elles sont surchargées par celles du matériau
    call diraidklv(for_discret%nomte, raide, klv)
    valre1(:) = 0.0
    valre1(1) = raide(1)
!   Caractéristiques du matériau
    call rcvala(icodma, ' ', 'DIS_CONTACT', 0, ' ', &
                [0.0d0], nbre1, nomre1, valre1, codre1, &
                0, nan='NON')
    rignor = valre1(1)
    rigtan = valre1(2)
    coulom = valre1(5)
!
    varpl(icalc) = EnVitesse
!   Élément avec 2 noeuds
    if (for_discret%nno .eq. 2) then
        dist12 = valre1(6)+valre1(7)
! Dans l'axe du discret
        utot = ulp(1+for_discret%nc)-ulp(1)
! Vitesse tangente
        vit2 = dvl(2+for_discret%nc)-dvl(2)
        vit3 = 0.0
        if (for_discret%ndim .eq. 3) then
            vit3 = dvl(3+for_discret%nc)-dvl(3)
        end if
! Longueur du discret
        do ii = 1, for_discret%ndim
            xd(ii) = xl(for_discret%ndim+ii)-xl(ii)
        end do
        b_n = to_blas_int(for_discret%ndim)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, dpe(1), b_incx, dirl, b_incy)
        b_n = to_blas_int(for_discret%ndim)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, dpe(1+for_discret%nc), b_incx, dirl(4), b_incy)
        depx = xd(1)-dist12+utot+dirl(4)-dirl(1)
        depx = depx-r8prem()
        depy = xd(2)+ulp(2+for_discret%nc)-ulp(2)+dirl(5)-dirl(2)
        depz = 0.0
        if (for_discret%ndim .eq. 3) then
            depz = xd(3)+ulp(3+for_discret%nc)-ulp(3)+dirl(6)-dirl(3)
        end if
        b_n = to_blas_int(for_discret%ndim)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, dve(1), b_incx, dirl, b_incy)
        b_n = to_blas_int(for_discret%ndim)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, dve(1+for_discret%nc), b_incx, dirl(4), b_incy)
! Vitesse tangente
        vity = vit2+dirl(5)-dirl(2)
        vitz = 0.0
        if (for_discret%ndim .eq. 3) then
            vitz = vit3+dirl(6)-dirl(3)
        end if
        if (depx .le. 0.0) then
            force(1) = rignor*depx
            if (force(1) .gt. 0.0) force(1) = 0.0
            psca = varmo(ify)*vity+varmo(ifz)*vitz
            if ((psca .ge. 0.0) .and. (nint(varmo(iidic)) .eq. EtatGliss)) then
                vitt = (vity**2+vitz**2)**0.5d0
                force(2) = 0.0
                force(3) = 0.0
                if (vitt .gt. 0.00) then
                    force(2) = -coulom*force(1)*vity/vitt
                    force(3) = -coulom*force(1)*vitz/vitt
                end if
                varpl(iidic) = EtatGliss
            else
                force(2) = rigtan*(depy-varmo(idepy))+varmo(ify)
                force(3) = rigtan*(depz-varmo(idepz))+varmo(ifz)
                varpl(iidic) = EtatAdher
                fort = (force(2)**2+force(3)**2)**0.5d0
                if (fort .gt. abs(coulom*force(1))) then
                    vitt = (vity**2+vitz**2)**0.5d0
                    force(2) = 0.0
                    force(3) = 0.0
                    if (vitt .gt. 0.0) then
                        force(2) = -coulom*force(1)*vity/vitt
                        force(3) = -coulom*force(1)*vitz/vitt
                        varpl(iidic) = EtatGliss
                    end if
                end if
            end if
            varpl(ifx) = force(1)
            varpl(ify) = force(2)
            varpl(ifz) = force(3)
            !
            if (abs(rigtan) .gt. r8prem()) then
                varpl(idepyp) = depy-varmo(ify)/rigtan
                varpl(idepzp) = depz-varmo(ifz)/rigtan
            end if
            !
            force(2) = force(2)+raide(2)*(ulp(2+for_discret%nc)-ulp(2))
            if (for_discret%ndim .eq. 3) then
                force(3) = force(3)+raide(3)*(ulp(3+for_discret%nc)-ulp(3))
            end if
! Actualisation de la matrice de raideur
            raide(1) = rignor
            call diklvraid(for_discret%nomte, klv, raide)
        else
            force(1:3) = 0.0
            varpl(ifx) = 0.0
            varpl(ify) = 0.0
            varpl(ifz) = 0.0
            varpl(iidic) = EtatDecol
            varpl(idepyp) = 0.0
            varpl(idepzp) = 0.0
            klv(1:78) = 0.0
        end if
!
!   Élément avec 1 noeud
    else
        dist12 = valre1(8)-valre1(6)
! Vitesse tangente
        vit2 = dvl(2)
        vit3 = 0.0
        if (for_discret%ndim .eq. 3) then
            vit3 = dvl(3)
        end if
        b_n = to_blas_int(for_discret%ndim)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, dpe(1), b_incx, dirl, b_incy)
        depx = ulp(1)+dist12+dirl(1)
        depy = ulp(2)+dirl(2)
        depz = 0.0
        if (for_discret%ndim .eq. 3) then
            depz = ulp(3)+dirl(3)
        end if
        b_n = to_blas_int(for_discret%ndim)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, dve(1), b_incx, dirl, b_incy)
! Vitesse tangente
        vity = vit2+dirl(2)
        vitz = 0.0
        if (for_discret%ndim .eq. 3) then
            vitz = vit3+dirl(3)
        end if
        if (depx .le. 0.0d0) then
            force(1) = rignor*depx
            if (force(1) .gt. 0.0) force(1) = 0.0
            psca = varmo(ify)*vity+varmo(ifz)*vitz
            if ((psca .ge. 0.0) .and. (nint(varmo(iidic)) .eq. EtatGliss)) then
                vitt = (vity**2+vitz**2)**0.5d0
                force(2) = 0.0
                force(3) = 0.0
                if (vitt .gt. 0.0) then
                    force(2) = -coulom*force(1)*vity/vitt
                    force(3) = -coulom*force(1)*vitz/vitt
                end if
                varpl(iidic) = EtatGliss
            else
                force(2) = rigtan*(depy-varmo(idepy))+varmo(ify)
                force(3) = rigtan*(depz-varmo(idepz))+varmo(ifz)
                varpl(iidic) = EtatAdher
                fort = (force(2)**2+force(3)**2)**0.5d0
                if (fort .gt. abs(coulom*force(1))) then
                    vitt = (vity**2+vitz**2)**0.5d0
                    force(2) = 0.0
                    force(3) = 0.0
                    if (vitt .gt. 0.0) then
                        force(2) = -coulom*force(1)*vity/vitt
                        force(3) = -coulom*force(1)*vitz/vitt
                        varpl(iidic) = EtatGliss
                    end if
                end if
            end if
            varpl(ifx) = force(1)
            varpl(ify) = force(2)
            varpl(ifz) = force(3)
            !
            if (abs(rigtan) .gt. r8prem()) then
                varpl(idepyp) = depy-varmo(ify)/rigtan
                varpl(idepzp) = depz-varmo(ifz)/rigtan
            end if
            !
            force(2) = force(2)+raide(2)*ulp(2)
            if (for_discret%ndim .eq. 3) then
                force(3) = force(3)+raide(3)*ulp(3)
            end if
! Actualisation de la matrice de raideur
            raide(1) = rignor
            call diklvraid(for_discret%nomte, klv, raide)
        else
            force(1:3) = 0.0
            varpl(ifx) = 0.0
            varpl(ify) = 0.0
            varpl(ifz) = 0.0
            varpl(iidic) = EtatDecol
            varpl(idepyp) = 0.0
            varpl(idepzp) = 0.0
            klv(1:78) = 0.0
        end if
    end if
    varpl(idepx) = depx
    varpl(idepy) = depy
    varpl(idepz) = depz
!
!   Prédiction en dynamique, on retourne les efforts précédents
!       Si on passe d'une formulation 'plastique' à une en 'vitesse'
!       On le fait à la fin, la raideur doit être mise comme il faut
    if (Predic .and. (nint(varmo(icalc)) .eq. EnPlasticite)) then
! Les efforts précédents
        force(1) = varmo(ifx)
        force(2) = varmo(ify)
        force(3) = varmo(ifz)
! On remet les varpl comme il faut. Elles ont peut-être été modifiées
        varpl(ifx) = varmo(ifx)
        varpl(ify) = varmo(ify)
        varpl(ifz) = varmo(ifz)
    end if
!
end subroutine
