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
subroutine dis_choc_frot_nosyme(DD, icodma, ulp, xg, klv, &
                                dpe, varmo, force, varpl)
!
    use te0047_type
    implicit none
#include "asterf_types.h"
#include "asterfort/rcvala.h"
#include "asterfort/utmess.h"
#include "asterfort/utpvgl.h"
#include "asterfort/ut2vgl.h"
!
    type(te0047_dscr), intent(in) :: DD
    integer(kind=8) :: icodma
    real(kind=8) :: ulp(*), klv(*), xg(*), varmo(*), varpl(*), force(*), dpe(*), klvp(144)
!
! --------------------------------------------------------------------------------------------------
!
!     RELATION DE COMPORTEMENT "DIS_CHOC"
!
! --------------------------------------------------------------------------------------------------
! in :
!       icodma : adresse du materiau code
!       ulp    : deplacement
!       dpe    : déplacement d'entrainement (en dynamique)
!       xg     : coordonnees des noeuds repere global
!       varmo  : variables internes (temps moins)
! in/out :
!       klv    : matrice de raideur symétrique initiale     (triangulaire supérieure)
!       klv    : matrice de raideur non-symétrique ou pas   (toujours pleine)
! out :
!       force  : efforts
!       varpl  : variables internes (temps plus)
! --------------------------------------------------------------------------------------------------
! person_in_charge: jean-luc.flejou at edf.fr
!
    integer(kind=8), parameter :: nbre1 = 9
    integer(kind=8) :: codre1(nbre1)
    real(kind=8) :: valre1(nbre1)
    character(len=8) :: nomre1(nbre1)
!   Index des variables internes
    integer(kind=8), parameter :: idepx = 1, idepy = 2, idepz = 3, iidic = 4, idepyp = 5, idepzp = 6
    integer(kind=8), parameter :: ifx = 7, ify = 8, ifz = 9, icalc = 10
!   État du discret : adhérent, glissant, décollé
    integer(kind=8), parameter :: EtatAdher = 0, EtatGliss = 1, EtatDecol = 2
    integer(kind=8), parameter :: EnPlasticite = 2
!
    integer(kind=8) :: indic, ii
    real(kind=8) :: xl(6), xd(3), rignor, rigtan
    real(kind=8) :: coulom, dist12, utotx, utoty, utotz, depx, depy, depz
    real(kind=8) :: lambda, fort, dist0, rtmp
!
    character(len=32) :: messak(3)
!
    data nomre1/'RIGI_NOR', 'RIGI_TAN', 'AMOR_NOR', 'AMOR_TAN', 'COULOMB', &
        'DIST_1', 'DIST_2', 'JEU', 'CONTACT'/
!
! --------------------------------------------------------------------------------------------------
!
!   Initialisation
    xl(:) = 0.d0; xd(:) = 0.d0; fort = 0.d0
!   Coordonnées dans le repere local
    if (DD%ndim .eq. 3) then
        call utpvgl(DD%nno, 3, DD%pgl, xg, xl)
    else
        call ut2vgl(DD%nno, 2, DD%pgl, xg, xl)
    end if
!   Raideurs du discret
!       ==> Elles sont surchargées par celles du matériau
    valre1(:) = 0.d0
    valre1(1) = klv(1); valre1(2) = klv(3)
!   Caractéristiques du matériau
    call rcvala(icodma, ' ', 'DIS_CONTACT', 0, ' ', &
                [0.0d0], nbre1, nomre1, valre1, codre1, &
                0, nan='NON')
    if (nint(valre1(9)) .ne. 0) then
        messak(1) = 'DIS_CONTACT'
        messak(2) = 'DIS_CHOC (cas non symétrique)'
        messak(3) = '"1D"'
        call utmess('F', 'DISCRETS_35', nk=3, valk=messak)
    end if
    rignor = abs(valre1(1))
    rigtan = abs(valre1(2))
    coulom = abs(valre1(5))
!
!   Élément avec 2 noeuds
    if (DD%nno .eq. 2) then
        dist12 = valre1(6)+valre1(7)
        utotx = ulp(1+DD%nc)-ulp(1)+dpe(4)-dpe(1)
        utoty = ulp(2+DD%nc)-ulp(2)+dpe(5)-dpe(2)
        if (DD%ndim .eq. 3) then
            utotz = ulp(3+DD%nc)-ulp(3)+dpe(6)-dpe(3)
        end if
        do ii = 1, DD%ndim
            xd(ii) = xl(ii+3)-xl(ii)
        end do
        depx = xd(1)-dist12+utotx
        depy = xd(2)+utoty
        if (DD%ndim .eq. 3) then
            depz = xd(3)+utotz
        else
            depz = 0.d0
        end if
!
!   Élément avec 1 noeud
    else
        dist12 = valre1(6)
        utotx = ulp(1)+dpe(1)
        utoty = ulp(2)+dpe(2)
        if (DD%ndim .eq. 3) then
            utotz = ulp(3)+dpe(3)
        end if
        dist0 = valre1(8)
        depx = utotx+dist0-dist12
        depy = ulp(2)
        if (DD%ndim .eq. 3) then
            depz = ulp(3)
        else
            depz = 0.d0
        end if
    end if
    varpl(idepx) = depx
    varpl(idepy) = depy
    varpl(idepz) = depz
    varpl(icalc) = EnPlasticite
!   Calcul des variables internes et efforts
    force(1:3) = 0.d0
    if (DD%lVari .or. DD%lVect) then
        if (depx .le. 0.d0) then
            force(1) = rignor*depx
            force(2) = rigtan*(depy-varmo(idepyp))
            if (DD%ndim .eq. 3) then
                force(3) = rigtan*(depz-varmo(idepzp))
            end if
            fort = (force(2)**2+force(3)**2)**0.5
            if (fort .gt. abs(coulom*force(1))) then
                lambda = 1.0-abs(coulom*force(1))/fort
                varpl(idepyp) = varmo(idepyp)+lambda*force(2)/rigtan
                force(2) = rigtan*(depy-varpl(idepyp))
                if (DD%ndim .eq. 3) then
                    varpl(idepzp) = varmo(idepzp)+lambda*force(3)/rigtan
                    force(3) = rigtan*(depz-varpl(idepzp))
                end if
                varpl(iidic) = EtatGliss
            else
                varpl(idepyp) = varmo(idepyp)
                if (DD%ndim .eq. 3) then
                    varpl(idepzp) = varmo(idepzp)
                end if
                varpl(iidic) = EtatAdher
            end if
        else
            force(1) = 0.0
            force(2) = 0.0
            force(3) = 0.0
            fort = 0.0
            varpl(idepyp) = depy
            varpl(idepzp) = depz
            varpl(iidic) = EtatDecol
        end if
    end if
!   Calcul de la matrice complete
    if (DD%lMatr) then
        if (DD%option(1:9) .eq. 'FULL_MECA') then
            indic = nint(varpl(iidic))
        else
            indic = nint(varmo(iidic))
            force(1) = rignor*depx
            force(2) = rigtan*(depy-varmo(idepyp))
            if (DD%ndim .eq. 3) then
                force(3) = rigtan*(depz-varmo(idepzp))
            end if
            fort = (force(2)**2+force(3)**2)**0.50
        end if
        klv(1:DD%neq*DD%neq) = 0.d0
        if (DD%nno .eq. 2) then
            ! Cas d'un élément à 2 noeuds
            if (indic .eq. EtatAdher) then
                klv(idx(1, 1)) = rignor
                klv(idx(4, 1)) = -rignor
                klv(idx(1, 4)) = -rignor
                klv(idx(4, 4)) = rignor
                klv(idx(2, 2)) = rigtan
                klv(idx(5, 2)) = -rigtan
                klv(idx(2, 5)) = -rigtan
                klv(idx(5, 5)) = rigtan
                if (DD%ndim .eq. 3) then
                    klv(idx(3, 3)) = rigtan
                    klv(idx(6, 3)) = -rigtan
                    klv(idx(3, 6)) = -rigtan
                    klv(idx(6, 6)) = rigtan
                end if
            else if (indic .eq. EtatGliss) then
                klv(idx(1, 1)) = rignor
                klv(idx(4, 1)) = -rignor
                klv(idx(1, 4)) = -rignor
                klv(idx(4, 4)) = rignor
                if ((coulom .gt. 0) .and. (fort .gt. 0)) then
                    rtmp = -rignor*rigtan*coulom/fort*(depy-varmo(idepyp))
                    klv(idx(2, 1)) = rtmp
                    klv(idx(5, 1)) = -rtmp
                    klv(idx(2, 4)) = -rtmp
                    klv(idx(5, 4)) = rtmp
                    rtmp = -coulom*force(1)*rigtan/fort*(1.0-rigtan**2*(depy-varmo(idepyp))**2 &
                                                         /fort**2)
                    klv(idx(2, 2)) = rtmp
                    klv(idx(2, 5)) = -rtmp
                    klv(idx(5, 2)) = -rtmp
                    klv(idx(5, 5)) = rtmp
                    if (DD%ndim .eq. 3) then
                        rtmp = -rignor*rigtan*coulom/fort*(depz-varmo(idepzp))
                        klv(idx(3, 1)) = rtmp
                        klv(idx(6, 1)) = -rtmp
                        klv(idx(3, 4)) = -rtmp
                        klv(idx(6, 4)) = rtmp
                        rtmp = -coulom*force(1)*rigtan/fort*(1.0-rigtan**2*(depz-varmo(idepzp))**2 &
                                                             /fort**2)
                        klv(idx(3, 3)) = rtmp
                        klv(idx(3, 6)) = -rtmp
                        klv(idx(6, 3)) = -rtmp
                        klv(idx(6, 6)) = rtmp
                        rtmp = rigtan**3*coulom*force(1)/fort**3*(depy-varmo(idepyp))* &
                               (depz-varmo(idepzp))
                        klv(idx(2, 3)) = rtmp
                        klv(idx(3, 2)) = rtmp
                        klv(idx(2, 6)) = -rtmp
                        klv(idx(6, 2)) = -rtmp
                        klv(idx(5, 3)) = -rtmp
                        klv(idx(3, 5)) = -rtmp
                        klv(idx(5, 6)) = rtmp
                        klv(idx(6, 5)) = rtmp
                    end if
                end if
            end if
        else
            ! Cas d'un élément à 1 noeud
            if (indic .eq. EtatAdher) then
                klv(idx(1, 1)) = rignor
                klv(idx(2, 2)) = rigtan
                klv(idx(3, 3)) = rigtan
            else if (indic .eq. EtatGliss) then
                klv(idx(1, 1)) = rignor
                if (coulom .gt. 0) then
                    rtmp = -rignor*rigtan*coulom/fort*(depy-varmo(idepyp))
                    klv(idx(2, 1)) = rtmp
                    rtmp = -coulom*force(1)*rigtan/fort*(1.0-rigtan**2*(depy-varmo(idepyp))**2 &
                                                         /fort**2)
                    klv(idx(2, 2)) = rtmp
                    if (DD%ndim .eq. 3) then
                        rtmp = -coulom*force(1)*rigtan/fort*(1.0-rigtan**2*(depz-varmo(idepzp))**2 &
                                                             /fort**2)
                        klv(idx(3, 3)) = rtmp
                        rtmp = -rignor*rigtan*coulom/fort*(depz-varmo(idepzp))
                        klv(idx(3, 1)) = rtmp
                        rtmp = rigtan**3*coulom*force(1)/fort**3*(depy-varmo(idepyp))* &
                               (depz-varmo(idepzp))
                        klv(idx(2, 3)) = rtmp
                        klv(idx(3, 2)) = rtmp
                    end if
                end if
            end if
        end if
    end if
    ! Elément élastique en parallèle
    varpl(ifx) = force(1)
    varpl(ify) = force(2)
    if (DD%ndim .eq. 3) then
        varpl(ifz) = force(3)
    end if

! Fonctions
contains
! --------------------------------------------------------------------------------------------------
! Indices de matrices stockées au format vecteur
!   12*12  : klv(ii) : klv(idx(i,j))   DIS_TR_L
!    6*6   : klv(ii) : klv(idx(i,j))   DIS_T_L ou DIS_TR_N
!    3*3   : klv(ii) : klv(idx(i,j))   DIS_T_N
    function idx(i, j)
        integer(kind=8) :: i, j
        integer(kind=8) :: i_offset, j_offset
        integer(kind=8) :: idx
        !
        if ((DD%nno .eq. 2) .and. (DD%nc .eq. 6)) then
            ! Cas d'un discret de type DIS_TR_L
            i_offset = merge(0, 3, i > 3)
            j_offset = merge(0, 2, j > 3)
        else
            i_offset = 0
            j_offset = 0
        end if
        !
        idx = (i+i_offset)+(j+j_offset-1)*DD%neq
    end function idx
! --------------------------------------------------------------------------------------------------
!
end subroutine
