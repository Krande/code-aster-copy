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
subroutine dis_choc_frot_nosyme(for_discret, icodma, ulp, xg, klv, &
                                dpe, varmo, force, varpl)
!
    use te0047_type
    implicit none
#include "asterf_types.h"
#include "asterfort/rcvala.h"
#include "asterfort/utmess.h"
#include "asterfort/utpvgl.h"
!
    type(te0047_dscr), intent(in) :: for_discret
    integer(kind=8) :: icodma
    real(kind=8) :: ulp(*), klv(*), xg(*), varmo(*), varpl(*), force(*), dpe(*), klvp(36)
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
    integer(kind=8), parameter :: nbre1 = 11
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
    real(kind=8) :: lambda, fort, dist0, rtmp, kp, kt
!
    character(len=32) :: messak(3)
!
    data nomre1/'RIGI_NOR', 'RIGI_TAN', 'AMOR_NOR', 'AMOR_TAN', 'COULOMB', &
        'DIST_1', 'DIST_2', 'JEU', 'CONTACT', 'KP', 'KT'/
!
! --------------------------------------------------------------------------------------------------
!
! Pour les index de matrice stockées en vecteur
!   6*6  : klv(idx) : klv( id6(i,j) )   T_L
!   3*3  : klv(idx) : klv( id3(i,j) )   T_N
#define id6(i,j) i+(j-1)*6
#define id3(i,j) i+(j-1)*3
!
!   Initialisation
    xl(:) = 0.d0; xd(:) = 0.d0; fort = 0.d0
!   Coordonnées dans le repere local
    call utpvgl(for_discret%nno, 3, for_discret%pgl, xg, xl)
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
    kp = valre1(10)
    kt = valre1(11)
!
!   Élément avec 2 noeuds
    if (for_discret%nno .eq. 2) then
        dist12 = valre1(6)+valre1(7)
        utotx = ulp(1+for_discret%nc)-ulp(1)+dpe(4)-dpe(1)
        utoty = ulp(2+for_discret%nc)-ulp(2)+dpe(5)-dpe(2)
        utotz = ulp(3+for_discret%nc)-ulp(3)+dpe(6)-dpe(3)
        do ii = 1, 3
            xd(ii) = xl(ii+3)-xl(ii)
        end do
        depx = xd(1)-dist12+utotx
        depy = xd(2)+utoty
        depz = xd(3)+utotz
!
!   Élément avec 1 noeud
    else
        dist12 = valre1(6)
        utotx = ulp(1)
        utoty = ulp(2)
        utotz = ulp(3)
        dist0 = valre1(8)
        depx = utotx+dist12-dist0
        depy = ulp(2)
        depz = ulp(3)
    end if
    varpl(idepx) = depx
    varpl(idepy) = depy
    varpl(idepz) = depz
    varpl(icalc) = EnPlasticite
!   Calcul des variables internes et efforts
    if (for_discret%lVari .or. for_discret%lVect) then
        if (depx .le. 0.d0) then
            force(1) = rignor*depx
            force(2) = rigtan*(depy-varmo(idepyp))
            force(3) = rigtan*(depz-varmo(idepzp))
            fort = (force(2)**2+force(3)**2)**0.5
            if (fort .gt. abs(coulom*force(1))) then
                lambda = 1.0-abs(coulom*force(1))/fort
                varpl(idepyp) = varmo(idepyp)+lambda*force(2)/rigtan
                varpl(idepzp) = varmo(idepzp)+lambda*force(3)/rigtan
                force(2) = rigtan*(depy-varpl(idepyp))
                force(3) = rigtan*(depz-varpl(idepzp))
                varpl(iidic) = EtatGliss
            else
                varpl(idepyp) = varmo(idepyp)
                varpl(idepzp) = varmo(idepzp)
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
    if (for_discret%lMatr) then
        if (for_discret%option(1:9) .eq. 'FULL_MECA') then
            indic = nint(varpl(iidic))
        else
            indic = nint(varmo(iidic))
            force(1) = rignor*depx
            force(2) = rigtan*(depy-varmo(idepyp))
            force(3) = rigtan*(depz-varmo(idepzp))
            fort = (force(2)**2+force(3)**2)**0.50
        end if
        if (for_discret%nno .eq. 2) then
            klv(1:36) = 0.d0
            if (indic .eq. EtatAdher) then
                klv(id6(1, 1)) = rignor
                klv(id6(4, 1)) = -rignor
                klv(id6(1, 4)) = -rignor
                klv(id6(4, 4)) = rignor
                klv(id6(2, 2)) = rigtan
                klv(id6(3, 3)) = rigtan
                klv(id6(5, 5)) = rigtan
                klv(id6(6, 6)) = rigtan
                klv(id6(5, 2)) = -rigtan
                klv(id6(2, 5)) = -rigtan
                klv(id6(6, 3)) = -rigtan
                klv(id6(3, 6)) = -rigtan
            else if (indic .eq. EtatGliss) then
                klv(id6(1, 1)) = rignor
                klv(id6(4, 1)) = -rignor
                klv(id6(1, 4)) = -rignor
                klv(id6(4, 4)) = rignor
                if ((coulom .gt. 0) .and. (fort .gt. 0)) then
                    rtmp = -rignor*rigtan*coulom/fort*(depy-varmo(idepyp))
                    klv(id6(2, 1)) = rtmp
                    klv(id6(5, 1)) = -rtmp
                    klv(id6(2, 4)) = -rtmp
                    klv(id6(5, 4)) = rtmp
                    rtmp = -rignor*rigtan*coulom/fort*(depz-varmo(idepzp))
                    klv(id6(3, 1)) = rtmp
                    klv(id6(6, 1)) = -rtmp
                    klv(id6(3, 4)) = -rtmp
                    klv(id6(6, 4)) = rtmp
                 rtmp = -coulom*force(1)*rigtan/fort*(1.0-rigtan**2*(depy-varmo(idepyp))**2/fort**2)
                    klv(id6(2, 2)) = rtmp
                    klv(id6(2, 5)) = -rtmp
                    klv(id6(5, 2)) = -rtmp
                    klv(id6(5, 5)) = rtmp
                 rtmp = -coulom*force(1)*rigtan/fort*(1.0-rigtan**2*(depz-varmo(idepzp))**2/fort**2)
                    klv(id6(3, 3)) = rtmp
                    klv(id6(3, 6)) = -rtmp
                    klv(id6(6, 3)) = -rtmp
                    klv(id6(6, 6)) = rtmp
                  rtmp = rigtan**3*coulom*force(1)/fort**3*(depy-varmo(idepyp))*(depz-varmo(idepzp))
                    klv(id6(2, 3)) = rtmp
                    klv(id6(3, 2)) = rtmp
                    klv(id6(2, 6)) = -rtmp
                    klv(id6(6, 2)) = -rtmp
                    klv(id6(5, 3)) = -rtmp
                    klv(id6(3, 5)) = -rtmp
                    klv(id6(5, 6)) = rtmp
                    klv(id6(6, 5)) = rtmp
                end if
            end if
            ! Elément élastique en parallèle
            klvp = 0.d0
            klvp(id6(1, 1)) = kp
            klvp(id6(4, 1)) = -kp
            klvp(id6(1, 4)) = -kp
            klvp(id6(4, 4)) = kp
            klvp(id6(2, 2)) = kt
            klvp(id6(3, 3)) = kt
            klvp(id6(5, 5)) = kt
            klvp(id6(6, 6)) = kt
            klvp(id6(5, 2)) = -kt
            klvp(id6(2, 5)) = -kt
            klvp(id6(6, 3)) = -kt
            klvp(id6(3, 6)) = -kt
            ! Cumul des raideurs
            klv(1:36) = klv(1:36)+klvp(1:36)
        else
            klv(1:9) = 0.d0
            if (indic .eq. EtatAdher) then
                klv(id3(1, 1)) = rignor
                klv(id3(2, 2)) = rigtan
                klv(id3(3, 3)) = rigtan
            else if (indic .eq. EtatGliss) then
                klv(id3(1, 1)) = rignor
                if (coulom .gt. 0) then
                    rtmp = -rignor*rigtan*coulom/fort*(depy-varmo(idepyp))
                    klv(id3(2, 1)) = rtmp
                    rtmp = -rignor*rigtan*coulom/fort*(depz-varmo(idepzp))
                    klv(id3(3, 1)) = rtmp
                 rtmp = -coulom*force(1)*rigtan/fort*(1.0-rigtan**2*(depy-varmo(idepyp))**2/fort**2)
                    klv(id3(2, 2)) = rtmp
                 rtmp = -coulom*force(1)*rigtan/fort*(1.0-rigtan**2*(depz-varmo(idepzp))**2/fort**2)
                    klv(id3(3, 3)) = rtmp
                  rtmp = rigtan**3*coulom*force(1)/fort**3*(depy-varmo(idepyp))*(depz-varmo(idepzp))
                    klv(id3(2, 3)) = rtmp
                    klv(id3(3, 2)) = rtmp
                end if
            end if
            ! Elément élastique en parallèle
            klvp = 0.d0
            klvp(id3(1, 1)) = rignor
            klvp(id3(2, 2)) = rigtan
            klvp(id3(3, 3)) = rigtan
            ! Cumul des raideurs
            klv(1:9) = klv(1:9)+klvp(1:9)

        end if
    end if
    ! Elément élastique en parallèle
    force(1) = force(1)+kp*utotx
    force(2) = force(2)+kt*utoty
    force(3) = force(3)+kt*utotz
    varpl(ifx) = force(1)
    varpl(ify) = force(2)
    varpl(ifz) = force(3)
!
end subroutine
