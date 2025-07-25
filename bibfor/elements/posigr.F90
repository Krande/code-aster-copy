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

subroutine posigr(nomte, efge, sigm)
!
!
! --------------------------------------------------------------------------------------------------
!
! Calcul du vecteur élémentaire contrainte réel pour les éléments de poutre EULER et TIMOSHENKO
!
! --------------------------------------------------------------------------------------------------
!
    implicit none
!
    character(len=*) :: nomte
    real(kind=8) :: sigm(*), efge(12)
!
#include "jeveux.h"
#include "asterfort/poutre_modloc.h"
#include "asterfort/utmess.h"
#include "asterfort/get_value_mode_local.h"
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8)      :: itsec
    real(kind=8) :: a, a2, hy1, hy2, hz1, hz2, r1, r2
    real(kind=8) :: zero, deux
    real(kind=8) :: smf1, smf2, smfy1, smfy2, smfz1, smfz2, sn1, sn2
    real(kind=8) :: xiy, xiy2, xiz, xiz2
    real(kind=8) :: xfly, xflz, xsiy, xsiz, xxy, xxz
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nb_cara = 6
    real(kind=8) :: vale_cara(nb_cara)
    character(len=8) :: noms_cara(nb_cara)
    data noms_cara/'A1', 'IY1', 'IZ1', 'A2', 'IY2', 'IZ2'/
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nb_cara1 = 7
    real(kind=8) :: vale_cara1(nb_cara1)
    character(len=8) :: noms_cara1(nb_cara1)
    data noms_cara1/'HY1', 'HZ1', 'HY2', 'HZ2', 'R1', 'R2', 'TSEC'/
!
    integer(kind=8)             :: retp(4), iret
    real(kind=8)        :: valr(4)
    character(len=8)    :: valp(4)
!
! --------------------------------------------------------------------------------------------------
!
    zero = 0.d0
    deux = 2.d0
!
!   Récupération des caractéristiques générales des sections
    call poutre_modloc('CAGNPO', noms_cara, nb_cara, lvaleur=vale_cara)
!   Section initiale
    a = vale_cara(1)
    xiy = vale_cara(2)
    xiz = vale_cara(3)
!   Section finale
    a2 = vale_cara(4)
    xiy2 = vale_cara(5)
    xiz2 = vale_cara(6)
!
    if (nomte .eq. 'MECA_POU_D_TG') then
        a2 = a
    else if (nomte .eq. 'MECA_POU_D_T') then
        valp(1:4) = ['C_FLEX_Y', 'C_FLEX_Z', 'I_SIGM_Y', 'I_SIGM_Z']
        call get_value_mode_local('PCAARPO', valp, valr, iret, retpara_=retp)
        xfly = 1.0; xflz = 1.0; xsiy = 1.0; xsiz = 1.0
        if (retp(1) .eq. 0) xfly = valr(1)
        if (retp(2) .eq. 0) xflz = valr(2)
        if (retp(3) .eq. 0) xsiy = valr(3)
        if (retp(4) .eq. 0) xsiz = valr(4)
!       prise en compte de l'indice de flexibilité
        xiy = xiy/xfly
        xiz = xiz/xflz
        xiy2 = xiy2/xfly
        xiz2 = xiz2/xflz
!       prise en compte de l'indice de contraintes
        xxy = xsiy/xfly
        xxz = xsiz/xflz
        xiy = xiy/xxy
        xiz = xiz/xxz
        xiy2 = xiy2/xxy
        xiz2 = xiz2/xxz
    end if
!
!   caractéristiques des sections cercle et rectangle
    call poutre_modloc('CAGEPO', noms_cara1, nb_cara1, lvaleur=vale_cara1)
    itsec = nint(vale_cara1(7))
!
!   sxx calculé à partir des 2 flexions et de l'effort normal
    sn1 = -efge(1)/a
    sn2 = efge(7)/a2
!
!   section rectangulaire: le max  et le min sont obtenus sur les coins
    if (itsec .eq. 1) then
        hy1 = vale_cara1(1)
        hz1 = vale_cara1(2)
        hy2 = vale_cara1(3)
        hz2 = vale_cara1(4)
        smfy1 = abs(efge(5)/xiy*hz1/deux)
        smfz1 = abs(efge(6)/xiz*hy1/deux)
        smfy2 = abs(efge(11)/xiy2*hz2/deux)
        smfz2 = abs(efge(12)/xiz2*hy2/deux)
        sigm(1) = sn1-smfy1-smfz1
        sigm(2) = sn1+smfy1+smfz1
        sigm(3) = sn2-smfy2-smfz2
        sigm(4) = sn2+smfy2+smfz2
!
!   section circulaire: xiy = xiz.
    else if (itsec .eq. 2) then
!       formule utilisee :  a cos(t) + b sin(t) = R cos(t-s)
!                         avec R= sqrt(a^2+b^2) et tan(s)= b/a
!       donc max de a cos(t) + b sin(t) = R
!       et   min de a cos(t) + b sin(t) = -R
        r1 = vale_cara1(5)
        r2 = vale_cara1(6)
        smf1 = (r1/xiy)*sqrt(efge(5)**2+efge(6)**2)
        smf2 = (r2/xiy2)*sqrt(efge(11)**2+efge(12)**2)
        sigm(1) = sn1-smf1
        sigm(2) = sn1+smf1
        sigm(3) = sn2-smf2
        sigm(4) = sn2+smf2
!
!   section generale: interdit
    else if (itsec .eq. 0) then
        call utmess('A', 'ELEMENTS4_4')
    end if
!
end subroutine
