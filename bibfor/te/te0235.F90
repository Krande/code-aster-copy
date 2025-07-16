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

subroutine te0235(option, nomte)
!
!
! --------------------------------------------------------------------------------------------------
!
!     CALCULE LA MATRICE DE RAIDEUR CENTRIFUGE DES ELEMENTS DE POUTRE
!     SANS/AVEC GAUCHISSEMENT, MULTIFIBRE OU NON
!
! --------------------------------------------------------------------------------------------------
!
!   IN
!       OPTION  : NOM DE L'OPTION A CALCULER
!           'RIGI_MECA_RO'      : CALCUL DE LA MATRICE DE RAIDEUR CENTRIFUGE
!       NOMTE   : NOM DU TYPE ELEMENT
!           'MECA_POU_D_E'  : POUTRE DROITE D'EULER (SANS GAUCHISSEMENT)
!           'MECA_POU_D_T'  : POUTRE DROITE DE TIMOSHENKO (SANS GAUCHISSEMENT)
!           'MECA_POU_D_TG' : POUTRE DROITE DE TIMOSHENKO (AVEC GAUCHISSEMENT)
!           'MECA_POU_D_TGM': POUTRE DROITE DE TIMOSHENKO MULTI-FIBRES (AVEC GAUCHISSEMENT)
! --------------------------------------------------------------------------------------------------
!
    implicit none
    character(len=*) :: option, nomte
!
#include "jeveux.h"
#include "asterfort/jevech.h"
#include "asterfort/get_elas_id.h"
#include "asterfort/lonele.h"
#include "asterfort/masstg.h"
#include "asterfort/matrot.h"
#include "asterfort/pmfinfo.h"
#include "asterfort/pmfitg.h"
#include "asterfort/pmfitx.h"
#include "asterfort/poriro.h"
#include "asterfort/poutre_modloc.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utpslg.h"
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i, lmater, j, lorien, lmat, nno, nc, kpg, spt
    integer(kind=8) :: itype, irota, jacf
    real(kind=8) :: omega(3), omegl(3), s, xl
    real(kind=8) :: e, g, xnu, rho, a1, a2, xiy1, xiy2, xiz1, xiz2, alfay1, alfay2, alfaz1, alfaz2
    real(kind=8) :: a, xiy, xiz
    real(kind=8) :: pgl(3, 3), mlv(105), matp1(78)
    real(kind=8) :: carsec(6), rbid, casrho(6), casece(6)
    character(len=8) :: fami, poum
!
    integer(kind=8) :: nbres
    parameter(nbres=6)
    integer(kind=8) :: codres(nbres)
    real(kind=8) :: valres(nbres)
    character(len=16) :: nomres(nbres), elas_keyword
    integer(kind=8) :: elas_id
    data nomres/'E', 'NU', 'RHO', 'PROF_RHO_F_INT', 'PROF_RHO_F_EXT', 'COEF_MASS_AJOU'/
!
    integer(kind=8) :: nbfibr, nbgrfi, tygrfi, nbcarm, nug(10)
!
! --------------------------------------------------------------------------------------------------
    integer(kind=8), parameter :: nb_cara = 25
    real(kind=8) :: vale_cara(nb_cara)
    character(len=8) :: noms_cara(nb_cara)
    data noms_cara/'A1      ', 'IY1     ', 'IZ1     ', 'AY1     ', 'AZ1     ', &
        'EY1     ', 'EZ1     ', 'JX1     ', 'RY1     ', 'RZ1     ', &
        'RT1     ', 'JG1     ', 'A2      ', 'IY2     ', 'IZ2     ', &
        'AY2     ', 'AZ2     ', 'EY2     ', 'EZ2     ', 'JX2     ', &
        'RY2     ', 'RZ2     ', 'RT2     ', 'JG2     ', 'TVAR    '/
! --------------------------------------------------------------------------------------------------
!
!   Caracteristiques des elements
    nno = 2
    if (nomte .eq. 'MECA_POU_D_E' .or. nomte .eq. 'MECA_POU_D_T') then
        nc = 6
    else if (nomte .eq. 'MECA_POU_D_TG' .or. nomte .eq. 'MECA_POU_D_TGM') then
        nc = 7
    end if
!   Caracteristiques generales des sections
    call poutre_modloc('CAGNPO', noms_cara, nb_cara, lvaleur=vale_cara)
    a1 = vale_cara(1)
    xiy1 = vale_cara(2)
    xiz1 = vale_cara(3)
    alfay1 = vale_cara(4)
    alfaz1 = vale_cara(5)
    a2 = vale_cara(13)
    xiy2 = vale_cara(14)
    xiz2 = vale_cara(15)
    alfay2 = vale_cara(16)
    alfaz2 = vale_cara(17)
    itype = nint(vale_cara(25))
!   Recuperation des caracteristiques materiaux
    call jevech('PMATERC', 'L', lmater)
!
    if (nomte .eq. 'MECA_POU_D_E' .or. nomte .eq. 'MECA_POU_D_T' &
        .or. nomte .eq. 'MECA_POU_D_TG') then
        valres(:) = 0.0d0
        fami = 'FPG1'
        kpg = 1
        spt = 1
        poum = '+'
        call get_elas_id(zi(lmater), elas_id, elas_keyword)
        call rcvalb(fami, kpg, spt, poum, zi(lmater), ' ', elas_keyword, 0, ' ', [0.0d0], &
                    3, nomres, valres, codres, 1)
        e = valres(1)
        xnu = valres(2)
        rho = valres(3)
        g = e/(2.0d0*(1.0d0+xnu))
!
    else if (nomte .eq. 'MECA_POU_D_TGM') then
!       calcul de E et G
        call pmfitx(zi(lmater), 1, casece, g)
!       calcul de RHO MOYEN
        call pmfitx(zi(lmater), 2, casrho, rbid)
!       Récupération des caractéristiques des fibres
        call pmfinfo(nbfibr, nbgrfi, tygrfi, nbcarm, nug, jacf=jacf)
!
        call pmfitg(tygrfi, nbfibr, nbcarm, zr(jacf), carsec)
        a = carsec(1)
        xiy = carsec(5)
        xiz = carsec(4)
        rho = casrho(1)/a
        e = casece(1)/a
!
    end if
!   Coordonnees des noeuds
    xl = lonele()
!   Récupération des orientations
    call jevech('PCAORIE', 'L', lorien)
!   Récupération du vecteur rotation
    call jevech('PROTATR', 'L', irota)
    omega(1) = zr(irota+1)*zr(irota)
    omega(2) = zr(irota+2)*zr(irota)
    omega(3) = zr(irota+3)*zr(irota)
    call matrot(zr(lorien), pgl)
    do i = 1, 3
        s = 0.d0
        do j = 1, 3
            s = s+pgl(i, j)*omega(j)
        end do
        omegl(i) = s
    end do
!   Calcul de la matrice de raideur centrifuge locale
    matp1(:) = 0.0d0
!
    if (nomte .eq. 'MECA_POU_D_E' .or. nomte .eq. 'MECA_POU_D_T' &
        .or. nomte .eq. 'MECA_POU_D_TG') then
        call poriro(itype, matp1, rho, omegl, e, a1, a2, xl, xiy1, xiy2, &
                    xiz1, xiz2, g, alfay1, alfay2, alfaz1, alfaz2)
    else if (nomte .eq. 'MECA_POU_D_TGM') then
        call poriro(itype, matp1, rho, omegl, e, a, a, xl, xiy, xiy, &
                    xiz, xiz, g, alfay1, alfay2, alfaz1, alfaz2)
    end if
!
    call jevech('PMATUUR', 'E', lmat)
!
    if (nomte .eq. "MECA_POU_D_E" .or. nomte .eq. "MECA_POU_D_T") then
        call utpslg(nno, nc, pgl, matp1, zr(lmat))
!
    else if (nomte .eq. "MECA_POU_D_TG" .or. nomte .eq. "MECA_POU_D_TGM") then
        call masstg(matp1, mlv)
        call utpslg(nno, nc, pgl, mlv, zr(lmat))
    end if
!
!
end subroutine
