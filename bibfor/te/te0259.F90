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

subroutine te0259(option, nomte)
! aslint: disable=
    implicit none
#include "jeveux.h"
#include "asterfort/jevech.h"
#include "asterfort/masstg.h"
#include "asterfort/matrot.h"
#include "asterfort/pogyro.h"
#include "asterfort/rcvalb.h"
#include "asterfort/upletr.h"
#include "asterfort/utmess.h"
#include "asterfort/utpalg.h"
!
    character(len=16) :: option, nomte
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
! IN  OPTION : K16 : NOM DE L'OPTION A CALCULER
!       'MECA_GYRO': CALCUL DE LA MATRICE D'AMORTISSEMENT GYROSCOPIQUE
! IN  NOMTE  : K16 : NOM DU TYPE ELEMENT
!       'MECA_POU_D_E'  : POUTRE DROITE D'EULER       (SECTION VARIABLE)
!       'MECA_POU_D_T'  : POUTRE DROITE DE TIMOSHENKO (SECTION VARIABLE)
!       'MECA_POU_D_EM' : POUTRE DROITE MULTIFIBRE D EULER (SECT. CONST)
!       'MECA_POU_D_TG' : POUTRE DROITE DE TIMOSHENKO (GAUCHISSEMENT)
!       'MECA_POU_D_TGM': POUTRE DROITE DE TIMOSHENKO (GAUCHISSEMENT)
!                         MULTI-FIBRES (SECTION CONSTANTE)
!
!
    integer(kind=8) :: nddl, nbres
    parameter(nbres=3)
    real(kind=8) :: valres(nbres)
    integer(kind=8) :: codres(nbres)
    character(len=8) :: nompar, fami, poum
    character(len=16) :: nomres(nbres)
    real(kind=8) :: pgl(3, 3), klv(78), klw(78), mlv(105)
    real(kind=8) :: e, rho
    real(kind=8) :: valpar, xnu, zero
    integer(kind=8) :: imate, lmat, lorien
    integer(kind=8) :: nbpar, nc, nno, kpg, spt
!     ------------------------------------------------------------------
    data nomres/'E', 'RHO', 'NU'/
!     ------------------------------------------------------------------
    zero = 0.d0
!     ------------------------------------------------------------------
!
!     --- CARACTERISTIQUES DES ELEMENTS
!
    if (nomte .eq. 'MECA_POU_D_E' .or. nomte .eq. 'MECA_POU_D_T' .or. nomte .eq. &
        'MECA_POU_D_EM') then
        nno = 2
        nc = 6
    elseif (nomte .eq. 'MECA_POU_D_TG' .or. nomte .eq. 'MECA_POU_D_TGM' &
            ) then
        nno = 2
        nc = 7
    else
        call utmess('F', 'ELEMENTS2_42', sk=nomte)
    end if
!
    nddl = nc*nno
!
!     --- RECUPERATION DES CARACTERISTIQUES MATERIAUX ---
!
    call jevech('PMATERC', 'L', imate)
!
    nbpar = 0
    nompar = ' '
    valpar = zero
    fami = 'FPG1'
    kpg = 1
    spt = 1
    poum = '+'
    call rcvalb(fami, kpg, spt, poum, zi(imate), &
                ' ', 'ELAS', nbpar, nompar, [valpar], &
                nbres, nomres, valres, codres, 1)
    e = valres(1)
    rho = valres(2)
    xnu = valres(3)
!
    call jevech('PMATUNS', 'E', lmat)
!
!     --- RECUPERATION DES ORIENTATIONS ---
!
    call jevech('PCAORIE', 'L', lorien)
!
!     --- CALCUL DE LA MATRICE GYROSCOPIQUE LOCALE ---
    call pogyro(nomte, rho, xnu, zi(imate), klv, &
                78)
!
    call matrot(zr(lorien), pgl)
!  CHANGEMENT DE BASE : LOCAL -> GLOBAL
    call utpalg(nno, 6, pgl, klv, klw)
!
    if (nomte .eq. 'MECA_POU_D_TG' .or. nomte .eq. 'MECA_POU_D_TGM') then
        call masstg(klw, mlv)
    end if
!
! CONSITUER UNE MATRICE PLEINE A PARTIR DE LA TRIANGULAIRE SUPERIEURE
!
    if (nomte .eq. 'MECA_POU_D_TG' .or. nomte .eq. 'MECA_POU_D_TGM') then
        call upletr(nddl, zr(lmat), mlv)
    else
        call upletr(nddl, zr(lmat), klw)
    end if
!
end subroutine
