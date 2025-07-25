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

subroutine porigy(nomte, rho, xnu, icdmat, klv, nl)
!
!
!
! --------------------------------------------------------------------------------------------------
!
!     CALCULE LA MATRICE GYROSCOPIQUE DES ELEMENTS DE POUTRE
!
! IN  NOMTE : NOM DU TYPE ELEMENT
!             'MECA_POU_D_E'  'MECA_POU_D_T'  'MECA_POU_D_TG'
!             'MECA_POU_D_EM' 'MECA_POU_D_TGM'
! --------------------------------------------------------------------------------------------------
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/lonele.h"
#include "asterfort/pmfitx.h"
#include "asterfort/poutre_modloc.h"
#include "asterfort/ptgy02.h"
#include "asterfort/utmess.h"
#include "asterfort/lteatt.h"
!
    integer(kind=8) :: icdmat
    character(len=*) :: nomte
    real(kind=8) :: rho, xnu, klv(*)
    integer(kind=8) :: nl
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: istruc, itype
    real(kind=8) :: rbid, casect(6)
    real(kind=8) :: ey, ez, xl
    real(kind=8) :: a, xiy, xiz, xjx, alfay, alfaz, alfinv
    real(kind=8) :: a2, xiy2, xiz2, xjx2, alfay2, alfaz2
    character(len=16) :: ch16
    aster_logical :: euler
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nb_cara = 17
    real(kind=8) :: vale_cara(nb_cara)
    character(len=8) :: noms_cara(nb_cara)
    data noms_cara/'A1', 'IY1', 'IZ1', 'AY1', 'AZ1', 'EY1', 'EZ1', 'JX1', &
        'A2', 'IY2', 'IZ2', 'AY2', 'AZ2', 'EY2', 'EZ2', 'JX2', 'TVAR'/
!
! --------------------------------------------------------------------------------------------------
!
    euler = lteatt('EULER', 'OUI')
!
!   recuperation des caracteristiques generales des sections
    xl = lonele()
    call poutre_modloc('CAGNPO', noms_cara, nb_cara, lvaleur=vale_cara)
!
    a = vale_cara(1)
    xiy = vale_cara(2)
    xiz = vale_cara(3)
    alfay = vale_cara(4)
    alfaz = vale_cara(5)
    xjx = vale_cara(8)
    a2 = vale_cara(9)
    xiy2 = vale_cara(10)
    xiz2 = vale_cara(11)
    alfay2 = vale_cara(12)
    alfaz2 = vale_cara(13)
    xjx2 = vale_cara(16)
    ey = (vale_cara(6)+vale_cara(14))/2.d0
    ez = (vale_cara(7)+vale_cara(15))/2.d0
    itype = nint(vale_cara(17))
!
    if (nomte .eq. 'MECA_POU_D_E') then
!       poutre droite d'euler a 6 ddl
        istruc = 1
        alfinv = 0.0d0
    else if (nomte .eq. 'MECA_POU_D_T') then
!       poutre droite de timoskenko a 6 ddl
        istruc = 1
        alfinv = 2.0d0/(alfay+alfaz)
    else if (nomte .eq. 'MECA_POU_D_EM' .or. nomte .eq. 'MECA_POU_D_TGM') then
!       poutre droite multifibre
        itype = 0
        istruc = 1
        alfinv = 0.0d0
!       on met rho=1
        rho = 1.d0
        call pmfitx(icdmat, 2, casect, rbid)
        a = casect(1)
        xiy = casect(5)
        xiz = casect(4)
        xjx = 0.d0
    else
        ch16 = nomte
        call utmess('F', 'ELEMENTS2_42', sk=ch16)
    end if
!
!
    if (itype .eq. 1 .or. itype .eq. 2) then
!       moyennage
        a = (a+a2)/2.0d0
        xiy = (xiy+xiy2)/2.0d0
        xiz = (xiz+xiz2)/2.0d0
        alfay = (alfay+alfay2)/2.0d0
        alfaz = (alfaz+alfaz2)/2.0d0
        xjx = (xjx+xjx2)/2.0d0
        if (euler) then
            alfinv = 0.0d0
        else
            alfinv = 2.0d0/(alfay+alfaz)
        end if
    end if
    call ptgy02(klv, nl, xnu, rho, a, xl, xiy, xiz, xjx, alfinv, ey, ez, istruc)
!
end subroutine
