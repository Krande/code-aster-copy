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

subroutine te0240(option, nomte)
!
!
! --------------------------------------------------------------------------------------------------
!
!     CALCULE LA MATRICE DE RIGIDITE ELEMENTAIRE DES ELEMENTS DE TUYAU
!     DROIT DE TIMOSHENKO
!
! --------------------------------------------------------------------------------------------------
!
! IN  OPTION : K16 : NOM DE L'OPTION A CALCULER
!        'RIGI_MECA '   : CALCUL DE LA MATRICE DE RIGIDITE
!
! IN  NOMTE  : K16 : NOM DU TYPE ELEMENT
!        'MEFS_POU_D_T' : POUTRE DROITE DE TIMOSHENKO
!                         (SECTION CONSTANTE OU NON)
!
! --------------------------------------------------------------------------------------------------
!
    implicit none
#include "jeveux.h"
#include "asterfort/jevech.h"
#include "asterfort/lonele.h"
#include "asterfort/matrot.h"
#include "asterfort/ptktfv.h"
#include "asterfort/ptktuf.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
#include "asterfort/utpslg.h"
    character(len=*) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: imate, itype, lmat, lorien, lsect, lsect2
    integer(kind=8) :: kpg, spt
    integer(kind=8) :: nbpar, nc, nno
!
    real(kind=8) :: rho, valpar
    real(kind=8) :: pgl(3, 3), mat(136)
    real(kind=8) :: e, nu, g, celer
    real(kind=8) :: a, ai, xiy, xiz, alfay, alfaz, xjx, xl
    real(kind=8) :: a2, ai2, xiy2, xiz2, alfay2, alfaz2, xjx2, ey, ez
!
    character(len=8) :: nompar, fami, poum
    character(len=16) :: ch16
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nbres
    parameter(nbres=2)
    integer(kind=8) :: codres(nbres)
    real(kind=8) :: valres(nbres)
    character(len=16) ::  nomres(nbres)
!
! --------------------------------------------------------------------------------------------------
!
!   recuperation des caracteristiques materiaux
    call jevech('PMATERC', 'L', imate)
    nbpar = 0
    nompar = '  '
    valpar = 0.0d0
    fami = 'FPG1'
    kpg = 1
    spt = 1
    poum = '+'
    valres(1:nbres) = 0.0d0
    nomres(1) = 'E'
    nomres(2) = 'NU'
!
    call rcvalb(fami, kpg, spt, poum, zi(imate), ' ', 'ELAS', nbpar, nompar, [valpar], &
                nbres, nomres, valres, codres, 1)
    e = valres(1)
    nu = valres(2)
    g = e/(2.0d0*(1.0d0+nu))
!
    valres(1) = 0.0d0
    valres(2) = 0.0d0
    nomres(1) = 'RHO'
    nomres(2) = 'CELE_R'
!
    call rcvalb(fami, kpg, spt, poum, zi(imate), ' ', 'FLUIDE', nbpar, nompar, [valpar], &
                nbres, nomres, valres, codres, 1)
    rho = valres(1)
    celer = valres(2)
!
!   recuperation des caracteristiques generales des sections
    call jevech('PCAGNPO', 'L', lsect)
    lsect = lsect-1
    itype = nint(zr(lsect+25))
    nno = 2
    nc = 8
!
!     --- SECTION INITIALE ---
    a = zr(lsect+1)
    xiy = zr(lsect+2)
    xiz = zr(lsect+3)
    alfay = zr(lsect+4)
    alfaz = zr(lsect+5)
    xjx = zr(lsect+8)
    ai = zr(lsect+12)
!
!   SECTION FINALE
    lsect2 = lsect+12
    a2 = zr(lsect2+1)
    xiy2 = zr(lsect2+2)
    xiz2 = zr(lsect2+3)
    alfay2 = zr(lsect2+4)
    alfaz2 = zr(lsect2+5)
    xjx2 = zr(lsect2+8)
    ai2 = zr(lsect2+12)
    ey = -(zr(lsect+6)+zr(lsect2+6))/2.0d0
    ez = -(zr(lsect+7)+zr(lsect2+7))/2.0d0
!
    if (nomte .ne. 'MEFS_POU_D_T') then
        ch16 = nomte
        call utmess('F', 'ELEMENTS2_42', sk=ch16)
    end if
!
!   Longueur de l'élément
    xl = lonele()
!   calcul des matrices elementaires
    mat(:) = 0.d0
    if (option .eq. 'RIGI_MECA') then
        if (itype .eq. 0) then
!           poutre droite a section constante
            call ptktuf(mat, e, rho, celer, a, &
                        ai, xl, xiy, xiz, xjx, &
                        g, alfay, alfaz, ey, ez)
        else if (itype .eq. 1 .or. itype .eq. 2) then
!           poutre droite a section variable (type 1 ou 2)
            call ptktfv(itype, mat, e, rho, celer, &
                        a, ai, a2, ai2, xl, &
                        xiy, xiy2, xiz, xiz2, xjx, &
                        xjx2, g, alfay, alfay2, alfaz, &
                        alfaz2, ey, ez)
        end if
!       recuperation des orientations alpha,beta,gamma
        call jevech('PCAORIE', 'L', lorien)
!       passage du repere local au reper global
        call matrot(zr(lorien), pgl)
        call jevech('PMATUUR', 'E', lmat)
        call utpslg(nno, nc, pgl, mat, zr(lmat))
    else
        ch16 = option
        call utmess('F', 'ELEMENTS2_47', sk=ch16)
    end if
!
end subroutine
