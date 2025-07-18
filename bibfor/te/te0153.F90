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

subroutine te0153(option, nomte)
!
!
    implicit none
    character(len=*) :: option, nomte
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "asterfort/lonele.h"
#include "asterfort/matrot.h"
#include "asterfort/pmavec.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
#include "asterfort/utpslg.h"
#include "asterfort/vecma.h"
!
! --------------------------------------------------------------------------------------------------
!
!     CALCULE LES MATRICES ELEMENTAIRES DES ELEMENTS DE BARRE
!
! --------------------------------------------------------------------------------------------------
!
! IN  OPTION : K16 : NOM DE L'OPTION A CALCULER
!        'RIGI_MECA'      : CALCUL DE LA MATRICE DE RAIDEUR
!        'MASS_MECA'      : CALCUL DE LA MATRICE DE MASSE
! IN  NOMTE  : K16 : NOM DU TYPE ELEMENT
!        'MECA_BARRE' : ELEMENT BARRE
!        'MECA_2D_BARRE' : ELEMENT BARRE
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: codres(1)
    integer(kind=8) :: i, iacce, imate, lmat, lorien, lsect
    integer(kind=8) :: lvec, nc, nno
    real(kind=8) :: e(1), rho(1), pgl(3, 3), mat(21), matr(21)
    real(kind=8) :: a, xl, xrig, xmas, matp(6, 6), mat2dm(4, 4), mat2dv(10)
    real(kind=8) :: r8b
    character(len=16) :: ch16
!
! --------------------------------------------------------------------------------------------------
!
    r8b = 0.d0
    call jevech('PCAGNBA', 'L', lsect)
    a = zr(lsect)
    nno = 2
    nc = 3
!
!   Longueur de l'élément
    if (nomte .eq. 'MECA_BARRE') then
        xl = lonele()
    else if (nomte .eq. 'MECA_2D_BARRE') then
        xl = lonele(dime=2)
    else
        xl = 0.0d0
        ASSERT(ASTER_FALSE)
    end if
!
!   recuperation des orientations alpha,beta,gamma
    call jevech('PCAORIE', 'L', lorien)
    if (option .eq. 'M_GAMMA') then
        call jevech('PVECTUR', 'E', lvec)
        call jevech('PACCELR', 'L', iacce)
    else
        call jevech('PMATUUR', 'E', lmat)
    end if
!
!   calcul des matrices elementaires
    mat(:) = 0.d0
    call jevech('PMATERC', 'L', imate)
    if (option .eq. 'RIGI_MECA') then
        call rcvalb('FPG1', 1, 1, '+', zi(imate), ' ', 'ELAS', 0, ' ', [r8b], &
                    1, 'E', e, codres, 1)
        xrig = e(1)*a/xl
        mat(1) = xrig
        mat(7) = -xrig
        mat(10) = xrig
!
    else if (option .eq. 'MASS_MECA' .or. option .eq. 'M_GAMMA') then
        call rcvalb('FPG1', 1, 1, '+', zi(imate), ' ', 'ELAS', 0, ' ', [r8b], &
                    1, 'RHO', rho, codres, 1)
        matr(:) = 0.d0
!
        xmas = rho(1)*a*xl/6.d0
        mat(1) = xmas*2.d0
        mat(3) = xmas*2.d0
        mat(6) = xmas*2.d0
        mat(10) = xmas*2.d0
        mat(15) = xmas*2.d0
        mat(21) = xmas*2.d0
!
        mat(7) = xmas
        mat(12) = xmas
        mat(18) = xmas
!
    else if ((option .eq. 'MASS_MECA_DIAG') .or. (option .eq. 'MASS_MECA_EXPLI')) then
        call rcvalb('FPG1', 1, 1, '+', zi(imate), ' ', 'ELAS', 0, ' ', [r8b], &
                    1, 'RHO', rho, codres, 1)
        xmas = rho(1)*a*xl/2.d0
        mat(1) = xmas
        mat(3) = xmas
        mat(6) = xmas
        mat(10) = xmas
        mat(15) = xmas
        mat(21) = xmas
!
    else
        ch16 = option
        call utmess('F', 'ELEMENTS2_47', sk=ch16)
    end if
!
!   passage du repere local au repere global
    call matrot(zr(lorien), pgl)
    call utpslg(nno, nc, pgl, mat, matr)
!
    if (option .eq. 'M_GAMMA') then
        if (nomte .eq. 'MECA_BARRE') then
            matp(:, :) = 0.d+0
            call vecma(matr, 21, matp, 6)
            call pmavec('ZERO', 6, matp, zr(iacce), zr(lvec))
        else
            mat2dv(1) = matr(1)
            mat2dv(2) = matr(2)
            mat2dv(3) = matr(3)
            mat2dv(4) = matr(7)
            mat2dv(5) = matr(8)
            mat2dv(6) = matr(10)
            mat2dv(7) = matr(11)
            mat2dv(8) = matr(12)
            mat2dv(9) = matr(14)
            mat2dv(10) = matr(15)
            mat2dm(:, :) = 0.d+0
            call vecma(mat2dv, 10, mat2dm, 4)
            call pmavec('ZERO', 4, mat2dm, zr(iacce), zr(lvec))
        end if
    else
!       ecriture dans le vecteur pmattur suivant l'element
        if (nomte .eq. 'MECA_BARRE') then
            do i = 1, 21
                zr(lmat+i-1) = matr(i)
            end do
        else if (nomte .eq. 'MECA_2D_BARRE') then
            zr(lmat) = matr(1)
            zr(lmat+1) = matr(2)
            zr(lmat+2) = matr(3)
            zr(lmat+3) = matr(7)
            zr(lmat+4) = matr(8)
            zr(lmat+5) = matr(10)
            zr(lmat+6) = matr(11)
            zr(lmat+7) = matr(12)
            zr(lmat+8) = matr(14)
            zr(lmat+9) = matr(15)
        end if
    end if
!
end subroutine
