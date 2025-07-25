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
subroutine vpnor1(norm, neq, nbmode, ddlexc, vecpro, &
                  isign, numddl, coef)
    implicit none
#include "asterfort/utmess.h"
    integer(kind=8) :: nbmode, neq, ddlexc(*), isign, numddl
    real(kind=8) :: vecpro(neq, *), coef(*)
    character(len=*) :: norm
    character(len=24) :: valk
!     NORMALISATION DE VECTEURS D'UN CONCEPT MODE_FLAMB
!     ------------------------------------------------------------------
! IN  NORM   : TYPE DE NORMALISATION
!          = 'AVEC_CMP'
!          = 'EUCL', 'EUCL_TRAN', ...
! IN  NEQ    : NOMBRE D'EQUATIONS
! IN  NBMODE : NOMBRE DE MODES
! IN  DDLEXC : TABLEAU DES DDL EXCLUS
!              = 0 SI EXCLUS
!              = 1 SI NON EXCLUS
! VAR VECPRO : TABLEAU DES VECTEURS PROPRES
! OUT COEF   : COEFFICIENTS
!     ------------------------------------------------------------------
    integer(kind=8) :: im, ie
    real(kind=8) :: xnorm, xx1
! DEB ------------------------------------------------------------------
!
    if (norm .eq. 'AVEC_CMP' .or. norm(1:4) .eq. 'EUCL') then
!
!     --- NORMALISATION SUR LES DDL NON EXCLUS
!
        do im = 1, nbmode
            xnorm = 0.0d0
            if (norm(1:4) .eq. 'EUCL') then
                do ie = 1, neq
                    xx1 = vecpro(ie, im)*ddlexc(ie)
                    xnorm = xnorm+xx1*xx1
                end do
                xnorm = sqrt(xnorm)
            else
                do ie = 1, neq
                    xx1 = vecpro(ie, im)*ddlexc(ie)
                    if (abs(xnorm) .lt. abs(xx1)) then
                        xnorm = xx1
                    end if
                end do
            end if
            xx1 = 1.0d0/xnorm
            coef(im) = xx1
            do ie = 1, neq
                vecpro(ie, im) = vecpro(ie, im)*xx1
            end do
        end do
!
    else
!
        valk = norm
        call utmess('F', 'ALGELINE4_77', sk=valk)
!
    end if
!
    if (isign .eq. 0) then
    else if (isign .eq. 1) then
        do im = 1, nbmode
            xx1 = vecpro(numddl, im)
            if (xx1 .lt. 0.0d0) then
                coef(im) = -coef(im)
                do ie = 1, neq
                    vecpro(ie, im) = -vecpro(ie, im)
                end do
            end if
        end do
    else if (isign .eq. -1) then
        do im = 1, nbmode
            xx1 = vecpro(numddl, im)
            if (xx1 .gt. 0.0d0) then
                coef(im) = -coef(im)
                do ie = 1, neq
                    vecpro(ie, im) = -vecpro(ie, im)
                end do
            end if
        end do
    end if
!
end subroutine
