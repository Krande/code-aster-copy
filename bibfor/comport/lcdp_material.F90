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

function lcdp_material(fami, kpg, ksp, imate, resi) result(mat)

    use lcdp_module, only: dp_material

    implicit none
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/rcvala.h"
#include "asterfort/rcvalb.h"

    aster_logical, intent(in) :: resi
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: imate, kpg, ksp
    type(dp_material) :: mat
! ----------------------------------------------------------------------
!  LOI DP - LECTURE DES PARAMETRES MATERIAUX
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
    integer(kind=8), parameter:: nbel = 2, nbdpli = 5, nbdpqua = 5, nbdpexp = 6
! ----------------------------------------------------------------------
    integer(kind=8) :: iok(nbel+nbdpexp)
    real(kind=8) :: valel(nbel), valdpli(nbdpli), valdpqua(nbdpqua), valdpexp(nbdpexp)
    real(kind=8) :: ltyped(1)
    character(len=16) :: nomel(nbel), nomdpli(nbdpli), nomdpqua(nbdpqua), nomdpexp(nbdpexp)
    character(len=1) :: poum
! ----------------------------------------------------------------------
    data nomel/'E', 'NU'/
    data nomdpli/'SY', 'H', 'P_ULTM', 'ALPHA', 'DILAT'/
    data nomdpqua/'SY', 'P_ULTM', 'SY_ULTM', 'ALPHA', 'DILAT'/
    data nomdpexp/'SY', 'P_C', 'SY_ULTM', 'DILAT_ULTM', 'ALPHA', 'DILAT'/
! ----------------------------------------------------------------------

    poum = merge('+', '-', resi)

!  Elasticity
    call rcvalb(fami, kpg, ksp, poum, imate, ' ', 'ELAS', 0, ' ', [0.d0], &
                nbel, nomel, valel, iok, 2)
    mat%lambda = valel(1)*valel(2)/((1+valel(2))*(1-2*valel(2)))
    mat%deuxmu = valel(1)/(1+valel(2))
    mat%troismu = 1.5d0*mat%deuxmu
    mat%troisk = valel(1)/(1.d0-2.d0*valel(2))
    mat%Eyoung = valel(1)

!  Drucker-Prager
!  Type d'ecrouissage
    call rcvala(imate, ' ', 'DRUCK_PRAGER', 0, ' ', [0.d0], 1, 'TYPE_DP', ltyped(1), iok, 1)
!  Lineaire
    if (nint(ltyped(1)) .eq. 1) then
        call rcvalb(fami, kpg, ksp, poum, imate, ' ', 'DRUCK_PRAGER', 0, ' ', [0.d0], &
                    nbdpli, nomdpli, valdpli, iok, 2)
        mat%sy = valdpli(1)
        mat%h = valdpli(2)
        mat%kau = valdpli(3)
        mat%a = valdpli(4)
        mat%b0 = valdpli(5)
        mat%troisa = 3*mat%a
        mat%syultm = 0.d0
        mat%bultm = 0.d0
        mat%kac = 0.d0
        mat%type_dp = nint(ltyped(1))
!  Parabolique
    else if (nint(ltyped(1)) .eq. 2) then
        call rcvalb(fami, kpg, ksp, poum, imate, ' ', 'DRUCK_PRAGER', 0, ' ', [0.d0], &
                    nbdpqua, nomdpqua, valdpqua, iok, 2)
        mat%sy = valdpqua(1)
        mat%h = 0.d0
        mat%kau = valdpqua(2)
        mat%a = valdpqua(4)
        mat%b0 = valdpqua(5)
        mat%troisa = 3*mat%a
        mat%syultm = valdpqua(3)
        mat%bultm = 0.d0
        mat%kac = 0.d0
        mat%type_dp = nint(ltyped(1))
! Exponentiel
    else if (nint(ltyped(1)) .eq. 3) then
        call rcvalb(fami, kpg, ksp, poum, imate, ' ', 'DRUCK_PRAGER', 0, ' ', [0.d0], &
                    nbdpexp, nomdpexp, valdpexp, iok, 2)
        mat%sy = valdpexp(1)
        mat%h = 0.d0
        mat%kau = 0.d0
        mat%a = valdpexp(5)
        mat%b0 = valdpexp(6)
        mat%troisa = 3*mat%a
        mat%syultm = valdpexp(3)
        mat%bultm = valdpexp(4)
        mat%kac = valdpexp(2)
        mat%type_dp = nint(ltyped(1))
    else
        ASSERT(.false.)
    end if

end function
