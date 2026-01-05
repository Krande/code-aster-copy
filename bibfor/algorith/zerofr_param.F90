! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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

subroutine zerofr_param(intini, algo, funcp, para, nb_para, x1, x2, &
                        tol, itmax, solu, iret, iter)
    implicit none
!
#include "asterfort/zerofr_param2.h"
    integer(kind=8) :: intini, itmax, iter, iret, nb_para
    character(len=*) :: algo
    real(kind=8) :: solu, tol, x1, x2, para(nb_para)
    interface
        function funcp(x, param)
            real(kind=8), intent(in) :: x
            real(kind=8), intent(in) :: param(*)
            real(kind=8) :: funcp
        end function
    end interface
!
! ----------------------------------------------------------------------
!     BUT : TROUVER LE ZERO D'UNE FONCTION SCALAIRE REELLE
!
!
! IN  INTINI : RECHERCHE DE L'INTERVALLE INITIAL
!              TEL QUE F CHANGE DE SIGNE
!              = 0 : PAS DE RECHERCHE
!              = 1 : RECHERCHE PAR BRACKETING CROISSANT
!              = 2 : RECHERCHE PAR BRACKETING CROISSANT A DROITE
! IN  ALGO   : ALGORITHME DE RECHERCHE DU ZERO : 'AUTO', 'SECANTE',
!              'DEKKER', 'DEKKER2', 'BRENT'
!              SI ALGO VAUT 'AUTO', ON PREND 'BRENT'
! IN  funcp   : FONCTION F
! IN  PARA   : Paramètre supplémentaire à passer à l'appel de funcp
! IN  X1, X2 : INTERVELLE DE RECHERCHE
! IN  TOL    : PRECISION ABSOLUE : LA SOLUTION X EST TELLE QUE F(X)<TOL
! IN  ITMAX  : NOMBRE D'ITERATIONS MAXIMUM
! OUT SOLU   : ZERO DE F
! OUT IRET   : CODE RETOUR : IRET = 0 : OK
!            :               SINON    : PROBLEME
! OUT ITER   : NOMBRE D'ITERATIONS EFFECTUEES
! ----------------------------------------------------------------------
!
    call zerofr_param2(intini, algo, func, funcp, para, nb_para, &
                       x1, x2, tol, itmax, solu, iret, iter)
!
contains
!
    function func(x)
        real(kind=8) :: x
        real(kind=8) :: func, x0
!
        x0 = x
        func = 0.0
    end function
end subroutine
