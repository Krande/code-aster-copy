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

subroutine zerofr_param2(intini, algo, func, funcp, para, nb_para, x1, x2, &
                         tol, itmax, solu, iret, iter)
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/encadr.h"
#include "asterfort/zerof2.h"
#include "asterfort/zerofb.h"
#include "asterfort/zerofc.h"
#include "asterfort/zerofo.h"
    integer(kind=8) :: intini, itmax, iter, iret
    character(len=*) :: algo
    integer(kind=8), intent(in) :: nb_para
    real(kind=8) :: solu, tol, x1, x2, para(nb_para)
    interface
        function funcp(x, param)
            real(kind=8), intent(in) :: x
            real(kind=8), intent(in) :: param(*)
            real(kind=8) :: funcp
        end function
        function func(x)
            real(kind=8) :: x
            real(kind=8) :: func
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
    character(len=8) :: algoz
    real(kind=8) :: a, b, fa, fb
!
    algoz = algo
    a = x1
    b = x2
!
!     ------------------------------------------------------------------
!     RECHERCHE DE L'INTERVALLE INITIAL [A,B]
!     ------------------------------------------------------------------
!
    ASSERT(intini .eq. 0 .or. intini .eq. 1 .or. intini .eq. 2)
!
    if (intini .eq. 1) then
!
!       BRACKETING CROISSANT A GAUCHE ET A DROITE
        call encadr(func, funcp, para, nb_para, a, b, fa, fb, &
                    itmax, 1.6d0, iret)
        if (iret .ne. 0) goto 999
!
    else if (intini .eq. 2) then
!
!       BRACKETING CROISSANT UNIQUEMNT A DROITE :
!       SOUVENT LE CAS POUR LES LOIS DE COMPORTEMENT
!       (SI F EST CROISSANTE ET F(A)<0, OU L'INVERSE),
!       CE QUI PERMET DE PRENDRE UN COEF MULT GRAND (10)
        call encadr(func, funcp, para, nb_para, a, b, fa, fb, &
                    itmax, 10.d0, iret)
        if (iret .ne. 0) goto 999
!
    end if
!
!     ------------------------------------------------------------------
!     RECHERCHE DU ZERO DE F ENTRE X1 ET X2
!     ------------------------------------------------------------------
!
    if (algoz .eq. 'AUTO') algoz = 'BRENT'
!
!
    if (algoz .eq. 'BRENT') then
!
        call zerofb(func, funcp, para, nb_para, a, b, tol, itmax, &
                    solu, iret, iter)
!
    else if (algoz .eq. 'SECANTE') then
!
        call zerofc(func, funcp, para, nb_para, a, b, tol, itmax, &
                    solu, iret, iter)
!
    else if (algoz .eq. 'DEKKER') then
!
        call zerofo(func, funcp, para, nb_para, a, b, tol, itmax, &
                    solu, iret, iter)
!
    else if (algoz .eq. 'DEKKER2') then
!
        call zerof2(func, funcp, para, nb_para, a, b, tol, itmax, &
                    solu, iret, iter)
!
    else
!
        ASSERT(.false.)
!
    end if
!
!
999 continue
end subroutine
