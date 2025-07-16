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
subroutine dlfext(nbVectAsse, nbLoad, temps, neq, liad, &
                  lifo, loadNameJv, loadInfoJv, loadFuncJv, model, &
                  materField, mateco, caraElem, numedd, f)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterc/r8vide.h"
#include "asterfort/asasve.h"
#include "asterfort/ascova.h"
#include "asterfort/fext.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/vechme.h"
#include "asterfort/vedime.h"
#include "blas/daxpy.h"
!
    integer(kind=8), intent(in) :: nbVectAsse, nbLoad, neq, liad(*)
    real(kind=8), intent(in) :: temps
    character(len=24), intent(in) :: lifo(*), loadInfoJv, loadFuncJv
    character(len=24), intent(in) :: model, caraElem, loadNameJv, materField, mateco, numedd
    real(kind=8), intent(out) :: f(*)
!
! --------------------------------------------------------------------------------------------------
!
!  CALCUL DU SECOND MEMBRE F* A PARTIR DE :
!      - VECT_ASSE
!      - CHARGE
!
! --------------------------------------------------------------------------------------------------
!
!  INPUT:
!        NVECA    : NOMBRE D'OCCURENCES DU MOT CLE VECT_ASSE
!        NCHAR    : NOMBRE D'OCCURENCES DU MOT CLE CHARGE
!        TEMPS    : INSTANT DE CALCUL
!        NEQ      : NOMBRE D'EQUATIONS (D.D.L. ACTIFS)
!        LIAD     : LISTE DES ADRESSES DES VECTEURS CHARGEMENT (NVECT)
!        LIFO     : LISTE DES NOMS DES FONCTIONS EVOLUTION (NVECT)
!        CHARGE   : LISTE DES CHARGES
!        INFOCH   : INFO SUR LES CHARGES
!        FOMULT   : LISTE DES FONC_MULT ASSOCIES A DES CHARGES
!        MODELE   : NOM DU MODELE
!        MATE     : NOM DU CHAMP DE MATERIAU
!        CARELE   : CARACTERISTIQUES DES POUTRES ET COQUES
!        NUMEDD   : NUME_DDL DE LA MATR_ASSE RIGID
!
!  OUTPUT:
!        F        : VECTEUR FORCE EXTERIEURE (NEQ)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=4), parameter :: typmat = "R", para = 'INST'
    integer(kind=8) :: iret, ieq, n1
    real(kind=8) :: partps(3)
    character(len=16) :: method
    character(len=24) :: vechmp, vachmp, cnchmp
    real(kind=8), pointer :: f1(:) => null()
    real(kind=8), pointer :: f2(:) => null()
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    partps(1:3) = [temps, r8vide(), r8vide()]
    vechmp = ' '
    vachmp = ' '
    cnchmp = ' '
    call getvtx('SCHEMA_TEMPS', 'SCHEMA', iocc=1, scal=method, nbret=n1)

! - Load from VECT_ASSE
    f(1:neq) = 0.d0
    if (nbVectAsse .ne. 0) then
        call fext(temps, neq, nbVectAsse, liad, lifo, f)
    end if
!
    if (nbLoad .ne. 0) then
! ----- Neumann loads
        call vechme('S', &
                    model, caraElem, materField, mateco, &
                    loadNameJv, loadInfoJv, &
                    partps, &
                    vechmp)
        call asasve(vechmp, numedd, typmat, vachmp)
        call ascova('D', vachmp, loadFuncJv, 'INST', temps, &
                    typmat, cnchmp)
        call jeveuo(cnchmp(1:19)//'.VALE', 'L', vr=f1)
!
        b_n = to_blas_int(neq)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0, f1, b_incx, f, b_incy)
!
! 2.2.2. ==> -- LES DIRICHLETS
!
        call vedime(model, loadNameJv, loadInfoJv, temps, typmat, &
                    vechmp)
        call asasve(vechmp, numedd, typmat, vachmp)
        call ascova('D', vachmp, loadFuncJv, para, temps, &
                    typmat, cnchmp)
        call jeveuo(cnchmp(1:19)//'.VALE', 'L', vr=f2)
!
! -- TEST DE PRESENCE DE CHARGEMENT DIRICHLET (DEPL IMPOSE NON NUL)
        iret = 0
        do ieq = 1, neq
            if (abs(f2(ieq)) .gt. r8prem()) iret = 1
        end do
        if ((iret .eq. 1) .and. (method .ne. 'NEWMARK')) then
            call utmess('F', 'DYNALINE1_20')
        end if
!
        b_n = to_blas_int(neq)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0, f2, b_incx, f, b_incy)
    end if
!
    call jedetr(cnchmp)
!
    call jedema()
!
end subroutine
