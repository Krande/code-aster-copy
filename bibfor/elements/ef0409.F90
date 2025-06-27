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
subroutine ef0409()
!     CALCUL DE EFGE_ELNO EN NON LINEAIRE
!     ------------------------------------------------------------------
    implicit none
!
#include "jeveux.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/ppgan2.h"
#include "asterfort/r8inir.h"
#include "asterfort/tecach.h"
#include "blas/dcopy.h"
! person_in_charge: sebastien.fayolle at edf.fr
!
    integer(kind=8) :: nnos, ipoids, ivf, idfdx, jgano, jtab(7), ipg, npg
    integer(kind=8) :: ichn, icontm, nno, igeom, ndim, iret
!
!     ---> POUR DKT MATELEM = 3 * 6 DDL = 171 TERMES STOCKAGE SYME
!     ---> POUR DKQ MATELEM = 4 * 6 DDL = 300 TERMES STOCKAGE SYME
!
!     ---> POUR DKT EFFINT = 24
!     ---> POUR DKQ EFFINT = 32
    real(kind=8) :: effint(32)
    blas_int :: b_incx, b_incy, b_n
!
! ---   RECUPERATION DES ADRESSES DANS ZR DES POIDS DES PG
!       DES FONCTIONS DE FORME DES VALEURS DES DERIVEES DES FONCTIONS
!       DE FORME ET DE LA MATRICE DE PASSAGE GAUSS -> NOEUDS
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfdx, jgano=jgano)
    call jevech('PEFFORR', 'E', ichn)
!
    call jevech('PGEOMER', 'L', igeom)
    call tecach('OOO', 'PCONTRR', 'L', iret, nval=7, &
                itab=jtab)
    call r8inir(32, 0.d0, effint, 1)
!
    do ipg = 1, npg
        icontm = jtab(1)+8*(ipg-1)
        b_n = to_blas_int(8)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(icontm), b_incx, effint(8*(ipg-1)+1), b_incy)
    end do
!
    call ppgan2(jgano, 1, 8, effint, zr(ichn))
!
end subroutine
