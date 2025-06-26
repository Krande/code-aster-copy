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
subroutine pcstru(n, in, ip, icpl, icpc, &
                  icpd, icpcx, icplx, niveau, complt, &
                  lca, imp, ier)
    implicit none
!
!----------------------------------------------------------------------
!   ENTREE
!   N          : TAILLE DE A
!   IN,IP      : MATRICE D'ENTREE FORMAT SYMETRIQUE
!   LCA        : LONGUEUR MAXI MATRICE FACTORISEE
!
!   SORTIE
!   ICPL,ICPC  : MATRICE APRES REMPLISSAGE FORMAT SYMETRIQUE
!   COMPLT    : FALSE OU TRUE
!   IER        : =0 TAILLE LCA SUFFISANTE
!              : =NN TAILLE LCA INSUFFISANTE IL FAUT NN
!
!   TRAVAIL
!   ICPD       : POINTEUR SUR LA DIAG DE LU
!   ICPCX      : IDEM ICPC
!   ICPLX      : IDEM ICPL
!----------------------------------------------------------------------
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/pcdiag.h"
#include "asterfort/pcfalu.h"
#include "asterfort/pcfull.h"
#include "asterfort/pcinfe.h"
!
    integer(kind=8) :: n, in(n)
    integer(kind=4) :: ip(*), icpc(*)
    integer(kind=8) :: icpl(0:n), icpd(n)
    integer(kind=8) :: icplx(0:n), icpcx(*)
!
    aster_logical :: complt
    integer(kind=8) :: i, ier, imp, k, k1, k2
    integer(kind=8) :: kk, lca, niv, niveau, nz
    integer(kind=8), pointer :: ind(:) => null()
!-----------------------------------------------------------------------
!
! IN-IP---> IPL-IPC
! =================
    AS_ALLOCATE(vi=ind, size=n)
    call pcfalu(n, in, ip, icpl, icpc, &
                ind, imp)
!
! INITIALISATION
    ier = 0
    complt = .false.
    call pcdiag(n, icpl, icpc, icpd)
!
! BOUCLE SUR LES NIVEAUX
! ======================
    do niv = 1, niveau
        nz = icpl(n)
        if (niv .lt. niveau) then
            call pcfull(n, icpl, icpc, icpd, icplx, &
                        icpcx, ind, lca, ier)
        else
            call pcinfe(n, icpl, icpc, icpd, icplx, &
                        icpcx, ind, lca, ier)
        end if
!
        if (ier .gt. 0) goto 50
!
        call pcdiag(n, icpl, icpc, icpd)
        if (icpl(n) .eq. nz) then
!         WRITE (6,4000) NIV
            complt = .true.
            goto 20
        end if
    end do
!
20  continue
!
! ICPL,ICPC FORMAT LU ---> FORMAT SYMETRIQUE
! ================================
    icpc(1) = 1
    kk = 1
    do i = 2, n
!                  ATTENTION ICPL(0:N)
        icpl(i-2) = kk
        k1 = icpl(i-1)+1
        k2 = icpd(i)
        do k = k1, k2
            kk = kk+1
            icpc(kk) = icpc(k)
        end do
!   TERME DIAG
        kk = kk+1
        icpc(kk) = int(i, 4)
    end do
    icpl(n-1) = kk
    goto 60
!
50  continue
60  continue
!
    AS_DEALLOCATE(vi=ind)
!
end subroutine
