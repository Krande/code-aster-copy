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

subroutine ndcent(igeom, ndim, lsn, nfiss, tx, txlsn, nnc)
    implicit none
!
#include "jeveux.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/elrfvf.h"
#include "asterfort/reerel.h"
    integer(kind=8) :: igeom, nnc, ndim, nfiss
    real(kind=8) :: tx(3, 7), lsn(*), txlsn(28)
!
!       CALCUL DES COORDONNEES ET DE LA LSN DES NOEUDS MILIEUX
!       CENTRAUX D UN ELEMENT QUADRATIQUE
!
!     ENTREE
!       IGEOM    : ADRESSE DES COORDONNÉES DES NOEUDS DE L'ELT PARENT
!       LSN      : LSN DES NOEUDS DE L'ELT PARENT
!       NFISS    : NOMBRE DE FISSURES
!     SORTIE
!       X        : COORDONNEES DES NOEUDS MILIEUX CENTRAUX
!       XLSN     : LSN AUX NOEUDS MILIEUX CENTRAUX
!       NNC      : NOMBRE DE NOEUDS MILIEUX CENTRAUX
!
    integer(kind=8) :: nbnomx
    parameter(nbnomx=20)
    integer(kind=8) :: i, j, nnop, ifiss
    real(kind=8) :: ff(nbnomx), xlsn, xe(3)
    character(len=8) :: elp
!
!
    call elref1(elp)
    call elrefe_info(elrefe=elp, fami='RIGI', nno=nnop)
    tx = 0.d0
    txlsn = 0.d0
!
!     INITIALIASATION PAR DEFAUT DU NOMBRE DE NOEUDS CENTRAUX A ZERO
!     (E.G. TRIANGLES ET TETRAHEDRES)
    nnc = 0
!
!     CALCUL DES COORDONNEES DU MILIEU DE [AB] DANS LE CAS 'H20'
    if (nnop .eq. 20) then
        nnc = 7
        tx(1, 1) = 0.d0
        tx(2, 1) = 0.d0
        tx(3, 1) = -1.0
        tx(1, 2) = 0.d0
        tx(2, 2) = -1.0
        tx(3, 2) = 0.d0
        tx(1, 3) = 1.0
        tx(2, 3) = 0.d0
        tx(3, 3) = 0.d0
        tx(1, 4) = 0.d0
        tx(2, 4) = 1.0
        tx(3, 4) = 0.d0
        tx(1, 5) = -1.0
        tx(2, 5) = 0.d0
        tx(3, 5) = 0.d0
        tx(1, 6) = 0.d0
        tx(2, 6) = 0.d0
        tx(3, 6) = 1.0
        tx(1, 7) = 0.d0
        tx(2, 7) = 0.d0
        tx(3, 7) = 0.d0
!
!     CALCUL DES COORDONNEES DU MILIEU DE [AB] DANS LE CAS 'P15'
    else if (nnop .eq. 15) then
        nnc = 3
        tx(1, 1) = 0.d0
        tx(2, 1) = 0.5
        tx(3, 1) = 0.5
        tx(1, 2) = 0.d0
        tx(2, 2) = 0.d0
        tx(3, 2) = 0.5
        tx(1, 3) = 0.d0
        tx(2, 3) = 0.5
        tx(3, 3) = 0.d0
!
!     CALCUL DES COORDONNEES DU MILIEU DE [AB] DANS LE CAS 'P13'
    else if (nnop .eq. 13) then
        nnc = 1
        tx(1, 1) = 0.d0
        tx(2, 1) = 0.d0
        tx(3, 1) = 0.d0
!     CALCUL DES COORDONNEES DU MILIEU DE [AB] DANS LE CAS 2D 'QU8'
    else if ((nnop .eq. 8) .and. (ndim .eq. 2)) then
        nnc = 1
        tx(1, 1) = 0.d0
        tx(2, 1) = 0.d0
!     CALCUL DES COORDONNEES DU MILIEU DE [AB] DANS LE CAS 3D 'QU8'
    else if ((nnop .eq. 8) .and. (ndim .eq. 3)) then
        nnc = 1
        tx(1, 1) = 0.d0
        tx(2, 1) = 0.d0
        tx(3, 1) = 0.d0
    end if
!
!.....................................................................
!     CALCUL DE LA LSN DU MILIEU
!
    do j = 1, nnc
        do i = 1, ndim
            xe(i) = tx(i, j)
        end do
        call elrfvf(elp, xe, ff)
        do ifiss = 1, nfiss
            xlsn = 0.d0
            do i = 1, nnop
                xlsn = xlsn+ff(i)*lsn((i-1)*nfiss+ifiss)
            end do
            txlsn((j-1)*nfiss+ifiss) = xlsn
        end do
    end do
!
!.....................................................................
!      CALCUL DES COORDONNES DANS L ELEMENT REEL
    do j = 1, nnc
        call reerel(elp, nnop, ndim, zr(igeom), tx(1:ndim, j), xe)
        do i = 1, ndim
            tx(i, j) = xe(i)
        end do
    end do
!
end subroutine
