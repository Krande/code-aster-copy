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
subroutine fonno7(noma, cnxinv, ndim, na, vecdir, &
                  hmax)
    implicit none
#include "jeveux.h"
#include "asterc/r8dgrd.h"
#include "asterc/r8prem.h"
#include "asterfort/conare.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "blas/ddot.h"
!
    integer(kind=8) :: na, ndim
    real(kind=8) :: vecdir(3), hmax
    character(len=8) :: noma
    character(len=19) :: cnxinv
!
!
!       ----------------------------------------------------------------
!       DETERMINATION DE LA TAILLE MAXIMALE DES MAILLES CONNECTEES AU
!       NOEUD DU FOND
!       ----------------------------------------------------------------
!    ENTREES
!       NOMA   : NOM DU MAILLAGE
!       CNXINV : CONNECTIVITE INVERSE
!       NDIM   : DIMENSION DU MODELE
!       NA     : NUMERO DU NOEUD SOMMET COURANT
!       VECDIR : VECTEUR TANGENT
!    SORTIE
!       HMAX  : TAILLE MAXIMALE DES MAILLES
!
!
    integer(kind=8) :: adra, ar(12, 3)
    integer(kind=8) :: iatyma, iar, ima, ino1, ino2, ityp
    integer(kind=8) :: jcncin, jconx2, jdrvlc, k
    integer(kind=8) :: nbar, nbmaca, ndime, nno, nno1, nno2, numac
    real(kind=8) :: coor(3), vect(3), p, cos70, cosinu, normv
    character(len=8) :: type
    integer(kind=8), pointer :: connex(:) => null()
    real(kind=8), pointer :: vale(:) => null()
    blas_int :: b_incx, b_incy, b_n
!     -----------------------------------------------------------------
!
    call jemarq()
!
!     RECUPERATION DES DONNEES SUR LE MAILLAGE
    call jeveuo(noma//'.CONNEX', 'L', vi=connex)
    call jeveuo(jexatr(noma//'.CONNEX', 'LONCUM'), 'L', jconx2)
    call jeveuo(noma//'.TYPMAIL', 'L', iatyma)
    call jeveuo(noma//'.COORDO    .VALE', 'L', vr=vale)
!
!     RECUPERATION DES COORDONNNES DE NA
    coor(1) = vale((na-1)*3+1)
    coor(2) = vale((na-1)*3+2)
    coor(3) = vale((na-1)*3+3)
!
!     RECUPERATION DE LA CONNECTIVITE INVERSE
    call jeveuo(jexatr(cnxinv, 'LONCUM'), 'L', jdrvlc)
    call jeveuo(jexnum(cnxinv, 1), 'L', jcncin)
!
!     MAILLES CONNECTEES A NA
    adra = zi(jdrvlc-1+na)
    nbmaca = zi(jdrvlc-1+na+1)-zi(jdrvlc-1+na)
!
    hmax = r8prem()
!
    do ima = 1, nbmaca
!       NUMERO DE LA MAILLE
        numac = zi(jcncin-1+adra+ima-1)
        ityp = iatyma-1+numac
        call jenuno(jexnum('&CATA.TM.NOMTM', zi(ityp)), type)
        call dismoi('DIM_TOPO', type, 'TYPE_MAILLE', repi=ndime)
!
!       ON ZAPPE LES MAILLES DE BORDS
        if (ndime .ne. ndim) goto 10
!
        call conare(type, ar, nbar)
!
!       BOUCLE SUR LE NOMBRE D'ARETES DE LA MAILLE NUMAC
        do iar = 1, nbar
!
            ino1 = ar(iar, 1)
            nno1 = connex(zi(jconx2+numac-1)+ino1-1)
            ino2 = ar(iar, 2)
            nno2 = connex(zi(jconx2+numac-1)+ino2-1)
!
            if (na .eq. nno1) then
                nno = nno2
            else if (na .eq. nno2) then
                nno = nno1
            else
                goto 100
            end if
!
!          VECTEUR REPRESENTANT L'ARETE NA-NNO
            do k = 1, ndim
                vect(k) = vale((nno-1)*3+k)-coor(k)
            end do
!
!          PROJECTION DE L'ARETE SUR LE VECTEUR TANGENT
            b_n = to_blas_int(ndim)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            p = ddot(b_n, vect, b_incx, vecdir, b_incy)
!
!          FILTRAGE DES ARETES A PRENDRE EN COMPTE:
!          L'ANGLE ENTRE VECT ET VECDIR DOIT ETRE <60
            b_n = to_blas_int(ndim)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            normv = sqrt(ddot(b_n, vect, b_incx, vect, b_incy))
!
            cosinu = p/normv
            cos70 = cos(70*r8dgrd())
            if (abs(cosinu) .lt. cos70) goto 100
!
!          ON PREND LE MAX DES PROJECTIONS
            p = abs(p)
            if (p .ge. hmax) hmax = p
!
100         continue
        end do
!
10      continue
    end do
!
    if (hmax .le. r8prem()) then
        call utmess('A', 'RUPTURE0_49')
    end if
!
    call jedema()
end subroutine
