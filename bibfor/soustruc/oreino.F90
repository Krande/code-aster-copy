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
subroutine oreino(noma, lnoeud, nbno, nori, next, &
                  coor, crit, prec, iera, ier)
    implicit none
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/utmess.h"
#include "asterfort/int_to_char8.h"
!
    integer(kind=8) :: lnoeud(*), nbno, nori, next, ier, iera
    real(kind=8) :: coor(*), prec
    character(len=8) :: noma
    character(len=*) :: crit
!     BUT : CLASSER DES NUMEROS DE NOEUDS SELON LEUR PROJECTION
!           SUR UN SEGMENT
!-----------------------------------------------------------------------
!     I/O : LNOEUD: NUMEROS DES NOEUDS
!     IN  : NBNO  : NOMBRE DE NOEUDS
!     IN  : IREA  : INDICATEUR DE COMMANDE : 1-DEFI_GROUP/SEGM_DROI_ORDO
!                                            2-POST_RELEVE_T/PRECISION
!                                            3-DEFI_FOND_FISS/PREC_NORM
!     IN  : NORI  : NUMERO DU NOEUD ORIGINE
!     IN  : NEXT  : NUMERO DU NOEUD EXTREMITE
!     IN  : COOR  : COORDONNEES DES NOEUDS
!     IN  : CRIT  : CRITERE
!     IN  : PREC  : PRECISION
!     IN  : IER   : CODE RETOUR,  = 0  OK
!-----------------------------------------------------------------------
    integer(kind=8) :: i, j, k, n, inoe, inod
    real(kind=8) :: xa, ya, za, xb, yb, zb, xab, yab, zab, ab2, xm, ym, zm, xam
    real(kind=8) :: yam, zam, c, c2, xv, yv, zv, v2, r8b, ecart, valr
    character(len=8) :: nomn
    character(len=24) :: valk(2)
    real(kind=8), pointer :: bary(:) => null()
! DEB ------------------------------------------------------------------
    call jemarq()
!
    ier = 0
!
    xa = coor(3*(nori-1)+1)
    ya = coor(3*(nori-1)+2)
    za = coor(3*(nori-1)+3)
!
    xb = coor(3*(next-1)+1)
    yb = coor(3*(next-1)+2)
    zb = coor(3*(next-1)+3)
!
    xab = xb-xa
    yab = yb-ya
    zab = zb-za
    ab2 = xab**2+yab**2+zab**2
    if (ab2 .eq. 0.0d0) then
        call utmess('A', 'SOUSTRUC_20')
        ier = ier+1
        goto 999
    end if
!
    AS_ALLOCATE(vr=bary, size=nbno)
!
!     --- CALCUL DE LA CORDONNEE BARYCENTRIQUE ---
!
    do inoe = 1, nbno
        inod = lnoeud(inoe)
        xm = coor(3*(inod-1)+1)
        ym = coor(3*(inod-1)+2)
        zm = coor(3*(inod-1)+3)
        xam = xm-xa
        yam = ym-ya
        zam = zm-za
        c = (xam*xab+yam*yab+zam*zab)/ab2
        c2 = xam**2+yam**2+zam**2
        xv = xam-c*xab
        yv = yam-c*yab
        zv = zam-c*zab
        v2 = xv**2+yv**2+zv**2
!        --- VERIFICATION QUE LA DISTANCE A L'AXE
!                         NE DEPASSE PAS LA TOLERANCE ---
        if (crit(1:4) .eq. 'ABSO') then
            r8b = v2
        else if (crit(1:4) .eq. 'RELA') then
            r8b = v2/ab2
        else
            call utmess('A', 'SOUSTRUC_21')
            ier = ier+1
            goto 999
        end if
        r8b = sqrt(r8b)
        if (r8b .gt. prec) then
            v2 = sqrt(v2)
            nomn = int_to_char8(inod)
            if (iera .eq. 3) then
                call utmess('A', 'SOUSTRUC_17', sk=nomn, sr=v2)
            else
                call utmess('A', 'SOUSTRUC_22', sk=nomn, sr=v2)
            end if
            ier = ier+1
        end if
!        --- VERIFICATION QUE LA PROJECTION EST BIEN
!                         SITUEE ENTRE LES POINTS A ET B ---
        ecart = (c2-ab2)/ab2
        if (c .lt. 0.0d0 .or. c2 .gt. ab2) then
            if (ecart .gt. r8prem()) then
                nomn = int_to_char8(inod)
                valk(1) = nomn
                valk(2) = nomn
                valr = c
                call utmess('A', 'SOUSTRUC_86', nk=2, valk=valk, sr=valr)
                ier = ier+1
            end if
        end if
        bary(inoe) = c
    end do
!
!     --- TRI PAR BUBBLE SORT ---
!
    do k = 1, nbno-1
        do i = nbno-1, k, -1
            j = i+1
            if (bary(i) .gt. bary(j)) then
                c = bary(j)
                bary(j) = bary(i)
                bary(i) = c
                n = lnoeud(j)
                lnoeud(j) = lnoeud(i)
                lnoeud(i) = n
            end if
        end do
    end do
!
!     --- VERIFICATION QUE DEUX NOEUDS CONSECUTIFS
!                          N'ONT PAS LA MEME PROJECTION ---
    do inoe = 1, nbno-1
        if (bary(inoe) .eq. bary(inoe+1)) then
            call utmess('A', 'SOUSTRUC_23')
            ier = ier+1
        end if
    end do
!
999 continue
!
    AS_DEALLOCATE(vr=bary)
!
    call jedema()
end subroutine
