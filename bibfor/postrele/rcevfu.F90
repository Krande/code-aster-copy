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
subroutine rcevfu(cnoc, cfat, fut)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
    real(kind=8) :: fut
    character(len=24) :: cnoc, cfat
!     OPERATEUR POST_RCCM, TYPE_RESU_MECA='EVOLUTION'
!     CALCUL DU FACTEUR D'USAGE TOTAL (FUT)
!
!     ------------------------------------------------------------------
!
    integer(kind=8) :: nbinst, jnocr, jfu, i1, i2, ind, noc1, noc2, i1m, i2m, noc1m
    integer(kind=8) :: noc2m, nbcycl
    integer(kind=8) :: indi, inds, k, l, ifm, niv
    real(kind=8) :: fum, fukl
    aster_logical :: encore
    real(kind=8), pointer :: matr_fu(:) => null()
    integer(kind=8), pointer :: nb_occ_k(:) => null()
    integer(kind=8), pointer :: nb_occ_l(:) => null()
! DEB ------------------------------------------------------------------
    call jemarq()
!
    call infniv(ifm, niv)
!
    call jelira(cnoc, 'LONMAX', nbinst)
    call jeveuo(cnoc, 'L', jnocr)
    call jeveuo(cfat, 'L', jfu)
!
    fut = 0.d0
!
    AS_ALLOCATE(vi=nb_occ_k, size=nbinst)
    AS_ALLOCATE(vi=nb_occ_l, size=nbinst)
    do i1 = 1, nbinst
        nb_occ_k(i1) = zi(jnocr+i1-1)
        nb_occ_l(i1) = zi(jnocr+i1-1)
    end do
!
    AS_ALLOCATE(vr=matr_fu, size=nbinst*nbinst)
    ind = 0
    do i1 = 1, nbinst
        indi = nbinst*(i1-1)+i1
        ind = ind+1
        matr_fu(indi) = zr(jfu-1+5*(ind-1)+4)
        do i2 = i1+1, nbinst
            inds = nbinst*(i1-1)+i2
            indi = nbinst*(i2-1)+i1
            ind = ind+1
            matr_fu(inds) = zr(jfu-1+5*(ind-1)+4)
            matr_fu(indi) = zr(jfu-1+5*(ind-1)+4)
        end do
    end do
!
    ifm = 6
    ind = 0
!
100 continue
    ind = ind+1
!
    if (niv .eq. 2) then
        if (ind .eq. 1) then
            write (ifm, *) 'MATRICE FACTEUR D''USAGE INITIALE'
        else
            write (ifm, *) 'MATRICE FACTEUR D''USAGE MODIFIEE'
        end if
        write (ifm, 1010) (nb_occ_l(l), l=1, nbinst)
        do k = 1, nbinst
            i1 = nbinst*(k-1)
            write (ifm, 1000) nb_occ_k(k), (matr_fu(i1+l), l=1, &
                                            nbinst)
        end do
    end if
!
    fum = 0.d0
!
    do i1 = 1, nbinst
        noc1 = nb_occ_k(i1)
        if (noc1 .eq. 0) goto 110
        k = nbinst*(i1-1)
!
        do i2 = 1, nbinst
            noc2 = nb_occ_l(i2)
            if (noc2 .eq. 0) goto 112
            l = i2
!
            fukl = matr_fu(k+l)
            if (fukl .lt. r8prem()) goto 112
            if (fukl .gt. fum) then
                noc1m = noc1
                noc2m = noc2
                i1m = i1
                i2m = i2
                fum = fukl
            end if
!
112         continue
        end do
!
110     continue
    end do
    nbcycl = min(noc1m, noc2m)
!
    if (fum .lt. r8prem()) goto 999
    if (niv .eq. 2) then
        write (ifm, 1020) '=> FACTEUR D''USAGE MAXI: ', fum, i1m, i2m
        write (ifm, 1030) '   NB_OCCUR = ', nbcycl
    end if
!
! --- ON CUMULE
!
    fut = fut+fum*dble(nbcycl)
!
! --- ON MET A ZERO LES FACTEURS D'USAGE INCRIMINES
!
    if (noc1m .eq. noc2m) then
        nb_occ_l(i1m) = 0
        nb_occ_l(i2m) = 0
        nb_occ_k(i1m) = 0
        nb_occ_k(i2m) = 0
        do k = 1, nbinst
            matr_fu((k-1)*nbinst+i2m) = 0.d0
            matr_fu((i2m-1)*nbinst+k) = 0.d0
            matr_fu((k-1)*nbinst+i1m) = 0.d0
            matr_fu((i1m-1)*nbinst+k) = 0.d0
        end do
    else if (noc1m .lt. noc2m) then
        nb_occ_l(i2m) = nb_occ_l(i2m)-noc1m
        nb_occ_k(i2m) = nb_occ_k(i2m)-noc1m
        nb_occ_k(i1m) = 0
        nb_occ_l(i1m) = 0
        do k = 1, nbinst
            matr_fu((i1m-1)*nbinst+k) = 0.d0
            matr_fu((k-1)*nbinst+i1m) = 0.d0
        end do
    else
        nb_occ_l(i2m) = 0
        nb_occ_k(i2m) = 0
        nb_occ_l(i1m) = nb_occ_l(i1m)-noc2m
        nb_occ_k(i1m) = nb_occ_k(i1m)-noc2m
        do k = 1, nbinst
            matr_fu((k-1)*nbinst+i2m) = 0.d0
            matr_fu((i2m-1)*nbinst+k) = 0.d0
        end do
    end if
!
! --- EXISTE-T-IL DES ETATS TELS QUE NB_OCCUR > 0
!
    encore = .false.
    do i1 = 1, nbinst
        if (nb_occ_k(i1) .gt. 0) then
            encore = .true.
        end if
    end do
    if (encore) goto 100
!
999 continue
!
    if (niv .eq. 2) write (ifm, *) '-->> FACTEUR D''USAGE CUMULE = ', fut
!
    AS_DEALLOCATE(vi=nb_occ_k)
    AS_DEALLOCATE(vi=nb_occ_l)
    AS_DEALLOCATE(vr=matr_fu)
!
1000 format(1p, i10, '|', 40(e10.3, '|'))
1010 format(1p, '   NB_OCCUR ', '|', 40(i10, '|'))
1020 format(1p, a28, e12.5, ', LIGNE:', i4, ', COLONNE:', i4)
1030 format(1p, a15, i8)
!
    call jedema()
end subroutine
