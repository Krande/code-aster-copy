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

subroutine amogen(mat19)
    implicit none
#include "jeveux.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=19) :: mat19
    character(len=8) :: masse, raid, listam
    character(len=16) :: nomcmd
    integer(kind=8) :: nbid, jamog, iamog, idiff
    integer(kind=8) :: vali(3)
    integer(kind=8) :: iamat, n, m, m2, i, iam, iak, j, nbamor, nlist
    integer(kind=8) :: iblo, iaconl, jrefa2, iadesc, n2, n1
    real(kind=8) :: kmin, valmin, kmax, rk
    integer(kind=8), pointer :: desc(:) => null()
    character(len=24), pointer :: refa(:) => null()
!
    call jemarq()
!
    nomcmd = 'CALC_AMOR_GENE'
    call getvid(nomcmd, 'MASS_GENE', iocc=1, scal=masse, nbret=nbid)
    call getvid(nomcmd, 'RIGI_GENE', iocc=1, scal=raid, nbret=nbid)
    call getvr8(nomcmd, 'AMOR_REDUIT', iocc=1, nbval=0, nbret=n1)
    call getvid(nomcmd, 'LIST_AMOR', iocc=1, nbval=0, nbret=n2)
    call jeveuo(masse//'           .DESC', 'E', vi=desc)
    n = desc(2)
    call jelira(masse//'           .VALM', 'LONMAX', m)
    call jelira(raid//'           .VALM', 'LONMAX', m2)
    if (m2 .ne. m) then
        vali(1) = m
        vali(2) = m2
        call utmess('F', 'ALGELINE5_28', ni=2, vali=vali)
    end if
!
    if (n1 .ne. 0) then
        nbamor = -n1
    else
        call getvid(nomcmd, 'LIST_AMOR', iocc=1, scal=listam, nbret=nlist)
        call jelira(listam//'           .VALE', 'LONMAX', nbamor)
    end if
!
    if (nbamor .gt. n) then
!
        vali(1) = n
        vali(2) = nbamor
        vali(3) = n
        call utmess('A', 'ALGELINE5_29', ni=3, vali=vali)
        call wkvect('&&AMORMA.AMORTI', 'V V R8', n, jamog)
        if (n1 .ne. 0) then
            call getvr8(nomcmd, 'AMOR_REDUIT', iocc=1, nbval=n, vect=zr(jamog), &
                        nbret=nbid)
        else
            call jeveuo(listam//'           .VALE', 'L', iamog)
            do i = 1, n
                zr(jamog+i-1) = zr(iamog+i-1)
            end do
        end if
    else
        call wkvect('&&AMORMA.AMORTI', 'V V R8', n, jamog)
        if (n1 .ne. 0) then
            call getvr8(nomcmd, 'AMOR_REDUIT', iocc=1, nbval=nbamor, vect=zr(jamog), &
                        nbret=nbid)
        else
            call jeveuo(listam//'           .VALE', 'L', iamog)
            do i = 1, nbamor
                zr(jamog+i-1) = zr(iamog+i-1)
            end do
        end if
        if (nbamor .lt. n) then
            do i = nbamor+1, n
                zr(jamog+i-1) = zr(jamog+nbamor-1)
            end do
!
            idiff = n-nbamor
            vali(1) = idiff
            vali(2) = n
            vali(3) = idiff
            call utmess('I', 'ALGELINE5_30', ni=3, vali=vali)
        end if
    end if
    iblo = 1
    call jeveuo(masse//'           .REFA', 'E', vk24=refa)
!
!   CREATION DES BASES DE DONNEES DE LA MATRICE A GENERER.
!   SUIVANT LE MODELE DE OP0071
!
    call jecrec(mat19//'.VALM', 'G V R', 'NU', 'DISPERSE', 'CONSTANT', &
                1)
    call jecroc(jexnum(mat19//'.VALM', iblo))
    call jeecra(mat19//'.VALM', 'LONMAX', m)
    call wkvect(mat19//'.CONL', 'G V R', n, iaconl)
    call wkvect(mat19//'.REFA', 'G V K24', 20, jrefa2)
    zk24(jrefa2-1+11) = 'MPI_COMPLET'
    zk24(jrefa2-1+1) = refa(1)
    zk24(jrefa2-1+2) = refa(2)
    zk24(jrefa2-1+9) = 'MS'
    zk24(jrefa2-1+10) = 'GENE'
!
! ----- CREATION DU .DESC
!
    call wkvect(mat19//'.DESC', 'G V I', 3, iadesc)
    zi(iadesc) = 2
    zi(iadesc+1) = n
    zi(iadesc+2) = 2
!
    do i = 1, n
        zr(iaconl+i-1) = 1.0d0
    end do
!
    iblo = 1
    call jeveuo(jexnum(masse//'           .VALM', iblo), 'L', iam)
    call jeveuo(jexnum(raid//'           .VALM', iblo), 'E', iak)
    call jeveuo(jexnum(mat19//'.VALM', iblo), 'E', iamat)
    do i = 1, m
        zr(iamat-1+i) = 0d0
    end do
    kmin = 0.d0
    kmax = 0.00001d0
    valmin = 1.d-4
    do i = 1, n
        if (m .eq. n*(n+1)/2) then
            j = i*(i+1)/2-1
        else if (m .eq. n) then
            j = i-1
        else
            goto 190
        end if
        if (zr(iak+j) .lt. (0.0d0)) then
            kmin = min(kmin, zr(iak+j))
            zr(iak+j) = 0.0d0
        else
            kmax = max(kmax, zr(iak+j))
        end if
        zr(iamat+j) = 2.0d0*zr(jamog+i-1)*sqrt(abs(zr(iam+j)*zr(iak+j)))
190     continue
    end do
    rk = kmin/kmax
    if (abs(rk) .ge. valmin) then
        call utmess('A', 'PREPOST4_20')
!         CALL UTMESS('F','PREPOST4_21')
    end if
    call jedetr('&&AMORMA.AMORTI')
    call jedema()
end subroutine
