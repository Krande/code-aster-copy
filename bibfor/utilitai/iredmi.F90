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

subroutine iredmi(macr)
    implicit none
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterc/r8pi.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/iredm1.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/rslipa.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=*) :: macr
!     INTERFACE ASTER - MISS3D : PROCEDURE  IMPR_MACR_ELEM
!     ------------------------------------------------------------------
    integer(kind=8) :: vali(2)
!
    character(len=8) :: k8b, mael, basemo, masse, noma, listam
    character(len=16) :: nomcmd
!     ------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, i2, iam, icamor, icmass, icrigi
    integer(kind=8) :: iret, isamor, ismass, isrigi, ival3
    integer(kind=8) :: j, j2, jamo2, jamor, jfreq, jmass, jordr
    integer(kind=8) :: jrefe, jrigi, k, lamor, n1, n2
    integer(kind=8) :: nbamor, nbmode, nbmods, nbmodt, ntriam, ntriar
    real(kind=8) :: petir8, pi
    real(kind=8), pointer :: mael_raid_vale(:) => null()
    real(kind=8), pointer :: mael_raid_vali(:) => null()
    real(kind=8), pointer :: mael_mass_vale(:) => null()
    real(kind=8), pointer :: mael_mass_vali(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
    mael = macr
    call getres(k8b, k8b, nomcmd)
    pi = r8pi()
    petir8 = 1.d-40
!
!     ----- RECUPERATION DES MODES -----
    call jeveuo(mael//'.MAEL_REFE', 'L', jrefe)
    basemo = zk24(jrefe) (1:8)
    noma = zk24(jrefe+1) (1:8)
    call jelira(basemo//'           .ORDR', 'LONMAX', nbmodt)
    call jeveuo(basemo//'           .ORDR', 'L', jordr)
!
    call dismoi('NB_MODES_DYN', basemo, 'RESULTAT', repi=nbmode)
    call dismoi('NB_MODES_STA', basemo, 'RESULTAT', repi=nbmods)
    nbmodt = nbmode+nbmods
!
    call jeveuo(mael//'.MAEL_MASS_REFE', 'L', jrefe)
    masse = zk24(jrefe+1)
!
!     ----- RECUPERATION DES FREQUENCES -----
    call rslipa(basemo, 'FREQ', '&&IREDMI.LIFREQ', jfreq, nbmodt)
!
!
!     ----- EXTRACTION DU MACRO-ELEMENT DYNAMIQUE -----
!
    if (nbmode .eq. 0) then
        call wkvect('&&IREDMI.DMASS', 'V V R', 1, jmass)
        call wkvect('&&IREDMI.DRIGI', 'V V R', 1, jrigi)
        call wkvect('&&IREDMI.DAMOR', 'V V R', 1, lamor)
    else
        call wkvect('&&IREDMI.DMASS', 'V V R', nbmode*nbmode, jmass)
        call wkvect('&&IREDMI.DRIGI', 'V V R', nbmode*nbmode, jrigi)
        call wkvect('&&IREDMI.DAMOR', 'V V R', nbmode*nbmode, lamor)
    end if
    if (nbmods .eq. 0) then
        call wkvect('&&IREDMI.SMASS', 'V V R', 1, ismass)
        call wkvect('&&IREDMI.SRIGI', 'V V R', 1, isrigi)
        call wkvect('&&IREDMI.SAMOR', 'V V R', 1, isamor)
    else
        call wkvect('&&IREDMI.SMASS', 'V V R', nbmods*nbmods, ismass)
        call wkvect('&&IREDMI.SRIGI', 'V V R', nbmods*nbmods, isrigi)
        call wkvect('&&IREDMI.SAMOR', 'V V R', nbmods*nbmods, isamor)
    end if
    if (nbmode .eq. 0 .or. nbmods .eq. 0) then
        call wkvect('&&IREDMI.CMASS', 'V V R', 1, icmass)
        call wkvect('&&IREDMI.CRIGI', 'V V R', 1, icrigi)
        call wkvect('&&IREDMI.CAMOR', 'V V R', 1, icamor)
    else
        call wkvect('&&IREDMI.CMASS', 'V V R', nbmode*nbmods, icmass)
        call wkvect('&&IREDMI.CRIGI', 'V V R', nbmode*nbmods, icrigi)
        call wkvect('&&IREDMI.CAMOR', 'V V R', nbmode*nbmods, icamor)
    end if
!
    call jelira(mael//'.MAEL_MASS_VALE', 'NMAXOC', ntriam)
    call jeveuo(jexnum(mael//'.MAEL_MASS_VALE', 1), 'L', vr=mael_mass_vale)
    if (ntriam .gt. 1) then
        call jeveuo(jexnum(mael//'.MAEL_MASS_VALE', 2), 'L', vr=mael_mass_vali)
    else
        call jeveuo(jexnum(mael//'.MAEL_MASS_VALE', 1), 'L', vr=mael_mass_vali)
    end if
!
    call jelira(mael//'.MAEL_RAID_VALE', 'NMAXOC', ntriar)
    call jeveuo(jexnum(mael//'.MAEL_RAID_VALE', 1), 'L', vr=mael_raid_vale)
    if (ntriar .gt. 1) then
        call jeveuo(jexnum(mael//'.MAEL_RAID_VALE', 2), 'L', vr=mael_raid_vali)
    else
        call jeveuo(jexnum(mael//'.MAEL_RAID_VALE', 1), 'L', vr=mael_raid_vali)
    end if
!
    do j = 1, nbmode
        do i = 1, j
            k = j*(j-1)/2+i
            zr(jmass+i-1+(j-1)*nbmode) = mael_mass_vale(k)+petir8
            zr(jmass+j-1+(i-1)*nbmode) = mael_mass_vale(k)+petir8
            zr(jrigi+i-1+(j-1)*nbmode) = mael_raid_vale(k)+petir8
            zr(jrigi+j-1+(i-1)*nbmode) = mael_raid_vale(k)+petir8
        end do
    end do
    do j = nbmode+1, nbmodt
        do i = 1, nbmode
            k = j*(j-1)/2+i
            j2 = j-nbmode
            zr(icmass+j2-1+(i-1)*nbmods) = mael_mass_vale(k)+petir8
            zr(icrigi+j2-1+(i-1)*nbmods) = mael_raid_vale(k)+petir8
        end do
        do i = nbmode+1, j
            k = j*(j-1)/2+i
            i2 = i-nbmode
            j2 = j-nbmode
            zr(ismass+i2-1+(j2-1)*nbmods) = mael_mass_vale(k)+petir8
            zr(ismass+j2-1+(i2-1)*nbmods) = mael_mass_vale(k)+petir8
            zr(isrigi+i2-1+(j2-1)*nbmods) = mael_raid_vale(k)+petir8
            zr(isrigi+j2-1+(i2-1)*nbmods) = mael_raid_vale(k)+petir8
        end do
    end do
!
    call jeexin(mael//'.MAEL_AMOR_VALE', iret)
    if (iret .ne. 0) then
        call jeveuo(mael//'.MAEL_AMOR_VALE', 'L', ival3)
        do j = 1, nbmode
            do i = 1, j
                k = j*(j-1)/2+i
                zr(lamor+i-1+(j-1)*nbmode) = zr(ival3+k-1)+petir8
                zr(lamor+j-1+(i-1)*nbmode) = zr(ival3+k-1)+petir8
            end do
        end do
        do j = nbmode+1, nbmodt
            do i = 1, nbmode
                k = j*(j-1)/2+i
                j2 = j-nbmode
                zr(icamor+j2-1+(i-1)*nbmods) = zr(ival3+k-1)+petir8
            end do
            do i = nbmode+1, j
                k = j*(j-1)/2+i
                i2 = i-nbmode
                j2 = j-nbmode
                zr(isamor+i2-1+(j2-1)*nbmods) = zr(ival3+k-1)+petir8
                zr(isamor+j2-1+(i2-1)*nbmods) = zr(ival3+k-1)+petir8
            end do
        end do
    else
        ival3 = 0
    end if
!
!     ----- RECUPERATION DES AMORTISSEMENTS -----
    call getvr8(' ', 'AMOR_REDUIT', nbval=0, nbret=n1)
    call getvid(' ', 'LIST_AMOR', nbval=0, nbret=n2)
    if (nbmode .eq. 0) then
        call wkvect('&&IREDMI.AMORTISSEMENT', 'V V R', 1, jamor)
    else
        call wkvect('&&IREDMI.AMORTISSEMENT', 'V V R', nbmode, jamor)
    end if
    if (n1 .ne. 0 .or. n2 .ne. 0) then
        if (n1 .ne. 0) then
            nbamor = -n1
            call getvr8(' ', 'AMOR_REDUIT', nbval=nbamor, vect=zr(jamor), nbret=n1)
        else
            call getvid(' ', 'LIST_AMOR', scal=listam, nbret=n2)
            call jelira(listam//'           .VALE', 'LONMAX', nbamor)
            call jeveuo(listam//'           .VALE', 'L', jamor)
        end if
        if (nbamor .gt. nbmode) then
            vali(1) = nbamor
            vali(2) = nbmode
            call utmess('F', 'UTILITAI6_44', ni=2, vali=vali)
        end if
        if (nbamor .lt. nbmode) then
            call wkvect('&&IREDMI.AMORTISSEMEN2', 'V V R', nbmode, jamo2)
            do iam = 1, nbamor
                zr(jamo2+iam-1) = zr(jamor+iam-1)
            end do
            do iam = nbamor+1, nbmode
                zr(jamo2+iam-1) = zr(jamor+nbamor-1)
            end do
            nbamor = nbmode
            jamor = jamo2
        end if
    else
        do k = 1, nbmode
            zr(jamor+k-1) = zr( &
                            lamor+(k-1)*(nbmode+1))/(4.d0*pi*zr(jfreq+k-1)*zr(jmass+(k-1)*(nbmo&
                            &de+1)) &
                            )
        end do
    end if
!
    call iredm1(masse, noma, basemo, nbmode, nbmods, &
                ival3, zr(jmass), zr(jrigi), zr(jamor), zr(jfreq), &
                zr(ismass), zr(isrigi), zr(isamor), zr(icmass), zr(icrigi), &
                zr(icamor))
!
!
! --- MENAGE
!
    call jedetr('&&IREDMI.LIFREQ')
    call jedetr('&&IREDMI.DMASS')
    call jedetr('&&IREDMI.DRIGI')
    call jedetr('&&IREDMI.DAMOR')
    call jedetr('&&IREDMI.SMASS')
    call jedetr('&&IREDMI.SRIGI')
    call jedetr('&&IREDMI.SAMOR')
    call jedetr('&&IREDMI.CMASS')
    call jedetr('&&IREDMI.CRIGI')
    call jedetr('&&IREDMI.CAMOR')
    call jedetr('&&IREDMI.AMORTISSEMENT')
    call jedetr('&&IREDMI.AMORTISSEMEN2')
!
    call jedema()
end subroutine
