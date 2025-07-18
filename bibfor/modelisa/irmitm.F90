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
subroutine irmitm(nbmode, ifmis, freq, tabrig, ibin)
    implicit none
!
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/wkvect.h"
!
    character(len=24) :: tabrig, tabfrq, tabri2, tabri0
    real(kind=8) :: a(3), nins2, nbm
!      INTEGER*8    LONG1,LONG2,LONG3
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, i1, i2, ic, ifmis, ifreq, jfrq
    integer(kind=8) :: jri0, jri2, jrig, nbmode, nfreq, ibin
    real(kind=8) :: freq, pas
!-----------------------------------------------------------------------
    call jemarq()
!
!
    tabri2 = '&&IRMITM.RIG2'
    tabri0 = '&&IRMITM.RIG0'
    tabfrq = '&&IRMITM.FREQ'
    call wkvect(tabri2, 'V V R', nbmode*nbmode, jri2)
    call wkvect(tabri0, 'V V R', nbmode*nbmode, jri0)
!
!   Lecture d'entiers INTEGER*8 en binaire venant de MISS3D
!   On convertit ensuite en INTEGER (*4 sur machine 32 bits, sinon *8).
!   Les reels ne posent pas de probleme : ce sont toujours des REAL*8
!
    if (ibin .eq. 0) then
        rewind ifmis
        read (ifmis, *) nins2, pas
    else
        open (unit=ifmis, form='unformatted', status='old', access='stream')
        read (ifmis) nins2, pas, nbm
    end if
!
    nfreq = int(nins2)
    write (6, *) 'ninst pas ', nfreq, pas
    call wkvect(tabfrq, 'V V R', nfreq, jfrq)
    ic = 1
    call jeveuo(tabrig, 'E', jrig)
    if (ibin .eq. 0) then
        do i = 1, nfreq
            zr(jfrq+i-1) = (i-1)*pas
        end do
    else
        read (ifmis) (zr(jfrq+i-1), i=1, nfreq)
    end if
    write (6, *) 'lfreq ', (zr(jfrq+i-1), i=1, nfreq)
    do i = 1, nfreq
        a(1) = zr(jfrq+i-1)
        if (freq .le. (a(1)+r8prem())) then
            ifreq = i
            if (i .gt. 1 .and. freq .lt. (a(1)-r8prem())) then
                ifreq = ifreq-1
            end if
            if (freq .le. r8prem()) ic = 2
            if (i .eq. 1 .and. nfreq .eq. 1) ic = 0
            if (i .eq. nfreq .and. freq .ge. (a(1)-r8prem())) then
                ic = 0
                ifreq = nfreq
            end if
            goto 7
        end if
    end do
    ifreq = nfreq
    ic = 0
7   continue
    if (ibin .eq. 0) then
        do i = 1, ifreq-1
            read (ifmis, *) a(1)
            read (ifmis, 100) ((zr(jri0+(i2-1)*nbmode+i1-1), i1=1, nbmode),&
        &i2=1, nbmode)
        end do
        read (ifmis, *) a(1)
        read (ifmis, 100) ((zr(jrig+(i2-1)*nbmode+i1-1),&
       &              i1=1, nbmode), i2=1, nbmode)
    else
        do i = 1, ifreq-1
            read (ifmis) ((zr(jri0+(i2-1)*nbmode+i1-1), i1=1, nbmode),&
        &i2=1, nbmode)
        end do
        read (ifmis) ((zr(jrig+(i2-1)*nbmode+i1-1),&
       &              i1=1, nbmode), i2=1, nbmode)
    end if
    if (ic .ge. 1) then
        if (ibin .eq. 0) then
            read (ifmis, *) a(1)
            read (ifmis, 100) ((zr(jri2+(i2-1)*nbmode+i1-1), i1=1, nbmode),&
          &i2=1, nbmode)
        else
            read (ifmis) ((zr(jri2+(i2-1)*nbmode+i1-1), i1=1, nbmode),&
          &i2=1, nbmode)
        end if
        do i1 = 1, nbmode
            do i2 = 1, nbmode
                zr(jrig+(i2-1)*nbmode+i1-1) = zr( &
                                              jrig+(i2-1)*nbmode+i1-1)+(freq-zr(jfrq+ifreq-1))&
                                              &/(zr(jfrq+ifreq)-zr(jfrq+ifreq-1))*(zr(jri2+(i2&
                                              &-1)*nbmode+i1-1)-zr(jrig+(i2-1)*nbmode+i1-1) &
                                              )
            end do
        end do
    end if
    if (ibin .ne. 0) then
        close (unit=ifmis)
    end if
!
    call jedetr(tabri0)
    call jedetr(tabri2)
    call jedetr(tabfrq)
!
100 format((6(1x, 1pe13.6)))
    call jedema()
end subroutine
