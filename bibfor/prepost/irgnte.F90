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
subroutine irgnte(ifi, nbordr, coord, connex, point, &
                  njvmai, nbmai, cnsv, partie, jtype, &
                  cnsd)
    implicit none
!
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jelibe.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
    integer(kind=8) :: ifi, nbordr, connex(*), point(*), cnsv(*), cnsd(*), jtype
    real(kind=8) :: coord(*)
    character(len=*) :: njvmai, partie
!
!     IMPRESSION D'UN CHAM_NO AU FORMAT GMSH :
!
!     ------------------------------------------------------------------
!
    integer(kind=8) :: imai, ima, ipoin, listno(8), j, jcnsv, jcnsd, ncmp, jmai, ior
    integer(kind=8) :: inoe, nbno, nbmai
    real(kind=8) :: zero
!     ------------------------------------------------------------------
!
    call jemarq()
!
    zero = 0.0d0
!
    call jeveuo(njvmai, 'L', jmai)
    if (njvmai(10:12) .eq. 'POI') then
        nbno = 1
    else if (njvmai(10:12) .eq. 'SEG') then
        nbno = 2
    else if (njvmai(10:12) .eq. 'TRI') then
        nbno = 3
    else if (njvmai(10:12) .eq. 'QUA') then
        nbno = 4
    else if (njvmai(10:12) .eq. 'TET') then
        nbno = 4
    else if (njvmai(10:12) .eq. 'PYR') then
        nbno = 5
    else if (njvmai(10:12) .eq. 'PRI') then
        nbno = 6
    else if (njvmai(10:12) .eq. 'HEX') then
        nbno = 8
    end if
!
    do imai = 1, nbmai
!
        ima = zi(jmai-1+imai)
!
        ipoin = point(ima)
!
        do inoe = 1, nbno
!
            listno(inoe) = connex(ipoin+inoe-1)
!
        end do
!
        do j = 1, 3
            write (ifi, 1000) (coord(3*(listno(inoe)-1)+j), inoe=1, nbno)
        end do
!
        do ior = 1, nbordr
!
            jcnsv = cnsv(ior)
            jcnsd = cnsd(ior)
            ncmp = zi(jcnsd-1+2)
            if (zk8(jtype-1+ior) .eq. 'R') then
!
                do inoe = 1, nbno
!
                    if (njvmai(10:12) .eq. 'SEG' .or. njvmai(10:12) .eq. 'TRI' .or. &
                        njvmai(10:12) .eq. 'QUA') then
!
                        write (ifi, 1000) zr(jcnsv-1+(listno(inoe)-1)* &
                                             ncmp+1), zr(jcnsv-1+(listno(inoe)-1)*ncmp+4), &
                            zero, zr(jcnsv-1+(listno(inoe)-1)*ncmp+4), &
                            zr(jcnsv-1+(listno(inoe)-1)*ncmp+2), zero, &
                            zero, zero, zr(jcnsv-1+(listno(inoe)-1)*ncmp+ &
                                           3)
!
                    elseif (njvmai(10:12) .eq. 'PYR' .or. njvmai(10:12) &
                            .eq. 'PRI' .or. njvmai(10:12) .eq. 'HEX') then
!
                        write (ifi, 1000) zr(jcnsv-1+(listno(inoe)-1)* &
                                             ncmp+1), zr(jcnsv-1+(listno(inoe)-1)*ncmp+4), &
                            zr(jcnsv-1+(listno(inoe)-1)*ncmp+5), zr(jcnsv- &
                                                          1+(listno(inoe)-1)*ncmp+4), zr(jcnsv-1+( &
                                                      listno(inoe)-1)*ncmp+2), zr(jcnsv-1+(listno( &
                                                    inoe)-1)*ncmp+6), zr(jcnsv-1+(listno(inoe)-1)* &
                                                     ncmp+5), zr(jcnsv-1+(listno(inoe)-1)*ncmp+6), &
                            zr(jcnsv-1+(listno(inoe)-1)*ncmp+3)
!
                    end if
!
                end do
!
            else if (zk8(jtype-1+ior) .eq. 'C') then
!
                if (partie .eq. 'REEL') then
!
                    do inoe = 1, nbno
!
                        if (njvmai(10:12) .eq. 'SEG' .or. njvmai(10:12) .eq. 'TRI' .or. &
                            njvmai(10:12) .eq. 'QUA') then
!
                            write (ifi, 1000) zr(jcnsv-1+(listno(inoe)- &
                                                          1)*ncmp+1), zr(jcnsv-1+(listno(inoe)-1)* &
                                                        ncmp+4), zero, zr(jcnsv-1+(listno(inoe)-1) &
                                                            *ncmp+4), zr(jcnsv-1+(listno(inoe)-1)* &
                                                           ncmp+2), zero, zero, zero, zr(jcnsv-1+( &
                                                                             listno(inoe)-1)*ncmp+3)
!
                        elseif (njvmai(10:12) .eq. 'PYR' .or. njvmai(10: &
                                                  12) .eq. 'PRI' .or. njvmai(10:12) .eq. 'HEX') then
!
                            write (ifi, 1000) zr(jcnsv-1+(listno(inoe)- &
                                                          1)*ncmp+1), zr(jcnsv-1+(listno(inoe)-1)* &
                                                        ncmp+4), zr(jcnsv-1+(listno(inoe)-1)*ncmp+ &
                                                          5), zr(jcnsv-1+(listno(inoe)-1)*ncmp+4), &
                                zr(jcnsv-1+(listno(inoe)-1)*ncmp+2), &
                                zr(jcnsv-1+(listno(inoe)-1)*ncmp+6), &
                                zr(jcnsv-1+(listno(inoe)-1)*ncmp+5), &
                                zr(jcnsv-1+(listno(inoe)-1)*ncmp+6), &
                                zr(jcnsv-1+(listno(inoe)-1)*ncmp+3)
!
                        end if
!
                    end do
                else if (partie .eq. 'IMAG') then
                    do inoe = 1, nbno
!
                        if (njvmai(10:12) .eq. 'SEG' .or. njvmai(10:12) .eq. 'TRI' .or. &
                            njvmai(10:12) .eq. 'QUA') then
!
                            write (ifi, 1000) zr(jcnsv-1+(listno(inoe)- &
                                                          1)*ncmp+1), zr(jcnsv-1+(listno(inoe)-1)* &
                                                        ncmp+4), zero, zr(jcnsv-1+(listno(inoe)-1) &
                                                            *ncmp+4), zr(jcnsv-1+(listno(inoe)-1)* &
                                                           ncmp+2), zero, zero, zero, zr(jcnsv-1+( &
                                                                             listno(inoe)-1)*ncmp+3)
!
                        elseif (njvmai(10:12) .eq. 'PYR' .or. njvmai(10: &
                                                  12) .eq. 'PRI' .or. njvmai(10:12) .eq. 'HEX') then
!
                            write (ifi, 1000) zr(jcnsv-1+(listno(inoe)- &
                                                          1)*ncmp+1), zr(jcnsv-1+(listno(inoe)-1)* &
                                                        ncmp+4), zr(jcnsv-1+(listno(inoe)-1)*ncmp+ &
                                                          5), zr(jcnsv-1+(listno(inoe)-1)*ncmp+4), &
                                zr(jcnsv-1+(listno(inoe)-1)*ncmp+2), &
                                zr(jcnsv-1+(listno(inoe)-1)*ncmp+6), &
                                zr(jcnsv-1+(listno(inoe)-1)*ncmp+5), &
                                zr(jcnsv-1+(listno(inoe)-1)*ncmp+6), &
                                zr(jcnsv-1+(listno(inoe)-1)*ncmp+3)
                        end if
                    end do
                end if
            end if
        end do
    end do
!
    call jelibe(njvmai)
!
    call jedema()
!
1000 format(1p, 9(e15.7e3, 1x))
!
end subroutine
