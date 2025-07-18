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
subroutine dlidef()
    implicit none
!
!     COMMANDE : DEFI_LIST_ENTI/OPERATION='DEFI'
!
!     ------------------------------------------------------------------
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterfort/getvis.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: idebut, irest, iii, ipdt
    integer(kind=8) :: vali(2)
    character(len=8) :: resu
    character(len=16) :: nomcmd, concep
    integer(kind=8) :: i, j, ndim, nbvale, nv, jval, jbor, nbocc, jnbp, n1, nbval
    integer(kind=8) :: kval, ico, np, jpas, iocc
!     ------------------------------------------------------------------
    call jemarq()
!
    nbval = 1
!
    call getres(resu, concep, nomcmd)
    call getvis(' ', 'VALE', nbval=0, nbret=nv)
    call getvis(' ', 'DEBUT', scal=idebut, nbret=n1)
    call getfac('INTERVALLE', nbocc)
!
    if (nv .ne. 0) then
!
    else
        call wkvect('&&DLIDEF.BORNE', 'V V I', nbocc+1, jbor)
        zi(jbor) = idebut
        do iocc = 1, nbocc
            call getvis('INTERVALLE', 'JUSQU_A', iocc=iocc, scal=zi(jbor+iocc), nbret=n1)
            iii = zi(jbor+iocc)-zi(jbor-1+iocc)
            if (iii .le. 0) then
                vali(1) = zi(jbor+iocc-1)
                vali(2) = zi(jbor+iocc)
                call utmess('F', 'ALGORITH13_78', ni=2, vali=vali)
            end if
            call getvis('INTERVALLE', 'PAS', iocc=iocc, nbval=0, nbret=np)
            if (np .ne. 0) then
                call getvis('INTERVALLE', 'PAS', iocc=iocc, scal=jpas, nbret=n1)
                jnbp = int(iii/jpas)
                irest = iii-jnbp*jpas
                if (irest .ne. 0) then
                    vali(1) = jpas
                    vali(2) = iocc
                    call utmess('F', 'ALGORITH13_79', ni=2, vali=vali)
                end if
!
            else
                call getvis('INTERVALLE', 'NOMBRE', iocc=iocc, scal=jnbp, nbret=n1)
                if (jnbp .gt. 0) then
                    ipdt = int(iii/jnbp)
                    irest = iii-jnbp*ipdt
                    if (irest .ne. 0) then
                        vali(1) = jnbp
                        vali(2) = iocc
                        call utmess('F', 'ALGORITH13_80', ni=2, vali=vali)
                    end if
                end if
            end if
        end do
    end if
!
!
    if (nv .ne. 0) then
        nbvale = -nv
        ndim = max(1, nbvale-1)
        call wkvect(resu//'           .LPAS', 'G V I', ndim, jpas)
        call wkvect(resu//'           .NBPA', 'G V I', ndim, jnbp)
        call wkvect(resu//'           .BINT', 'G V I', nbvale, jbor)
        call wkvect(resu//'           .VALE', 'G V I', nbvale, jval)
        call wkvect('&&DLIDEF.VALE', 'V V I', nbvale, kval)
        call getvis(' ', 'VALE', nbval=nbvale, vect=zi(kval), nbret=nv)
        do i = 1, nbvale-1
            if (zi(kval+i-1) .ge. zi(kval+i)) then
                vali(1) = zi(kval+i-1)
                vali(2) = zi(kval+i)
                call utmess('F', 'ALGORITH13_81', ni=2, vali=vali)
            end if
            zi(jpas+i-1) = zi(kval+i)-zi(kval+i-1)
            zi(jnbp+i-1) = 1
            zi(jbor+i-1) = zi(kval+i-1)
            zi(jval+i-1) = zi(kval+i-1)
        end do
        zi(jbor+nbvale-1) = zi(kval+nbvale-1)
        zi(jval+nbvale-1) = zi(kval+nbvale-1)
!
    else
!
        call wkvect(resu//'           .LPAS', 'G V I', max(1, nbocc), jpas)
        call wkvect(resu//'           .NBPA', 'G V I', max(1, nbocc), jnbp)
        call wkvect(resu//'           .BINT', 'G V I', nbocc+1, jbor)
!
        zi(jbor) = idebut
        do iocc = 1, nbocc
            call getvis('INTERVALLE', 'JUSQU_A', iocc=iocc, scal=zi(jbor+iocc), nbret=n1)
            iii = zi(jbor+iocc)-zi(jbor-1+iocc)
            call getvis('INTERVALLE', 'PAS', iocc=iocc, nbval=0, nbret=np)
            if (np .ne. 0) then
                call getvis('INTERVALLE', 'PAS', iocc=iocc, scal=zi(jpas+iocc-1), nbret=n1)
                zi(jnbp+iocc-1) = iii/zi(jpas+iocc-1)
!
            else
                call getvis('INTERVALLE', 'NOMBRE', iocc=iocc, scal=zi(jnbp+iocc-1), nbret=n1)
                zi(jpas+iocc-1) = iii/zi(jnbp+iocc-1)
            end if
            nbval = nbval+zi(jnbp+iocc-1)
        end do
!
!        --- ALLOCATION DE .VALE ET REMPLISSAGE DE CE DERNIER ---
        call wkvect(resu//'           .VALE', 'G V I', nbval, jval)
        zi(jval) = zi(jbor)
        ico = 0
        do i = 1, nbocc
            ipdt = zi(jpas-1+i)
            do j = 1, zi(jnbp-1+i)-1
                ico = ico+1
                zi(jval+ico) = zi(jval+ico-1)+ipdt
            end do
            ico = ico+1
            zi(jval+ico) = zi(jbor+i)
        end do
    end if
!
    call jedema()
end subroutine
