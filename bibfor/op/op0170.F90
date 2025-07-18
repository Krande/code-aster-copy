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
subroutine op0170()
    implicit none
!
!     CALCUL FATIGUE ALEATOIRE
!
!     -----------------------------------------------------------------
!
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterc/r8vide.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/infmaj.h"
#include "asterfort/jeveuo.h"
#include "asterfort/pdadom.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
#include "asterfort/tbexp2.h"
#include "asterfort/tbexve.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: ibid, nbtab, nbmom, n1, nbpfat, ivmom, i, ilign, nbl0, inbl0
    integer(kind=8) :: nbl2, inbl2, nbl4, inbl4
    parameter(nbpfat=5)
    real(kind=8) :: xm0, xm2, xm4, rundf, dom, rduree, valer(5)
    complex(kind=8) :: c16b
    character(len=8) :: k8b, nomres, table, typfat(nbpfat)
    character(len=16) :: nomcmd, concep, nopfat(nbpfat), nopfa2(4)
    character(len=24) :: nomob1, nomob2, nomob3
!     -----------------------------------------------------------------
    data nopfat/'MOMENT_SPEC_0', 'MOMENT_SPEC_2',&
     &               'MOMENT_SPEC_4',&
     &               'DUREE', 'DOMMAGE'/
    data nopfa2/'MOMENT_SPEC_0', 'MOMENT_SPEC_2',&
     &               'DUREE', 'DOMMAGE'/
    data typfat/'R', 'R', 'R', 'R', 'R'/
!     -----------------------------------------------------------------
!
    call infmaj()
    c16b = (0.d0, 0.d0)
    ibid = 0
    rundf = r8vide()
    xm4 = rundf
    ivmom = 0
!
    call getres(nomres, concep, nomcmd)
!
    call getvr8(' ', 'DUREE', scal=rduree, nbret=n1)
!
    call getvid(' ', 'TABL_POST_ALEA', nbval=0, nbret=nbtab)
!
    if (nbtab .ne. 0) then
        call getvid(' ', 'TABL_POST_ALEA', scal=table, nbret=n1)
!        CALL TBEXP2(TABLE,'GRANDEUR')
        call tbexp2(table, 'LAMBDA_00')
        call tbexp2(table, 'LAMBDA_02')
        call tbexp2(table, 'LAMBDA_04')
        nomob1 = '&&OP0170.LAMBDA_0'
        call tbexve(table, 'LAMBDA_00', nomob1, 'V', nbl0, &
                    k8b)
        call jeveuo(nomob1, 'L', inbl0)
        nomob2 = '&&OP0170.LAMBDA_2'
        call tbexve(table, 'LAMBDA_02', nomob2, 'V', nbl2, &
                    k8b)
        if (nbl2 .ne. nbl0) then
            call utmess('F', 'MODELISA2_89')
        end if
        call jeveuo(nomob2, 'L', inbl2)
        nomob3 = '&&OP0170.LAMBDA_4'
        call tbexve(table, 'LAMBDA_04', nomob3, 'V', nbl4, &
                    k8b)
        if (nbl4 .ne. nbl0) then
            call utmess('F', 'ALGELINE_7')
        end if
        call jeveuo(nomob3, 'L', inbl4)
        nbmom = nbl0
        call wkvect('&&OP0170.MOMENT', 'V V R', 3*nbmom, ivmom)
        do i = 1, nbl0
            zr(ivmom+(i-1)*3) = zr(inbl0+i-1)
            zr(ivmom+(i-1)*3+1) = zr(inbl2+i-1)
            zr(ivmom+(i-1)*3+2) = zr(inbl4+i-1)
        end do
!
    else
!
        call getvr8(' ', 'MOMENT_SPEC_0', scal=xm0, nbret=n1)
        call getvr8(' ', 'MOMENT_SPEC_2', scal=xm2, nbret=n1)
        call getvr8(' ', 'MOMENT_SPEC_4', scal=xm4, nbret=n1)
        nbmom = 1
        call wkvect('&&OP0170.MOMENT', 'V V R', 3*nbmom, ivmom)
        zr(ivmom) = xm0
        zr(ivmom+1) = xm2
        zr(ivmom+2) = xm4
!
    end if
!
    if (nbmom .eq. 0) then
        call utmess('A', 'PREPOST4_17')
    end if
!
    call tbcrsd(nomres, 'G')
    call tbajpa(nomres, nbpfat, nopfat, typfat)
!
    ilign = 0
    do i = 1, nbmom
        xm0 = zr(ivmom+(i-1)*3)
        xm2 = zr(ivmom+(i-1)*3+1)
        xm4 = zr(ivmom+(i-1)*3+2)
        call pdadom(xm0, xm2, xm4, dom)
        dom = dom*rduree
        if (xm4 .eq. rundf) then
            valer(1) = xm0
            valer(2) = xm2
            valer(3) = rduree
            valer(4) = dom
            call tbajli(nomres, 4, nopfa2, [ibid], valer, &
                        [c16b], k8b, ilign)
        else
            valer(1) = xm0
            valer(2) = xm2
            valer(3) = xm4
            valer(4) = rduree
            valer(5) = dom
            call tbajli(nomres, nbpfat, nopfat, [ibid], valer, &
                        [c16b], k8b, ilign)
        end if
    end do
!
    call titre()
!
!
end subroutine
