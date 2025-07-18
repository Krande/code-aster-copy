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
subroutine pecag3(ndim, nsymx, nsymy, noma, motcle, &
                  nbmail, noment, valpar)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8dgrd.h"
#include "asterc/r8maem.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/char8_to_int.h"
!
    integer(kind=8) :: ndim, nbmail
    real(kind=8) :: valpar(*)
    character(len=*) :: noment(*), noma
    character(len=*) :: motcle
    aster_logical :: nsymx, nsymy
!     OPERATEUR   POST_ELEM
!     TRAITEMENT DU MOT CLE-FACTEUR "CARA_GEOM"
!     ------------------------------------------------------------------
!
!
    character(len=8) :: noma8
    character(len=24) :: mlggma, mlgval, mlgcox
!     ------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ibid, ig, im, in, jcoor, jdes
    integer(kind=8) :: jgro, nbma, nbno, nbnoeu, numail, nuno
    real(kind=8) :: alpha, cdx, cdy, cosa, r
    real(kind=8) :: rmax, rx, ry, sina, tamp, x, x0
    real(kind=8) :: xmax, xmin, y, y0, ymax, ymin, zmax
    real(kind=8) :: zmin
!-----------------------------------------------------------------------
    call jemarq()
    noma8 = noma
    mlggma = noma8//'.GROUPEMA'
    mlgcox = noma8//'.CONNEX'
    mlgval = noma8//'.COORDO    .VALE'
    call jeveuo(mlgval, 'L', jcoor)
    call jelira(mlgval, 'LONMAX', nbnoeu)
    nbnoeu = nbnoeu/3
!
    if (ndim .eq. 2) then
        cdx = valpar(13)
        cdy = valpar(14)
        alpha = r8dgrd()*valpar(20)
        cosa = cos(alpha)
        sina = sin(alpha)
    else
        call utmess('F', 'UTILITAI3_48')
        cdx = valpar(19)
        cdy = valpar(20)
    end if
    xmax = -r8maem()
    xmin = r8maem()
    ymax = -r8maem()
    ymin = r8maem()
    zmax = -r8maem()
    zmin = r8maem()
    rmax = -r8maem()
!
!
    if (motcle(1:4) .eq. 'TOUT') then
        do i = 1, nbnoeu
            x0 = zr(jcoor-1+3*(i-1)+1)-cdx
            y0 = zr(jcoor-1+3*(i-1)+2)-cdy
            x = x0*cosa+y0*sina
            y = y0*cosa-x0*sina
            r = sqrt(x*x+y*y)
            xmax = max(xmax, x)
            xmin = min(xmin, x)
            ymax = max(ymax, y)
            ymin = min(ymin, y)
            zmax = 0.d0
            zmin = 0.d0
            rmax = max(rmax, r)
        end do
!
    else if (motcle(1:6) .eq. 'MAILLE') then
        do im = 1, nbmail
            ibid = char8_to_int(noment(im))
            call jeveuo(jexnum(mlgcox, ibid), 'L', jdes)
            call jelira(jexnum(mlgcox, ibid), 'LONMAX', nbno)
            do in = 1, nbno
                nuno = zi(jdes+in-1)
                x0 = zr(jcoor-1+3*(nuno-1)+1)-cdx
                y0 = zr(jcoor-1+3*(nuno-1)+2)-cdy
                x = x0*cosa+y0*sina
                y = y0*cosa-x0*sina
                r = sqrt(x*x+y*y)
                xmax = max(xmax, x)
                xmin = min(xmin, x)
                ymax = max(ymax, y)
                ymin = min(ymin, y)
                zmax = 0.d0
                zmin = 0.d0
                rmax = max(rmax, r)
            end do
        end do
!
    else if (motcle(1:8) .eq. 'GROUP_MA') then
        do ig = 1, nbmail
            call jeveuo(jexnom(mlggma, noment(ig)), 'L', jgro)
            call jelira(jexnom(mlggma, noment(ig)), 'LONUTI', nbma)
            do im = 1, nbma
                numail = zi(jgro+im-1)
                call jeveuo(jexnum(mlgcox, numail), 'L', jdes)
                call jelira(jexnum(mlgcox, numail), 'LONMAX', nbno)
                do in = 1, nbno
                    nuno = zi(jdes+in-1)
                    x0 = zr(jcoor-1+3*(nuno-1)+1)-cdx
                    y0 = zr(jcoor-1+3*(nuno-1)+2)-cdy
                    x = x0*cosa+y0*sina
                    y = y0*cosa-x0*sina
                    r = sqrt(x*x+y*y)
                    xmax = max(xmax, x)
                    xmin = min(xmin, x)
                    ymax = max(ymax, y)
                    ymin = min(ymin, y)
                    zmax = 0.d0
                    zmin = 0.d0
                    rmax = max(rmax, r)
                end do
            end do
        end do
    end if
    rx = max(abs(xmax), (abs(xmin)))
    ry = max(abs(ymax), (abs(ymin)))
!
    if (nsymx) then
        x0 = 1.d0
        y0 = 0.d0
        x = x0*cosa+y0*sina
        y = y0*cosa-x0*sina
        if (abs(abs(x)-1.d0) .le. 1.d-5) then
            tamp = max(abs(ymin), abs(ymax))
            ymin = -tamp
            ymax = tamp
        else
            tamp = max(abs(xmin), abs(xmax))
            xmin = -tamp
            xmax = tamp
        end if
    end if
    if (nsymy) then
        x0 = 0.d0
        y0 = 1.d0
        x = x0*cosa+y0*sina
        y = y0*cosa-x0*sina
        if (abs(abs(y)-1.d0) .le. 1.d-5) then
            tamp = max(abs(xmin), abs(xmax))
            xmin = -tamp
            xmax = tamp
        else
            tamp = max(abs(ymin), abs(ymax))
            ymin = -tamp
            ymax = tamp
        end if
    end if
    if (ndim .eq. 2) then
        valpar(7) = xmax
        valpar(8) = ymax
        valpar(9) = xmin
        valpar(10) = ymin
        valpar(11) = rmax
        valpar(42) = rx
        valpar(43) = ry
    else
        valpar(11) = xmax
        valpar(12) = ymax
        valpar(13) = zmax
        valpar(14) = xmin
        valpar(15) = ymin
        valpar(16) = zmin
        valpar(17) = rmax
    end if
!
    call jedema()
!
end subroutine
