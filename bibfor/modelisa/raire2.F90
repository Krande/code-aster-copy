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
subroutine raire2(noma, rigi, nbgr, ligrma, nbnoeu, &
                  nbno, tabnoe, rignoe)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/compma.h"
#include "asterfort/fointe.h"
#include "asterfort/getvem.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/iunifi.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/provec.h"
#include "asterfort/utmess.h"
#include "asterfort/char8_to_int.h"
#include "blas/ddot.h"
!
    integer(kind=8) :: nbgr, nbno, nbnoeu, tabnoe(nbnoeu)
    character(len=8) :: noma
    character(len=24) :: ligrma(nbgr)
    real(kind=8) :: rignoe(6*nbnoeu)
!
    character(len=8) :: k8b
    character(len=8) :: nomnoe
    character(len=24) :: nomgr, magrno, magrma, manoma
    real(kind=8) :: zero, x(8), y(8), z(8), rigi(6)
    real(kind=8) :: a(3), b(3), c(3), u(3)
    aster_logical :: lfonc
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ifr, ii
    integer(kind=8) :: ij, im, in, inoe, iret
    integer(kind=8) :: ldgm, ldgn, ldnm, nb, nbma, ncf
    integer(kind=8) :: ncg, nfg, ngn, nm, nn, nno, noemax
!
    real(kind=8) :: coef, dist, hc, r1, r2, r3
    real(kind=8) :: r4, r5, r6, rig4, rig45, rig46, rig5
    real(kind=8) :: rig56, rig6, surf, surtot, xc, xg, xx
    real(kind=8) :: yc, yg, yy, zg, zz
    real(kind=8), pointer :: coegro(:) => null()
    real(kind=8), pointer :: coeno(:) => null()
    character(len=8), pointer :: fongro(:) => null()
    integer(kind=8), pointer :: parno(:) => null()
    real(kind=8), pointer :: surmai(:) => null()
    real(kind=8), pointer :: vale(:) => null()
    blas_int :: b_incx, b_incy, b_n
!-----------------------------------------------------------------------
    call jemarq()
    zero = 0.d0
    ifr = iunifi('RESULTAT')
    lfonc = .false.
!
!
!     --- ON RECUPERE LES POINTS D'ANCRAGE ---
!
!
!
!        --- ON ECLATE LE GROUP_NO EN NOEUDS ---
    call compma(noma, nbgr, ligrma, nbma)
    magrno = noma//'.GROUPENO'
    magrma = noma//'.GROUPEMA'
    manoma = noma//'.CONNEX'
    noemax = 0
!
!     --- DESCRIPTION NOEUDS STRUCTURE ---
    call jeveuo(noma//'.COORDO    .VALE', 'L', vr=vale)
!
!       RECUPERATION DU CENTRE
!
    xg = zero
    yg = zero
    zg = zero
    call getvr8('ENER_SOL', 'COOR_CENTRE', iocc=1, nbval=0, nbret=ncg)
    call getvem(noma, 'NOEUD', 'ENER_SOL', 'NOEUD_CENTRE', 1, &
                0, k8b, nno)
    call getvem(noma, 'GROUP_NO', 'ENER_SOL', 'GROUP_NO_CENTRE', 1, &
                0, k8b, ngn)
    if (ncg .ne. 0) then
        call getvr8('ENER_SOL', 'COOR_CENTRE', iocc=1, nbval=3, vect=c, &
                    nbret=ncg)
        xg = c(1)
        yg = c(2)
        zg = c(3)
    else if (nno .ne. 0) then
        call getvem(noma, 'NOEUD', 'ENER_SOL', 'NOEUD_CENTRE', 1, &
                    1, nomnoe, nno)
        inoe = char8_to_int(nomnoe)
        xg = vale(1+3*(inoe-1)+1-1)
        yg = vale(1+3*(inoe-1)+2-1)
        zg = vale(1+3*(inoe-1)+3-1)
    else if (ngn .ne. 0) then
        call getvem(noma, 'GROUP_NO', 'ENER_SOL', 'GROUP_NO_CENTRE', 1, &
                    1, nomgr, ngn)
        call jeveuo(jexnom(magrno, nomgr), 'L', ldgn)
        inoe = zi(ldgn)
        xg = vale(1+3*(inoe-1)+1-1)
        yg = vale(1+3*(inoe-1)+2-1)
        zg = vale(1+3*(inoe-1)+3-1)
    end if
!
!       RECUPERATION DES COEFS OU FONCTIONS DE GROUPE
!
    call getvr8('ENER_SOL', 'COEF_GROUP', iocc=1, nbval=0, nbret=ncg)
    if (ncg .ne. 0) then
        AS_ALLOCATE(vr=coegro, size=nbgr)
        call getvr8('ENER_SOL', 'COEF_GROUP', iocc=1, nbval=nbgr, vect=coegro, &
                    nbret=ncg)
    else
        call getvid('ENER_SOL', 'FONC_GROUP', iocc=1, nbval=0, nbret=ncf)
        if (ncf .eq. 0) then
            call utmess('F', 'MODELISA6_33')
        end if
        AS_ALLOCATE(vk8=fongro, size=nbgr)
        lfonc = .true.
        call getvid('ENER_SOL', 'FONC_GROUP', iocc=1, nbval=nbgr, vect=fongro, &
                    nbret=nfg)
    end if
!
    do i = 1, nbgr
        call jelira(jexnom(magrma, ligrma(i)), 'LONUTI', nb)
        call jeveuo(jexnom(magrma, ligrma(i)), 'L', ldgm)
        do in = 0, nb-1
            call jelira(jexnum(manoma, zi(ldgm+in)), 'LONMAX', nm)
            call jeveuo(jexnum(manoma, zi(ldgm+in)), 'L', ldnm)
            do nn = 1, nm
                inoe = zi(ldnm+nn-1)
                noemax = max(noemax, inoe)
            end do
        end do
    end do
    AS_ALLOCATE(vr=coeno, size=noemax)
!
!        TABLEAU DE PARTICIPATION DES NOEUDS DE L INTERFACE
!
    AS_ALLOCATE(vi=parno, size=noemax)
!
!
!     CALCUL DES SURFACES ELEMENTAIRES ET DE LA SURFACE TOTALE
!
    AS_ALLOCATE(vr=surmai, size=nbma)
    im = 0
    surtot = zero
    do i = 1, nbgr
        call jelira(jexnom(magrma, ligrma(i)), 'LONUTI', nb)
        call jeveuo(jexnom(magrma, ligrma(i)), 'L', ldgm)
        do in = 0, nb-1
            im = im+1
            call jelira(jexnum(manoma, zi(ldgm+in)), 'LONMAX', nm)
            call jeveuo(jexnum(manoma, zi(ldgm+in)), 'L', ldnm)
            xc = zero
            yc = zero
            hc = zero
            do nn = 1, nm
                inoe = zi(ldnm+nn-1)
                parno(inoe) = parno(inoe)+1
                x(nn) = vale(1+3*(inoe-1)+1-1)
                y(nn) = vale(1+3*(inoe-1)+2-1)
                z(nn) = vale(1+3*(inoe-1)+3-1)
                xc = xc+x(nn)
                yc = yc+y(nn)
                hc = hc+z(nn)
            end do
            xc = xc/nm
            yc = yc/nm
            hc = hc/nm
            a(1) = x(3)-x(1)
            a(2) = y(3)-y(1)
            a(3) = z(3)-z(1)
            if (nm .eq. 3 .or. nm .eq. 6) then
                b(1) = x(2)-x(1)
                b(2) = y(2)-y(1)
                b(3) = z(2)-z(1)
            else if (nm .eq. 4 .or. nm .eq. 8) then
                b(1) = x(4)-x(2)
                b(2) = y(4)-y(2)
                b(3) = z(4)-z(2)
            else
                call utmess('F', 'MODELISA6_34')
            end if
            call provec(a, b, c)
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            surf = ddot(b_n, c, b_incx, c, b_incy)
            surmai(im) = sqrt(surf)*0.5d0
            if (lfonc) then
                u(1) = xg-xc
                u(2) = yg-yc
                u(3) = zg-hc
                b_n = to_blas_int(3)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                dist = ddot(b_n, u, b_incx, u, b_incy)
                dist = sqrt(dist)
                call fointe('F ', fongro(i), 1, ['X'], [dist], &
                            coef, iret)
                surmai(im) = surmai(im)*coef
            else
                surmai(im) = surmai(im)*coegro(i)
            end if
            surtot = surtot+surmai(im)
            surmai(im) = surmai(im)/nm
        end do
    end do
!
!     CALCUL DES PONDERATIONS ELEMENTAIRES
!
    im = 0
    do i = 1, nbgr
        call jelira(jexnom(magrma, ligrma(i)), 'LONUTI', nb)
        call jeveuo(jexnom(magrma, ligrma(i)), 'L', ldgm)
        do in = 0, nb-1
            im = im+1
            call jelira(jexnum(manoma, zi(ldgm+in)), 'LONMAX', nm)
            call jeveuo(jexnum(manoma, zi(ldgm+in)), 'L', ldnm)
            do nn = 1, nm
                do ij = 1, noemax
                    if (parno(ij) .eq. 0) goto 37
                    if (zi(ldnm+nn-1) .eq. ij) then
                        coeno(ij) = coeno(ij)+surmai(im)/surtot
                    end if
37                  continue
                end do
            end do
        end do
    end do
    nbma = im
!
!     CALCUL DES RAIDEURS DE TORSION
!
    ii = 0
    rig4 = zero
    rig5 = zero
    rig6 = zero
    rig45 = zero
    rig46 = zero
    rig56 = zero
    do ij = 1, noemax
        if (parno(ij) .eq. 0) goto 50
        ii = ii+1
        xx = vale(1+3*(ij-1)+1-1)-xg
        yy = vale(1+3*(ij-1)+2-1)-yg
        zz = vale(1+3*(ij-1)+3-1)-zg
        rig4 = rig4+(rigi(2)*zz**2+rigi(3)*yy**2)*coeno(ij)
        rig5 = rig5+(rigi(1)*zz**2+rigi(3)*xx**2)*coeno(ij)
        rig6 = rig6+(rigi(2)*xx**2+rigi(1)*yy**2)*coeno(ij)
        rig45 = rig45-rigi(3)*xx*yy*coeno(ij)
        rig46 = rig46-rigi(2)*xx*zz*coeno(ij)
        rig56 = rig56-rigi(1)*yy*zz*coeno(ij)
50      continue
    end do
    nbno = ii
    rig4 = rigi(4)-rig4
    rig5 = rigi(5)-rig5
    rig6 = rigi(6)-rig6
    write (ifr, 1001) rig4, rig5, rig6
!
    ii = 0
    do ij = 1, noemax
        if (parno(ij) .eq. 0) goto 51
        ii = ii+1
        r1 = rigi(1)*coeno(ij)
        r2 = rigi(2)*coeno(ij)
        r3 = rigi(3)*coeno(ij)
        r4 = rig4*coeno(ij)
        r5 = rig5*coeno(ij)
        r6 = rig6*coeno(ij)
        rignoe(6*(ii-1)+1) = r1
        rignoe(6*(ii-1)+2) = r2
        rignoe(6*(ii-1)+3) = r3
        rignoe(6*(ii-1)+4) = r4
        rignoe(6*(ii-1)+5) = r5
        rignoe(6*(ii-1)+6) = r6
        tabnoe(ii) = ij
51      continue
    end do
!
1001 format(1x, 'RAIDEURS DE ROTATION A REPARTIR:', /&
&      1x, ' KRX: ', 1x, 1pe12.5, ' KRY: ', 1x, 1pe12.5,&
&      ' KRZ: ', 1x, 1pe12.5)
    AS_DEALLOCATE(vr=coegro)
    AS_DEALLOCATE(vk8=fongro)
    AS_DEALLOCATE(vr=coeno)
    AS_DEALLOCATE(vi=parno)
    AS_DEALLOCATE(vr=surmai)
!
    call jedema()
end subroutine
