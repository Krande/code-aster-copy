! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
subroutine rairep(noma, ioc, km, rigi, nbgr, &
                  ligrma, nbno, tabnoe, rignoe, rigto, &
                  amoto, rirot, ndim)
!
!
    implicit none
    integer :: ioc, nbgr, nbno, ndim
    character(len=8) :: noma, tabnoe(*), km
    character(len=24) :: ligrma(nbgr)
    real(kind=8) :: rignoe(*), rigto(*), amoto(*), rirot(3)
!
! --------------------------------------------------------------------------------------------------
! person_in_charge: jean-luc.flejou at edf.fr
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/compma.h"
#include "asterfort/dismoi.h"
#include "asterfort/fointe.h"
#include "asterfort/getvem.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/provec.h"
#include "asterfort/utmess.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "blas/ddot.h"
!
! --------------------------------------------------------------------------------------------------
    character(len=8) :: k8b
    character(len=8) :: nomnoe, typm
    character(len=24) :: nomgr, magrno, manono, magrma, manoma, matyma
    real(kind=8) :: zero, x(9), y(9), z(9), rigi(6)
    real(kind=8) :: a(3), b(3), c(3), u(3)
    aster_logical :: lfonc, trans
    integer :: appui
! --------------------------------------------------------------------------------------------------
    integer :: i, ii
    integer :: ij, im, in, inoe, iret
    integer :: ldgm, ldgn, ldnm, ltyp, nb, nbma, ncg
    integer :: nfg, ngn, nm, nn, noemax, ntopo
    integer :: numa
    real(kind=8) :: coef, dist, hc, r1, r2, r3
    real(kind=8) :: r4, r5, r6, rig3, rig4, rig45, rig46
    real(kind=8) :: rig5, rig56, rig6, surf, surtot, xc, xg
    real(kind=8) :: xx, yc, yg, yy, zg, zz
    real(kind=8), pointer :: coegro(:) => null()
    real(kind=8), pointer :: coeno(:) => null()
    character(len=8), pointer :: fongro(:) => null()
    integer, pointer :: parno(:) => null()
    real(kind=8), pointer :: surmai(:) => null()
    real(kind=8), pointer :: vale(:) => null()
    blas_int :: b_incx, b_incy, b_n
! --------------------------------------------------------------------------------------------------
    call jemarq()
    zero = 0.d0
    lfonc = .false.
!
!   Récupere les points d'ancrage ---
!   eclate le GROUP_NO en noeuds
    call compma(noma, nbgr, ligrma, nbma)
    magrno = noma//'.GROUPENO'
    manono = noma//'.NOMNOE'
    magrma = noma//'.GROUPEMA'
    manoma = noma//'.CONNEX'
    matyma = noma//'.TYPMAIL'
!
    noemax = 0
!   Description noeuds structure ---
    call jeveuo(noma//'.COORDO    .VALE', 'L', vr=vale)
!
!   Récupération du centre
    call getvr8('RIGI_PARASOL', 'COOR_CENTRE', iocc=ioc, nbval=0, nbret=ncg)
    call getvem(noma, 'GROUP_NO', 'RIGI_PARASOL', 'GROUP_NO_CENTRE', ioc, &
                0, k8b, ngn)
    xg = 0.0
    yg = 0.0
    zg = 0.0
    if (ncg .ne. 0) then
        call getvr8('RIGI_PARASOL', 'COOR_CENTRE', iocc=ioc, nbval=3, vect=c, &
                    nbret=ncg)
        xg = c(1)
        yg = c(2)
        zg = c(3)
    else if (ngn .ne. 0) then
        call getvem(noma, 'GROUP_NO', 'RIGI_PARASOL', 'GROUP_NO_CENTRE', ioc, &
                    1, nomgr, ngn)
        call jeveuo(jexnom(magrno, nomgr), 'L', ldgn)
        inoe = zi(ldgn)
        call jenuno(jexnum(manono, inoe), nomnoe)
        xg = vale(1+3*(inoe-1)+1-1)
        yg = vale(1+3*(inoe-1)+2-1)
        zg = vale(1+3*(inoe-1)+3-1)
    else
        ASSERT(ASTER_FALSE)
    end if
!
!   Récuperation des coefs ou fonctions de groupe
    call getvr8('RIGI_PARASOL', 'COEF_GROUP', iocc=ioc, nbval=0, nbret=ncg)
    if (ncg .ne. 0) then
        AS_ALLOCATE(vr=coegro, size=nbgr)
        call getvr8('RIGI_PARASOL', 'COEF_GROUP', iocc=ioc, nbval=nbgr, vect=coegro, &
                    nbret=ncg)
    else
        AS_ALLOCATE(vk8=fongro, size=nbgr)
        lfonc = .true.
        call getvid('RIGI_PARASOL', 'FONC_GROUP', iocc=ioc, nbval=nbgr, vect=fongro, &
                    nbret=nfg)
    end if
!
    if (ndim .eq. 2) then
        appui = 1
    else
!       la dimension de l'appui n'est pas encore determinée
        appui = -1
    end if
!
    call jeveuo(matyma, 'L', ltyp)
    do i = 1, nbgr
        call jelira(jexnom(magrma, ligrma(i)), 'LONUTI', nb)
        call jeveuo(jexnom(magrma, ligrma(i)), 'L', ldgm)
        do in = 0, nb-1
            numa = zi(ldgm+in)
            ASSERT(numa .gt. 0)
            call jelira(jexnum(manoma, numa), 'LONMAX', nm)
            call jeveuo(jexnum(manoma, numa), 'L', ldnm)
            call jenuno(jexnum('&CATA.TM.NOMTM', zi(ltyp-1+numa)), typm)
            call dismoi('DIM_TOPO', typm, 'TYPE_MAILLE', repi=ntopo)
!
            if (appui .eq. -1) then
!               la dimension de la première maille définit l'appui
                appui = ntopo
            else if ((appui .eq. 1) .or. (appui .eq. 2)) then
                if (appui .ne. ntopo) then
                    call utmess('F', 'MODELISA6_35')
                end if
            else
                call utmess('F', 'MODELISA6_29')
            end if
            do nn = 1, nm
                inoe = zi(ldnm+nn-1)
                noemax = max(noemax, inoe)
            end do
        end do
    end do
    ASSERT(appui .ne. -1)
!
    AS_ALLOCATE(vr=coeno, size=noemax)
!   Tableau de participation des noeuds de l interface
    AS_ALLOCATE(vi=parno, size=noemax)
!   Calcul des surfaces élémentaires et de la surface totale
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
!
            if (appui .eq. 1) then
                a(1) = x(2)-x(1)
                a(2) = y(2)-y(1)
                a(3) = z(2)-z(1)
                b_n = to_blas_int(2)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                surf = ddot(b_n, a, b_incx, a, b_incy)
                surmai(im) = sqrt(surf)
            else if (appui .eq. 2) then
                a(1) = x(3)-x(1)
                a(2) = y(3)-y(1)
                a(3) = z(3)-z(1)
                if (nm .eq. 3 .or. nm .eq. 6 .or. nm .eq. 7) then
                    b(1) = x(2)-x(1)
                    b(2) = y(2)-y(1)
                    b(3) = z(2)-z(1)
                else if (nm .eq. 4 .or. nm .eq. 8 .or. nm .eq. 9) then
                    b(1) = x(4)-x(2)
                    b(2) = y(4)-y(2)
                    b(3) = z(4)-z(2)
                else
                    ASSERT(.false.)
                end if
                call provec(a, b, c)
                b_n = to_blas_int(3)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                surf = ddot(b_n, c, b_incx, c, b_incy)
                surmai(im) = sqrt(surf)*0.5d0
            else
                ASSERT(.false.)
            end if
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
!   Calcul des pondérations élémentaires
    im = 0
    do i = 1, nbgr
        call jelira(jexnom(magrma, ligrma(i)), 'LONUTI', nb)
        call jeveuo(jexnom(magrma, ligrma(i)), 'L', ldgm)
        do in = 0, nb-1
            im = im+1
            call jelira(jexnum(manoma, zi(ldgm+in)), 'LONMAX', nm)
            call jeveuo(jexnum(manoma, zi(ldgm+in)), 'L', ldnm)
            do nn = 1, nm
                cij1: do ij = 1, noemax
                    if (parno(ij) .eq. 0) cycle cij1
                    if (zi(ldnm+nn-1) .eq. ij) then
                        coeno(ij) = coeno(ij)+surmai(im)/surtot
                    end if
                end do cij1
            end do
        end do
    end do
    nbma = im
!
!   Calcul des raideurs de rotation
    ii = 0
    rig4 = zero
    rig5 = zero
    rig6 = zero
    rig45 = zero
    rig46 = zero
    rig56 = zero
    rig3 = 0.d0
    cij2: do ij = 1, noemax
        if (parno(ij) .eq. 0) cycle cij2
        ii = ii+1
        xx = vale(1+3*(ij-1)+1-1)-xg
        yy = vale(1+3*(ij-1)+2-1)-yg
        zz = vale(1+3*(ij-1)+3-1)-zg
        if (ndim .eq. 3) then
            rig4 = rig4+(rigi(2)*zz**2+rigi(3)*yy**2)*coeno(ij)
            rig5 = rig5+(rigi(1)*zz**2+rigi(3)*xx**2)*coeno(ij)
            rig6 = rig6+(rigi(2)*xx**2+rigi(1)*yy**2)*coeno(ij)
            rig45 = rig45-rigi(3)*xx*yy*coeno(ij)
            rig46 = rig46-rigi(2)*xx*zz*coeno(ij)
            rig56 = rig56-rigi(1)*yy*zz*coeno(ij)
            rig3 = 0.d0
        else
            rig3 = rig3+(rigi(2)*xx**2+rigi(1)*yy**2)*coeno(ij)
        end if
    end do cij2
    nbno = ii
!
    trans = (km(1:7) .eq. 'K_T_D_N') .or. (km(1:7) .eq. 'K_T_D_L') .or. (km(1:7) .eq. 'A_T_D_N') &
            .or. (km(1:7) .eq. 'A_T_D_L')
!
    if (trans) then
!       pas de raideur en rotation sur les discrets
        if (ndim .eq. 2) then
            rigi(3) = zero
            rig3 = zero
        end if
        rigi(4) = zero
        rigi(5) = zero
        rigi(6) = zero
        rig4 = zero
        rig5 = zero
        rig6 = zero
        rirot(1) = zero
        rirot(2) = zero
        rirot(3) = zero
    else
        rig3 = rigi(3)-rig3
        rig4 = rigi(4)-rig4
        rig5 = rigi(5)-rig5
        rig6 = rigi(6)-rig6
        if (ndim .eq. 3) then
            rirot(1) = rig4
            rirot(2) = rig5
            rirot(3) = rig6
        else
            rirot(1) = rig3
        end if
    end if
!
    ii = 0
    cij3: do ij = 1, noemax
        if (parno(ij) .eq. 0) cycle cij3
        ii = ii+1
        r1 = rigi(1)*coeno(ij)
        r2 = rigi(2)*coeno(ij)
        if (ndim .eq. 3) then
            r3 = rigi(3)*coeno(ij)
            r4 = rig4*coeno(ij)
            r5 = rig5*coeno(ij)
            r6 = rig6*coeno(ij)
        else
            r3 = rig3*coeno(ij)
            r4 = zero
            r5 = zero
            r6 = zero
        end if
        call jenuno(jexnum(manono, ij), nomnoe)
        if (km(1:1) .eq. 'K') then
            rigto(6*(ij-1)+1) = r1+rigto(6*(ij-1)+1)
            rigto(6*(ij-1)+2) = r2+rigto(6*(ij-1)+2)
            rigto(6*(ij-1)+3) = r3+rigto(6*(ij-1)+3)
            rigto(6*(ij-1)+4) = r4+rigto(6*(ij-1)+4)
            rigto(6*(ij-1)+5) = r5+rigto(6*(ij-1)+5)
            rigto(6*(ij-1)+6) = r6+rigto(6*(ij-1)+6)
            r1 = rigto(6*(ij-1)+1)
            r2 = rigto(6*(ij-1)+2)
            r3 = rigto(6*(ij-1)+3)
            r4 = rigto(6*(ij-1)+4)
            r5 = rigto(6*(ij-1)+5)
            r6 = rigto(6*(ij-1)+6)
        else if (km(1:1) .eq. 'A') then
            amoto(6*(ij-1)+1) = r1+amoto(6*(ij-1)+1)
            amoto(6*(ij-1)+2) = r2+amoto(6*(ij-1)+2)
            amoto(6*(ij-1)+3) = r3+amoto(6*(ij-1)+3)
            amoto(6*(ij-1)+4) = r4+amoto(6*(ij-1)+4)
            amoto(6*(ij-1)+5) = r5+amoto(6*(ij-1)+5)
            amoto(6*(ij-1)+6) = r6+amoto(6*(ij-1)+6)
            r1 = amoto(6*(ij-1)+1)
            r2 = amoto(6*(ij-1)+2)
            r3 = amoto(6*(ij-1)+3)
            r4 = amoto(6*(ij-1)+4)
            r5 = amoto(6*(ij-1)+5)
            r6 = amoto(6*(ij-1)+6)
        end if
        rignoe(6*(ii-1)+1) = r1
        rignoe(6*(ii-1)+2) = r2
        rignoe(6*(ii-1)+3) = r3
        rignoe(6*(ii-1)+4) = r4
        rignoe(6*(ii-1)+5) = r5
        rignoe(6*(ii-1)+6) = r6
        tabnoe(ii) = nomnoe
    end do cij3
!
    AS_DEALLOCATE(vr=coegro)
    AS_DEALLOCATE(vk8=fongro)
    AS_DEALLOCATE(vr=coeno)
    AS_DEALLOCATE(vi=parno)
    AS_DEALLOCATE(vr=surmai)
!
    call jedema()
end subroutine
