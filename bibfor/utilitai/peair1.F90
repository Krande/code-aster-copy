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
subroutine peair1(mesh, nbma, lisma, aire, long)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/normev.h"
#include "asterfort/provec.h"
#include "asterfort/utmess.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "blas/ddot.h"
#include "asterfort/int_to_char8.h"
!
    integer(kind=8) :: nbma, lisma(*)
    real(kind=8) :: aire, long
    character(len=8), intent(in) :: mesh
!
!
!     CALCUL DE L'AIRE_INTERNE A UN CONTOUR
!     IN : LISMA (NBMA) : LISTE DES NUMEROS DE MAILLES DU CONTOUR FERME
!     OUT : AIRE : AIRE DELIMITEE PAR LE CONTOUR
!     OUT : LONG : LONGUEUR DU CONTOUR
!
!
    integer(kind=8) :: jma, ifm, niv, ima, numa
    integer(kind=8) :: nutyma, nbel, jdno, nbext1, nbext2, iext1
    integer(kind=8) :: iext2, ni1, ni2, nj1, nj2, nbe, nj3, nj0, jdco, i
    real(kind=8) :: orig(3), zero, vgn1(3), vn1n2(3), aire1, aire2, vgn3(3)
    real(kind=8) :: vn1n3(3)
    real(kind=8) :: xx1(3), xx2(3), xx3(3), xn(3), pv(3), xnorm, vn3n2(3)
    real(kind=8) :: vgn2(3)
    real(kind=8) :: x1, y1, z1, x2, y2, z2, xxl
    character(len=8) :: nomail, typel
    character(len=24) :: mlgcnx, mlgcoo
    character(len=24) :: valk(2)
    integer(kind=8), pointer :: mailles(:) => null()
    integer(kind=8), pointer :: noeud1(:) => null()
    integer(kind=8), pointer :: noeud2(:) => null()
    real(kind=8), pointer :: vale(:) => null()
    integer(kind=8), pointer :: typmail(:) => null()
    blas_int :: b_incx, b_incy, b_n
!
    call jemarq()
!
    call infniv(ifm, niv)
!
    zero = 0.0d0
    orig(1) = zero
    orig(2) = zero
    orig(3) = zero
!
    mlgcnx = mesh//'.CONNEX'
!
    call jeveuo(mesh//'.TYPMAIL', 'L', vi=typmail)
    call jeveuo(mesh//'.COORDO    .VALE', 'L', vr=vale)
!
    AS_ALLOCATE(vi=noeud1, size=nbma*3)
    AS_ALLOCATE(vi=noeud2, size=nbma*3)
    AS_ALLOCATE(vi=mailles, size=nbma)
!
!     VERIFICATION DU TYPE DES MAILLES ET STOCKAGE DES CONNECTIVITES
!
    long = 0.d0
    nbel = 0
    do ima = 1, nbma
        numa = lisma(ima)
        nomail = int_to_char8(numa)
!
!        TYPE DE LA MAILLE COURANTE :
!
        nutyma = typmail(numa)
        call jenuno(jexnum('&CATA.TM.NOMTM', nutyma), typel)
!
        if (typel(1:3) .ne. 'SEG') then
            nomail = int_to_char8(numa)
            valk(1) = nomail
            valk(2) = typel
            call utmess('F', 'UTILITY_1', nk=2, valk=valk)
        end if
        nbel = nbel+1
        call jeveuo(jexnum(mlgcnx, numa), 'L', jdno)
        noeud1(3*nbel-2) = zi(jdno)
        noeud1(3*nbel-1) = zi(jdno+1)
        if (typel(1:4) .eq. 'SEG3') then
            noeud1(3*nbel) = zi(jdno+2)
            x1 = vale(1+3*(zi(jdno)-1)+1-1)
            y1 = vale(1+3*(zi(jdno)-1)+2-1)
            z1 = vale(1+3*(zi(jdno)-1)+3-1)
            x2 = vale(1+3*(zi(jdno+2)-1)+1-1)
            y2 = vale(1+3*(zi(jdno+2)-1)+2-1)
            z2 = vale(1+3*(zi(jdno+2)-1)+3-1)
            xxl = (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2)
            long = long+sqrt(xxl)
            x1 = vale(1+3*(zi(jdno+1)-1)+1-1)
            y1 = vale(1+3*(zi(jdno+1)-1)+2-1)
            z1 = vale(1+3*(zi(jdno+1)-1)+3-1)
            xxl = (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2)
            long = long+sqrt(xxl)
        else
            x1 = vale(1+3*(zi(jdno)-1)+1-1)
            y1 = vale(1+3*(zi(jdno)-1)+2-1)
            z1 = vale(1+3*(zi(jdno)-1)+3-1)
            x2 = vale(1+3*(zi(jdno+1)-1)+1-1)
            y2 = vale(1+3*(zi(jdno+1)-1)+2-1)
            z2 = vale(1+3*(zi(jdno+1)-1)+3-1)
            xxl = (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2)
            long = long+sqrt(xxl)
        end if
    end do
    ASSERT(nbma .eq. nbel)
!
!     VERIFICATION QUE LE CONTOUR EST FERME
!
    nbext1 = 0
    nbext2 = 0
    do ima = 1, nbel
        iext1 = 0
        iext2 = 0
        ni1 = noeud1(3*ima-2)
        ni2 = noeud1(3*ima-1)
        mailles(ima) = ima
        do jma = 1, nbel
            if (jma .ne. ima) then
                nj1 = noeud1(3*jma-2)
                nj2 = noeud1(3*jma-1)
                if ((ni1 .eq. nj2) .or. (ni1 .eq. nj1)) iext1 = 1
                if ((ni2 .eq. nj1) .or. (ni2 .eq. nj2)) iext2 = 1
            end if
        end do
        if (iext1 .eq. 0) nbext1 = nbext1+1
        if (iext2 .eq. 0) nbext2 = nbext2+1
    end do
    if ((nbext1 .ne. 0) .and. (nbext2 .ne. 0)) then
        call utmess('F', 'UTILITY_2')
    end if
!
!     VERIFICATION QUE LE CONTOUR EST CONTINU ET REORIENTATION
!
    nbe = 1
    mailles(1) = 0
    noeud2(1) = noeud1(1)
    noeud2(1+1) = noeud1(1+1)
    noeud2(1+2) = noeud1(1+2)
41  continue
    ni1 = noeud2(3*nbe-2)
    ni2 = noeud2(3*nbe-1)
    do jma = 1, nbel
        if ((mailles(jma) .ne. 0)) then
            nj1 = noeud1(3*jma-2)
            nj2 = noeud1(3*jma-1)
            nj3 = noeud1(3*jma)
            if (ni2 .eq. nj1) then
                nbe = nbe+1
                noeud2(3*nbe-2) = nj1
                noeud2(3*nbe-1) = nj2
                if (nj3 .ne. 0) noeud2(3*nbe) = nj3
                goto 43
            else if (ni2 .eq. nj2) then
                nbe = nbe+1
                noeud2(3*nbe-2) = nj2
                noeud2(3*nbe-1) = nj1
                if (nj3 .ne. 0) noeud2(3*nbe) = nj3
                goto 43
            end if
        end if
    end do
    call utmess('F', 'UTILITY_2')
43  continue
    mailles(jma) = 0
    if (nbe .ge. nbma) then
        goto 11
    else
        goto 41
    end if
11  continue
    ASSERT(nbma .eq. nbe)
    nj2 = noeud2(3*nbe-1)
    nj0 = noeud2(1)
    ASSERT(nj2 .eq. nj0)
!
!     CALCUL DU CDG APPROXIMATIF
!
    mlgcoo = mesh//'.COORDO    .VALE'
    call jeveuo(mlgcoo, 'L', jdco)
    do ima = 1, nbma
        nj1 = noeud2(3*ima-2)
        orig(1) = orig(1)+zr(jdco-1+3*nj1-2)
        orig(2) = orig(2)+zr(jdco-1+3*nj1-1)
        orig(3) = orig(3)+zr(jdco-1+3*nj1)
    end do
    orig(1) = orig(1)/nbma
    orig(2) = orig(2)/nbma
    orig(3) = orig(3)/nbma
!
!     CALCUL DE L'AIRE GM.VECT.DL
!
    nj1 = noeud2(1)
    nj2 = noeud2(2)
!
!     CALCUL DE LA NORMALE A LA COURBE SUPPOSEE PLANE
!
    xx1(1) = zr(jdco-1+3*nj1-2)
    xx1(2) = zr(jdco-1+3*nj1-1)
    xx1(3) = zr(jdco-1+3*nj1)
    xx2(1) = zr(jdco-1+3*nj2-2)
    xx2(2) = zr(jdco-1+3*nj2-1)
    xx2(3) = zr(jdco-1+3*nj2)
    do i = 1, 3
        vgn1(i) = xx1(i)-orig(i)
        vgn2(i) = xx2(i)-orig(i)
    end do
    call provec(vgn1, vgn2, xn)
    call normev(xn, xnorm)
    aire = 0.d0
    do ima = 1, nbma
        nj1 = noeud2(3*ima-2)
        nj2 = noeud2(3*ima-1)
        nj3 = noeud2(3*ima)
        if (nj3 .eq. 0) then
            xx1(1) = zr(jdco-1+3*nj1-2)
            xx1(2) = zr(jdco-1+3*nj1-1)
            xx1(3) = zr(jdco-1+3*nj1)
            xx2(1) = zr(jdco-1+3*nj2-2)
            xx2(2) = zr(jdco-1+3*nj2-1)
            xx2(3) = zr(jdco-1+3*nj2)
            do i = 1, 3
                vgn1(i) = xx1(i)-orig(i)
                vn1n2(i) = xx2(i)-xx1(i)
            end do
            call provec(vgn1, vn1n2, pv)
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            aire1 = ddot(b_n, pv, b_incx, xn, b_incy)
            aire = aire+aire1/2.d0
        else
            xx1(1) = zr(jdco-1+3*nj1-2)
            xx1(2) = zr(jdco-1+3*nj1-1)
            xx1(3) = zr(jdco-1+3*nj1)
            xx2(1) = zr(jdco-1+3*nj2-2)
            xx2(2) = zr(jdco-1+3*nj2-1)
            xx2(3) = zr(jdco-1+3*nj2)
            xx3(1) = zr(jdco-1+3*nj3-2)
            xx3(2) = zr(jdco-1+3*nj3-1)
            xx3(3) = zr(jdco-1+3*nj3)
            do i = 1, 3
                vgn1(i) = xx1(i)-orig(i)
                vn1n3(i) = xx3(i)-xx1(i)
            end do
            call provec(vgn1, vn1n3, pv)
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            aire1 = ddot(b_n, pv, b_incx, xn, b_incy)
            aire = aire+aire1/2.d0
            do i = 1, 3
                vgn3(i) = xx3(i)-orig(i)
                vn3n2(i) = xx2(i)-xx3(i)
            end do
            call provec(vgn3, vn3n2, pv)
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            aire2 = ddot(b_n, pv, b_incx, xn, b_incy)
            aire = aire+aire2/2.d0
        end if
    end do
    AS_DEALLOCATE(vi=noeud1)
    AS_DEALLOCATE(vi=noeud2)
    AS_DEALLOCATE(vi=mailles)
    call jedema()
end subroutine
