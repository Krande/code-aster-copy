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
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/jeexin.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/normev.h"
#include "asterfort/provec.h"
#include "asterfort/utmess.h"
#include "asterfort/vector_ghosts_comm.h"
#include "asterfort/wkvect.h"
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
    integer(kind=8) :: jma, ifm, niv, ima, numa, ino
    integer(kind=8) :: nutyma, nbel, jdno, nbext1, nbext2, iext1
    integer(kind=8) :: iext2, ni1, ni2, nj1, nj2, nbe, nj3, nj0, jdco, i
    integer(kind=8) :: rang, nbproc, nbma_int, no_1, no_2, ma_ext, iret
    integer(kind=8) :: nbma_tot, node_nb, nbma_ext
    real(kind=8) :: orig(3), zero, vgn1(3), vn1n2(3), aire1, aire2, vgn3(3)
    real(kind=8) :: vn1n3(3)
    real(kind=8) :: xx1(3), xx2(3), xx3(3), xn(3), pv(3), xnorm, vn3n2(3)
    real(kind=8) :: vgn2(3)
    real(kind=8) :: x1, y1, z1, x2, y2, z2, xxl
    character(len=8) :: nomail, typel
    character(len=24) :: mlgcnx, mlgcoo
    character(len=24) :: valk(2), vec_name
    integer(kind=8), pointer :: mailles(:) => null()
    integer(kind=8), pointer :: noeud1(:) => null()
    integer(kind=8), pointer :: noeud2(:) => null()
    real(kind=8), pointer :: vale(:) => null()
    integer(kind=8), pointer :: typmail(:) => null()
    integer(kind=8), pointer :: maex(:) => null()
    integer(kind=8), pointer :: noex(:) => null()
    integer(kind=8), pointer :: lisma_int(:) => null()
    integer(kind=8), pointer :: node_line(:) => null()
    blas_int :: b_incx, b_incy, b_n
    mpi_int :: mrank, msize
    aster_logical :: l_par
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
    call jeexin(mesh//'.MAEX', iret)
    l_par = .false.
    if (iret .ne. 0) then
        call jeveuo(mesh//'.MAEX', 'L', vi=maex)
        call jeveuo(mesh//'.NOEX', 'L', vi=noex)
        call dismoi('NB_NO_MAILLA', mesh, 'MAILLAGE', repi=node_nb)
        vec_name = '&&PEAIR.NODE_LINE'
        call wkvect(vec_name, 'V V I', node_nb, vi=node_line)
        l_par = .true.
    end if
    AS_ALLOCATE(vi=lisma_int, size=nbma)
    call asmpi_info(rank=mrank, size=msize)
    rang = to_aster_int(mrank)
    nbproc = to_aster_int(msize)
    nbma_int = 0
    do ima = 1, nbma
        numa = lisma(ima)
        if (l_par) then
            if (maex(numa) .eq. rang) then
                nbma_int = nbma_int+1
                lisma_int(nbma_int) = numa
            end if
        else
            nbma_int = nbma_int+1
            lisma_int(nbma_int) = numa
        end if
    end do
!
    AS_ALLOCATE(vi=noeud1, size=nbma_int*3)
    AS_ALLOCATE(vi=noeud2, size=nbma_int*3)
    AS_ALLOCATE(vi=mailles, size=nbma_int)
!
!     VERIFICATION DU TYPE DES MAILLES ET STOCKAGE DES CONNECTIVITES
!
    long = 0.d0
    nbel = 0
    do ima = 1, nbma_int
        numa = lisma_int(ima)
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
    ASSERT(nbma_int .eq. nbel)
!
!     VERIFICATION QUE LE CONTOUR EST FERME
!
    ma_ext = 0
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
        if (iext1 .eq. 0) then
            if (l_par) then
                if (noex(ni1) .ne. rang) then
                    node_line(ni1) = 1
                end if
            end if
            ma_ext = ima
            no_1 = 2
            no_2 = 1
            nbext1 = nbext1+1
        end if
        if (iext2 .eq. 0) then
            if (l_par) then
                if (noex(ni2) .ne. rang) then
                    node_line(ni2) = 1
                end if
            end if
            ma_ext = ima
            no_1 = 1
            no_2 = 2
            nbext2 = nbext2+1
        end if
    end do
    if (ma_ext .eq. 0) then
        ma_ext = 1
        no_1 = 2
        no_2 = 1
    end if
    if ((nbext1 .ne. 0) .and. (nbext2 .ne. 0)) then
!       Verification que le contour est ferme en parallele
!       On regarde que tous les noeuds ghosts d'un sous-domaine sont
!       relies a une autre partie du contour d'un autre sous-domaine
        if (l_par) then
            call vector_ghosts_comm(vec_name, mesh)
            do ino = 1, node_nb
                if (node_line(ino) .ne. 0) then
                    call utmess('F', 'UTILITY_2')
                end if
            end do
        else
            call utmess('F', 'UTILITY_2')
        end if
    end if
!
!     VERIFICATION QUE LE CONTOUR EST CONTINU ET REORIENTATION
!
    nbma_ext = 0
    nbe = 1
    mailles(ma_ext) = 0
    noeud2(1) = noeud1(3*ma_ext-no_1)
    noeud2(1+1) = noeud1(3*ma_ext-no_2)
    noeud2(1+2) = noeud1(3*ma_ext)
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
    if (l_par) then
        nbma_ext = nbma_ext+1
    else
        call utmess('F', 'UTILITY_2')
    end if
43  continue
    mailles(jma) = 0
    if (nbe .ge. nbma_int) then
        goto 11
    else
        goto 41
    end if
11  continue
    if (l_par) then
!       En parallele, si on nbma_ext > 0 et qu'on n'a pas parcouru
!       encore toutes les mailles, c'est qu'on n'a pas reussi à
!       construire un ruban continu unique sur ce processeur
!       Si on continuait, il manquerait des noeuds et la mesure
!       du contour serait fausse
        if (nbma_ext .gt. 0 .and. nbe .ne. nbma_int) then
            call utmess('F', 'UTILITY_2')
        end if
    end if
    ASSERT(nbma_int .eq. nbe)
    nj2 = noeud2(3*nbe-1)
    nj0 = noeud2(1)
!
!     CALCUL DU CDG APPROXIMATIF
!
    mlgcoo = mesh//'.COORDO    .VALE'
    call jeveuo(mlgcoo, 'L', jdco)
    do ima = 1, nbma_int
        nj1 = noeud2(3*ima-no_1)
        orig(1) = orig(1)+zr(jdco-1+3*nj1-2)
        orig(2) = orig(2)+zr(jdco-1+3*nj1-1)
        orig(3) = orig(3)+zr(jdco-1+3*nj1)
    end do
    nbma_tot = nbma_int
    if (l_par) then
        call jedetr(vec_name)
        call asmpi_comm_vect('MPI_SUM', 'R', nbval=3, vr=orig)
        call asmpi_comm_vect('MPI_SUM', 'I', sci=nbma_tot)
    end if
    orig(1) = orig(1)/nbma_tot
    orig(2) = orig(2)/nbma_tot
    orig(3) = orig(3)/nbma_tot
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
    do ima = 1, nbma_int
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
    AS_DEALLOCATE(vi=lisma_int)
    if (l_par) then
        call asmpi_comm_vect('MPI_SUM', 'R', scr=aire)
        call asmpi_comm_vect('MPI_SUM', 'R', scr=long)
    end if
    call jedema()
end subroutine
