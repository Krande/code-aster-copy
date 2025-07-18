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
subroutine iredm1(masse, noma, basemo, nbmode, nbmods, &
                  iamor, mass, rigi, amored, freq, &
                  smass, srigi, samor, cmass, crigi, &
                  camor)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/copmod.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvem.h"
#include "asterfort/getvis.h"
#include "asterfort/getvtx.h"
#include "asterfort/iunifi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/nbec.h"
#include "asterfort/ulopen.h"
#include "asterfort/utmess.h"
#include "asterfort/vtcreb.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "blas/dcopy.h"
!
    character(len=8) :: masse, noma, basemo
    real(kind=8) :: mass(*), rigi(*), smass(*), srigi(*), samor(*), cmass(*)
    real(kind=8) :: crigi(*), camor(*), amored(*), freq(*)
!              INTERFACE ASTER - MISS3D : PROCEDURE  IMPR_MACR_ELEM
!     ------------------------------------------------------------------
!
    integer(kind=8) :: aprno, gd, tabl(8), tab2(8)
    character(len=8) :: k8b, typi, impmod, impmec, interf, formim
    character(len=14) :: nume
    character(len=16) :: nomcmd
    character(len=24) :: magrma, manoma, nprno
    character(len=24) :: nomch0
    character(len=80) :: titre
    aster_logical :: lamor
!     ------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iamor, ibid, ic, idbase
    integer(kind=8) :: iddl, iddl0, idgm2, idgm3, idgm4, idgm5
    integer(kind=8) :: ifmis, ii, ij, imess, in
    integer(kind=8) :: ino, inoe, j, j2
    integer(kind=8) :: k, l, ldgm, ldnm, nb
    integer(kind=8) :: nbgr, nbgr2, nbgr3, nbgr4, nbgr5, nbma, nbma2
    integer(kind=8) :: nbma3, nbma4, nbma5, nbmode, nbmods, nbmodt, nbno
    integer(kind=8) :: nbnoeu, nbv, ncmp, nec, neq, nf, ni
    integer(kind=8) :: nm, nn, nti, nu
    real(kind=8) :: zero
    complex(kind=8) :: cbid
    character(len=24), pointer :: group_solstru(:) => null()
    integer(kind=8), pointer :: noeud(:) => null()
    integer(kind=8), pointer :: parno(:) => null()
    real(kind=8), pointer :: vect1(:) => null()
    real(kind=8), pointer :: vect2(:) => null()
    character(len=8), pointer :: idc_type(:) => null()
    real(kind=8), pointer :: vale(:) => null()
    blas_int :: b_incx, b_incy, b_n
    cbid = dcmplx(0.d0, 0.d0)
!-----------------------------------------------------------------------
    call jemarq()
    imess = iunifi('MESSAGE')
    zero = 0.d0
    nomch0 = '&&IREDM1.CHAMNO'
    magrma = noma//'.GROUPEMA'
    manoma = noma//'.CONNEX'
    lamor = iamor .ne. 0
    call getres(k8b, k8b, nomcmd)
    call getvis(' ', 'UNITE', scal=ifmis, nbret=nu)
    call ulopen(ifmis, ' ', ' ', 'NEW', 'O')
    call getvtx(' ', 'IMPR_MODE_STAT', scal=impmod, nbret=ni)
    call getvtx(' ', 'IMPR_MODE_MECA', scal=impmec, nbret=ni)
    call getvtx(' ', 'FORMAT_R', scal=formim, nbret=nf)
    call getvtx(' ', 'SOUS_TITRE', scal=titre, nbret=nti)
!
!
!     --- ON RECUPERE LE TYPE D'INTERFACE ---
!
    call dismoi('REF_INTD_PREM', basemo, 'RESU_DYNA', repk=interf, arret='C')
    if (interf .ne. ' ') then
        call jeveuo(interf//'.IDC_TYPE', 'L', vk8=idc_type)
        typi = idc_type(1)
    else
        typi = 'CRAIGB'
    end if
!
    write (imess, '(1X,I6,1X,''MODES DYNAMIQUES'',1X,A8)') nbmode, typi
    write (imess, '(1X,I6,1X,''MODES STATIQUES'' ,2X,A8)') nbmods, typi
!
    call dismoi('NOM_NUME_DDL', masse, 'MATR_ASSE', repk=nume)
    call dismoi('NB_EQUA', masse, 'MATR_ASSE', repi=neq)
    call dismoi('NB_NO_MAILLA', noma, 'MAILLAGE', repi=nbnoeu)
    call dismoi('NUM_GD_SI', nume, 'NUME_DDL', repi=gd)
    if (interf .eq. ' ') call vtcreb(nomch0, 'V', 'R', nume_ddlz=nume, nb_equa_outz=neq)
!CC
!     ----- DEBUT DES IMPRESSIONS DE MISS3D -----
!CC
!
    write (ifmis, 120) 'DYNA', nbmode, typi
    write (ifmis, 120) 'STAT', nbmods, typi
    nbmodt = nbmode+nbmods
!
    if (nti .ne. 0) then
        write (ifmis, '(''TITRE'',/A80)') titre
        write (imess, '(A80)') titre
    end if
    call getvem(noma, 'GROUP_MA', ' ', 'GROUP_MA_INTERF', 1, &
                0, k8b, nbgr)
    nbgr = -nbgr
    AS_ALLOCATE(vk24=group_solstru, size=nbgr)
    call getvem(noma, 'GROUP_MA', ' ', 'GROUP_MA_INTERF', 1, &
                nbgr, group_solstru, nbv)
    call getvem(noma, 'GROUP_MA', ' ', 'GROUP_MA_FLU_STR', 1, &
                0, k8b, nbgr2)
    nbgr2 = -nbgr2
    if (nbgr2 .eq. 0) then
        call wkvect('&&IREDM1.GROUP_FLUSTRU', 'V V K24', 1, idgm2)
    else
        call wkvect('&&IREDM1.GROUP_FLUSTRU', 'V V K24', nbgr2, idgm2)
    end if
    call getvem(noma, 'GROUP_MA', ' ', 'GROUP_MA_FLU_STR', 1, &
                nbgr2, zk24(idgm2), nbv)
    call getvem(noma, 'GROUP_MA', ' ', 'GROUP_MA_FLU_SOL', 1, &
                0, k8b, nbgr3)
    nbgr3 = -nbgr3
    if (nbgr3 .eq. 0) then
        call wkvect('&&IREDM1.GROUP_FLUSOL', 'V V K24', 1, idgm3)
    else
        call wkvect('&&IREDM1.GROUP_FLUSOL', 'V V K24', nbgr3, idgm3)
    end if
    call getvem(noma, 'GROUP_MA', ' ', 'GROUP_MA_FLU_SOL', 1, &
                nbgr3, zk24(idgm3), nbv)
    call getvem(noma, 'GROUP_MA', ' ', 'GROUP_MA_SOL_SOL', 1, &
                0, k8b, nbgr4)
    nbgr4 = -nbgr4
    if (nbgr4 .eq. 0) then
        call wkvect('&&IREDM1.GROUP_LIBRE', 'V V K24', 1, idgm4)
    else
        call wkvect('&&IREDM1.GROUP_LIBRE', 'V V K24', nbgr4, idgm4)
    end if
    call getvem(noma, 'GROUP_MA', ' ', 'GROUP_MA_SOL_SOL', 1, &
                nbgr4, zk24(idgm4), nbv)
    call getvem(noma, 'GROUP_MA', ' ', 'GROUP_MA_CONTROL', 1, &
                0, k8b, nbgr5)
    nbgr5 = -nbgr5
    if (nbgr5 .eq. 0) then
        call wkvect('&&IREDM1.GROUP_CONTROL', 'V V K24', 1, idgm5)
    else
        call wkvect('&&IREDM1.GROUP_CONTROL', 'V V K24', nbgr5, idgm5)
    end if
    call getvem(noma, 'GROUP_MA', ' ', 'GROUP_MA_CONTROL', 1, &
                nbgr5, zk24(idgm5), nbv)
!
!
!        TABLEAU DE PARTICIPATION DES NOEUDS DE L INTERFACE
!
    AS_ALLOCATE(vi=parno, size=nbnoeu)
!
    nbma = 0
    nbma2 = 0
    nbma3 = 0
    nbma4 = 0
    nbma5 = 0
    do i = 1, nbgr
        call jelira(jexnom(magrma, group_solstru(i)), 'LONUTI', nb)
        call jeveuo(jexnom(magrma, group_solstru(i)), 'L', ldgm)
        nbma = nbma+nb
        do in = 0, nb-1
            call jelira(jexnum(manoma, zi(ldgm+in)), 'LONMAX', nm)
            if (nm .ne. 3 .and. nm .ne. 4 .and. nm .ne. 6 .and. nm .ne. 8) then
                call utmess('F', 'UTILITAI2_36')
            end if
            call jeveuo(jexnum(manoma, zi(ldgm+in)), 'L', ldnm)
            do nn = 1, nm
                inoe = zi(ldnm+nn-1)
                parno(inoe) = parno(inoe)+1
            end do
        end do
    end do
    do i = 1, nbgr2
        call jelira(jexnom(magrma, zk24(idgm2+i-1)), 'LONUTI', nb)
        call jeveuo(jexnom(magrma, zk24(idgm2+i-1)), 'L', ldgm)
        nbma2 = nbma2+nb
        do in = 0, nb-1
            call jelira(jexnum(manoma, zi(ldgm+in)), 'LONMAX', nm)
            if (nm .ne. 3 .and. nm .ne. 4 .and. nm .ne. 6 .and. nm .ne. 8) then
                call utmess('F', 'UTILITAI2_37')
            end if
            call jeveuo(jexnum(manoma, zi(ldgm+in)), 'L', ldnm)
            do nn = 1, nm
                inoe = zi(ldnm+nn-1)
                parno(inoe) = parno(inoe)+1
            end do
        end do
    end do
    do i = 1, nbgr3
        call jelira(jexnom(magrma, zk24(idgm3+i-1)), 'LONUTI', nb)
        call jeveuo(jexnom(magrma, zk24(idgm3+i-1)), 'L', ldgm)
        nbma3 = nbma3+nb
        do in = 0, nb-1
            call jelira(jexnum(manoma, zi(ldgm+in)), 'LONMAX', nm)
            if (nm .ne. 3 .and. nm .ne. 4 .and. nm .ne. 6 .and. nm .ne. 8) then
                call utmess('F', 'UTILITAI2_38')
            end if
            call jeveuo(jexnum(manoma, zi(ldgm+in)), 'L', ldnm)
            do nn = 1, nm
                inoe = zi(ldnm+nn-1)
                parno(inoe) = parno(inoe)+1
            end do
        end do
    end do
    do i = 1, nbgr4
        call jelira(jexnom(magrma, zk24(idgm4+i-1)), 'LONUTI', nb)
        call jeveuo(jexnom(magrma, zk24(idgm4+i-1)), 'L', ldgm)
        nbma4 = nbma4+nb
        do in = 0, nb-1
            call jelira(jexnum(manoma, zi(ldgm+in)), 'LONMAX', nm)
            if (nm .ne. 3 .and. nm .ne. 4 .and. nm .ne. 6 .and. nm .ne. 8) then
                call utmess('F', 'UTILITAI2_39')
            end if
            call jeveuo(jexnum(manoma, zi(ldgm+in)), 'L', ldnm)
            do nn = 1, nm
                inoe = zi(ldnm+nn-1)
                parno(inoe) = parno(inoe)+1
            end do
        end do
    end do
    do i = 1, nbgr5
        call jelira(jexnom(magrma, zk24(idgm5+i-1)), 'LONUTI', nb)
        call jeveuo(jexnom(magrma, zk24(idgm5+i-1)), 'L', ldgm)
        nbma5 = nbma5+nb
        do in = 0, nb-1
            call jeveuo(jexnum(manoma, zi(ldgm+in)), 'L', ldnm)
            inoe = zi(ldnm)
            parno(inoe) = parno(inoe)+1
        end do
    end do
!
    nbno = 0
    do ij = 1, nbnoeu
        if (parno(ij) .ne. 0) then
            nbno = nbno+1
        end if
    end do
!
    AS_ALLOCATE(vi=noeud, size=nbno)
    ii = 0
    do ij = 1, nbnoeu
        if (parno(ij) .ne. 0) then
            ii = ii+1
            noeud(ii) = ij
        end if
    end do
!
!
!     --- ECRITURE DESCRIPTION NOEUDS STRUCTURE ---
    call jeveuo(noma//'.COORDO    .VALE', 'L', vr=vale)
    nprno = nume//'.NUME.PRNO'
    call jenonu(jexnom(nprno(1:19)//'.LILI', '&MAILLA'), ibid)
    call jeveuo(jexnum(nprno, ibid), 'L', aprno)
    nec = nbec(gd)
    write (imess, '(1X,I6,1X,''NOEUDS'')') nbno
    write (ifmis, '(''NOEU'',1X,I6)') nbno
    do ino = 1, nbno
        inoe = noeud(ino)
        ncmp = zi(aprno+(nec+2)*(inoe-1)+2-1)
        write (ifmis, '(3(1X,1PE12.5))') (vale(1+3*(inoe-1)+in-1), &
                                          in=1, 3)
    end do
    write (imess, '(1X,I6,1X,''ELEMENTS SOLSTRU'')') nbma
    write (ifmis, '(''ELEM'',1X,I6)') nbma
    do i = 1, nbgr
        call jelira(jexnom(magrma, group_solstru(i)), 'LONUTI', nb)
        call jeveuo(jexnom(magrma, group_solstru(i)), 'L', ldgm)
        do in = 0, nb-1
            do k = 1, 8
                tabl(k) = 0
                tab2(k) = 0
            end do
            call jelira(jexnum(manoma, zi(ldgm+in)), 'LONMAX', nm)
            call jeveuo(jexnum(manoma, zi(ldgm+in)), 'L', ldnm)
            do nn = 1, nm
                do ij = 1, nbno
                    if (zi(ldnm+nn-1) .eq. noeud(ij)) tab2(nn) = ij
                end do
                if (nm .le. 4) tabl(2*nn-1) = tab2(nn)
                if (nm .eq. 6 .and. nn .le. 3) tabl(2*nn-1) = tab2(nn)
                if (nm .eq. 6 .and. nn .gt. 3) tabl(2*nn-nm) = tab2(nn)
                if (nm .eq. 8 .and. nn .le. 4) tabl(2*nn-1) = tab2(nn)
                if (nm .eq. 8 .and. nn .gt. 4) tabl(2*nn-nm) = tab2(nn)
            end do
            write (ifmis, '(8(1X,I6))') (tabl(k), k=1, 8)
        end do
    end do
    write (imess, '(1X,I6,1X,''ELEMENTS FLUSTRU'')') nbma2
    if (nbma2 .ne. 0) write (ifmis, '(''ELEM'',1X,I6)') nbma2
    do i = 1, nbgr2
        call jelira(jexnom(magrma, zk24(idgm2+i-1)), 'LONUTI', nb)
        call jeveuo(jexnom(magrma, zk24(idgm2+i-1)), 'L', ldgm)
        do in = 0, nb-1
            do k = 1, 8
                tabl(k) = 0
                tab2(k) = 0
            end do
            call jelira(jexnum(manoma, zi(ldgm+in)), 'LONMAX', nm)
            call jeveuo(jexnum(manoma, zi(ldgm+in)), 'L', ldnm)
            do nn = 1, nm
                do ij = 1, nbno
                    if (zi(ldnm+nn-1) .eq. noeud(ij)) tab2(nn) = ij
                end do
                if (nm .le. 4) tabl(2*nn-1) = tab2(nn)
                if (nm .eq. 6 .and. nn .le. 3) tabl(2*nn-1) = tab2(nn)
                if (nm .eq. 6 .and. nn .gt. 3) tabl(2*nn-nm) = tab2(nn)
                if (nm .eq. 8 .and. nn .le. 4) tabl(2*nn-1) = tab2(nn)
                if (nm .eq. 8 .and. nn .gt. 4) tabl(2*nn-nm) = tab2(nn)
            end do
            write (ifmis, '(8(1X,I6))') (tabl(k), k=1, 8)
        end do
    end do
    write (imess, '(1X,I6,1X,''ELEMENTS FLUSOL'')') nbma3
    if (nbma3 .ne. 0) write (ifmis, '(''ELEM'',1X,I6)') nbma3
    do i = 1, nbgr3
        call jelira(jexnom(magrma, zk24(idgm3+i-1)), 'LONUTI', nb)
        call jeveuo(jexnom(magrma, zk24(idgm3+i-1)), 'L', ldgm)
        do in = 0, nb-1
            do k = 1, 8
                tabl(k) = 0
                tab2(k) = 0
            end do
            call jelira(jexnum(manoma, zi(ldgm+in)), 'LONMAX', nm)
            call jeveuo(jexnum(manoma, zi(ldgm+in)), 'L', ldnm)
            do nn = 1, nm
                do ij = 1, nbno
                    if (zi(ldnm+nn-1) .eq. noeud(ij)) tab2(nn) = ij
                end do
                if (nm .le. 4) tabl(2*nn-1) = tab2(nn)
                if (nm .eq. 6 .and. nn .le. 3) tabl(2*nn-1) = tab2(nn)
                if (nm .eq. 6 .and. nn .gt. 3) tabl(2*nn-nm) = tab2(nn)
                if (nm .eq. 8 .and. nn .le. 4) tabl(2*nn-1) = tab2(nn)
                if (nm .eq. 8 .and. nn .gt. 4) tabl(2*nn-nm) = tab2(nn)
            end do
            write (ifmis, '(8(1X,I6))') (tabl(k), k=1, 8)
        end do
    end do
    write (imess, '(1X,I6,1X,''ELEMENTS LIBRE'')') nbma4
    if (nbma4 .ne. 0) write (ifmis, '(''ELEM'',1X,I6)') nbma4
    do i = 1, nbgr4
        call jelira(jexnom(magrma, zk24(idgm4+i-1)), 'LONUTI', nb)
        call jeveuo(jexnom(magrma, zk24(idgm4+i-1)), 'L', ldgm)
        do in = 0, nb-1
            do k = 1, 8
                tabl(k) = 0
                tab2(k) = 0
            end do
            call jelira(jexnum(manoma, zi(ldgm+in)), 'LONMAX', nm)
            call jeveuo(jexnum(manoma, zi(ldgm+in)), 'L', ldnm)
            do nn = 1, nm
                do ij = 1, nbno
                    if (zi(ldnm+nn-1) .eq. noeud(ij)) tab2(nn) = ij
                end do
                if (nm .le. 4) tabl(2*nn-1) = tab2(nn)
                if (nm .eq. 6 .and. nn .le. 3) tabl(2*nn-1) = tab2(nn)
                if (nm .eq. 6 .and. nn .gt. 3) tabl(2*nn-nm) = tab2(nn)
                if (nm .eq. 8 .and. nn .le. 4) tabl(2*nn-1) = tab2(nn)
                if (nm .eq. 8 .and. nn .gt. 4) tabl(2*nn-nm) = tab2(nn)
            end do
            write (ifmis, '(8(1X,I6))') (tabl(k), k=1, 8)
        end do
    end do
    write (imess, '(1X,I6,1X,''POINTS CONTROLE'')') nbma5
    if (nbma5 .ne. 0) write (ifmis, '(''POINT'',1X,I6)') nbma5
    do i = 1, nbgr5
        call jelira(jexnom(magrma, zk24(idgm5+i-1)), 'LONUTI', nb)
        call jeveuo(jexnom(magrma, zk24(idgm5+i-1)), 'L', ldgm)
        do in = 0, nb-1
            call jeveuo(jexnum(manoma, zi(ldgm+in)), 'L', ldnm)
            do ij = 1, nbno
                if (zi(ldnm) .eq. noeud(ij)) inoe = ij
            end do
            write (ifmis, '(1X,I6)') inoe
        end do
    end do
!
    call wkvect('&&IREDM1.BASEMO', 'V V R', nbmodt*neq, idbase)
    call copmod(basemo, numer=nume, bmodr=zr(idbase))
!
! --- ALLOCATION VECTEUR DE TRAVAIL
!
    AS_ALLOCATE(vr=vect1, size=neq)
    AS_ALLOCATE(vr=vect2, size=neq)
!
    if (typi(1:5) .ne. 'CRAIG' .or. impmec .eq. 'OUI') then
        do j = 1, nbmode
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(idbase+(j-1)*neq), b_incx, vect1, b_incy)
            write (ifmis, '(''MODE DYNA INTER'',1X,I6)') j
            do ino = 1, nbno
                inoe = noeud(ino)
                iddl = zi(aprno+(nec+2)*(inoe-1)+1-1)-1
                ncmp = zi(aprno+(nec+2)*(inoe-1)+2-1)
                ncmp = min(6, ncmp)
                iddl0 = iddl+1
                if (iddl0 .eq. 0) then
                    write (ifmis, 110) zero, zero, zero, zero, zero, zero
                else
                    write (ifmis, 110) (vect1(1+iddl+ic-1), ic=1, &
                                        ncmp)
                end if
            end do
        end do
    end if
!
    if (formim .eq. '1PE16.9') then
        write (ifmis, 100) 'DYNA FREQ', (freq(k), k=1, nbmode)
        write (ifmis, 100) 'DYNA AMOR', (amored(k), k=1, nbmode)
        write (ifmis, 100) 'DYNA MASS', (mass(k+(k-1)*nbmode), k=1, &
                                         nbmode)
        write (ifmis, 100) 'DYNA RIGI', (rigi(k+(k-1)*nbmode), k=1, &
                                         nbmode)
    else
        write (ifmis, 130) 'DYNA FREQ', (freq(k), k=1, nbmode)
        write (ifmis, 130) 'DYNA AMOR', (amored(k), k=1, nbmode)
        write (ifmis, 130) 'DYNA MASS', (mass(k+(k-1)*nbmode), k=1, &
                                         nbmode)
        write (ifmis, 130) 'DYNA RIGI', (rigi(k+(k-1)*nbmode), k=1, &
                                         nbmode)
    end if
!
    if (typi(1:5) .ne. 'CRAIG' .or. impmod .eq. 'OUI') then
        do j = 1, nbmods
            j2 = j+nbmode
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(idbase+(j2-1)*neq), b_incx, vect2, b_incy)
            write (ifmis, '(''MODE STAT INTER'',1X,I6)') j
            do ino = 1, nbno
                inoe = noeud(ino)
                iddl = zi(aprno+(nec+2)*(inoe-1)+1-1)-1
                ncmp = zi(aprno+(nec+2)*(inoe-1)+2-1)
                ncmp = min(6, ncmp)
                iddl0 = iddl+1
                if (iddl0 .eq. 0) then
                    write (ifmis, 110) zero, zero, zero, zero, zero, zero
                else
                    write (ifmis, 110) (vect2(1+iddl+ic-1), ic=1, &
                                        ncmp)
                end if
            end do
        end do
    end if
    if (formim .eq. '1PE16.9') then
        write (ifmis, 100) 'STAT MASS', ((smass(k+(l-1)*nbmods), k=1, &
                                          nbmods), l=1, nbmods)
        write (ifmis, 100) 'STAT RIGI', ((srigi(k+(l-1)*nbmods), k=1, &
                                          nbmods), l=1, nbmods)
        if (lamor) write (ifmis, 100) 'STAT AMOR', &
            ((samor(k+(l-1)*nbmods), k=1, nbmods), l=1, nbmods)
        write (ifmis, '(''COUPL'',2(1X,I6))') nbmode, nbmods
        write (ifmis, 100) 'COUPL MASS', ((cmass(k+(l-1)*nbmods), k=1, &
                                           nbmods), l=1, nbmode)
        write (ifmis, 100) 'COUPL RIGI', ((crigi(k+(l-1)*nbmods), k=1, &
                                           nbmods), l=1, nbmode)
        if (lamor) write (ifmis, 100) 'COUPL AMOR', &
            ((camor(k+(l-1)*nbmods), k=1, nbmods), l=1, nbmode)
    else
        write (ifmis, 130) 'STAT MASS', ((smass(k+(l-1)*nbmods), k=1, &
                                          nbmods), l=1, nbmods)
        write (ifmis, 130) 'STAT RIGI', ((srigi(k+(l-1)*nbmods), k=1, &
                                          nbmods), l=1, nbmods)
        if (lamor) write (ifmis, 130) 'STAT AMOR', &
            ((samor(k+(l-1)*nbmods), k=1, nbmods), l=1, nbmods)
        write (ifmis, '(''COUPL'',2(1X,I6))') nbmode, nbmods
        write (ifmis, 130) 'COUPL MASS', ((cmass(k+(l-1)*nbmods), k=1, &
                                           nbmods), l=1, nbmode)
        write (ifmis, 130) 'COUPL RIGI', ((crigi(k+(l-1)*nbmods), k=1, &
                                           nbmods), l=1, nbmode)
        if (lamor) write (ifmis, 130) 'COUPL AMOR', &
            ((camor(k+(l-1)*nbmods), k=1, nbmods), l=1, nbmode)
    end if
!
    if (formim .eq. '1PE16.9') then
        write (ifmis, '(A)') 'FORMAT_REAL_LONG'
    else
        write (ifmis, '(A)') 'FORMAT_REAL_COURT'
    end if
!
!CC
!     ----- FIN DES IMPRESSIONS DE MISS3D -----
!CC
!
100 format(a, /, 4(2x, 1p, d16.9))
110 format(6(1x, 1p, d12.5))
120 format(a4, 1x, i6, 1x, a8)
130 format(a, /, 6(1x, 1p, d12.5))
!
!
! --- MENAGE
!
    call detrsd('CHAM_NO', '&&IREDM1.CHAMNO')
    AS_DEALLOCATE(vk24=group_solstru)
    call jedetr('&&IREDM1.GROUP_FLUSTRU')
    call jedetr('&&IREDM1.GROUP_FLUSOL')
    call jedetr('&&IREDM1.GROUP_LIBRE')
    call jedetr('&&IREDM1.GROUP_CONTROL')
    AS_DEALLOCATE(vi=parno)
    AS_DEALLOCATE(vi=noeud)
    call jedetr('&&IREDM1.BASEMO')
    AS_DEALLOCATE(vr=vect1)
    AS_DEALLOCATE(vr=vect2)
!
    call jedema()
end subroutine
