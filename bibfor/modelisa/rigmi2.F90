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
subroutine rigmi2(noma, nogr, ifreq, nfreq, ifmis, &
                  rigma, rigma2)
    implicit none
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/irmiim.h"
#include "asterfort/iunifi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/r8inir.h"
#include "asterfort/wkvect.h"
#include "asterfort/int_to_char8.h"
!
    integer(kind=8) :: ifmis
    integer(kind=8) :: ifreq, nfreq
    character(len=8) :: noma
    character(len=24) :: nogr
    real(kind=8) :: rigma(*), rigma2(*)
!
! ------------------------------------------------------------------
!
    character(len=8) :: nommai
    character(len=24) :: magrma, manoma, tabrig
!
! -----------------------------------------------------------------------
    integer(kind=8) :: i1, i2, ifr, ii, ij, im
    integer(kind=8) :: in, inoe, iret, isoto
    integer(kind=8) :: jrig, ldgm, ldnm, nb, nbmode, nbno, noemax
!
    real(kind=8) :: r1, r2, r3
    integer(kind=8), pointer :: noeud(:) => null()
    integer(kind=8), pointer :: parno(:) => null()
    real(kind=8), pointer :: sompar(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
    ifr = iunifi('RESULTAT')
!
    magrma = noma//'.GROUPEMA'
    manoma = noma//'.CONNEX'
    noemax = 0
!
!
    call jelira(jexnom(magrma, nogr), 'LONUTI', nb)
    call jeveuo(jexnom(magrma, nogr), 'L', ldgm)
    do in = 0, nb-1
        call jeveuo(jexnum(manoma, zi(ldgm+in)), 'L', ldnm)
        inoe = zi(ldnm)
        noemax = max(noemax, inoe)
        inoe = zi(ldnm+1)
        noemax = max(noemax, inoe)
    end do
!
!        TABLEAU DE PARTICIPATION DES NOEUDS DE L INTERFACE
!
    AS_ALLOCATE(vi=parno, size=noemax)
!
    do in = 0, nb-1
        call jeveuo(jexnum(manoma, zi(ldgm+in)), 'L', ldnm)
        inoe = zi(ldnm)
        parno(inoe) = parno(inoe)+1
        inoe = zi(ldnm+1)
        parno(inoe) = parno(inoe)+1
    end do
!
    nbno = 0
    do ij = 1, noemax
        if (parno(ij) .eq. 0) goto 25
        nbno = nbno+1
25      continue
    end do
!
    AS_ALLOCATE(vi=noeud, size=nbno)
    ii = 0
    do ij = 1, noemax
        if (parno(ij) .eq. 0) goto 26
        ii = ii+1
        noeud(ii) = ij
26      continue
    end do
!
!     LECTURE DES RIGIDITES ELEMENTAIRES
!
    tabrig = '&&ACEARM.RIGM'
    call jeexin(tabrig, iret)
    if (iret .eq. 0) call irmiim(ifmis, ifreq, nfreq, nbno, tabrig)
    call jeveuo(tabrig, 'L', jrig)
    nbmode = 3*nbno
    call wkvect('&&RIGMI2.SOMTOT', 'V V R', nbmode, isoto)
    AS_ALLOCATE(vr=sompar, size=nbmode)
    do i1 = 1, nbno
        do i2 = 1, nbno
            if (i1 .ne. i2) then
                zr(isoto+3*i1-3) = zr(isoto+3*i1-3)+zr(jrig+(3*i2-3)*nbmode+3*i1-3)
                zr(isoto+3*i1-2) = zr(isoto+3*i1-2)+zr(jrig+(3*i2-2)*nbmode+3*i1-2)
                zr(isoto+3*i1-1) = zr(isoto+3*i1-1)+zr(jrig+(3*i2-1)*nbmode+3*i1-1)
            end if
        end do
    end do
!
    i1 = 0
    i2 = 0
    do in = 0, nb-1
        im = zi(ldgm+in)
        call jeveuo(jexnum(manoma, zi(ldgm+in)), 'L', ldnm)
        do ii = 1, nbno
            if (zi(ldnm) .eq. noeud(ii)) i1 = ii
            if (zi(ldnm+1) .eq. noeud(ii)) i2 = ii
        end do
        sompar(1+3*i1-3) = sompar(1+3*i1-3)+zr(jrig+(3*i2-3)*nbmode+3*i1-3)
        sompar(1+3*i2-3) = sompar(1+3*i2-3)+zr(jrig+(3*i2-3)*nbmode+3*i1-3)
        sompar(1+3*i1-2) = sompar(1+3*i1-2)+zr(jrig+(3*i2-2)*nbmode+3*i1-2)
        sompar(1+3*i2-2) = sompar(1+3*i2-2)+zr(jrig+(3*i2-2)*nbmode+3*i1-2)
        sompar(1+3*i1-1) = sompar(1+3*i1-1)+zr(jrig+(3*i2-1)*nbmode+3*i1-1)
        sompar(1+3*i2-1) = sompar(1+3*i2-1)+zr(jrig+(3*i2-1)*nbmode+3*i1-1)
    end do
!
    do in = 0, nb-1
        im = zi(ldgm+in)
        call jeveuo(jexnum(manoma, zi(ldgm+in)), 'L', ldnm)
        do ii = 1, nbno
            if (zi(ldnm) .eq. noeud(ii)) i1 = ii
            if (zi(ldnm+1) .eq. noeud(ii)) i2 = ii
        end do
        rigma(3*in+1) = 0.5d0*zr(jrig+(3*i2-3)*nbmode+3*i1-3)* &
                        (zr(isoto+3*i1-3)/sompar(1+3*i1-3)+zr(isoto+3*i2-3)/sompar(1+3*i2-3))
        rigma(3*in+2) = 0.5d0*zr(jrig+(3*i2-2)*nbmode+3*i1-2)* &
                        (zr(isoto+3*i1-2)/sompar(1+3*i1-2)+zr(isoto+3*i2-2)/sompar(1+3*i2-2))
        rigma(3*in+3) = 0.5d0*zr(jrig+(3*i2-1)*nbmode+3*i1-1)* &
                        (zr(isoto+3*i1-1)/sompar(1+3*i1-1)+zr(isoto+3*i2-1)/sompar(1+3*i2-1))
    end do
!
    call r8inir(3*nbno, 0.d0, rigma2, 1)
    do in = 0, nb-1
        im = zi(ldgm+in)
        call jeveuo(jexnum(manoma, zi(ldgm+in)), 'L', ldnm)
        do ii = 1, nbno
            if (zi(ldnm) .eq. noeud(ii)) i1 = ii
            if (zi(ldnm+1) .eq. noeud(ii)) i2 = ii
        end do
        r1 = rigma(3*in+1)
        r2 = rigma(3*in+2)
        r3 = rigma(3*in+3)
!
        rigma(3*in+1) = r1
        rigma(3*in+2) = r2
        rigma(3*in+3) = r3
        rigma2(3*(i1-1)+1) = r1+rigma2(3*(i1-1)+1)
        rigma2(3*(i1-1)+2) = r2+rigma2(3*(i1-1)+2)
        rigma2(3*(i1-1)+3) = r3+rigma2(3*(i1-1)+3)
        rigma2(3*(i2-1)+1) = r1+rigma2(3*(i2-1)+1)
        rigma2(3*(i2-1)+2) = r2+rigma2(3*(i2-1)+2)
        rigma2(3*(i2-1)+3) = r3+rigma2(3*(i2-1)+3)
        nommai = int_to_char8(im)
        write (ifr, 100) nommai, -r1, -r2, -r3
    end do
!
100 format(2x, '_F ( MAILLE=''', a8, ''',', 1x, 'CARA= ''K_T_D_L'' , ', &
           /7x, 'VALE=(', 1x, 3(1x, 1pe12.5, ','), 1x, '),', &
           /'   ),')
!
    AS_DEALLOCATE(vi=parno)
    AS_DEALLOCATE(vi=noeud)
    call jedetr('&&RIGMI2.SOMTOT')
    AS_DEALLOCATE(vr=sompar)
!
    call jedema()
end subroutine
