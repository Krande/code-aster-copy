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
subroutine tbtrtb(tabin, basout, tabout, npara, lipara, &
                  lcrit, prec, crit)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/jecreo.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/tbtr01.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: npara
    real(kind=8) :: prec
    character(len=8) :: crit
    character(len=*) :: tabin, basout, tabout, lipara(*), lcrit(*)
!     TRI DE LA TABLE.
!     LE TRI NE PORTE QUE SUR LES TYPES I , R  ET  K
!     ARRET EN FATAL SUR LES COMPLEXES
! ----------------------------------------------------------------------
! IN  : TABIN  : NOM DE LA TABLE DONT ON VEUT TRIER DES LIGNES
! IN  : BASOUT : BASE DE CREATION DE "TABOUT"
! OUT : TABOUT : NOM DE LA TABLE QUI CONTIENDRA LES LIGNES TRIEES
! IN  : NPARA  : NOMBRE DE PARAMETRES IMPLIQUES DANS LE TRI
! IN  : LIPARA : LISTE DES PARAMETRES A TRIER
! IN  : LCRIT  : TYPES DE CRITERES: CR  CROISSANT
!                                   DE  DECROISSANT
! IN  : PREC   : PRECISION
! IN  : CRIT   : RELATIF / ABSOLU
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
    integer(kind=8) :: iret, nbpara, nblign, jnume, ii, jj
    integer(kind=8) :: i, j, k, n, m, ideb, ifin, nbuti, ndim
    integer(kind=8) :: jvall, kvall, jvale, kvale, ktbnp, ktbba
    character(len=1) :: base
    character(len=4) :: type, knume
    character(len=19) :: nomtab, nomta2
    character(len=24) :: nomjv, nojv2, nomjvl, nojvl2, inpar, jnpar
    character(len=24) :: valk
    aster_logical :: lok
    integer(kind=8), pointer :: tri2(:) => null()
    integer(kind=8), pointer :: tbnp(:) => null()
    character(len=24), pointer :: tblp(:) => null()
    character(len=24), pointer :: nktblp(:) => null()
! ----------------------------------------------------------------------
!
    call jemarq()
!
    nomtab = tabin
    nomta2 = tabout
    base = basout(1:1)
!
!     --- VERIFICATION DE LA BASE ---
!
    ASSERT(base .eq. 'V' .or. base .eq. 'G')
!
!     --- VERIFICATION DE LA TABLE ---
!
    call jeexin(nomtab//'.TBBA', iret)
    if (iret .eq. 0) then
        call utmess('F', 'UTILITAI4_64')
    end if
!
    call jeveuo(nomtab//'.TBNP', 'E', vi=tbnp)
    nbpara = tbnp(1)
    nblign = tbnp(2)
    if (nbpara .eq. 0) then
        call utmess('F', 'UTILITAI4_65')
    end if
    if (nblign .eq. 0) then
        call utmess('F', 'UTILITAI4_66')
    end if
!
!     --- VERIFICATION QUE LES PARAMETRES EXISTENT DANS LA TABLE ---
!
    call jeveuo(nomtab//'.TBLP', 'L', vk24=tblp)
    do i = 1, npara
        inpar = lipara(i)
        do j = 1, nbpara
            jnpar = tblp(1+4*(j-1))
            if (inpar .eq. jnpar) then
                type = tblp(1+4*(j-1)+1) (1:4)
                if (type(1:1) .eq. 'C') then
                    valk = inpar
                    call utmess('F', 'UTILITAI7_2', sk=valk)
                end if
                goto 10
            end if
        end do
        valk = inpar
        call utmess('F', 'UTILITAI6_89', sk=valk)
10      continue
    end do
!
    call wkvect('&&TBTRTB.TRI', 'V V I', nblign, jnume)
    AS_ALLOCATE(vi=tri2, size=nblign)
    do i = 1, nblign
        zi(jnume+i-1) = i
    end do
!
    call tbtr01(tabin, nbpara, lipara(1), nblign, zi(jnume))
!
    if (lcrit(1) (1:2) .eq. 'DE') then
        do i = 1, nblign
            tri2(i) = zi(jnume-1+i)
        end do
        do i = 1, nblign
            zi(jnume-1+i) = tri2(nblign-i+1)
        end do
    end if
!
    if (npara .eq. 1) then
        goto 104
    else if (npara .eq. 2) then
    else if (npara .gt. 2) then
        call utmess('F', 'UTILITAI4_87')
    end if
    i = 1
    inpar = lipara(i)
    do j = 1, nbpara
        jnpar = tblp(1+4*(j-1))
        if (inpar .eq. jnpar) then
            type = tblp(1+4*(j-1)+1) (1:4)
            nomjv = tblp(1+4*(j-1)+2)
            nomjvl = tblp(1+4*(j-1)+3)
            call jeveuo(nomjv, 'L', jvale)
            call jeveuo(nomjvl, 'L', jvall)
            ii = 0
30          continue
            ii = ii+1
            if (ii .ge. nblign) goto 32
            ideb = ii
            ifin = ii
            n = zi(jnume+ii-1)
            do jj = ii+1, nblign
                ifin = jj
                m = zi(jnume+jj-1)
                if (type(1:1) .eq. 'I') then
                    if (zi(jvale+n-1) .ne. zi(jvale+m-1)) then
                        ifin = ifin-1
                        goto 36
                    end if
                else if (type(1:1) .eq. 'R') then
!              IF ( ZR(JVALE+N-1) .NE. ZR(JVALE+M-1) ) THEN
!             TEST D'EGALITE A PREC PRES
                    if (crit .eq. 'ABSOLU  ') then
                        lok = (abs(zr(jvale+n-1)-zr(jvale+m-1)) .le. prec*abs(zr(jvale+m-1)))
                    else
                        lok = (abs(zr(jvale+n-1)-zr(jvale+m-1)) .le. prec)
                    end if
                    if (.not. lok) then
                        ifin = ifin-1
                        goto 36
                    end if
                else if (type(1:1) .eq. 'C') then
                    if (zc(jvale+n-1) .ne. zc(jvale+m-1)) then
                        ifin = ifin-1
                        goto 36
                    end if
                else if (type(1:3) .eq. 'K80') then
                    if (zk80(jvale+n-1) .ne. zk80(jvale+m-1)) then
                        ifin = ifin-1
                        goto 36
                    end if
                else if (type(1:3) .eq. 'K32') then
                    if (zk32(jvale+n-1) .ne. zk32(jvale+m-1)) then
                        ifin = ifin-1
                        goto 36
                    end if
                else if (type(1:3) .eq. 'K24') then
                    if (zk24(jvale+n-1) .ne. zk24(jvale+m-1)) then
                        ifin = ifin-1
                        goto 36
                    end if
                else if (type(1:3) .eq. 'K16') then
                    if (zk16(jvale+n-1) .ne. zk16(jvale+m-1)) then
                        ifin = ifin-1
                        goto 36
                    end if
                else if (type(1:2) .eq. 'K8') then
                    if (zk8(jvale+n-1) .ne. zk8(jvale+m-1)) then
                        ifin = ifin-1
                        goto 36
                    end if
                end if
            end do
36          continue
            nbuti = ifin-ideb+1
            if (nbuti .gt. 1) then
!           --- ON TESTE AVEC LE PARAMETRE SUIVANT ---
                call tbtr01(tabin, nbpara, lipara(i+1), nbuti, zi(jnume-1+ideb))
                if (lcrit(i+1) (1:2) .eq. 'DE') then
                    do k = 1, nbuti
                        tri2(k) = zi(jnume-1+ideb+k-1)
                    end do
                    do k = 1, nbuti
                        zi(jnume-1+ideb+k-1) = tri2(nbuti-k+1)
                    end do
                end if
            end if
            ii = ifin
            goto 30
32          continue
            goto 104
        end if
    end do
104 continue
!
!     --- ON DUPLIQUE LA TABLE ---
!
!     -- .TBBA :
    call wkvect(nomta2//'.TBBA', base//' V K8', 1, ktbba)
    zk8(ktbba) = base
!
!     -- .TBNP :
    call wkvect(nomta2//'.TBNP', base//' V I', 2, ktbnp)
    zi(ktbnp) = nbpara
    zi(ktbnp+1) = nblign
!
!     -- .TBLP :
    ndim = 4*nbpara
    call jecreo(nomta2//'.TBLP', base//' V K24')
    call jeecra(nomta2//'.TBLP', 'LONMAX', ndim)
    call jeecra(nomta2//'.TBLP', 'LONUTI', ndim)
    call jeveuo(nomta2//'.TBLP', 'E', vk24=nktblp)
    do i = 1, nbpara
        nktblp(1+4*(i-1)) = tblp(1+4*(i-1))
        nktblp(1+4*(i-1)+1) = tblp(1+4*(i-1)+1)
!
        call codent(i, 'D0', knume)
        nomjv = nomta2//'.'//knume
        nktblp(1+4*(i-1)+2) = nomjv
        type = tblp(1+4*(i-1)+1) (1:4)
        call jecreo(nomjv, base//' V '//type)
        call jeecra(nomjv, 'LONMAX', nblign)
        call jeecra(nomjv, 'LONUTI', nblign)
        call jeveuo(nomjv, 'E', kvale)
!
        nomjv = nomta2(1:17)//'LG.'//knume
        nktblp(1+4*(i-1)+3) = nomjv
        call jecreo(nomjv, base//' V I')
        call jeecra(nomjv, 'LONMAX', nblign)
        call jeecra(nomjv, 'LONUTI', nblign)
        call jeveuo(nomjv, 'E', kvall)
!
        nojv2 = tblp(1+4*(i-1)+2)
        nojvl2 = tblp(1+4*(i-1)+3)
        call jeveuo(nojv2, 'L', jvale)
        call jeveuo(nojvl2, 'L', jvall)
!
        do j = 1, nblign
            zi(kvall+j-1) = zi(jvall+zi(jnume+j-1)-1)
            if (type(1:1) .eq. 'I') then
                zi(kvale+j-1) = zi(jvale+zi(jnume+j-1)-1)
            else if (type(1:1) .eq. 'R') then
                zr(kvale+j-1) = zr(jvale+zi(jnume+j-1)-1)
            else if (type(1:1) .eq. 'C') then
                zc(kvale+j-1) = zc(jvale+zi(jnume+j-1)-1)
            else if (type(1:3) .eq. 'K80') then
                zk80(kvale+j-1) = zk80(jvale+zi(jnume+j-1)-1)
            else if (type(1:3) .eq. 'K32') then
                zk32(kvale+j-1) = zk32(jvale+zi(jnume+j-1)-1)
            else if (type(1:3) .eq. 'K24') then
                zk24(kvale+j-1) = zk24(jvale+zi(jnume+j-1)-1)
            else if (type(1:3) .eq. 'K16') then
                zk16(kvale+j-1) = zk16(jvale+zi(jnume+j-1)-1)
            else if (type(1:3) .eq. 'K8') then
                zk8(kvale+j-1) = zk8(jvale+zi(jnume+j-1)-1)
            end if
        end do
!
    end do
!
    call jedetr('&&TBTRTB.TRI')
    AS_DEALLOCATE(vi=tri2)
!
    call jedema()
end subroutine
