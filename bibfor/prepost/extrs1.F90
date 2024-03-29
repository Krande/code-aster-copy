! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
subroutine extrs1(resu0, nbrang, nuordr, nbpara, nompar, &
                  nbarch, nuarch, nbexcl, chexcl, nbnosy)
    implicit none
#include "jeveux.h"
#include "asterc/isnnem.h"
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/extrs3.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsmena.h"
#include "asterfort/rsutch.h"
!
    integer :: nbrang, nuordr(*), nbarch, nbpara, nuarch(*), nbexcl, nbnosy
    character(len=16) :: nompar(*), chexcl(*)
    character(len=*) :: resu0
!     ------------------------------------------------------------------
! person_in_charge: jacques.pellet at edf.fr
!     EXTR_RESU / ARCHIVAGE + REUSE
! ----------------------------------------------------------------------
!
!
! 0.3. ==> VARIABLES LOCALES
!
!
    integer :: irang, i, j, k, jtach, iadin, iadou, ire1
    integer :: iundf, iordr
    real(kind=8) :: rundf
    character(len=3) :: type
    character(len=16) :: nomsym
    character(len=16) :: nopara
    character(len=19) :: nomsdr
    character(len=19) :: chamin, nomch1, nomch2
    integer, pointer :: ordr(:) => null()
!     ------------------------------------------------------------------
!
    call jemarq()
    rundf = r8vide()
    iundf = isnnem()
!
    nomsdr = resu0
!
!
!     1. -- ON COMPACTE LES CHAMPS ARCHIVES :
!     -------------------------------------------------------
    do i = 1, nbnosy
        call jenuno(jexnum(nomsdr//'.DESC', i), nomsym)
        call jeveuo(jexnum(nomsdr//'.TACH', i), 'E', jtach)
        do j = 1, nbexcl
            if (chexcl(j) .eq. nomsym) then
                do k = 1, nbrang
                    if (zk24(jtach+k-1) (1:1) .eq. ' ') goto 10
                    call rsexch('F', nomsdr, nomsym, nuordr(k), chamin, &
                                ire1)
                    call detrsd('CHAMP_GD', chamin)
                    zk24(jtach+k-1) = ' '
10                  continue
                end do
                goto 50
!
            end if
        end do
!
        irang = 0
        do j = 1, nbrang
            if (zk24(jtach+j-1) (1:1) .eq. ' ') goto 30
            call rsexch('F', nomsdr, nomsym, nuordr(j), chamin, &
                        ire1)
            if (nuarch(j) .eq. 0) then
                call detrsd('CHAMP_GD', chamin)
            else
                irang = irang+1
                zk24(jtach+irang-1) = chamin
            end if
30          continue
        end do
!
        do k = irang+1, nbrang
            zk24(jtach+k-1) = ' '
        end do
50      continue
    end do
!
!
!     2. -- ON COMPACTE LES PARAMETRES ARCHIVES :
!     -------------------------------------------
    irang = 0
    do i = 1, nbrang
        if (nuarch(i) .eq. 0) goto 70
        irang = irang+1
        do j = 1, nbpara
            nopara = nompar(j)
            call rsadpa(nomsdr, 'L', 1, nopara, nuordr(i), &
                        1, sjv=iadin, styp=type, istop=0)
            call extrs3(nomsdr, nopara, irang, 'E', 1, &
                        type, iadou)
            if (type(1:1) .eq. 'I') then
                zi(iadou) = zi(iadin)
            else if (type(1:1) .eq. 'R') then
                zr(iadou) = zr(iadin)
            else if (type(1:1) .eq. 'C') then
                zc(iadou) = zc(iadin)
            else if (type(1:3) .eq. 'K80') then
                zk80(iadou) = zk80(iadin)
            else if (type(1:3) .eq. 'K32') then
                zk32(iadou) = zk32(iadin)
            else if (type(1:3) .eq. 'K24') then
                zk24(iadou) = zk24(iadin)
            else if (type(1:3) .eq. 'K16') then
                zk16(iadou) = zk16(iadin)
            else if (type(1:2) .eq. 'K8') then
                zk8(iadou) = zk8(iadin)
            end if
        end do
70      continue
    end do
    ASSERT(irang .eq. nbarch)
!
!
!     3. -- ON COMPACTE LES NUME_ORDRE ARCHIVES :
!     -------------------------------------------
    call jeecra(nomsdr//'.ORDR', 'LONUTI', nbarch)
    call jeveuo(nomsdr//'.ORDR', 'E', vi=ordr)
    irang = 0
    do i = 1, nbrang
        if (nuarch(i) .eq. 0) goto 80
        irang = irang+1
        ordr(irang) = nuordr(i)
80      continue
    end do
    ASSERT(irang .eq. nbarch)
!
!
!     4. -- ON MET A "ZERO" LES IRANG INUTILISES :
!     -------------------------------------------------
    do irang = nbarch+1, nbrang
        ordr(irang) = iundf
        do j = 1, nbpara
            nopara = nompar(j)
            call extrs3(nomsdr, nopara, irang, 'E', 1, &
                        type, iadou)
            if (type(1:1) .eq. 'I') then
                zi(iadou) = iundf
            else if (type(1:1) .eq. 'R') then
                zr(iadou) = rundf
            else if (type(1:1) .eq. 'C') then
                zc(iadou) = dcmplx(rundf, rundf)
            else if (type(1:3) .eq. 'K80') then
                zk80(iadou) = ' '
            else if (type(1:3) .eq. 'K32') then
                zk32(iadou) = ' '
            else if (type(1:3) .eq. 'K24') then
                zk24(iadou) = ' '
            else if (type(1:3) .eq. 'K16') then
                zk16(iadou) = ' '
            else if (type(1:2) .eq. 'K8') then
                zk8(iadou) = ' '
            end if
        end do
    end do
!
!
!     5. -- IL FAUT RENOMMER LES CHAMPS POUR QU'ILS RESPECTENT
!           LA REGLE DE NOMMAGE DE RSUTCH.F :
!     ---------------------------------------------------------
    do i = 1, nbnosy
        call jenuno(jexnum(nomsdr//'.DESC', i), nomsym)
        call jeveuo(jexnum(nomsdr//'.TACH', i), 'E', jtach)
        do j = 1, nbarch
            iordr = ordr(j)
            nomch1 = zk24(jtach-1+j)
            if (nomch1 .eq. ' ') goto 41
            call rsutch(nomsdr, nomsym, iordr, nomch2, .false._1)
            if (nomch1 .ne. nomch2) then
                call copisd('CHAMP', 'G', nomch1, nomch2)
                call detrsd('CHAMP', nomch1)
                zk24(jtach-1+j) = nomch2
            end if
41          continue
        end do
    end do
!
!
!     6. -- IL FAUT ENCORE DETRUIRE LES SCORIES INUTILES :
!     ----------------------------------------------------
    call rsmena(nomsdr)
!
    call jedema()
end subroutine
