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
subroutine tbcarapou(nomta, nomsec, nbcarac, nomcarac, valcarac, okcarac)
    implicit none
!
! --------------------------------------------------------------------------------------------------
!   Recherche dans la table des caractéristiques d'une section
!
!           nomta    : nom de la table issue de macr_cara_poutre
!           nomsec   : nom de la section ==> paramètre : LIEU
!           nbcarac  : nombre de paramètres à trouver
!           nomcarac : nom des paramètres à trouver
!           valcarac : valeur des paramètres
!           okcarac  : si 0 le paramètre est trouvé
!                         1 le paramètre n'est pas trouvé ==> valcarac = NaN
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), intent(in)             :: nbcarac
    character(len=*), intent(in)    :: nomta, nomsec, nomcarac(nbcarac)
    integer(kind=8), intent(out)            :: okcarac(nbcarac)
    real(kind=8), intent(out)       :: valcarac(nbcarac)
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8nnem.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
!
! --------------------------------------------------------------------------------------------------
    integer(kind=8) :: nbcolo, nblign, collieu, ii, iret, jj, iisec, itabl
    character(len=19) :: nomtab
    character(len=24) :: kpara, typca

    integer(kind=8), pointer :: colpara(:) => null()
    integer(kind=8), pointer :: tbnp(:) => null()
    character(len=24), pointer :: tblp(:) => null()
! --------------------------------------------------------------------------------------------------
    call jemarq()
!
    nomtab = nomta
    call jeexin(nomtab//'.TBBA', iret)
    if (iret .eq. 0) then
        call utmess('F', 'UTILITAI4_64')
    end if
!
    call jeveuo(nomtab//'.TBNP', 'E', vi=tbnp)
    nbcolo = tbnp(1)
    nblign = tbnp(2)
    if (nbcolo .eq. 0) then
        call utmess('F', 'UTILITAI4_65')
    end if
    if (nblign .eq. 0) then
        call utmess('F', 'UTILITAI4_66')
    end if
!
    call jeveuo(nomtab//'.TBLP', 'L', vk24=tblp)
!   On recherche :
!       - le paramètre LIEU
!       - les paramètres nomcarac et on vérifie qu'ils sont R
    collieu = 0
    AS_ALLOCATE(size=nbcarac, vi=colpara)
    colpara(1:nbcarac) = 0
    okcarac(:) = 1
    do ii = 1, nbcolo
        kpara = tblp((ii-1)*4+1)
        typca = tblp((ii-1)*4+2)
        if (kpara(1:4) .eq. 'LIEU') then
            collieu = ii
        end if
        jja: do jj = 1, nbcarac
            if (kpara .eq. nomcarac(jj)) then
                if (typca .ne. 'R') exit jja
                colpara(jj) = ii
                okcarac(jj) = 0
                exit jja
            end if
        end do jja
    end do
!   Si on n'a pas trouve LIEU
    if (collieu .eq. 0) then
        call utmess('F', 'UTILITAI4_63', sk='LIEU')
    end if
!   On recherche la ligne avec LIEU='nomsec'
    call jeveuo(tblp((collieu-1)*4+3), 'L', itabl)
    typca = tblp((collieu-1)*4+2)
    iisec = 0
    if (typca(1:2) .eq. 'K8') then
        ii8: do ii = 1, nblign
            if (zk8(itabl-1+ii) .eq. nomsec) then
                iisec = ii
                exit ii8
            end if
        end do ii8
    else if (typca(1:3) .eq. 'K24') then
        ii24: do ii = 1, nblign
            if (zk24(itabl-1+ii) (1:8) .eq. nomsec) then
                iisec = ii
                exit ii24
            end if
        end do ii24
    end if
    if (iisec .eq. 0) then
        call utmess('F', 'UTILITAI4_63', sk=nomsec)
    end if
!   Le paramètre jj est sur la ligne 'iisec' et sur la colonne colpara(jj)
    jj1: do jj = 1, nbcarac
        if (okcarac(jj) .eq. 0) then
            call jeveuo(tblp((colpara(jj)-1)*4+3), 'L', itabl)
            valcarac(jj) = zr(itabl-1+iisec)
!             write(*,*) 'DBG : tbcarapou >> ', nomcarac(jj),  valcarac(jj)
        else
            valcarac(jj) = r8nnem()
        end if
    end do jj1
!
    AS_DEALLOCATE(vi=colpara)
    call jedema()
end subroutine
