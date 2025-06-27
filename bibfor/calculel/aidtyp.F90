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
subroutine aidtyp(impr)
    implicit none
! person_in_charge: jacques.pellet at edf.fr
! ----------------------------------------------------------------------
!    BUT:
!       ECRIRE SUR LE FICHIER "IMPR"
!       LES COUPLES (OPTION, TYPE_ELEMENT) POSSIBLES DANS LES CATALOGUES
!      (POUR VERIFIER LA COMPLETUDE)
! ----------------------------------------------------------------------
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
!
    character(len=16) :: nophen, note, noop, nomodl
    character(len=80) :: ligne
    integer(kind=8) :: impr, nbtm, nbphen, nbte, nbop
    integer(kind=8) :: iop, iphen, nbmodl, imodl, iamodl, itm, ite, ioptte
    integer(kind=8) :: iaopmo, nucalc
    integer(kind=8), pointer :: optte(:) => null()
    integer(kind=8), pointer :: vnbop(:) => null()
    integer(kind=8), pointer :: vnbte(:) => null()
    character(len=16), pointer :: nop2(:) => null()
    character(len=80), pointer :: not2(:) => null()
!
!
    call jemarq()
!
!
!
    ligne(1:40) = '========================================'
    ligne(41:80) = '========================================'
!
    call jelira('&CATA.TM.NOMTM', 'NOMMAX', nbtm)
    call jelira('&CATA.PHENOMENE', 'NOMUTI', nbphen)
    call jeveuo('&CATA.TE.OPTTE', 'L', vi=optte)
    call jelira('&CATA.TE.NOMTE', 'NOMUTI', nbte)
    call jelira('&CATA.OP.NOMOPT', 'NOMUTI', nbop)
!
    AS_ALLOCATE(vi=vnbop, size=nbop)
    AS_ALLOCATE(vi=vnbte, size=nbte)
    AS_ALLOCATE(vk80=not2, size=nbte)
    AS_ALLOCATE(vk16=nop2, size=nbop)
!
!
!     -- REMPLISSAGE DE .NOP2:
!     ------------------------
    do iop = 1, nbop
        call jenuno(jexnum('&CATA.OP.NOMOPT', iop), noop)
        nop2(iop) = noop
    end do
!
!
!     -- REMPLISSAGE DE .NOT2:
!     ------------------------
    do iphen = 1, nbphen
        call jenuno(jexnum('&CATA.PHENOMENE', iphen), nophen)
        call jelira('&CATA.'//nophen, 'NUTIOC', nbmodl)
        do imodl = 1, nbmodl
            call jeveuo(jexnum('&CATA.'//nophen, imodl), 'L', iamodl)
            call jenuno(jexnum('&CATA.'//nophen(1:13)//'.MODL', imodl), nomodl)
            do itm = 1, nbtm
                ite = zi(iamodl-1+itm)
                if (ite .eq. 0) goto 3
                call jenuno(jexnum('&CATA.TE.NOMTE', ite), note)
                not2(ite) = nophen//' '//nomodl//' '//note
3               continue
            end do
        end do
    end do
!
!     ON COMPLETE .NOT2 AVEC LES ELEMENTS N'APPARTENANT A AUCUNE
!        MODELISATION NI PHENOMENE:
    do ite = 1, nbte
        if (not2(ite) (1:1) .eq. ' ') then
            call jenuno(jexnum('&CATA.TE.NOMTE', ite), note)
            not2(ite) (35:50) = note
        end if
    end do
!
!
!     -- ECRITURE DES COUPLES (TE,OPT)
!     --------------------------------
    write (impr, '(A80)') ligne
    write (impr, *) ' NOMBRE D''OPTION        : ', nbop
    write (impr, *) ' NOMBRE DE TYPE_ELEMENT : ', nbte
    write (impr, '(A80)') ligne
    do ite = 1, nbte
        do iop = 1, nbop
            ioptte = optte(nbop*(ite-1)+iop)
            if (ioptte .eq. 0) goto 101
            call jeveuo(jexnum('&CATA.TE.OPTMOD', ioptte), 'L', iaopmo)
            nucalc = zi(iaopmo)
            if (nucalc .eq. 0) goto 101
            vnbte(ite) = vnbte(ite)+1
            vnbop(iop) = vnbop(iop)+1
            write (impr, 1001) not2(ite) (1:50), nop2( &
                iop), nucalc
101         continue
        end do
    end do
!
!
!     -- ECRITURE RESUME TYPE_ELEMENT:
!     --------------------------------
    write (impr, '(A80)') ligne
    write (impr, *) ' RESUME TYPE_ELEMENTS : '
    do ite = 1, nbte
        write (impr, 1001) not2(ite) (1:50), ' NB_OPT_CALC: ', &
            vnbte(ite)
    end do
!
!
!     -- ECRITURE RESUME OPTIONS:
!     ---------------------------
    write (impr, '(A80)') ligne
    write (impr, *) ' RESUME OPTIONS : '
    do iop = 1, nbop
        write (impr, *) nop2(iop), ' NB_TYP_CALC: ', vnbop( &
            iop)
    end do
    write (impr, '(A80)') ligne
!
!
! --- MENAGE
!
    AS_DEALLOCATE(vi=vnbop)
    AS_DEALLOCATE(vi=vnbte)
    AS_DEALLOCATE(vk80=not2)
    AS_DEALLOCATE(vk16=nop2)
!
!
1001 format(a50, 1x, a16, 1x, i5)
    call jedema()
end subroutine
