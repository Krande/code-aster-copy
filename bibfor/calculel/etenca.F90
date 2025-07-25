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

subroutine etenca(chinz, ligrlz, iret)
    implicit none
!
! person_in_charge: jacques.pellet at edf.fr
!     ARGUMENTS:
!     ----------
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/mailla.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=19) :: ligrel, chin
    character(len=*) :: ligrlz, chinz
    integer(kind=8) :: iret
! ----------------------------------------------------------------------
!     ENTREES:
!     CHINZ  : NOM DE LA CARTE A ETENDRE
!     LIGRLZ : NOM DU LIGREL SUR LEQUEL ETENDRE LA CARTE
!
!     SORTIES:
!     ON CREE LES OBJETS CHIN.PTMA ET CHIN.PTMS SUR LA VOLATILE
!
!        (SI LES OBJETS .PTMA ET/OU .PTMS EXISTENT DEJA, C'EST QUE LA
!         CARTE EST DEJA ETENDUE. CAS OU LA CARTE EXISTE PLUSIEURS FOIS
!         DANS LA LISTE LCHIN(*). ON SUPPOSE QU'ELLE EST ETENDUE SUR LE
!         BON LIGREL PUISQUE CALCUL CES OBJETS EN FIN D'ALGORITHME.)
!
!
!     IRET : CODE RETOUR   0 --> OK
!                          1 --> PROBLEME
! ----------------------------------------------------------------------
!
!     FONCTIONS EXTERNES:
!     -------------------
!
!     VARIABLES LOCALES:
!     ------------------
    integer(kind=8) :: nma, nms, nbedit, igd, code, ient, i, ii, nb
    integer(kind=8) :: desc, grpma, lima, ialima, illima, jmalut
    integer(kind=8) :: ptma, ptms, noli, iexi
    aster_logical :: bonlig, lalloc
    character(len=8) :: ma
    character(len=24) :: ligri
    integer(kind=8) :: vali(3)
    character(len=24) :: valk
!
!
    call jemarq()
    ligrel = ligrlz
    chin = chinz
!
    iret = 0
    call jeveuo(chin//'.NOLI', 'L', noli)
!
!
!     ----ALLOCATION DES OBJETS PTMA ET PTMS:
    ma = mailla(ligrel)
    call dismoi('NB_MA_MAILLA', ma, 'MAILLAGE', repi=nma)
    call dismoi('NB_MA_SUP', ligrel, 'LIGREL', repi=nms)
!
    call jeexin(ma//'.GROUPEMA', iexi)
    if (iexi .gt. 0) then
        call jeveuo(jexatr(ma//'.GROUPEMA', 'LONUTI'), 'L', jmalut)
    else
        jmalut = 0
    end if
!
    lalloc = .false.
    if (nma .gt. 0) then
        call jeexin(chin//'.PTMA', iexi)
        if (iexi .eq. 0) then
            lalloc = .true.
            call wkvect(chin//'.PTMA', 'V V I', nma, ptma)
        end if
    end if
!
    if (nms .gt. 0) then
        call jeexin(chin//'.PTMS', iexi)
        if (iexi .eq. 0) then
            lalloc = .true.
            call wkvect(chin//'.PTMS', 'V V I', nms, ptms)
        end if
    end if
!
!       -- LA CARTE EST DEJA ETENDUE:
    if (.not. lalloc) goto 999
!
!
!     ----MISE EN MEMOIRE DE LA COLLECTION .LIMA :
    call jeexin(chin//'.LIMA', iexi)
    if (iexi .gt. 0) then
        call jeveuo(chin//'.LIMA', 'L', ialima)
        call jeveuo(jexatr(chin//'.LIMA', 'LONCUM'), 'L', illima)
    end if
!
!
!     ----REMPLISSAGE DES OBJETS PTMA ET PTMS:
    call jeveuo(chin//'.DESC', 'L', desc)
    nbedit = zi(desc-1+3)
    do igd = 1, nbedit
        code = zi(desc-1+3+2*igd-1)
        ient = zi(desc-1+3+2*igd)
!
!        ------- ON NOTE SI LE LIGREL ASSOCIE A IGD EST LE MEME
!        QUE CELUI SUR LEQUEL ON ETEND:
        ligri = zk24(noli-1+igd)
        if (ligri(1:19) .eq. ligrel) then
            bonlig = .true.
        else
            bonlig = .false.
        end if
!
!        ------ GROUPE PREDEFINI "TOUT":
        if (code .eq. 1) then
            do i = 1, nma
                zi(ptma-1+i) = igd
            end do
            goto 60
        end if
        if ((code .eq. -1) .and. bonlig) then
            do i = 1, nms
                zi(ptms-1+i) = igd
            end do
            goto 60
        end if
!
!        ------- GROUPE DE MAILLES DU MAILLAGE:
        if (code .eq. 2) then
            ASSERT(jmalut .ne. 0)
            nb = zi(jmalut-1+ient)
            call jeveuo(jexnum(ma//'.GROUPEMA', ient), 'L', grpma)
            do i = 1, nb
                ii = zi(grpma-1+i)
                zi(ptma-1+ii) = igd
            end do
            goto 60
        end if
!
!        ------- LISTE TARDIVE DE MAILLES ASSOCIEE A LA CARTE:
        if (abs(code) .eq. 3) then
            nb = zi(illima+ient)-zi(illima+ient-1)
            lima = ialima+zi(illima-1+ient)-1
!
            if (code .gt. 0) then
                do i = 1, nb
                    ii = zi(lima-1+i)
                    if (ii .le. 0) then
                        valk = chin
                        vali(1) = ient
                        vali(2) = i
                        vali(3) = ii
                        call utmess('F', 'CALCULEL5_85', sk=valk, ni=3, vali=vali)
                    end if
                    zi(ptma-1+ii) = igd
                end do
            else
                if (bonlig) then
                    do i = 1, nb
                        ii = zi(lima-1+i)
                        ASSERT(ii .lt. 0)
                        zi(ptms-1-ii) = igd
                    end do
                end if
            end if
            goto 60
        end if
60      continue
    end do
999 continue
    call jedema()
end subroutine
