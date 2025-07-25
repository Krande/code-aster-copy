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
subroutine rc36cm(iocc, etat, nbma, listma, nbchar, &
                  lichar, chmome)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/cesfus.h"
#include "asterfort/cesqua.h"
#include "asterfort/cesred.h"
#include "asterfort/detrsd.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: iocc, nbma, listma(*), nbchar, lichar(*)
    character(len=1) :: etat
    character(len=24) :: chmome
!
!     OPERATEUR POST_RCCM, TRAITEMENT DE FATIGUE_B3600
!     CALCUL DU TORSEUR PAR SOMMATION ALGEBRIQUE DES TORSEURS
!     CORRESPONDANT AUX DIFFERENTS CAS DE CHARGE DE LA SITUATION
!
! IN  : IOCC   : NUMERO D'OCCURRENCE DE SITUATION
! IN  : ETAT   : ETAT STABILISE A OU B POUR LE MESSAGE D'ERREUR
! IN  : NBMA   : NOMBRE DE MAILLES D'ANALYSE
! IN  : LISTMA : LISTE DES MAILLES D'ANALYSE
! IN  : NBCHAR : NOMBRE DE CAS DE CHARGE POUR UN ETAT STABILISE
! IN  : LICHAR : LISTE DES CAS DE CHARGE POUR UN ETAT STABILISE
! OUT : CHNOME : TORSEUR RESULTAT
!     ------------------------------------------------------------------
!
    integer(kind=8) :: nbresu, nbcmp, icha, ir
    integer(kind=8) :: vali(2)
    aster_logical :: seisme, autre
    character(len=8) :: nocmp(3)
    character(len=24) :: chams0
    complex(kind=8) :: cbid
    character(len=24), pointer :: lich(:) => null()
    aster_logical, pointer :: licm(:) => null()
    real(kind=8), pointer :: licr(:) => null()
    character(len=24), pointer :: champ(:) => null()
    integer(kind=8), pointer :: nume_char(:) => null()
! DEB ------------------------------------------------------------------
    call jemarq()
    cbid = (0.d0, 0.d0)
!
    call jeveuo('&&RC3600.NUME_CHAR', 'L', vi=nume_char)
    call jeveuo('&&RC3600.CHAMP', 'L', vk24=champ)
    call jelira('&&RC3600.NUME_CHAR', 'LONMAX', nbresu)
!
    nbcmp = 3
    nocmp(1) = 'MT'
    nocmp(2) = 'MFY'
    nocmp(3) = 'MFZ'
!
    seisme = .false.
    autre = .false.
!
    AS_ALLOCATE(vk24=lich, size=nbchar)
    AS_ALLOCATE(vl=licm, size=nbchar)
    AS_ALLOCATE(vr=licr, size=nbchar)
!
    do icha = 1, nbchar, 1
        do ir = 1, nbresu, 1
            if (lichar(icha) .eq. nume_char(ir)) goto 114
        end do
        vali(1) = iocc
        vali(2) = lichar(icha)
        call utmess('F', 'POSTRCCM_28', ni=2, vali=vali)
114     continue
        if (etat .eq. 'S') then
            seisme = .true.
        else
            autre = .true.
        end if
        lich(icha) = champ(ir)
        licm(icha) = .true.
        licr(icha) = 1.d0
    end do
!
    if (seisme .and. autre) then
        call utmess('F', 'POSTRCCM_29', si=iocc)
    end if
!
    if (nbchar .eq. 1) then
        chams0 = lich(1)
        call cesred(chams0, nbma, listma, nbcmp, nocmp, &
                    'V', chmome)
    else
!
        chams0 = '&&RC36CM.CHAMS0'
        if (autre) then
            call cesfus(nbchar, lich, licm, licr, [cbid], &
                        .false._1, 'V', chams0)
        else
            call cesqua(nbchar, lich, licm, 'V', chams0)
        end if
        call cesred(chams0, nbma, listma, nbcmp, nocmp, &
                    'V', chmome)
        call detrsd('CHAM_ELEM_S', chams0)
    end if
!
    AS_DEALLOCATE(vk24=lich)
    AS_DEALLOCATE(vl=licm)
    AS_DEALLOCATE(vr=licr)
!
    call jedema()
end subroutine
