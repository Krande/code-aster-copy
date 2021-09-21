! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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

subroutine compStressField(result, model, mater, mateco, cara_elem, list_load, &
                           l_sief_elga, l_strx_elga, rank, time)
!
implicit none
!
#include "asterc/r8vide.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/cochre.h"
#include "asterfort/compStress.h"
#include "asterfort/compStrx.h"
#include "asterfort/dismoi.h"
#include "asterfort/fointe.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/detrsd.h"
#include "asterfort/mecham.h"
#include "asterfort/mechti.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnoch.h"
#include "asterfort/utmess.h"
#include "asterfort/vrcins.h"
#include "asterfort/vrcref.h"
#include "jeveux.h"
!
    character(len=*), intent(in) :: model, cara_elem, list_load, result, mater, mateco
    aster_logical, intent(in) :: l_sief_elga, l_strx_elga
    integer, intent(in) :: rank
    real(kind=8), intent(in) :: time

!
! --------------------------------------------------------------------------------------------------
!
! Compute SIEF_ELGA and STRX_ELGA
!
! --------------------------------------------------------------------------------------------------
!
    integer :: nchar, jchar, nh, jfonc
    integer :: iocc, nfon, iret, nbchre
    real(kind=8) :: alpha, rundf
    character(len=1) :: base, typcoe
    character(len=2) :: codret
    character(len=8) :: k8bla, nomode, kstr
    character(len=19) :: ligrel
    character(len=24) :: charge, fomult, nomfon, charep
    character(len=24) :: chtime, chamgd, chamel, chstrx
    character(len=24) :: chgeom, chcara(18), chharm
    character(len=24) :: chvarc, chvref, compor
    aster_logical :: exipou
    complex(kind=8) :: calpha
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    if( .not.l_sief_elga .and. .not.l_strx_elga ) then
        goto 999
    end if
    rundf=r8vide()
!
! -- INITIALISATIONS
!
    base ='G'
    chtime = ' '
    charge = ' '
    nh = 0
    typcoe = ' '
    charep = ' '
    k8bla = ' '
    chstrx = ' '
    alpha = 0.d0
    calpha = (0.d0 , 0.d0)
    nfon = 0
    chvarc='&&CSF.VARC'
    chvref='&&CSF.VREF'
    exipou = ASTER_FALSE

!
    nomode = model(1:8)
    ligrel = nomode//'.MODELE'
    compor = mater(1:8)//'.COMPOR'
    charge = list_load//'.LCHA'
    fomult = list_load//'.FCHA'
!
!
!   A-t-on des POU_D_EM qui utilisent le champ STRX_ELGA en lineaire
    call dismoi('EXI_STR2', nomode, 'MODELE', repk=kstr)
!
!   A-t-on des VARC
    call dismoi('EXI_VARC', mater, 'CHAM_MATER', repk=k8bla)
!   On interdit provisoirement les POU_D_EM avec les VARC
    if ((k8bla(1:3).eq.'OUI') .and. (kstr(1:3).eq.'OUI')) then
        call utmess('F', 'MECASTATIQUE_1')
    endif
!
! - Compute field
!
    call rsexch(' ', result, 'DEPL', rank, chamgd, iret)
    if (iret .gt. 0) goto 999
!
    call mecham(' ', nomode, cara_elem, nh, chgeom,&
        chcara, chharm, iret)
    if (iret .ne. 0) goto 999
!
    call mechti(chgeom(1:8), time, rundf, rundf, chtime)
    call vrcins(model, mater, cara_elem, time, chvarc(1:19), codret)
    call vrcref(model(1:8), mater(1:8), cara_elem(1:8), chvref(1:19))
!
    if ( l_strx_elga ) then
        call rsexch(' ', result, 'STRX_ELGA', rank, chstrx, iret)
!         -- SI LE CHAMP A DEJA ETE CALCULE :
        if (iret .ne. 0) then
!
            call dismoi('EXI_POUX', model, 'MODELE', repk=k8bla)
            if (k8bla(1:3) .eq. 'OUI') exipou = .true.
!
            if (exipou) then
                call jeveuo(charge, 'L', jchar)
                call jelira(charge, 'LONMAX', nchar)
                call cochre(zk24(jchar), nchar, nbchre, iocc)
                if (nbchre .gt. 1) then
                    call utmess('F', 'MECASTATIQUE_2')
                endif
!
                typcoe = 'R'
                alpha = 1.d0
                if (iocc .gt. 0) then
                    charep = zk24(jchar-1+iocc)
                    call jelira(fomult, 'LONMAX', nfon)
                    if( nfon > 0) then
                        ASSERT(nfon >= iocc)
                        call jeveuo(fomult, 'L', jfonc)
                        nomfon = zk24(jfonc-1+iocc)
                        call fointe('F ', nomfon, 1, ['INST'], [time], alpha, iret)
                    end if
                endif
            endif
!
            call compStrx(nomode, ligrel, compor,&
                        chamgd, chgeom, mateco  , chcara ,&
                        chvarc, chvref, &
                        base  , chstrx, iret  ,&
                        exipou, charep, typcoe, alpha, calpha)
!
            call rsnoch(result, 'STRX_ELGA', rank)
        end if
    endif
    if ( l_sief_elga ) then
        call rsexch(' ', result, 'SIEF_ELGA', rank, chamel, iret)
!           -- SI LE CHAMP A DEJE ETE CALCULE :
        if (iret .ne. 0) then
            call compStress(nomode, ligrel, compor,&
                            chamgd, chgeom, mateco  ,&
                            chcara, chtime, chharm,&
                            chvarc, chvref, chstrx,&
                            base  , chamel, iret  )
            call rsnoch(result, 'SIEF_ELGA', rank)
        end if
    endif
!
    call detrsd("CHAM_ELEM", chvarc)
    call detrsd("CHAM_ELEM", chvref)
!
999 continue
!
    call jedema()
!
end subroutine
