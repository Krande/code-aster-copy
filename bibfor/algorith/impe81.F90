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

subroutine impe81(nomres, impe, basemo)
    implicit none
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterc/r8pi.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/infmaj.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=8) :: nomres, basemo
    character(len=19) :: impe
!
!     BUT:
!       REMPLIR
!
!
!     ARGUMENTS:
!     ----------
!
!      ENTREE :
!-------------
! IN   NOMRES    : NOM DE LA SD_RESULTAT
! IN   IMPE      : NOM DE LA MATRICE D'IMPEDANCE
! IN   BASEMO    : NOM DE LA BASE MODALE DE PROJECTION
!
!      SORTIE :
!-------------
!
! ......................................................................
!
!
!
!
    integer(kind=8) :: i, j, nbmode
    integer(kind=8) :: ldblo, ldbloi, lddesa, lddesm, lddesr, ldrefa, ldrefm
    integer(kind=8) :: ldrefr, ldresa, ldresm, ldresr, ldresi, ldrefi
    integer(kind=8) :: nbdef, nbmodd, nbmods, nfr, nim, ntail
    integer(kind=8) :: ntail1, ntail2, ntail3
    integer(kind=8) :: nk, nc, nm, ldblok, ldbloc, ldblom, nma, nb
!
    real(kind=8) :: partr, parti, partr0, parti0
    real(kind=8) :: amso, dpi, freq
!
    character(len=8) :: blanc
    character(len=16) :: typres, nomcom
    character(len=19) :: impini
    character(len=19) :: impk, impm, impc
!
    data blanc/'        '/
!
!-----------------------------------------------------------------------
!
    call jemarq()
    call infmaj()
!
    call getres(nomres, typres, nomcom)
    dpi = 2.d0*r8pi()
! --- MACR_ELEM_DYNA EXISTANT
    call getvid(' ', 'MACR_ELEM_DYNA', nbval=0, nbret=nma)
!
! --- RECUPERATION DES ARGUMENTS DE LA COMMANDE
!
    call getvr8(' ', 'FREQ_EXTR', scal=freq, nbret=nfr)
    call getvr8(' ', 'AMOR_SOL', scal=amso, nbret=nfr)
    call getvid(' ', 'MATR_IMPE_INIT', scal=impini, nbret=nim)
    call getvid(' ', 'MATR_IMPE_RIGI', scal=impk, nbret=nk)
    call getvid(' ', 'MATR_IMPE_MASS', scal=impm, nbret=nm)
    call getvid(' ', 'MATR_IMPE_AMOR', scal=impc, nbret=nc)
    if (nfr .ne. 0) amso = 2.d0*amso

    if (nma .ne. 0) then
        call getvid(' ', 'BASE_MODALE', scal=basemo, nbret=nb)
        goto 10
    end if
!
    call wkvect(nomres//'.MAEL_RAID_REFE', 'G V K24', 2, ldrefr)
    zk24(ldrefr) = basemo
    zk24(ldrefr+1) = blanc
!
    call wkvect(nomres//'.MAEL_MASS_REFE', 'G V K24', 2, ldrefm)
    zk24(ldrefm) = basemo
    zk24(ldrefm+1) = blanc
!
    call wkvect(nomres//'.MAEL_AMOR_REFE', 'G V K24', 2, ldrefa)
    zk24(ldrefa) = basemo
    zk24(ldrefa+1) = blanc
!
!
10  continue

    call dismoi('NB_MODES_DYN', basemo, 'RESULTAT', repi=nbmodd)
    call dismoi('NB_MODES_STA', basemo, 'RESULTAT', repi=nbmods)
    nbmode = nbmodd+nbmods
!
! --- RECUPERATION DES DIMENSIONS DE LA BASE MODALE
!
    nbdef = nbmode
!
! --- ALLOCATION DE LA MATRICE RESULTAT
!
    ntail = nbdef*(nbdef+1)/2
    if (nma .ne. 0) then
        call jelira(nomres//'.MAEL_RAID_VALE', 'LONMAX', ntail1)
        if (ntail .ne. ntail1) call utmess('F', 'ALGORITH_52')
        goto 11
    end if
    call jecrec(nomres//'.MAEL_RAID_VALE', 'G V R', 'NU', 'DISPERSE', &
                'CONSTANT', 1)
    call jeecra(nomres//'.MAEL_RAID_VALE', 'LONMAX', ntail)
    call jecroc(jexnum(nomres//'.MAEL_RAID_VALE', 1))
11  continue
    call jeveuo(jexnum(nomres//'.MAEL_RAID_VALE', 1), 'E', ldresr)
!
    if (nma .ne. 0) then
        call jelira(nomres//'.MAEL_MASS_VALE', 'LONMAX', ntail2)
        if (ntail .ne. ntail2) call utmess('F', 'ALGORITH_52')
        goto 12
    end if
    call jecrec(nomres//'.MAEL_MASS_VALE', 'G V R', 'NU', 'DISPERSE', &
                'CONSTANT', 1)
    call jeecra(nomres//'.MAEL_MASS_VALE', 'LONMAX', ntail)
    call jecroc(jexnum(nomres//'.MAEL_MASS_VALE', 1))
12  continue
    call jeveuo(jexnum(nomres//'.MAEL_MASS_VALE', 1), 'E', ldresm)
!
    if (nma .ne. 0) then
        call jelira(nomres//'.MAEL_AMOR_VALE', 'LONMAX', ntail3)
        if (ntail .ne. ntail3) call utmess('F', 'ALGORITH_52')
        goto 13
    end if
    call jecrec(nomres//'.MAEL_AMOR_VALE', 'G V R', 'NU', 'DISPERSE', &
                'CONSTANT', 1)
    call jeecra(nomres//'.MAEL_AMOR_VALE', 'LONMAX', ntail)
    call jecroc(jexnum(nomres//'.MAEL_AMOR_VALE', 1))
13  continue
    call jeveuo(jexnum(nomres//'.MAEL_AMOR_VALE', 1), 'E', ldresa)
!
!
!
!        BOUCLE SUR LES COLONNES DE LA MATRICE ASSEMBLEE
!
!
    call jeveuo(jexnum(impe//'.VALM', 1), 'L', ldblo)
    if (nim .ne. 0) call jeveuo(jexnum(impini//'.VALM', 1), 'L', ldbloi)
    if (nk .ne. 0) call jeveuo(jexnum(impk//'.VALM', 1), 'L', ldblok)
    if (nm .ne. 0) call jeveuo(jexnum(impm//'.VALM', 1), 'L', ldblom)
    if (nc .ne. 0) call jeveuo(jexnum(impc//'.VALM', 1), 'L', ldbloc)
    do i = 1, nbmode
!
! --------- BOUCLE SUR LES INDICES VALIDES DE LA COLONNE I
!
        do j = 1, i
!
!
            zr(ldresr+i*(i-1)/2+j-1) = 0.d0
            zr(ldresa+i*(i-1)/2+j-1) = 0.d0
            zr(ldresm+i*(i-1)/2+j-1) = 0.d0
            if (i .gt. nbmodd .and. j .gt. nbmodd) then
!
! ----------- STOCKAGE DANS LE .UALF A LA BONNE PLACE (1 BLOC)
!
                if ((nk+nm+nc) .eq. 0) then
                    partr = dble(zc(ldblo+i*(i-1)/2+j-1))
                    parti = dimag(zc(ldblo+i*(i-1)/2+j-1))
                    zr(ldresr+i*(i-1)/2+j-1) = partr
                    zr(ldresa+i*(i-1)/2+j-1) = (parti-amso*partr)/(dpi* &
                                                                   freq)
                    if (nim .ne. 0) then
                        partr0 = dble(zc(ldbloi+i*(i-1)/2+j-1))
                        parti0 = dimag(zc(ldbloi+i*(i-1)/2+j-1))
                        zr(ldresr+i*(i-1)/2+j-1) = partr0
                        zr(ldresa+i*(i-1)/2+j-1) = (parti-parti0)/(dpi*freq)
                        zr(ldresm+i*(i-1)/2+j-1) = (partr0-partr)/(dpi*freq)**2
                    end if
                else
                    if (nk .ne. 0) zr(ldresr+i*(i-1)/2+j-1) = dble(zc(ldblok+i*(i-1)/2+j-1))
                    if (nm .ne. 0) zr(ldresm+i*(i-1)/2+j-1) = dble(zc(ldblom+i*(i-1)/2+j-1))
                    if (nc .ne. 0) zr(ldresa+i*(i-1)/2+j-1) = dble(zc(ldbloc+i*(i-1)/2+j-1))
                end if
            end if
!
        end do
    end do
!
! --- CREATION DU .DESC
!
    if (nma .ne. 0) then
        goto 14
    end if
    call wkvect(nomres//'.MAEL_RAID_DESC', 'G V I', 3, lddesr)
    zi(lddesr) = 2
    zi(lddesr+1) = nbdef
    zi(lddesr+2) = 2
    call wkvect(nomres//'.MAEL_MASS_DESC', 'G V I', 3, lddesm)
    zi(lddesm) = 2
    zi(lddesm+1) = nbdef
    zi(lddesm+2) = 2
    call wkvect(nomres//'.MAEL_AMOR_DESC', 'G V I', 3, lddesa)
    zi(lddesa) = 2
    zi(lddesa+1) = nbdef
    zi(lddesa+2) = 2
!     INER
    call wkvect(nomres//'.MAEL_INER_REFE', 'G V K24', 2, ldrefi)
    zk24(ldrefr) = basemo
    zk24(ldrefr+1) = blanc
    call wkvect(nomres//'.MAEL_INER_VALE', 'G V R', 3*nbdef, ldresi)
14  continue
!
    call jedema()
end subroutine
