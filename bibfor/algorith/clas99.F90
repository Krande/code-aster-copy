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

subroutine clas99(nomres)
    implicit none
!  P. RICHARD   DATE 09/07/91
!-----------------------------------------------------------------------
!  BUT : ROUTINE DE CREATION D'UNE BASE MODALE CLASSIQUE
!        BASE MODALE DE TYPE MIXTE CRAIG-BAMPTON, MAC-NEAL OU AUCUN
!
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/camoat.h"
#include "asterfort/camoch.h"
#include "asterfort/camoco.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/moco99.h"
#include "asterfort/rscrsd.h"
#include "asterfort/rsorac.h"
#include "asterfort/wkvect.h"
!
!
!
    character(len=8) :: nomres, intf, kbid
    character(len=19):: numddl, raid, mass, raidlt
    complex(kind=8)  :: cbid
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! --- RECUPERATION DES CONCEPTS AMONT
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ibid, ii, inor, lrang
    integer(kind=8) :: ltmome, ltnbmo, nbid, nbnmaxmode, nbmod, nbmodo(1), nbmoma
    integer(kind=8) :: nbmome, nbmout, nbsdd, nmaxmode
    real(kind=8) :: bid, ebid
    integer(kind=8), pointer :: idc_desc(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
!
    call dismoi('REF_RIGI_PREM', nomres, 'RESU_DYNA', repk=raid, arret='C')
    call dismoi('REF_MASS_PREM', nomres, 'RESU_DYNA', repk=mass, arret='C')
    call dismoi('NUME_DDL', nomres, 'RESU_DYNA', repk=numddl)
    call dismoi('REF_INTD_PREM', nomres, 'RESU_DYNA', repk=intf)
!
!----ON AJOUT .NUME POUR OBTENIR LE NUME_EQUA
    numddl(15:19) = '.NUME'
!
! --- RECUPERATION DU NOMBRE DE MODE_MECA A PRENDRE EN COMPTE
!
    call getvid('CLASSIQUE', 'MODE_MECA', iocc=1, nbval=0, nbret=nbmome)
    nbmome = -nbmome
!
! --- RECUPERATION DU NOMBRE DE VALEURS DE LA LISTE NMAX_MODE
!
    call getvis('CLASSIQUE', 'NMAX_MODE', iocc=1, nbval=0, nbret=nbnmaxmode)
    nbnmaxmode = -nbnmaxmode
!
! --- CREATION DES OBJETS TEMPORAIRES
!
    call wkvect('&&CLAS99.LIST.MODE_MECA', 'V V K8', nbmome, ltmome)
    call wkvect('&&CLAS99.LIST.NBMOD', 'V V I', nbmome, ltnbmo)
!
    call getvid('CLASSIQUE', 'MODE_MECA', iocc=1, nbval=nbmome, vect=zk8(ltmome), &
                nbret=ibid)
!
    if (nbnmaxmode .ge. 1) then
!      length of NMAX_MODE list will be equal to length of MODE_MECA list
        call getvis('CLASSIQUE', 'NMAX_MODE', iocc=1, nbval=nbmome, vect=zi(ltnbmo), &
                    nbret=nbid)
        if (nbnmaxmode .eq. 1) then
            nmaxmode = 0
            call getvis('CLASSIQUE', 'NMAX_MODE', iocc=1, scal=nmaxmode, nbret=ibid)
            do i = 1, nbmome
                zi(ltnbmo+i-1) = nmaxmode
            end do
        end if
    end if
!
! --- DETERMINATION DU NOMBRE TOTAL DE MODES PROPRES DE LA BASE
!
    nbmod = 0
    nbmoma = 0
!
    do i = 1, nbmome
        call rsorac(zk8(ltmome-1+i), 'LONUTI', ibid, bid, kbid, &
                    cbid, ebid, 'ABSOLU', nbmodo, 1, &
                    nbid)
!
!       if NMAX_MODE is set by the user, one takes it into account
!       otherwise one takes all modes in each MODE_MECA
        if (nbnmaxmode .ge. 1) then
            nbmout = zi(ltnbmo-1+i)
            if (nbmout .lt. nbmodo(1)) nbmodo(1) = nbmout
        end if
!
        zi(ltnbmo+i-1) = nbmodo(1)
        nbmoma = max(nbmoma, nbmodo(1))
        nbmod = nbmod+nbmodo(1)
    end do
!
    call wkvect('&&CLAS99.NUME.RANG', 'V V I', nbmoma, lrang)
    do ii = 1, nbmoma
        zi(lrang+ii-1) = ii
    end do
!
!
! --- DETERMINATION NOMBRE TOTAL DE MODES ET DEFORMEES
!
    ASSERT(intf(1:8) .ne. ' ')
    call jeveuo(intf//'.IDC_DESC', 'L', vi=idc_desc)
    nbsdd = nbmod+idc_desc(5)
!      NBSDD1=ZI(LLDESC+4)
!
!
!
! --- NOMBRE DE DEFORMEES STATIQUES A CALCULER
!
!
! --- ALLOCATION DE LA STRUCTURE DE DONNEES MODE_MECA
!
    call rscrsd('G', nomres, 'MODE_MECA', nbsdd)
    raidlt = ' '
!
! --- COPIE DES MODES DYNAMIQUES
!
    inor = 1
    do i = 1, nbmome
        call moco99(nomres, zk8(ltmome+i-1), zi(ltnbmo+i-1), zi(lrang), inor, &
                    .true._1)
    end do
    if (nbmoma .gt. 0) call jedetr('&&CLAS99.NUME.ORD')
    if (nbmome .gt. 0) call jedetr('&&CLAS99.LIST.MODE_MECA')
    if (nbmome .gt. 0) call jedetr('&&CLAS99.LIST.NBMOD')
!      CALL UTIMSD(6,2,.TRUE.,.TRUE.,NOMRES(1:8),1,'G')
!
! --- CALCUL DES MODES D'ATTACHE
    call camoat(nomres, numddl, intf, raid, raidlt, &
                inor)
!
! --- CALCUL DES MODES CONTRAINTS
    call camoco(nomres, numddl, intf, raid, raidlt, &
                inor)
!
! --- CALCUL DES MODES CONTRAINTS HARMONIQUES
    call camoch(nomres, numddl, intf, raid, mass, &
                raidlt, inor)
!
! --- DESTRUCTION MATRICE FACTORISEE
!
    if (raidlt(1:1) .ne. ' ') then
        call detrsd('MATR_ASSE', raidlt)
        raidlt = ' '
    end if
!
    call jedema()
end subroutine
