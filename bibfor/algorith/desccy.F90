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

subroutine desccy(nomres)
!    P. RICHARD     DATE 07/03/91
!-----------------------------------------------------------------------
!  BUT:  CREATION DE LA NUMEROTATION GENERALISEE POUR LE PROBLEME
    implicit none
!        CYCLIQUE
!-----------------------------------------------------------------------
!
! NOMRES   /I/: NOM UTILISATEUR DU RESULTAT
! BASMOD   /I/: NOM UTILISATEUR DE L'EVENTUELLE BASE MODALE (OU BLANC)
! RESCYC   /I/: NOM UTILISATEUR EVENTUEL CONCEPT MODE CYCLIQUE(OU BLANC)
! NUMCYC   /O/: NOM K24 DE LA NUMEROTATION RESULTAT
!
!-----------------------------------------------------------------------
!
#include "jeveux.h"
#include "asterfort/bmnodi.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvis.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: vali(3)
!
!
!
    character(len=8) :: intf, kbid, basmod, nomres
    character(len=24) :: noeint
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: ibid(1), ldnoea, ldnoed, ldnoeg, ldnumg
    integer(kind=8) :: lldesc, nba, nbd, nbda, nbdd, nbdg
    integer(kind=8) :: nbg, nbmcal, nbmod, nbmod1, nbmod2, nbnot, nboc
    integer(kind=8) :: nbtemp, numa, numd, numg, n1
    integer(kind=8), pointer :: cycl_nuin(:) => null()
    character(len=24), pointer :: cycl_refe(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
    kbid = ' '
!-----------------------------------------------------------------------
!
!------------------RECUPERATION CONCEPT AMONT---------------------------
!
    call jeveuo(nomres//'.CYCL_REFE', 'L', vk24=cycl_refe)
    intf = cycl_refe(2)
    basmod = cycl_refe(3)
!
!-----RECUPERATION NUMEROS INTERFACES DROITE GAUCHE ET AXE--------------
!
    call jeveuo(nomres//'.CYCL_NUIN', 'L', vi=cycl_nuin)
    numd = cycl_nuin(1)
    numg = cycl_nuin(2)
    numa = cycl_nuin(3)
!
!----------RECUPERATION DU DESCRIPTEUR DES DEFORMEES STATIQUES----------
!
    call jeveuo(intf//'.IDC_DEFO', 'L', lldesc)
    call jelira(intf//'.IDC_DEFO', 'LONMAX', nbnot)
    nbnot = nbnot/3
!
!--------RECUPERATION DES DEFINITIONS DES INTERFACES DROITE ET GAUCHE---
!
    noeint = intf//'.IDC_LINO'
!
    call jeveuo(jexnum(noeint, numd), 'L', ldnoed)
    call jelira(jexnum(noeint, numd), 'LONMAX', nbd)
!
    call jeveuo(jexnum(noeint, numg), 'L', ldnoeg)
    call jelira(jexnum(noeint, numg), 'LONMAX', nbg)
!
    if (numa .gt. 0) then
        call jeveuo(jexnum(noeint, numa), 'L', ldnoea)
        call jelira(jexnum(noeint, numa), 'LONMAX', nba)
    end if
!
    if (nbg .ne. nbd) then
        vali(1) = nbd
        vali(2) = nbg
        call utmess('F', 'ALGORITH12_79', ni=2, vali=vali)
    end if
!
!------COMPTAGE DEFORMEES STATIQUES INTERFACE DROITE GAUCHE-------------
!
    call bmnodi(basmod, kbid, '        ', numd, 0, &
                ibid, nbdd)
    kbid = ' '
    call bmnodi(basmod, kbid, '        ', numg, 0, &
                ibid, nbdg)
!
!--------------TEST SUR NOMBRE DE DDL AUX INTERFACES--------------------
!
    if (nbdd .ne. nbdg) then
        vali(1) = nbdd
        vali(2) = nbdg
        call utmess('F', 'ALGORITH12_80', ni=2, vali=vali)
    end if
!
!-----COMPTAGE NOMBRE DEFORMEES STATIQUE SUR EVENTUELLE INTERFACE AXE---
!
    nbda = 0
    if (numa .gt. 0) then
        kbid = ' '
        call bmnodi(basmod, kbid, '        ', numa, 0, &
                    ibid, nbda)
    else
        nbda = 0
    end if
!
!--------DETERMINATION DU NOMBRE DE MODES PROPRES DE LA BASE------------
!
!  NOMBRE DE MODES DEMANDES
!
    call getvis('   ', 'NB_MODE', iocc=1, scal=nbmod1, nbret=n1)
!
!  NOMBRE DE MODES EXISTANTS
    call dismoi('NB_MODES_DYN', basmod, 'RESULTAT', repi=nbmod2)
!
!  TEST
!
    if (n1 .eq. 0) then
        nbmod1 = nbmod2
    end if

    nbmod = min(nbmod1, nbmod2)
    if (nbmod .eq. 0) then
        call utmess('F', 'ALGORITH12_81')
    end if
!
!---------DETERMINATION DU NOMBRE DE MODES PROPRES A CALCULER-----------
!
    call getvis('CALCUL', 'NMAX_FREQ', iocc=1, nbval=0, nbret=nboc)
!
    if (nboc .eq. 0) then
        nbmcal = nbmod
    else
        call getvis('CALCUL', 'NMAX_FREQ', iocc=1, scal=nbmcal, nbret=ibid(1))
    end if
!
    if (nbmcal .gt. nbmod) then
        nbtemp = nbmcal-nbmod
        vali(1) = nbmcal
        vali(2) = nbmod
        vali(3) = nbtemp
        call utmess('A', 'ALGORITH12_82', ni=3, vali=vali)
    end if
!
!----------------ALLOCATION DE L'OBJET .DESC----------------------------
!
    call wkvect(nomres//'.CYCL_DESC', 'G V I', 4, ldnumg)
!
!------------------CREATION DE LA NUMEROTATION--------------------------
!
    zi(ldnumg) = nbmod
    zi(ldnumg+1) = nbdd
    zi(ldnumg+2) = nbda
    zi(ldnumg+3) = nbmcal
!
    call jedema()
end subroutine
