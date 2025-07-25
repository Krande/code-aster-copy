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

subroutine op0188()
!
! person_in_charge: samuel.geniaut at edf.fr
!
    implicit none
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/alcart.h"
#include "asterfort/assert.h"
#include "asterfort/cncinv.h"
#include "asterfort/cnocns.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/nocart.h"
#include "asterfort/tecart.h"
#include "asterfort/wkvect.h"
!
!
! ----------------------------------------------------------------------
!
! OPERATEUR RAFF_XFEM_ZONE
!
! CALCUL D'UN INDICATEUR BINAIRE (APPELEE PAR RAFF_XFEM)
!
!
    integer(kind=8) :: ibid, iret, i, j, ino, nuno, numa, nbnozo
    integer(kind=8) ::   ncmp, jmafon, nmafon, jma, nbma, nbno, nbmali
    integer(kind=8) ::   jnoeu, nbmac, jadr, adrvlc, acncin
    integer(kind=8) :: idlima, nbmazo
    real(kind=8) :: rayon, dist
    character(len=8) :: fiss, ma, chout
    character(len=16) :: typdis, k16bid
    character(len=19) :: carte, cnslt, cnsln
    character(len=24) :: mafond, listma, cnxinv, lisnoz, lismaz
    real(kind=8), pointer :: valv(:) => null()
    character(len=8), pointer :: vncmp(:) => null()
    real(kind=8), pointer :: lsn(:) => null()
    real(kind=8), pointer :: lst(:) => null()
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
!     ------------------------------------
!     1) INITIALISATIONS
!     ------------------------------------
!
!     NOM DU CONCEPT EN SORTIE : CHOUT
    call getres(chout, k16bid, k16bid)
!
!     NOM DU CONCEPT FISSURE
    call getvid(' ', 'FISSURE', scal=fiss, nbret=ibid)
!
!     RECUP DU RAYON DE LA ZONE
    call getvr8(' ', 'RAYON', scal=rayon, nbret=ibid)
!
!     TYPE DE SD_FISS_XFEM EN ENTREE (FISSURE/INTERFACE)
    call dismoi('TYPE_DISCONTINUITE', fiss, 'FISS_XFEM', repk=typdis)
!
!     MAILLAGE ASSOCIE A LA FISSURE/INTERFACE
    call dismoi('NOM_MAILLA', fiss, 'FISS_XFEM', repk=ma)
    call dismoi('NB_MA_MAILLA', ma, 'MAILLAGE', repi=nbma)
    call dismoi('NB_NO_MAILLA', ma, 'MAILLAGE', repi=nbno)
!
!     INITIALISATION DE LA CARTE AVEC LA VALEUR 0
    carte = chout
    call alcart('G', carte, ma, 'NEUT_R')
    call jeveuo(carte//'.NCMP', 'E', vk8=vncmp)
    call jeveuo(carte//'.VALV', 'E', vr=valv)
    ncmp = 1
    vncmp(1) = 'X1'
    valv(1) = 0.d0
    call nocart(carte, 1, ncmp)
!
!     CREATION DE LA LISTE DES MAILLES QUI AURONT LA VALEUR 1
!     on surdimensionne a 2 fois NBMA car au pire on a NMAFOND = NBMA
!     et NBMAZO = NBMA
    listma = '&&OP0188.LISTMA'
    call wkvect(listma, 'V V I', 2*nbma, jma)
!
!     -------------------------------------------------------------
!     2) REMPLISSAGE DE LA LISTE AVEC LES MAILLES CONTENANT LE FOND
!        OU L'INTERFACE (ON PARLERA DE 'FOND'  DANS LES 2 CAS)
!     --------------------------------------------------------------
!
    if (typdis .eq. 'FISSURE') then
        mafond = fiss//'.MAILFISS.MAFOND'
    else if (typdis .eq. 'INTERFACE') then
        mafond = fiss//'.MAILFISS.HEAV'
    else
        ASSERT(.false.)
    end if
!
    call jeexin(mafond, iret)
    if (iret .eq. 0) then
        nmafon = 0
    else
        call jeveuo(mafond, 'L', jmafon)
        call jelira(mafond, 'LONMAX', nmafon)
    end if
!
    do i = 1, nmafon
        zi(jma-1+i) = zi(jmafon-1+i)
    end do
!
!     ------------------------------------------------------------------
!     3) REMPLISSAGE DE LA LISTE AVEC LES MAILLES DONT UN NOEUD EST
!        DANS LA ZONE
!     ------------------------------------------------------------------
!
!     ON CREE D'ABORD LA LISTE DES NOEUDS QUI SONT LA ZONE
    lisnoz = '&&OP0188.NOEU'
    call wkvect(lisnoz, 'V V I', nbno, jnoeu)
!
!     RECUP DES LEVEL SETS
    cnslt = '&&OP0188.CNSLT'
    cnsln = '&&OP0188.CNSLN'
    call cnocns(fiss//'.LNNO', 'V', cnsln)
    call jeveuo(cnsln//'.CNSV', 'L', vr=lsn)
    if (typdis .eq. 'FISSURE') then
        call cnocns(fiss//'.LTNO', 'V', cnslt)
        call jeveuo(cnslt//'.CNSV', 'L', vr=lst)
    end if
!
!     REMPLISSAGE DE LA LISTE DES NOEUDS QUI SONT LA ZONE
    nbnozo = 0
    do ino = 1, nbno
        if (typdis .eq. 'FISSURE') then
            dist = sqrt(lst(ino)**2+lsn(ino)**2)
        else if (typdis .eq. 'INTERFACE') then
            dist = sqrt(lsn(ino)**2)
        end if
        if (dist .le. rayon) then
            nbnozo = nbnozo+1
            zi(jnoeu-1+nbnozo) = ino
        end if
    end do
!
!     EMULATION DE DEFI_GROUP/CREA_GROUP_MA/OPTION='APPUI'
    lismaz = '&&OP0188.LISMAZ'
    call wkvect(lismaz, 'V V I', nbma, idlima)
!
!     CONNECTIVITE INVERSE
    cnxinv = '&&OP0188.CNXINV'
    call jeexin(cnxinv, iret)
    if (iret .eq. 0) then
!       ON AIMERAIT LA STOCKER DANS LA BASE GLOBALE AU CAS OU ON EN AIT
!       ENCORE BESOIN (POUR LA FISSURE SUIVANTE) MAIS ON A PAS LE DROIT
        call cncinv(ma, [ibid], 0, 'V', cnxinv)
    end if
    call jeveuo(jexatr(cnxinv, 'LONCUM'), 'L', adrvlc)
    call jeveuo(jexnum(cnxinv, 1), 'L', acncin)
!
    do i = 1, nbnozo
        nuno = zi(jnoeu+i-1)
        nbmac = zi(adrvlc+nuno+1-1)-zi(adrvlc+nuno-1)
        jadr = zi(adrvlc+nuno-1)
        do j = 1, nbmac
            numa = zi(acncin+jadr-1+j-1)
            zi(idlima+numa-1) = 1
        end do
    end do
!
!     REMPLISSAGE DE LA LISTE A LA SUITE
    nbmazo = 0
    do i = 1, nbma
        if (zi(idlima+i-1) .eq. 1) then
            nbmazo = nbmazo+1
            zi(jma-1+nmafon+nbmazo) = i
        end if
    end do
!
!     NB DE MAILLES DANS LA LISTE
    nbmali = nmafon+nbmazo
!
!     -------------------------------------------------------------
!     4) STOCKAGE DANS LA CARTE ET TRANSFORMATION EN CHAM_ELEM
!     --------------------------------------------------------------
!
!     STOCKAGE DANS LA CARTE
    vncmp(1) = 'X1'
    valv(1) = 1.d0
    call nocart(carte, 3, ncmp, mode='NUM', nma=nbmali, &
                limanu=zi(jma))
    call tecart(carte)
!
    call jedema()
end subroutine
