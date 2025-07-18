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

subroutine xstama(noma, nbma, nmafis, jmafis, &
                  ncouch, lisnoe, stano, cnslt, cnsln, &
                  jmafon, jmaen1, jmaen2, jmaen3, nmafon, &
                  nmaen1, nmaen2, nmaen3, typdis)
!
! person_in_charge: samuel.geniaut at edf.fr
!
    implicit none
#include "jeveux.h"
#include "asterc/r8maem.h"
#include "asterfort/dismoi.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jerazo.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/loncar.h"
#include "asterfort/xstam1.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
!
    integer(kind=8) :: nmafis, nmafon, nmaen1, nmaen2, nmaen3, nbma, jmafis
    integer(kind=8) :: ncouch, stano(*), jmafon, jmaen1, jmaen2, jmaen3
    character(len=8) :: noma
    character(len=16) :: typdis
    character(len=19) :: cnslt, cnsln
    character(len=24) :: lisnoe
!
! ----------------------------------------------------------------------
!
! ROUTINE XFEM
!
! CALCUL DU STATUT DES MAILLES
!   + PRISE EN COMPTE SI NECESSAIRE DE L'ENRICHISSEMENT A NB_COUCHES
!   -> MAJ DE STANO
!
! ----------------------------------------------------------------------
!
!
! IN  NOMA   : NOM DU MAILLAGE
! IN  NBMA   : NOMBRE DE MAILLES DU MAILLAGE
! IN  NMAFIS : NOMBRE DE MAILLES DE LA ZONE FISSURE
! IN  JMAFIS : ADRESSE DES MAILLES DE LA ZONE FISSURE
! IN  NCOUCH : NOMBRE DE COUCHES D'ENRICHISSEMENT GEOMETRIQUE
! IN  LISNOE : NOM DE LA LISTE DES NOEUDS DE GROUP_ENRI
! IN  CNSLT  : LEVEL-SET TANGENTE (TRACE DE LA FISSURE)
! IN  CNSLN  : LEVEL-SET NORMALE  (PLAN DE LA FISSURE)
!
! OUT  NMAFON : NOMBRE DE MAILLES CONTENANT LE FOND DE FISSURE
! OUT  NMAEN1 : NOMBRE DE MAILLES 'HEAVISIDE'
! OUT  NMAEN2 : NOMBRE DE MAILLES 'CRACKTIP'
! OUT  NMAEN3 : NOMBRE DE MAILLES 'HEAVISIDE-CRACKTIP'
! OUT  JMAFON : POINTEUR SUR MAILLES 'CONTENANT LE FOND DE FISSURE
! OUT  JMAEN1 : POINTEUR SUR MAILLES 'HEAVISIDE'
! OUT  JMAEN2 : POINTEUR SUR MAILLES 'CRACKTIP'
! OUT  JMAEN3 : POINTEUR SUR MAILLES 'HEAVISIDE-CRACKTIP'
!
! IN/OUT  STANO  : VECTEUR STATUT DES NOEUDS
!
!
!
!
!
    integer(kind=8) :: jma, jconx2
    integer(kind=8) :: ima, itypma, j, idim, ndim
    integer(kind=8) :: nuno, ifm, niv
    integer(kind=8) :: nbnoe, ino, nabs, jdlino, nbnoma
    real(kind=8) :: hff, diam, lsn, lst, rayon
    character(len=8) :: typma
    character(len=19) :: mai
    real(kind=8), pointer :: macoord(:) => null()
    real(kind=8), pointer :: lnsv(:) => null()
    real(kind=8), pointer :: ltsv(:) => null()
    integer(kind=8), pointer :: connex(:) => null()
    real(kind=8), pointer :: vale(:) => null()
!
!
    call jemarq()
    call infdbg('XFEM', ifm, niv)
!
    mai = noma//'.TYPMAIL'
    call jeveuo(mai, 'L', jma)
    call jeveuo(noma//'.COORDO    .VALE', 'L', vr=vale)
    call jeveuo(noma//'.CONNEX', 'L', vi=connex)
    call jeveuo(jexatr(noma//'.CONNEX', 'LONCUM'), 'L', jconx2)
!
!     1) STATUT DES MAILLES SANS TENIR COMPTE DE NB_COUCHES
!     --------------------------------------------------
    call xstam1(noma, nbma, nmafis, zi(jmafis), &
                stano, zi(jmafon), zi(jmaen1), zi(jmaen2), zi(jmaen3), &
                nmafon, nmaen1, nmaen2, nmaen3, typdis, cnslt)
!
!     S'IL N'Y A PAS DE MAILLES DE FOND, ON SORT
    if (nmafon .eq. 0) goto 999
!
!     SI NB_COUCH N'EST PAS DEFINI, ON SORT
    if (ncouch .eq. 0) goto 999
!
!     2) POUR TENIR COMPTE DE L'ENRICHISSEMENT GEOMETRIQUE A NB_COUCH
!     ------------------------------------------------------------
    if (ncouch .gt. 0) then
!
!       DIMENSINO DU MAILLAGE
        call dismoi('DIM_GEOM', noma, 'MAILLAGE', repi=ndim)
!
!       LEVEL SETS
        call jeveuo(cnslt//'.CNSV', 'L', vr=ltsv)
        call jeveuo(cnsln//'.CNSV', 'L', vr=lnsv)
!
!       CALCUL DE HFF : TAILLE MINIMALE D'UNE MAILLE DE MAFON
        hff = r8maem()
        do j = 1, nmafon
            ima = zi(jmafon-1+j)
            itypma = zi(jma-1+ima)
            call jenuno(jexnum('&CATA.TM.NOMTM', itypma), typma)
!
!         CONSTRUCTION DES COORDONNES DE LA MAILLE
            nbnoma = zi(jconx2+ima)-zi(jconx2+ima-1)
            AS_ALLOCATE(vr=macoord, size=ndim*nbnoma)
            do ino = 1, nbnoma
                nuno = connex(zi(jconx2+ima-1)+ino-1)
                do idim = 1, ndim
                    macoord(ndim*(ino-1)+idim) = vale(3*(nuno- &
                                                         1)+idim)
                end do
            end do
!
            call loncar(ndim, typma, macoord, diam)
            AS_DEALLOCATE(vr=macoord)
            hff = min(hff, diam)
        end do
!
        rayon = hff*ncouch
        write (ifm, *) 'LE RAYON D ENRICHISSEMENT EQUIVALENT EST ', rayon
!
!       ON MODIFIE L'ENRICHISSEMENT DES NOEUDS (MAJ STANO)
!       SI ANCIEN STANO = 0 -> 2
!       SI ANCIEN STANO = 1 -> 3
        call jelira(lisnoe, 'LONMAX', nbnoe)
        call jeveuo(lisnoe, 'L', jdlino)
        do ino = 1, nbnoe
            nabs = zi(jdlino-1+(ino-1)+1)
            if (abs(stano(nabs)) .le. 1) then
                lsn = lnsv((nabs-1)+1)
                lst = ltsv((nabs-1)+1)
                if (sqrt(lsn**2+lst**2) .le. rayon) then
                    stano(nabs) = stano(nabs)+2
                end if
            end if
        end do
!
        call jerazo('&&XENRCH.MAFOND', nmafis, 1)
        call jerazo('&&XENRCH.MAENR1', nbma, 1)
        call jerazo('&&XENRCH.MAENR2', nbma, 1)
        call jerazo('&&XENRCH.MAENR3', nbma, 1)
!
!       ON RECOMMENCE L'ENRICHISSEMENT DES MAILLES AVEC LE NOUVEAU STANO
        call xstam1(noma, nbma, nmafis, zi(jmafis), &
                    stano, zi(jmafon), zi(jmaen1), zi(jmaen2), zi(jmaen3), &
                    nmafon, nmaen1, nmaen2, nmaen3, typdis, cnslt)
!
    end if
!
999 continue
!
    call jedema()
end subroutine
