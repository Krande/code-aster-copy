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
subroutine orilma(noma, ndim, listCellNume, nbCell, norien, &
                  ntrait, reorie, nbmavo, mailvo)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/oriema.h"
#include "asterfort/utmasu.h"
#include "asterfort/utmess.h"
#include "asterfort/int_to_char8.h"
!
    integer(kind=8) :: ndim, nbCell, norien, ntrait, nbmavo, mailvo(*)
    integer(kind=8), pointer :: listCellNume(:)
    character(len=8) :: noma
    aster_logical :: reorie
!.======================================================================
!
!     ORILMA  --  LE BUT EST DE REORIENTER, SI C'EST NECESSAIRE,
!                 LES MAILLES DE PEAU D'UNE LISTE DE MAILLES
!                 LA NORMALE A LA MAILLE DE PEAU DOIT ETRE
!                 EXTERIEURE AU VOLUME.
!                 DANS LE CAS OU REORIE EST FAUX, L'ORIENTATION
!                 GEOMETRIQUE N'EST PAS UTILISEE, CECI PERMET DE
!                 TESTER UNE SURFACE POUR UNE CONDITION AUX LIMITES
!                 DE PRESSION
!
!   ARGUMENT        E/S  TYPE         ROLE
!    NOMA           IN    K8      NOM DU MAILLAGE
!    NDIM           IN    I       DIMENSION DU PROBLEME
!    LISTMA         IN    I       LISTE DES NUMEROS DE MAILLE
!                                 A REORIENTER
!    NBMAIL         IN    I       NB DE MAILLES DE LA LISTE
!    NORIEN        VAR            NOMBRE DE MAILLES REORIENTEES
!    REORIE         IN    L       INDIQUE SI L'ON DOIT APPELER ORIEMA
!    MAILVO         IN    I       SI ORIE_PEAU ("GROUP_MA_INTERNE")
!                                   = LISTE DES MAILLES INTERNES
!                                     UTILES A LA REORIENTATION
!                                 SINON: MAILVO N'EST PAS UTILISE
!    NBMAVO         IN    I       NB DE MAILLES DE MAILVO
!.========================= DEBUT DES DECLARATIONS ====================
! -----  VARIABLES LOCALES
    integer(kind=8) :: ifm, niv, iCell, cellNume, cellTypeNume, nbnmai, numa3d, noriem, norieg
    integer(kind=8) :: p1, p2, jm3d, jdesm, jdes3d, ier
    aster_logical :: hasSkin1D, hasSkin2D, lcolle
    character(len=2) :: kdim
    character(len=8) :: cellTypeName, nomail, typ3d
    character(len=24) :: nomob1
    character(len=24) :: valk(2)
    integer(kind=8), pointer :: ori1(:) => null()
    integer(kind=8), pointer :: ori2(:) => null()
    character(len=8), pointer :: ori3(:) => null()
    character(len=8), pointer :: ori4(:) => null()
    integer(kind=8), pointer :: typmail(:) => null()
    real(kind=8), pointer :: vale(:) => null()
!
!.========================= DEBUT DU CODE EXECUTABLE ==================
!
    call jemarq()
    if (nbCell .eq. 0) goto 999
!
! --- INITIALISATIONS :
!     ---------------
    call infniv(ifm, niv)
!
! --- VECTEUR DU TYPE DES MAILLES DU MAILLAGE :
!     ---------------------------------------
    call jeveuo(noma//'.TYPMAIL        ', 'L', vi=typmail)
!
! --- COORDONNEES DES NOEUDS DU MAILLAGE :
!     ----------------------------------
    call jeveuo(noma//'.COORDO    .VALE', 'L', vr=vale)
!
! --- RECUPERATION DE LA CONNECTIVITE DES MAILLES :
!     -------------------------------------------
    call jeveuo(jexatr(noma//'.CONNEX', 'LONCUM'), 'L', p2)
    call jeveuo(noma//'.CONNEX', 'E', p1)
!
!     ALLOCATIONS :
!     -----------
    AS_ALLOCATE(vi=ori1, size=nbCell)
    AS_ALLOCATE(vi=ori2, size=nbCell)
    AS_ALLOCATE(vk8=ori3, size=nbCell)
    AS_ALLOCATE(vk8=ori4, size=nbCell)
!
! --- VERIFICATION DU TYPE DES MAILLES
! --- (ON DOIT AVOIR DES MAILLES DE PEAU) :
!     -----------------------------------
    hasSkin1D = .false.
    hasSkin2D = .false.
    lcolle = .false.
    call jeexin(noma//'.NOMMAI', ier)
    if (ier .ne. 0) then
        lcolle = .true.
    end if
    do iCell = 1, nbCell
        cellNume = listCellNume(iCell)
        nomail = int_to_char8(cellNume, lcolle, noma, 'MAILLE')
        ori3(iCell) = nomail
        ori1(iCell) = zi(p2+cellNume+1-1)-zi(p2+cellNume-1)
        ori2(iCell) = zi(p2+cellNume-1)
!
! ---   TYPE DE LA MAILLE COURANTE :
!       --------------------------
        cellTypeNume = typmail(cellNume)
        call jenuno(jexnum('&CATA.TM.NOMTM', cellTypeNume), cellTypeName)
        ori4(iCell) = cellTypeName
!
        if (cellTypeName(1:4) .eq. 'QUAD') then
            hasSkin2D = .true.
        else if (cellTypeName(1:4) .eq. 'TRIA') then
            hasSkin2D = .true.
        else if (cellTypeName(1:3) .eq. 'SEG') then
            hasSkin1D = .true.
        else
            valk(1) = nomail
            valk(2) = cellTypeName
            call utmess('F', 'MODELISA5_94', nk=2, valk=valk)
        end if
        if (hasSkin1D .and. hasSkin2D) then
            call utmess('F', 'MODELISA5_98')
        end if
!
    end do
!
! --- RECHERCHE DES MAILLES SUPPORTS
!
    kdim = ' '
    if (hasSkin1D) kdim = '2D'
    if (hasSkin2D) kdim = '3D'
    ASSERT(kdim .ne. ' ')
    nomob1 = '&&ORILMA.MAILLE_3D'
    call utmasu(noma, kdim, nbCell, listCellNume, nomob1, &
                vale, nbmavo, mailvo, .false._1)
    call jeveuo(nomob1, 'L', jm3d)
!
    norieg = 0
    ntrait = 0
!
    do iCell = 1, nbCell
!
        nomail = ori3(iCell)
        cellTypeName = ori4(iCell)
        nbnmai = ori1(iCell)
        jdesm = ori2(iCell)
        numa3d = zi(jm3d-1+iCell)
        if (numa3d .eq. 0) then
            ntrait = ntrait+1
            goto 100
        end if
        jdes3d = zi(p2+numa3d-1)
        call jenuno(jexnum('&CATA.TM.NOMTM', typmail(numa3d)), typ3d)
!
        call oriema(nomail, cellTypeName, nbnmai, zi(p1+jdesm-1), typ3d, &
                    zi(p1+jdes3d-1), ndim, vale, reorie, noriem, &
                    ifm, niv)
!
        norieg = norieg+noriem
!
100     continue
    end do
!
    norien = norien+norieg
!
    AS_DEALLOCATE(vi=ori1)
    AS_DEALLOCATE(vi=ori2)
    AS_DEALLOCATE(vk8=ori3)
    AS_DEALLOCATE(vk8=ori4)
    call jedetr(nomob1)
!
999 continue
    call jedema()
end subroutine
