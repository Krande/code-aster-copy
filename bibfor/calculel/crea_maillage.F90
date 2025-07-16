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

subroutine crea_maillage(noma, noma2, base, nbno, lino)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cargeo.h"
#include "asterfort/dismoi.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecreo.h"
#include "asterfort/jedema.h"
#include "asterfort/jedupo.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/wkvect.h"
!
    character(len=8) :: noma, noma2
    character(len=1) :: base
    integer(kind=8) :: nbno, lino(*)

! person_in_charge: jacques.pellet at edf.fr

! ======================================================================
!  But: creer un "petit" maillage contenant une liste de noeuds

!  noma   : in  : maillage de depart
!  noma2  : out : maillage a creer
!  base   : in  : 'G' ou 'V' : base pour la creation de noma2
!  nbno   : in  : nombre de noeuds de la liste lino
!  lima   : in  : liste des numeros des noeuds
! ======================================================================

    integer(kind=8) ::  nbnoin, ino, jdim, jcorou, iad, ntgeo, nbnoou
    integer(kind=8) ::  ino2, typpoi, jadou, itypou, k
    character(len=4) :: docu
    character(len=24) ::  cooval, coodsc
    character(len=24) ::  dimin, dimou, connex, typmai
    real(kind=8), pointer :: vale(:) => null()
!------------------------------------------------------------------------------

    call jemarq()

    ASSERT(noma .ne. noma2)
    ASSERT(base .eq. 'V' .or. base .eq. 'G')

! -1- PRELIMINAIRES
!     ============
    call dismoi('NB_NO_MAILLA', noma, 'MAILLAGE', repi=nbnoin)
    nbnoou = nbno

! -2- CREATION DU NOUVEAU MAILLAGE
!     ============================
    cooval = noma2//'.COORDO    .VALE'
    coodsc = noma2//'.COORDO    .DESC'

! --- OBJET .DIME
    dimin = noma//'.DIME'
    dimou = noma2//'.DIME'
    call jedupo(dimin, base, dimou, .false._1)
    call jeveuo(dimou, 'E', jdim)
    zi(jdim-1+1) = nbnoou
    zi(jdim-1+3) = 0

! --- OBJET .COORDO.VALE
    call wkvect(cooval, base//' V R', nbnoou*3, jcorou)
    call jelira(noma//'.COORDO    .VALE', 'DOCU', cval=docu)
    call jeecra(cooval, 'DOCU', cval=docu)
    call jeveuo(noma//'.COORDO    .VALE', 'L', vr=vale)
    do ino2 = 1, nbnoou
        ino = lino(ino2)
        zr(jcorou+3*(ino2-1)) = vale(1+3*(ino-1))
        zr(jcorou+3*(ino2-1)+1) = vale(1+3*(ino-1)+1)
        zr(jcorou+3*(ino2-1)+2) = vale(1+3*(ino-1)+2)
    end do

! --- OBJET COORDO.DESC
    call jecreo(coodsc, base//' V I')
    call jeecra(coodsc, 'LONMAX', 3)
    call jeecra(coodsc, 'DOCU', cval='CHGO')
    call jeveuo(coodsc, 'E', iad)
    call jenonu(jexnom('&CATA.GD.NOMGD', 'GEOM_R'), ntgeo)
    zi(iad) = ntgeo
    zi(iad+1) = -3
    zi(iad+2) = 14

! --- Pour qu'on puisse voire le maillage de noeuds avec salome,
!     on ajoute des POI1 sur tous les noeuds
!----------------------------------------------------
    connex = noma2//'.CONNEX'
    typmai = noma2//'.TYPMAIL'

! --- OBJET .TYPMAIL
    call wkvect(typmai, base//' V I', nbnoou, itypou)
    call jenonu(jexnom('&CATA.TM.NOMTM', 'POI1'), typpoi)
    do k = 1, nbnoou
        zi(itypou-1+k) = typpoi
    end do

! --- OBJET .CONNEX
    call jecrec(connex, base//' V I', 'NU', 'CONTIG', 'VARIABLE', nbnoou)
    call jeecra(connex, 'LONT', nbnoou, ' ')
    do k = 1, nbnoou
        call jeecra(jexnum(connex, k), 'LONMAX', 1)
        call jeveuo(jexnum(connex, k), 'E', jadou)
        zi(jadou+1-1) = k
    end do

    call cargeo(noma2)

    call jedema()

end subroutine
