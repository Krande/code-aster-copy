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
subroutine aceaca(nomu, noma, lmax, nbocc)
    implicit none
#include "jeveux.h"
#include "asterfort/alcart.h"
#include "asterfort/getvem.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nocart.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: lmax, nbocc
    character(len=8) :: nomu, noma
!     AFFE_CARA_ELEM
!     AFFECTATION DES CARACTERISTIQUES POUR L'ELEMENT CABLE
! ----------------------------------------------------------------------
! IN  : NOMU   : NOM UTILISATEUR DE LA COMMANDE
! IN  : NOMA   : NOM DU MAILLAGE
! IN  : LMAX   : NOMBRE MAX DE MAILLES OU GROUPE DE MAILLES
! IN  : NBOCC  : NOMBRE D'OCCURENCE DE NOT CLE CABLE
! ----------------------------------------------------------------------
    real(kind=8) :: sct, tens
    character(len=8) :: fcx
    character(len=19) :: cartca, cartcf
    character(len=24) :: tmpnca, tmpvca, tmpncf, tmpvcf
!     ------------------------------------------------------------------
!
! --- CONSTRUCTION DES CARTES ET ALLOCATION
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ioc, jdcc, jdccf, jdls, jdvc, jdvcf
    integer(kind=8) :: nfcx, ng, nt, nv
!-----------------------------------------------------------------------
    call jemarq()
    cartca = nomu//'.CARCABLE'
    tmpnca = cartca//'.NCMP'
    tmpvca = cartca//'.VALV'
    call alcart('G', cartca, noma, 'CACABL_R')
    call jeveuo(tmpnca, 'E', jdcc)
    call jeveuo(tmpvca, 'E', jdvc)
!
    cartcf = nomu//'.CVENTCXF'
    tmpncf = cartcf//'.NCMP'
    tmpvcf = cartcf//'.VALV'
    call jeveuo(tmpncf, 'E', jdccf)
    call jeveuo(tmpvcf, 'E', jdvcf)
!
    call wkvect('&&TMPCABLE', 'V V K24', lmax, jdls)
!
    zk8(jdcc) = 'SECT'
    zk8(jdcc+1) = 'TENS'
!
    zk8(jdccf) = 'FCXP'
!
!
! --- LECTURE DES VALEURS ET AFFECTATION DANS LA CARTE CARTCA
    do ioc = 1, nbocc
        sct = 0.d0
        call getvem(noma, 'GROUP_MA', 'CABLE', 'GROUP_MA', ioc, &
                    lmax, zk24(jdls), ng)
!
        call getvr8('CABLE', 'SECTION', iocc=ioc, scal=sct, nbret=nv)
        if (nv .eq. 0) then
            call getvr8('CABLE', 'A', iocc=ioc, scal=sct, nbret=nv)
        end if
        zr(jdvc) = sct
        call getvr8('CABLE', 'N_INIT', iocc=ioc, scal=tens, nbret=nt)
        zr(jdvc+1) = tens
!
        fcx = '.'
        call getvid('CABLE', 'FCX', iocc=ioc, scal=fcx, nbret=nfcx)
        zk8(jdvcf) = fcx
!
! ---    "GROUP_MA" = TOUTES LES MAILLES DE LA LISTE DE GROUPES MAILLES
        if (ng .gt. 0) then
            do i = 1, ng
                call nocart(cartca, 2, 2, groupma=zk24(jdls+i-1))
                call nocart(cartcf, 2, 1, groupma=zk24(jdls+i-1))
            end do
        end if
!
    end do
!
    call jedetr('&&TMPCABLE')
    call jedetr(tmpnca)
    call jedetr(tmpvca)
!
    call jedema()
end subroutine
