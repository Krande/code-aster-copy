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

subroutine pjngco(corres, noma1, noma2, method, cnref, &
                  base)
    implicit none
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/dismoi.h"
#include "asterfort/jacopo.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utnuav.h"
#include "asterfort/wkvect.h"
    character(len=*) :: corres, noma1, noma2, method, cnref, base
! person_in_charge: jacques.pellet at edf.fr
!----------------------------------------------------------------------
! BUT : CREER LA SD CORRESP_2_MAILLA SI METHODE='NUAGE_DEG_0/1'
!
! ARGUMENTS :
!  CORRES : IN/JXOUT  SD_CORRESP_2_MAILLA A CONSTRUIRE
!  NOMA1  : IN/JXIN   SD_MAILLAGE "1"
!  NOMA2  : IN/JXIN   SD_MAILLAGE "2"
!  METHOD : IN        METHODE CHOISIE POUR LA PROJECTION
!  CNREF  : IN/JXIN   SD_CHAM_NO "MODELE"
!  BASE   : IN   G/V  BASE DE CREATION POUR CORRES
!----------------------------------------------------------------------
    character(len=8) :: noma3
    character(len=16) :: corr16
    character(len=24) :: lno1, lno2
    character(len=4) :: num
    character(len=1) :: bas1
    integer(kind=8) :: nbocc, jngi1, jngi2, ioc, jnb12, lon1, lon2, nb1, nb2, idec1
    integer(kind=8) :: idec2, jlno1, jlno2, jxxk1
!----------------------------------------------------------------------
!
    call jemarq()
    ASSERT(method .eq. 'NUAGE_DEG_0' .or. method .eq. 'NUAGE_DEG_1')
    corr16 = corres
    bas1 = base
    ASSERT(bas1 .eq. 'G' .or. bas1 .eq. 'V')
!
    call wkvect(corr16//'.PJXX_K1', bas1//' V K24', 5, jxxk1)
    zk24(jxxk1-1+1) = noma1
    zk24(jxxk1-1+2) = noma2
    zk24(jxxk1-1+3) = method
    zk24(jxxk1-1+4) = cnref
!
    call getfac('VIS_A_VIS', nbocc)
    ASSERT(nbocc .le. 999)
!
    call dismoi('NOM_MAILLA', cnref, 'CHAMP', repk=noma3)
    ASSERT(noma3 .eq. noma2)
!
!
    if (nbocc .eq. 0) then
        call wkvect(corr16//'.PJNG_I1', bas1//' V I', 1, jngi1)
        zi(jngi1-1+1) = 0
        goto 50
!
    end if
!
!      -- ON CALCULE LES LISTES DE NOEUDS DANS DES LISTES TEMPORAIRES :
    do ioc = 1, nbocc
        call codent(ioc, 'D0', num)
        lno1 = '&&PJNGCO.'//num//'.LNO1'
        call utnuav(noma1, 1, ioc, lno1)
        lno2 = '&&PJNGCO.'//num//'.LNO2'
        call utnuav(noma2, 2, ioc, lno2)
    end do
!
!
!      -- ON RECOPIE LES LISTES DANS .PJNG_I1/2 :
    call wkvect('&&PJNGCO.NB1NB2', 'V V I', 2*nbocc, jnb12)
    lon1 = 0
    lon2 = 0
    do ioc = 1, nbocc
        call codent(ioc, 'D0', num)
        lno1 = '&&PJNGCO.'//num//'.LNO1'
        lno2 = '&&PJNGCO.'//num//'.LNO2'
        call jelira(lno1, 'LONMAX', nb1)
        call jelira(lno2, 'LONMAX', nb2)
        zi(jnb12-1+2*(ioc-1)+1) = nb1
        zi(jnb12-1+2*(ioc-1)+2) = nb2
        lon1 = lon1+nb1
        lon2 = lon2+nb2
    end do
!
    call wkvect(corr16//'.PJNG_I1', bas1//' V I', 1+nbocc+lon1, jngi1)
    zi(jngi1-1+1) = nbocc
    call wkvect(corr16//'.PJNG_I2', bas1//' V I', 1+nbocc+lon2, jngi2)
    zi(jngi2-1+1) = nbocc
    idec1 = 1+nbocc
    idec2 = 1+nbocc
    do ioc = 1, nbocc
        call codent(ioc, 'D0', num)
        lno1 = '&&PJNGCO.'//num//'.LNO1'
        lno2 = '&&PJNGCO.'//num//'.LNO2'
        nb1 = zi(jnb12-1+2*(ioc-1)+1)
        nb2 = zi(jnb12-1+2*(ioc-1)+2)
        call jeveuo(lno1, 'L', jlno1)
        call jeveuo(lno2, 'L', jlno2)
        call jacopo(nb1, 'I', jlno1, jngi1+idec1)
        call jacopo(nb2, 'I', jlno2, jngi2+idec2)
        idec1 = idec1+nb1
        idec2 = idec2+nb2
        zi(jngi1-1+1+ioc) = nb1
        zi(jngi2-1+1+ioc) = nb2
    end do
!
!
!
!        -- MENAGE :
    do ioc = 1, nbocc
        call codent(ioc, 'D0', num)
        lno1 = '&&PJNGCO.'//num//'.LNO1'
        lno2 = '&&PJNGCO.'//num//'.LNO2'
        call jedetr(lno1)
        call jedetr(lno2)
    end do
!
50  continue
    call jedema()
end subroutine
