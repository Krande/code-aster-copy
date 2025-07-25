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

subroutine temess(typ)
    use calcul_module, only: ca_option_, ca_nomte_, ca_icaelk_, ca_ialiel_, &
                             ca_illiel_, ca_igr_, ca_iel_, ca_nomtm_, ca_iamaco_, ca_ilmaco_

!

! But : Cette routine complete le message d'erreur emis par utmess
!       pendant un calcul elementaire :
!         * nom de l'option calculee
!         * nom de la maille
!         * coordonnees des noeuds
!         * nom du maillage
!         * noms de GROUP_MA contenant la maille
!         ...
!-----------------------------------------------------------------------------------
    implicit none

#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/utmess_core.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/jeexin.h"
#include "asterfort/jenuno.h"
#include "asterfort/jemarq.h"
#include "asterfort/jedema.h"

    character(len=1), intent(in) :: typ
    integer(kind=8) :: ima, iexi, k, n1, k1, jgrma, nbgrma, nbgrmt
    integer(kind=8) :: nno, ino, nuno, jcoor
    character(len=24) :: valkc(9), ligrma(4), nomgrm, grpmav
    character(len=8) :: ma
    character(len=256) :: ufname
    real(kind=8) :: valrc(3)
!-----------------------------------------------------------------------------------
    call jemarq()

    ma = zk24(ca_icaelk_-1+1) (1:8)

    nbgrma = 0
    nno = 0
    ufname = ' '
    ima = zi(ca_ialiel_-1+zi(ca_illiel_+ca_igr_-1)+ca_iel_-1)
    if (ima .gt. 0) then
        nno = zi(ca_ilmaco_-1+ima+1)-zi(ca_ilmaco_-1+ima)
    end if

!   -- recherche des 4 premiers GROUP_MA qui contiennent ima :
    if (ima .gt. 0) then
        grpmav = ma//'.GROUPEMA'
        call jeexin(grpmav, iexi)
        if (iexi .gt. 0) then
            call jelira(grpmav, 'NMAXOC', nbgrmt)
            do k = 1, nbgrmt
                call jeexin(jexnum(grpmav, k), iexi)
                if (iexi .eq. 0) cycle
                call jenuno(jexnum(grpmav, k), nomgrm)
                call jelira(jexnum(grpmav, k), 'LONUTI', n1)
                call jeveuo(jexnum(grpmav, k), 'L', jgrma)
                do k1 = 1, n1
                    if (zi(jgrma-1+k1) .eq. ima) then
                        nbgrma = nbgrma+1
                        ligrma(nbgrma) = nomgrm
                        if (nbgrma .eq. 4) goto 100
                    end if
                end do
            end do
        end if
100     continue

!       -- calcul du centre de gravite de la maille :
        valrc(:) = 0.d0
        call jeveuo(ma//'.COORDO    .VALE', 'L', jcoor)
        do ino = 1, nno
            nuno = zi(ca_iamaco_-1+zi(ca_ilmaco_+ima-1)+ino-1)
            valrc(1) = valrc(1)+zr(jcoor-1+3*(nuno-1)+1)/nno
            valrc(2) = valrc(2)+zr(jcoor-1+3*(nuno-1)+2)/nno
            valrc(3) = valrc(3)+zr(jcoor-1+3*(nuno-1)+3)/nno
        end do
    end if

    valkc(:) = ' '
    valkc(1) = ca_option_
    valkc(2) = ca_nomte_
    valkc(3) = ma
    valkc(4) = ca_nomtm_
    do k = 1, nbgrma
        valkc(4+k) = ligrma(k)
    end do
    if (nbgrma .eq. 4) valkc(4+4) = '...'

    call utmess_core(typ, 'CALCUL_49', 9, valkc, 1, &
                     [ima], 3, valrc, 0, ufname)

    call jedema()
end subroutine temess
