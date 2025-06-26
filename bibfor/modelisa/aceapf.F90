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
subroutine aceapf(nomu, noma, lmax, nbocc)
    implicit none
#include "jeveux.h"
#include "asterfort/alcart.h"
#include "asterfort/getvem.h"
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
!     AFFECTATION DES CARACTERISTIQUES POUR L'ELEMENT POUTRE_FLUIDE
! ----------------------------------------------------------------------
! IN  : NOMU   : NOM UTILISATEUR DE LA COMMANDE
! IN  : NOMA   : NOM DU MAILLAGE
! IN  : LMAX   : LONGUEUR
! IN  : NBOCC  : NOMBRE D'OCCURENCES DU MOT CLE POUFL ------------------
    real(kind=8) :: b(3), afl, ace, rapp
    character(len=19) :: cartpf
    character(len=24) :: tmpnpf, tmpvpf
!     ------------------------------------------------------------------
!
! --- CONSTRUCTION DES CARTES ET ALLOCATION
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ioc, jdcc, jdls, jdvc, nace, nafl
    integer(kind=8) :: nb1, nb2, nb3, ng, nm, nr
!-----------------------------------------------------------------------
    call jemarq()
    cartpf = nomu//'.CARPOUFL'
    tmpnpf = cartpf//'.NCMP'
    tmpvpf = cartpf//'.VALV'
    call alcart('G', cartpf, noma, 'CAPOUF_R')
    call jeveuo(tmpnpf, 'E', jdcc)
    call jeveuo(tmpvpf, 'E', jdvc)
!
    call wkvect('&&TMPPOUFL', 'V V K24', lmax, jdls)
!
    zk8(jdcc) = 'B_T'
    zk8(jdcc+1) = 'B_N'
    zk8(jdcc+2) = 'B_TN'
    zk8(jdcc+3) = 'A_FLUI'
    zk8(jdcc+4) = 'A_CELL'
    zk8(jdcc+5) = 'COEF_ECH'
!
! --- LECTURE DES VALEURS ET AFFECTATION DANS LA CARTE CARTPF
    do ioc = 1, nbocc
        b(1) = 0.d0
        b(2) = 0.d0
        b(3) = 0.d0
        afl = 0.d0
        ace = 0.d0
        rapp = 0.d0
        call getvem(noma, 'GROUP_MA', 'POUTRE_FLUI', 'GROUP_MA', ioc, &
                    lmax, zk24(jdls), ng)
        call getvr8('POUTRE_FLUI', 'B_T', iocc=ioc, scal=b(1), nbret=nb1)
        call getvr8('POUTRE_FLUI', 'B_N', iocc=ioc, scal=b(2), nbret=nb2)
        call getvr8('POUTRE_FLUI', 'B_TN', iocc=ioc, scal=b(3), nbret=nb3)
        call getvr8('POUTRE_FLUI', 'A_FLUI', iocc=ioc, scal=afl, nbret=nafl)
        call getvr8('POUTRE_FLUI', 'A_CELL', iocc=ioc, scal=ace, nbret=nace)
        call getvr8('POUTRE_FLUI', 'COEF_ECHELLE', iocc=ioc, scal=rapp, nbret=nr)
!
        if (nb2 .eq. 0) b(2) = b(1)
        zr(jdvc) = b(1)
        zr(jdvc+1) = b(2)
        zr(jdvc+2) = b(3)
        zr(jdvc+3) = afl
        zr(jdvc+4) = ace
        zr(jdvc+5) = rapp
!
! ---    "GROUP_MA" = TOUTES LES MAILLES DE LA LISTE DE GROUPES MAILLES
        if (ng .gt. 0) then
            do i = 1, ng
                call nocart(cartpf, 2, 6, groupma=zk24(jdls+i-1))
            end do
        end if
!
    end do
!
    call jedetr('&&TMPPOUFL')
    call jedetr(tmpnpf)
    call jedetr(tmpvpf)
!
    call jedema()
end subroutine
