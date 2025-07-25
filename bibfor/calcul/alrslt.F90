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

subroutine alrslt(nout, lchout, lpaout, base)
    use calcul_module, only: ca_iachoi_, ca_iachok_, ca_iaobtr_, ca_nbobtr_, &
                             ca_ligrel_, ca_option_, ca_ldist_, ca_nuop_
    implicit none

! person_in_charge: jacques.pellet at edf.fr

#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/alchml.h"
#include "asterfort/alresl.h"
#include "asterfort/assert.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/grdeur.h"
#include "asterfort/jeexin.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/wkvect.h"

    integer(kind=8) :: nout
    character(len=*) :: base, lchout(*)
    character(len=8) :: lpaout(*)
! ----------------------------------------------------------------------
!     entrees:
!       lchout : liste des noms des champs de sortie
!       lpaout : liste des parametres associes aux champs de sortie
!       nout   : nombre de champs de sortie
!       base   : 'G', 'V'
!     sorties:
!       creation des champs globaux resultats
! ----------------------------------------------------------------------
    integer(kind=8) :: gd, descgd, code, i, iret1, iret2, iret
    character(len=19) :: nochou, dcel
    character(len=8) :: nompar
    character(len=8) :: nomgd, tsca, tych
    character(len=24), pointer :: noli(:) => null()
    character(len=24), pointer :: celk(:) => null()
! ----------------------------------------------------------------------

!   -- Allocation des champs resultats :
    do i = 1, nout
        nompar = lpaout(i)
        nochou = lchout(i)
        gd = grdeur(nompar)
        call jeveuo(jexnum('&CATA.GD.DESCRIGD', gd), 'L', descgd)
        code = zi(descgd-1+1)

!        -- si gd est 1 grandeur_simple --> cham_elem
        if (code .eq. 1) then
            call detrsd('CHAM_ELEM', nochou)
            call exisd('CHAM_ELEM_S', nochou, iret)
            if (iret .gt. 0) then
                dcel = nochou
            else
                dcel = ' '
            end if
            call alchml(ca_ligrel_, ca_option_, nompar, base, nochou, &
                        iret, dcel)
!           -- les cham_elems sont incomplets si ca_ldist_
            if (ca_ldist_) then
                call jeveuo(nochou//'.CELK', 'E', vk24=celk)
                celk(7) = 'MPI_INCOMPLET'
            end if

        else
!           -- sinon --> resuelem
            call detrsd('RESUELEM', nochou)
            ASSERT((code .ge. 3) .and. (code .le. 5))
            call alresl(ca_nuop_, ca_ligrel_, nochou, nompar, base)
!           -- les resu_elems sont incomplets en ca_ldist_
            if (ca_ldist_) then
                call jeveuo(nochou//'.NOLI', 'E', vk24=noli)
                noli(3) = 'MPI_INCOMPLET'
            end if
        end if
    end do

    call wkvect('&&CALCUL.LCHOU_I', 'V V I', max(2*nout, 2), ca_iachoi_)
    ca_nbobtr_ = ca_nbobtr_+1
    zk24(ca_iaobtr_-1+ca_nbobtr_) = '&&CALCUL.LCHOU_I'
    call wkvect('&&CALCUL.LCHOU_K8', 'V V K8', max(2*nout, 2), ca_iachok_)
    ca_nbobtr_ = ca_nbobtr_+1
    zk24(ca_iaobtr_-1+ca_nbobtr_) = '&&CALCUL.LCHOU_K8'

    do i = 1, nout
        nompar = lpaout(i)
        nochou = lchout(i)
        gd = grdeur(nompar)
        call jeveuo(jexnum('&CATA.GD.DESCRIGD', gd), 'L', descgd)
        code = zi(descgd-1+1)
        call jeexin(nochou//'.DESC', iret1)
        call jeexin(nochou//'.CELD', iret2)
        if ((iret1+iret2) .eq. 0) goto 30

        call dismoi('NOM_GD', nochou, 'CHAMP', repk=nomgd)
        call dismoi('TYPE_SCA', nomgd, 'GRANDEUR', repk=tsca)
        zk8(ca_iachok_-1+2*(i-1)+2) = tsca
        call dismoi('TYPE_CHAMP', nochou, 'CHAMP', repk=tych)
        if (tych(1:2) .eq. 'EL') then
            call jeveuo(nochou//'.CELD', 'E', zi(ca_iachoi_-1+2*(i-1)+1))
            call jeveuo(nochou//'.CELV', 'E', zi(ca_iachoi_-1+2*(i-1)+2))
            zk8(ca_iachok_-1+2*(i-1)+1) = 'CHML'
        else
            call jeveuo(nochou//'.DESC', 'E', zi(ca_iachoi_-1+2*(i-1)+1))
            zk8(ca_iachok_-1+2*(i-1)+1) = 'RESL'
        end if
30      continue
    end do

end subroutine
