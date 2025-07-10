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

subroutine matr_asse_compute_kinematic_rhs(matasz, vcinez, vrhsz)
! person_in_charge: jacques.pellet at edf.fr
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/csmbgg.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jelira.h"
#include "asterfort/mtdscr.h"

    character(len=*), intent(in) :: matasz, vcinez, vrhsz
!-----------------------------------------------------------------------
! But : Construire le second membre induit par des charges cinématiques
!       et l'ajouter au vecteur en entrée
!
!  in/jxvar  k* matasz : nom de la sd_matr_asse
!  in/jxvar  k* vcinez : nom du vecteur de charge cinématique
! out/jxvar  k* vrhsz  : nom du vecteur second membre
!
!  ATTENTION : vrhs est modifié - on lui ajoute la contribution
!              des charges cinématiques
!-----------------------------------------------------------------------

    character(len=19) :: matas, vcine, vrhs
    character(len=1) :: rouc
    integer(kind=8) :: lmat, neq0, neq1, neq2
    real(kind=8), pointer :: cine(:) => null()
    real(kind=8), pointer :: rhs(:) => null()
    complex(kind=8) :: cbid
!-------------------------------------------------------------------
    cbid = dcmplx(0.d0, 0.d0)
    call jemarq()
    matas = matasz
    vcine = vcinez
    vrhs = vrhsz

    call dismoi('NB_EQUA', matas, 'MATR_ASSE', repi=neq0)
    call dismoi('NB_EQUA', vcine, 'CHAM_NO', repi=neq1)
    call dismoi('NB_EQUA', vrhs, 'CHAM_NO', repi=neq2)
    ASSERT(neq0 .eq. neq1)
    ASSERT(neq2 .eq. neq1)

    call mtdscr(matas)
    call jeveuo(matas//'.&INT', 'L', lmat)

    call jeveuo(vcine//'.VALE', 'L', vr=cine)
    call jeveuo(vrhs//'.VALE', 'E', vr=rhs)

    call jelira(vcine//'.VALE', 'TYPE', cval=rouc)
    ASSERT(rouc .eq. 'R')
    call csmbgg(lmat, rhs, cine, [cbid], [cbid], 'R')

    call jedema()
end subroutine
