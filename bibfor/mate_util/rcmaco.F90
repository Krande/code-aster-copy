! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
! person_in_charge: j-pierre.lefebvre at edf.fr
!
subroutine rcmaco(chmat, chmatgrp, indmat, nbmat, imate, l_ther, basename, base_)
!
implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jeveut.h"
#include "asterfort/jexnum.h"
#include "asterfort/matcod.h"
#include "asterfort/utmess.h"
!
character(len=8) :: chmat, basename
character(len=24) :: chmatgrp
integer :: indmat, nbmat, imate
aster_logical, intent(in) :: l_ther
character(len=1), intent(in), optional :: base_
!
! ----------------------------------------------------------------------
!
!     BUT: CREER L'OBJET BASENAME//'      .CODI' ,LE REMPLIR ET RENVOYER
!          SON ADRESSE PAR RAPPORT A ZI
!
! ----------------------------------------------------------------------
!
    integer :: nbcmp,  igrp
    character(len=8) :: nomgd
    character(len=19) :: codi
    integer, pointer :: desc(:) => null()
    character(len=1) :: base
!
! ----------------------------------------------------------------------
!
    call jemarq()
    if( present(base_) ) then
        base = base_
    else
        base = 'V'
    endif
!
    call jeveut(chmatgrp, 'L', igrp)
    call jeveuo(chmat(1:8)//'.CHAMP_MAT .DESC', 'L', vi=desc)
    call jenuno(jexnum('&CATA.GD.NOMCMP', desc(1)), nomgd)
    call dismoi('NB_CMP_MAX', nomgd, 'GRANDEUR', repi=nbcmp)
    if (imate .gt. 9999) then
        call utmess('F', 'CALCULEL6_11')
    endif
!
    call matcod(chmat, indmat, nbmat, imate, igrp,&
                    basename, codi, l_ther, base)
!
    call jedema()
!
end subroutine
