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
function idensd(dsTypeZ, dsName1Z, dsName2Z)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/idenob.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "jeveux.h"
!
    aster_logical :: idensd
    character(len=*), intent(in) :: dsName1Z, dsName2Z, dsTypeZ
!
! --------------------------------------------------------------------------------------------------
!
!  BUT : DETERMINER L'IDENTITE DE 2 SD D'ASTER.
!  IN   TYPESD : TYPE DE  SD1 ET SD2
!       SD1   : NOM DE LA 1ERE SD
!       SD2   : NOM DE LA 2EME SD
!
!     RESULTAT:
!       IDENSD : .TRUE.    SI SD1 == SD2
!                .FALSE.   SINON
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: iden
    character(len=16) :: dsType
    character(len=19) :: pchn1, pchn2, ligrel1, ligrel2
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    dsType = dsTypeZ
    idensd = ASTER_TRUE
!
    if (dsName1Z .eq. dsName2Z) goto 999
!
    if (dsType .eq. 'NUME_EQUA') then
        pchn1 = dsName1Z
        pchn2 = dsName2Z
        iden = idenob(pchn1//'.LILI', pchn2//'.LILI')
        if (.not. iden) goto 998
        iden = idenob(pchn1//'.PRNO', pchn2//'.PRNO')
        if (.not. iden) goto 998
        iden = idenob(pchn1//'.DEEQ', pchn2//'.DEEQ')
        if (.not. iden) goto 998
        iden = idenob(pchn1//'.NUEQ', pchn2//'.NUEQ')
        if (.not. iden) goto 998
        iden = idenob(pchn1//'.NEQU', pchn2//'.NEQU')
        if (.not. iden) goto 998
        iden = idenob(pchn1//'.DELG', pchn2//'.DELG')
        if (.not. iden) goto 998
        iden = idenob(pchn1//'.REFN', pchn2//'.REFN')
        if (.not. iden) goto 998

    elseif (dsType .eq. 'LIGREL') then
        ligrel1 = dsName1Z
        ligrel2 = dsName2Z
        iden = idenob(ligrel1//'.LGRF', ligrel2//'.LGRF')
        if (.not. iden) goto 998
        iden = idenob(ligrel1//'.NBNO', ligrel2//'.NBNO')
        if (.not. iden) goto 998
        iden = idenob(ligrel1//'.PRNM', ligrel2//'.PRNM')
        if (.not. iden) goto 998
        iden = idenob(ligrel1//'.LIEL', ligrel2//'.LIEL')
        if (.not. iden) goto 998
        iden = idenob(ligrel1//'.REPE', ligrel2//'.REPE')
        if (.not. iden) goto 998
        iden = idenob(ligrel1//'.NVGE', ligrel2//'.NVGE')
        if (.not. iden) goto 998
        iden = idenob(ligrel1//'.SSSA', ligrel2//'.SSSA')
        if (.not. iden) goto 998
        iden = idenob(ligrel1//'.NEMA', ligrel2//'.NEMA')
        if (.not. iden) goto 998
        iden = idenob(ligrel1//'.PRNS', ligrel2//'.LGNS')
        if (.not. iden) goto 998
        iden = idenob(ligrel1//'.LGNS', ligrel2//'.LGRF')
        if (.not. iden) goto 998

    else
        ASSERT(ASTER_FALSE)
    end if
!
    goto 999
998 continue
    idensd = ASTER_FALSE
!
999 continue
    call jedema()
end function
