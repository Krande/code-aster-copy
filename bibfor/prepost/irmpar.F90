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
subroutine irmpar(nomcon, ifichi, paraListNb, paraListName)
!
    use as_med_module, only: as_med_open
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/as_mficom.h"
#include "asterfort/as_mficlo.h"
#include "asterfort/as_mprcre.h"
#include "asterfort/codent.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/ulisog.h"
#include "asterfort/utmess.h"
!
    character(len=*) :: nomcon
    integer(kind=8) :: ifichi
    integer(kind=8), intent(in) :: paraListNb
    character(len=16), pointer :: paraListName(:)
!
! --------------------------------------------------------------------------------------------------
!
!     CREATION D'UN PARAMETRE DANS UN FICHIER MED
!
! --------------------------------------------------------------------------------------------------
!
! IN  NOMCON : K8  : NOM DU CONCEPT A IMPRIMER
! IN  IFICHI : IS  : UNITE LOGIQUE D'ECRITURE
! In  paraListNb       : length of list of parameter names
! Ptr paraListName     : pointer to the list of parameter names
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: edleaj, codret, hdfok, medok
    med_idt :: idfimd
    integer(kind=8) :: iPara
    integer(kind=8) :: typflo
    parameter(typflo=6)
    character(len=1) :: saux01
    character(len=8) :: saux08
    character(len=16) :: saux16, paraName
    character(len=64) :: saux64
    character(len=200) :: nofimd, nopar2
    character(len=255) :: kfic
    aster_logical :: ficexi
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    call ulisog(ifichi, kfic, saux01)
    if (kfic(1:1) .eq. ' ') then
        call codent(ifichi, 'G', saux08)
        nofimd = 'fort.'//saux08
    else
        nofimd = kfic(1:200)
    end if
    inquire (file=nofimd, exist=ficexi)
    if (ficexi) then
        call as_mficom(nofimd, hdfok, medok, codret)
        if (medok .eq. 0 .or. hdfok .eq. 0 .or. codret .ne. 0) then
            edleaj = 3
            call as_med_open(idfimd, nofimd, edleaj, codret)
        else
            edleaj = 1
            call as_med_open(idfimd, nofimd, edleaj, codret)
        end if
    else
        edleaj = 3
        call as_med_open(idfimd, nofimd, edleaj, codret)
    end if
    if (codret .ne. 0) then
        saux08 = 'mfiope'
        call utmess('F', 'DVP_97', sk=saux08, si=codret)
    end if
!
    do iPara = 1, paraListNb
        paraName = paraListName(iPara)
        saux64 = nomcon//paraName
        nopar2 = paraName
        saux16 = ' '
!       TROISIEME ARGUMENT : 6 TYPE FLOTTANT DANS MED
        call as_mprcre(idfimd, saux64, typflo, nopar2, saux16, codret)
        if (codret .ne. 0) then
            saux08 = 'mprcre'
            call utmess('F', 'DVP_97', sk=saux08, si=codret)
        end if
    end do
!
    call as_mficlo(idfimd, codret)
    if (codret .ne. 0) then
        saux08 = 'mficlo'
        call utmess('F', 'DVP_97', sk=saux08, si=codret)
    end if
!
    call jedema()
end subroutine
