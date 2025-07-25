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
subroutine lccata(iunit)
    implicit none
#include "asterfort/caver1.h"
#include "asterfort/lctel2.h"
#include "asterfort/lctel3.h"
#include "asterfort/lecojb.h"
    integer(kind=8) :: iunit, mxobj, iret, i
    parameter(mxobj=50)
!
    character(len=24) :: nomobj
!
!
!     LECTURE DU FICHIER CONTENANT LES OBJETS JEVEUX AU FORMAT OJB :
!     --------------------------------------------------------------
    rewind (iunit)
    do i = 1, mxobj
        call lecojb(nomobj, iunit, 'G', iret)
        if (iret .gt. 0) goto 2
        write (6, *) ' OBJET LU :', nomobj
    end do
2   continue
    write (6, *) ' NB_OBJETS LUS :', i
!
!
!     CREATION D'OBJETS SUPPLEMENTAIRES
!     ----------------------------------------------------
    call lctel2()
    call lctel3()
!
!
!     VERIFICATION DE COHERENCE DES CATALOGUES
!     ----------------------------------------------------
    call caver1()
!
!
!
end subroutine
