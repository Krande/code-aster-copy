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
subroutine cnomax(cnoz, ncmp, licmp, rmax, numno)
    implicit none
#include "jeveux.h"
#include "asterfort/cnocns.h"
#include "asterfort/cnsred.h"
#include "asterfort/detrsd.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
    character(len=*) :: cnoz
    integer(kind=8) :: ncmp
    character(len=8) :: licmp(ncmp)
    real(kind=8) :: rmax
    integer(kind=8) :: numno
!
! ======================================================================
! ROUTINE APPELEE PAR : CVGCNT
! ======================================================================
!
! CALCULER LE MAX DE LA NORME DU DEPL. (DX DY DZ) DE CNO
!
! IN  CNO    : SD CHAM_NO
! IN  NCMP   : NOMBRE DE COMPOSANTES DANS LICMP
! IN  LICMP  : COMPOSANTES SUR LESQUELLES LE MAX EST CALCULE
! OUT RMAX   : MAX DE LA NORME DU DEPL.
! OUT NUMNO  : NUMERO DU NOEUD REALISANT LE MAX DE DEPL.
!
!
!
!
    integer(kind=8) :: jcnsl
    integer(kind=8) :: nbno, k, ino
    character(len=19) :: cno, cns1, cns
    real(kind=8) :: norme
    integer(kind=8), pointer :: cnsd(:) => null()
    real(kind=8), pointer :: cnsv(:) => null()
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
    cno = cnoz
    rmax = 0.d0
    cns1 = '&&CNOMAX.CNS1'
    cns = '&&CNOMAX.CNS'
    numno = 0
    call cnocns(cno, 'V', cns1)
    call cnsred(cns1, 0, [0], ncmp, licmp, &
                'V', cns)
!
    call jeveuo(cns//'.CNSD', 'L', vi=cnsd)
    call jeveuo(cns//'.CNSV', 'L', vr=cnsv)
    call jeveuo(cns//'.CNSL', 'L', jcnsl)
!
    nbno = cnsd(1)
    do ino = 1, nbno
        norme = 0.d0
        do k = 1, ncmp
            if (zl(jcnsl-1+(ino-1)*ncmp+k)) then
                norme = norme+cnsv((ino-1)*ncmp+k)**2
            end if
        end do
        if (sqrt(norme) .ge. rmax) then
            rmax = sqrt(norme)
            numno = ino
        end if
    end do
!
!
    call detrsd('CHAM_NO_S', cns1)
    call detrsd('CHAM_NO_S', cns)
    call jedema()
end subroutine
