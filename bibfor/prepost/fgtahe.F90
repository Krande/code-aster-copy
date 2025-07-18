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

subroutine fgtahe(kdomm, nbcycl, epsmin, epsmax, dom)
    implicit none
#include "jeveux.h"
#include "asterfort/fgtaep.h"
#include "asterfort/fgtaes.h"
#include "asterfort/getvid.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/rccome.h"
#include "asterfort/rcpare.h"
#include "asterfort/utmess.h"
    character(len=*) :: kdomm
    real(kind=8) :: epsmin(*), epsmax(*)
    real(kind=8) :: dom(*)
    integer(kind=8) :: nbcycl
!     ROUTINE CHAPEAU POUR LE CALCUL DU DOMMAGE PAR LOIS DE TAHERI
!     ------------------------------------------------------------------
! IN  KDOMM  : K   : LOI DE DOMMAGE TAHERI_MANSON/TAHERI_MIXTE
! IN  NBCYCL : I   : NOMBRE DE CYCLES
! IN  EPSMIN : R   : DEFORMATIONS MINIMALES DES CYCLES
! IN  EPSMAX : R   : DEFORMATIONS MAXIMALES DES CYCLES
! OUT DOM    : R   : VALEURS DES DOMMAGES ELEMENTAIRES
!     ------------------------------------------------------------------
!
!
    integer(kind=8) :: icodwo, icodma, icodba, icodhs, icodre(3)
    character(len=8) :: nommat, nomfo1, nomnap
    character(len=16) :: cara
    character(len=32) :: pheno
!
!-----------------------------------------------------------------------
    integer(kind=8) :: nbval
!-----------------------------------------------------------------------
    call jemarq()
!
! --- CALCUL DU DOMMAGE ELEMENTAIRE DE TAHERI_MANSON_COFFIN
!
    if (kdomm(1:13) .eq. 'TAHERI_MANSON') then
        call getvid(' ', 'MATER', nbval=0, nbret=nbval)
        if (nbval .eq. 0) then
            call utmess('F', 'FATIGUE1_8')
        end if
        call getvid(' ', 'TAHERI_FONC', nbval=0, nbret=nbval)
        if (nbval .eq. 0) then
            call utmess('F', 'FATIGUE1_9')
        end if
        call getvid(' ', 'TAHERI_NAPPE', nbval=0, nbret=nbval)
        if (nbval .eq. 0) then
            call utmess('F', 'FATIGUE1_10')
        end if
        call getvid(' ', 'MATER', scal=nommat, nbret=nbval)
        pheno = 'FATIGUE'
        call rccome(nommat, pheno, icodre(1))
        if (icodre(1) .eq. 1) then
            call utmess('F', 'FATIGUE1_24')
        end if
        cara = 'MANSON_COFFIN'
        call rcpare(nommat, pheno, cara, icodma)
        if (icodma .ne. 0) then
            call utmess('F', 'FATIGUE1_11')
        end if
        call getvid(' ', 'TAHERI_FONC', scal=nomfo1, nbret=nbval)
        call getvid(' ', 'TAHERI_NAPPE', scal=nomnap, nbret=nbval)
        call fgtaep(nommat, nomfo1, nomnap, nbcycl, epsmin, &
                    epsmax, dom)
    end if
!
! --- CALCUL DU DOMMAGE ELEMENTAIRE DE TAHERI_MIXTE
!
    if (kdomm(1:14) .eq. 'TAHERI_MIXTE') then
        call getvid(' ', 'MATER', nbval=0, nbret=nbval)
        if (nbval .eq. 0) then
            call utmess('F', 'FATIGUE1_12')
        end if
        call getvid(' ', 'TAHERI_NAPPE', nbval=0, nbret=nbval)
        if (nbval .eq. 0) then
            call utmess('F', 'FATIGUE1_10')
        end if
        call getvid(' ', 'MATER', scal=nommat, nbret=nbval)
        pheno = 'FATIGUE'
        call rccome(nommat, pheno, icodre(1))
        if (icodre(1) .eq. 1) then
            call utmess('F', 'FATIGUE1_24')
        end if
        cara = 'MANSON_COFFIN'
        call rcpare(nommat, pheno, cara, icodma)
        cara = 'WOHLER'
        call rcpare(nommat, pheno, cara, icodwo)
        cara = 'A_BASQUIN'
        call rcpare(nommat, pheno, cara, icodba)
        cara = 'A0'
        call rcpare(nommat, pheno, cara, icodhs)
        if (icodma .ne. 0) then
            call utmess('F', 'FATIGUE1_13')
        end if
        if (icodwo .ne. 0 .and. icodba .ne. 0 .and. icodhs .ne. 0) then
            call utmess('F', 'FATIGUE1_14')
        end if
        call getvid(' ', 'TAHERI_NAPPE', scal=nomnap, nbret=nbval)
        call fgtaes(nommat, nomnap, nbcycl, epsmin, epsmax, &
                    dom)
    end if
!
    call jedema()
end subroutine
