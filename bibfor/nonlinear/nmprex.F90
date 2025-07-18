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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine nmprex(numedd, depmoi, solalg, sddisc, numins, &
                  incest, depest)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/copisd.h"
#include "asterfort/diinst.h"
#include "asterfort/dismoi.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmchex.h"
#include "asterfort/r8inir.h"
#include "asterfort/utmess.h"
#include "blas/daxpy.h"
!
    character(len=24) :: numedd
    character(len=19) :: sddisc, incest
    character(len=19) :: solalg(*), depmoi, depest
    integer(kind=8) :: numins
!
! ----------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME - PREDICTION)
!
! PREDICTION PAR EXTRAPOLATION
!
! ----------------------------------------------------------------------
!
!
! IN  NUMEDD : NUME_DDL
! IN  NUMINS : NUMERO INSTANT COURANT
! IN  SDDISC : SD DISC_INST
! IN  DEPMOI : DEPL. EN T-
! IN  SOLALG : VARIABLE CHAPEAU POUR INCREMENTS SOLUTIONS
! OUT INCEST : INCREMENT DE DEPLACEMENT EN PREDICTION
! OUT DEPEST : DEPLACEMENT ESTIME
!
!
!
!
    integer(kind=8) :: neq
    character(len=19) :: depold
    integer(kind=8) :: ifm, niv
    real(kind=8) :: instam, instap, instaa, coef
    real(kind=8), pointer :: depes(:) => null()
    real(kind=8), pointer :: inces(:) => null()
    real(kind=8), pointer :: old(:) => null()
    blas_int :: b_incx, b_incy, b_n
!
! ----------------------------------------------------------------------
!
    call jemarq()
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<MECANONLINE> ... PAR EXTRAPOLATION'
    end if
!
! --- INITIALISATIONS
!
    call dismoi('NB_EQUA', numedd, 'NUME_DDL', repi=neq)
    instam = diinst(sddisc, numins-1)
    instap = diinst(sddisc, numins)
!
! --- DECOMPACTION DES VARIABLES CHAPEAUX
!
    call nmchex(solalg, 'SOLALG', 'DEPOLD', depold)
!
! --- INITIALISATION DU DEPLACEMENT ESTIME
!
    call copisd('CHAMP_GD', 'V', depmoi, depest)
    call jeveuo(depest(1:19)//'.VALE', 'E', vr=depes)
!
! --- INITIALISATION DE L'INCREMENT
!
    call jeveuo(incest(1:19)//'.VALE', 'E', vr=inces)
    call r8inir(neq, 0.d0, inces, 1)
!
! --- EXTRAPOLATION DES DEPLACEMENTS S'IL EXISTE UN PAS PRECEDENT
!
    if (numins .ge. 2) then
        instaa = diinst(sddisc, numins-2)
        if (instaa .eq. instam) then
            call utmess('F', 'ALGORITH8_28')
        end if
        coef = (instap-instam)/(instam-instaa)
        call jeveuo(depold(1:19)//'.VALE', 'L', vr=old)
        b_n = to_blas_int(neq)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, coef, old, b_incx, depes, &
                   b_incy)
        b_n = to_blas_int(neq)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, coef, old, b_incx, inces, &
                   b_incy)
    end if
!
    call jedema()
end subroutine
