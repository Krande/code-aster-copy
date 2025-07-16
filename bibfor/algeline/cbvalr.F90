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
subroutine cbvalr(rouc, neq, smhc, smdi, idlexc, &
                  coefr, coefc, valmi, valmr, valmc)
    implicit none
!     BUT : ACCUMULTATION DANS VALMR (OU VALMC) DE COEF*VALMI
!     ROUC=
!        /'RR' : ON UTILISE VALMR ET COEFR
!        /'RC' : ON UTILISE VALMR ET COEFC
!        /'CR' : ON UTILISE VALMC ET COEFR
!        /'CC' : ON UTILISE VALMC ET COEFC
!-------------------------------------------------------------------
#include "asterfort/assert.h"
    character(len=2) :: rouc
    integer(kind=4) :: smhc(*)
    integer(kind=8) :: neq, smdi(*), idlexc(*)
    integer(kind=8) :: kin, idebli, ilig, ifinli, ind, jcol
    real(kind=8) :: coefr, valmi(*), valmr(*)
    complex(kind=8) :: coefc, valmc(*)
!     ------------------------------------------------------------------
    kin = 0
    idebli = 1
!
!
    if (rouc .eq. 'RR') then
!     -------------------------------
        do ilig = 1, neq
            ifinli = smdi(ilig)
            do ind = idebli, ifinli
                kin = kin+1
                jcol = smhc(ind)
                valmr(kin) = valmr(kin)+coefr*valmi(kin)*(1-idlexc(jcol))*(1-idlexc(ilig))
            end do
            idebli = smdi(ilig)+1
        end do
!
!
    else if (rouc .eq. 'RC') then
!     -------------------------------
        do ilig = 1, neq
            ifinli = smdi(ilig)
            do ind = idebli, ifinli
                kin = kin+1
                jcol = smhc(ind)
                valmr(kin) = valmr(kin)+dble(coefc*valmi(kin)*(1-idlexc(jcol))*(1-idlexc(ilig&
                             &)))
            end do
            idebli = smdi(ilig)+1
        end do
!
!
    else if (rouc .eq. 'CR') then
!     -------------------------------
        do ilig = 1, neq
            ifinli = smdi(ilig)
            do ind = idebli, ifinli
                kin = kin+1
                jcol = smhc(ind)
                valmc(kin) = valmc(kin)+coefr*valmi(kin)*(1-idlexc(jcol))*(1-idlexc(ilig))
            end do
            idebli = smdi(ilig)+1
        end do
!
!
    else if (rouc .eq. 'CC') then
!     -------------------------------
        do ilig = 1, neq
            ifinli = smdi(ilig)
            do ind = idebli, ifinli
                kin = kin+1
                jcol = smhc(ind)
                valmc(kin) = valmc(kin)+coefc*valmi(kin)*(1-idlexc(jcol))*(1-idlexc(ilig))
            end do
            idebli = smdi(ilig)+1
        end do
!
!
    else
        ASSERT(.false.)
    end if
!
end subroutine
