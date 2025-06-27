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
subroutine prelog(ndim, lgpg, vim, gn, lamb, &
                  logl, fPrev, fCurr, epslPrev, epslIncr, &
                  tlogPrev, lCorr, iret)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/deflog.h"
#include "asterfort/lcdetf.h"
#include "asterc/r8prem.h"
!
    integer(kind=8), intent(in) :: ndim, lgpg
    real(kind=8), intent(in) :: vim(lgpg)
    real(kind=8), intent(in) :: fPrev(3, 3), fCurr(3, 3)
    real(kind=8), intent(out) :: epslPrev(6), epslIncr(6)
    real(kind=8), intent(out) :: tlogPrev(6)
    real(kind=8), intent(out) :: gn(3, 3), lamb(3), logl(3)
    aster_logical, intent(in) :: lCorr
    integer(kind=8), intent(out) :: iret
!
! --------------------------------------------------------------------------------------------------
!
!  BUT:  CALCUL DES GRANDES DEFORMATIONS  LOG 2D (D_PLAN ET AXI) ET 3D
!     SUIVANT ARTICLE MIEHE APEL LAMBRECHT CMAME 2002
!
! --------------------------------------------------------------------------------------------------
!
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  LGPG    : DIMENSION DU VECTEUR DES VAR. INTERNES POUR 1 PT GAUSS
! IN  VIM     : VARIABLES INTERNES EN T-
! OUT GN      : TERMES UTILES AU CALCUL DE TL DANS POSLOG
! OUT LAMB    : TERMES UTILES AU CALCUL DE TL DANS POSLOG
! OUT LOGL    : TERMES UTILES AU CALCUL DE TL DANS POSLOG
! IN FM       : GRADIENT TRANSFORMATION EN T-
! IN FP       : GRADIENT TRANSFORMATION EN T+
! OUT EPSML   : DEFORAMTIONS LOGARITHMIQUES EN T-
! OUT DEPS    : ACCROISSEEMENT DE DEFORMATIONS LOGARITHMIQUES
! OUT TN      : CONTRAINTES ASSOCIEES AUX DEF. LOGARITHMIQUES EN T-
! OUT IRET    : 0=OK, 1=vp(Ft.F) trop petites (compression infinie)
! IN  RESI    : .TRUE. SI FULL_MECA/RAPH_MECA .FALSE. SI RIGI_MECA_TANG
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: epslCurr(6), detf
    real(kind=8), dimension(6), parameter  :: vrac2 = (/1.d0, 1.d0, 1.d0 &
                                                        , sqrt(2.d0), sqrt(2.d0), sqrt(2.d0)/)
!
! --------------------------------------------------------------------------------------------------
!
    epslPrev = 0.d0
    epslIncr = 0.d0
    tlogPrev = 0.d0
    gn = 0.d0
    lamb = 0.d0
    logl = 0.d0
    iret = 0

! - Compute kinematic at begin of step
    call lcdetf(ndim, fPrev, detf)
    if (detf .le. r8prem()) then
        iret = 1
        goto 999
    end if
    call deflog(ndim, fPrev, epslPrev, gn, lamb, logl, iret)
    if (iret .ne. 0) then
        goto 999
    end if

! - Compute kinematic at end of step
    if (lCorr) then
        call lcdetf(ndim, fCurr, detf)
        if (detf .le. r8prem()) then
            iret = 1
            goto 999
        end if
        call deflog(ndim, fCurr, epslCurr, gn, lamb, logl, iret)
        if (iret .ne. 0) then
            goto 999
        end if
        epslIncr(1:6) = epslCurr(1:6)-epslPrev(1:6)
    end if

! - Get previous stress from internal state variables
    tlogPrev(1:2*ndim) = vim(lgpg-6+1:lgpg-6+1+2*ndim)*vrac2
!
999 continue
!
end subroutine
