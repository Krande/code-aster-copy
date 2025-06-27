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
! aslint: disable=W0413,W0104
! comparaison aver r8gaem et -r8gaem uniquement
! passage de sigm comme argunement si des lois devaient un jour en avoir besoin (COMP_INCR)

subroutine pil000(typilo, compor, neps, dtau, mat, &
                  vim, sigm, epsm, epsd_cste, epsd_pilo, &
                  typmod, etamin, etamax, copilo)
!
!
    implicit none
#include "asterc/r8gaem.h"
#include "asterfort/assert.h"
#include "asterfort/lcmfbo.h"
#include "asterfort/lcmfga.h"
#include "asterfort/lcmfma.h"
#include "asterfort/lcqubo.h"
#include "asterfort/lcquga.h"
#include "asterfort/lcquma.h"
#include "asterfort/pidegv.h"
#include "asterfort/pieigv.h"
#include "asterfort/piesgv.h"
#include "asterfort/pipeab.h"
#include "asterfort/pipeex.h"
#include "asterfort/pipeou.h"
#include "asterfort/pipetc.h"
#include "asterfort/utmess.h"
    character(len=8), intent(in) :: typmod(*)
    character(len=16), intent(in) :: compor(*), typilo
    integer(kind=8), intent(in) :: neps, mat
    real(kind=8), intent(in) :: dtau, epsm(neps), epsd_pilo(neps), epsd_cste(neps), etamin, etamax
    real(kind=8), intent(in) :: vim(:), sigm(neps)
    real(kind=8), intent(out) :: copilo(5)
!
!---------------------------------------------------------------------------------------------------
!     PILOTAGE PRED_ELAS : BRANCHEMENT SELON COMPORTEMENT
!---------------------------------------------------------------------------------------------------
! IN  TYPILO  TYPE DE PILOTAGE : 'PRED_ELAS' OU 'DEFORMATION'
! IN  NEPS    DIMENSION DES DEFORMATIONS
! IN  TAU     INCREMENT DE PILOTAGE
! IN  MAT     NATURE DU MATERIAU                             (PRED_ELAS)
! IN  VIM     VARIABLES INTERNES EN T-                       (PRED_ELAS)
! IN  SIGM    CONTRAINTES EN T- (SI NECESSAIRE)              (PRED_ELAS)
! IN  EPSM    CHAMP DE DEFORMATION EN T-
! IN  EPSD_CSTE    INCREMENT FIXE
! IN  EPSD_PILO    INCREMENT PILOTE
! IN  ETAMIN  BORNE INF DU PILOTAGE (SI UTILE)               (PRED_ELAS)
! IN  ETAMAX  BORNE SUP DU PILOTAGE (SI UTILE)               (PRED_ELAS)
! OUT COPILO  COEFFICIENT DE PILOTAGE : F := A0+A1*ETA = TAU
!---------------------------------------------------------------------------------------------------
    integer(kind=8) :: ndim, nsol, sgn(2)
    real(kind=8):: mu_cste(3), su_cste(3), mu_pilo(3), su_pilo(3)
    real(kind=8):: sol(2)
!---------------------------------------------------------------------------------------------------

! EN ATTENTE D'UNE HARMONISATION DU TRAITEMENT (S'INSPIRER DE CZM_LAB_MIX ET DE LC0000)

!---------------------------------------------------------------------------------------------------
! MODELISATION A GRADIENT DE VARIABLES INTERNES
!---------------------------------------------------------------------------------------------------

    if (typmod(2) .eq. 'GRADVARI') then

        ! PILOTAGE 'DEFORMATION'
        if (typilo .eq. 'DEFORMATION') then
            call pidegv(neps, dtau, epsm, epsd_cste, epsd_pilo, copilo)

            ! PILOTAGE 'PRED_ELAS'
        else
            if (etamin .eq. -r8gaem() .or. etamax .eq. r8gaem()) &
                call utmess('F', 'MECANONLINE_60', sk=compor(1))

            if (compor(1) .eq. 'ENDO_SCALAIRE') then
                call piesgv(neps, dtau, mat, lcquma, vim, epsm, epsd_cste, epsd_pilo, typmod, &
                            lcquga, etamin, etamax, lcqubo, copilo)

            else if (compor(1) .eq. 'ENDO_FISS_EXP') then
                call piesgv(neps, dtau, mat, lcmfma, vim, epsm, epsd_cste, epsd_pilo, typmod, &
                            lcmfga, etamin, etamax, lcmfbo, copilo)

            else if (compor(1) .eq. 'ENDO_ISOT_BETON') then
                call pieigv(neps, dtau, mat, vim, epsm, epsd_cste, epsd_pilo, typmod, &
                            etamin, etamax, copilo)

            else
                call utmess('F', 'MECANONLINE_59')
            end if
        end if

!---------------------------------------------------------------------------------------------------
! MODELISATION CZM INTERFACE
!---------------------------------------------------------------------------------------------------

    else if (typmod(2) .eq. 'INTERFAC') then

        ASSERT(typilo .ne. 'DEFORMATION')
        ndim = neps/2
        su_cste = 0
        su_pilo = 0
        mu_cste = 0
        mu_pilo = 0
        su_cste(1:ndim) = epsm(1:ndim)+epsd_cste(1:ndim)
        su_pilo(1:ndim) = epsd_pilo(1:ndim)
        mu_cste(1:ndim) = epsm(ndim+1:2*ndim)+epsd_cste(ndim+1:2*ndim)
        mu_pilo(1:ndim) = epsd_pilo(ndim+1:2*ndim)

        if (compor(1) .eq. 'CZM_TAC_MIX') then
            call pipetc(mat, su_cste, su_pilo, mu_cste, mu_pilo, &
                        vim, dtau, copilo)
        else if (compor(1) .eq. 'CZM_OUV_MIX') then
            call pipeou(mat, su_cste, su_pilo, mu_cste, mu_pilo, &
                        vim, dtau, copilo)
        else if (compor(1) .eq. 'CZM_EXP_MIX') then
            call pipeex(mat, su_cste, su_pilo, mu_cste, mu_pilo, &
                        vim, dtau, copilo)
        else if (compor(1) .eq. 'CZM_LAB_MIX') then
            call pipeab(mat, dtau, vim(:), su_cste, su_pilo, mu_cste, mu_pilo, nsol, sol, sgn)

            if (nsol .eq. 0) then
                copilo(5) = 0.d0
            else if (nsol .eq. 1) then
                copilo(1) = dtau-sgn(1)*sol(1)
                copilo(2) = sgn(1)
            else if (nsol .eq. 2) then
                copilo(1) = dtau-sgn(1)*sol(1)
                copilo(2) = sgn(1)
                copilo(3) = dtau-sgn(2)*sol(2)
                copilo(4) = sgn(2)
            end if

        else
            call utmess('F', 'MECANONLINE_59')
        end if

    else
        call utmess('F', 'MECANONLINE_61', sk=typmod(2))
    end if

end subroutine
