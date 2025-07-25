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
subroutine nmelas_elas(BEHinteg, &
                       fami, kpg, ksp, typmod, &
                       imate, eps, option, sig, &
                       vi, dsidep)
!
    use Behaviour_type
    use tenseur_dime_module, only: sph_norm, deviator, kron, voigt, proten, identity
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/verift.h"
#include "asterfort/verifh.h"
#include "asterfort/verifs.h"
#include "asterfort/verifepsa.h"
#include "asterfort/get_elas_para.h"
!
    type(Behaviour_Integ), intent(in) :: BEHinteg
    character(len=*), intent(in)      :: fami
    character(len=8), intent(in)      :: typmod(*)
    character(len=16), intent(in)     :: option
    integer(kind=8), intent(in)               :: imate, kpg, ksp
    real(kind=8), intent(in)          :: eps(:)
    real(kind=8), intent(out)         :: sig(:), vi(1), dsidep(:, :)
! --------------------------------------------------------------------------------------------------
!     REALISE LA LOI DE VON MISES ISOTROPE ET ELASTIQUE POUR LES
!     ELEMENTS ISOPARAMETRIQUES EN PETITES DEFORMATIONS
! --------------------------------------------------------------------------------------------------
! IN  KPG,KSP : NUMERO DU (SOUS)POINT DE GAUSS
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  TYPMOD  : TYPE DE MODELISATION
! IN  IMATE   : ADRESSE DU MATERIAU CODE
! IN  COMPOR  : COMPORTEMENT : RELCOM ET DEFORM
! IN  CRIT    : CRITERES DE CONVERGENCE LOCAUX
! IN  EPS     : DEFORMATION A L'INSTANT ACTUEL
! IN  OPTION  : OPTION DEMANDEE : RIGI_MECA_TANG , FULL_MECA , RAPH_MECA
! OUT SIG     : CONTRAINTES A L'INSTANT ACTUEL
! OUT VI      : VARIABLES INTERNES A L'INSTANT ACTUEL
! OUT DSIDEP  : MATRICE CARREE (INUTILISE POUR RAPH_MECA)
!
!               ATTENTION LES TENSEURS ET MATRICES SONT RANGES DANS
!               L'ORDRE :  XX,YY,ZZ,SQRT(2)*XY,SQRT(2)*XZ,SQRT(2)*YZ
! --------------------------------------------------------------------------------------------------
    integer(kind=8), parameter :: elas_id = 1
    character(len=16), parameter :: elas_keyword = 'ELAS'
! --------------------------------------------------------------------------------------------------
    aster_logical :: cplan, resi, rigi
    integer(kind=8)       :: ndimsi, k, l
    real(kind=8)  :: kr(size(eps))
    real(kind=8)  :: e, nu, lambda, deuxmu, troisk
    real(kind=8)  :: eps_th, eps_hy, eps_se, eps_an(6), eps_vc(size(eps))
    real(kind=8)  :: eps_me(size(eps))
    character(len=1):: poum
! --------------------------------------------------------------------------------------------------

    ! Initialisation
    ndimsi = size(eps)
    kr = kron(ndimsi)
    cplan = typmod(1) .eq. 'C_PLAN'
    rigi = option(1:10) .eq. 'RIGI_MECA_' .or. option(1:9) .eq. 'FULL_MECA'
    resi = option(1:9) .eq. 'RAPH_MECA' .or. option(1:9) .eq. 'FULL_MECA' .or. &
           option .eq. 'RIGI_MECA_IMPLEX'
    poum = merge('+', '-', resi)

    ! Caracteristiques elastiques et contraintes initiales
    call get_elas_para(fami, imate, poum, kpg, ksp, elas_id, elas_keyword, &
                       e_=e, nu_=nu, BEHinteg=BEHinteg)
    lambda = e*nu/((1-2*nu)*(1+nu))
    deuxmu = e/(1+nu)
    troisk = e/(1-2*nu)

    ! Variables de commande
    call verift(fami, kpg, ksp, poum, imate, epsth_=eps_th)
    call verifh(fami, kpg, ksp, poum, imate, eps_hy)
    call verifs(fami, kpg, ksp, poum, imate, eps_se)
    call verifepsa(fami, kpg, ksp, poum, eps_an)
    eps_vc = (eps_th+eps_hy+eps_se)*kr+eps_an(1:ndimsi)*voigt(ndimsi)

    ! Calcul de la contrainte

    if (resi) then
        eps_me = eps-eps_vc
        if (cplan) eps_me(3) = -lambda*(eps_me(1)+eps_me(2))/(lambda+deuxmu)

        sig = lambda*dot_product(kr, eps_me)*kr+deuxmu*eps_me
        vi(1) = 0.d0
    end if

    ! Matrice tangente (ancienne formulation de la prediction)

    if (rigi) then
        dsidep = lambda*proten(kr, kr)+deuxmu*identity(ndimsi)

        if (cplan) then
            do k = 1, ndimsi
                if (k .ne. 3) then
                    do l = 1, ndimsi
                        if (l .ne. 3) then
                            dsidep(k, l) = dsidep(k, l)-1.d0/dsidep(3, 3)*dsidep(k, 3)*dsidep(3, l)
                        end if
                    end do
                end if
            end do
        end if
    end if

end subroutine
