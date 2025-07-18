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

subroutine nmelru(fami, kpg, ksp, &
                  imate, compor, epseq, p_arg, divu, gonf, inco, &
                  nonlin, ener)
!
! FONCTION REALISEE:
!
!     REALISE LE CALCUL DE L'ENERGIE LIBRE, DE LA DERIVEE DE L'ENERGIE
!             LIBRE PAR RAPPORT A LA TEMPERATURE (POUR LE CALCUL DE G)
!     ET DE SA DERIVEE PAR RAPPORT A UNE VARIATION DE DOMAINE (EN DP
!     ELASTIQUE ISOTROPE LINEAIRE).
!
! IN  IMATE   : NATURE DU MATERIAU
! IN  COMPOR  : COMPORTEMENT
! IN  EPSEQ   : DEFORMATION EQUIVALENTE
! IN  P       : DEFORMATION ELASTIQUE CUMULEE
! IN  DIVU    : TRACE DES DEFORMATIONS
! IN  NONLIN  : NON LINEARITE DU MATERIAU
! OUT ENER(1) : ENERGIE LIBRE (POUR LE CALCUL DE G)
! OUT ENER(2) : DERIVEE DE L'ENERGIE LIBRE / TEMPERATURE
!
! ----------------------------------------------------------------------

    use Behaviour_type
    implicit none
!
! DECLARATION PARAMETRES D'APPELS
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/rcfonc.h"
#include "asterfort/rctrac.h"
#include "asterfort/rctype.h"
#include "asterfort/rcvad2.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/utmess.h"
#include "asterfort/Behaviour_type.h"
!
    character(len=*), intent(in)  :: fami
    integer(kind=8), intent(in)          :: kpg
    integer(kind=8), intent(in)          :: ksp
    integer(kind=8), intent(in)          :: imate
    character(len=16), intent(in) :: compor(*)
    real(kind=8), intent(in)     :: epseq
    real(kind=8), intent(in)     :: p_arg
    real(kind=8), intent(in)     :: divu
    real(kind=8), intent(in)     :: gonf
    aster_logical, intent(in)    :: inco
    aster_logical, intent(in)    :: nonlin
    real(kind=8), intent(out)    :: ener(2)
! --------------------------------------------------------------------------------------------------
    character(len=1), parameter:: poum = '+'
    integer(kind=8) :: icodre(3)
    integer(kind=8) :: jprol, jvale, nbvale, iret1, iret2
    real(kind=8) :: temp, tref, p
    real(kind=8) :: e, nu, demu, k, k3, alpha, para_vale
    real(kind=8) :: de, dnu, demudt, dk, dalpha
    real(kind=8) :: dsde, sigy, rprim, rp, airep, coco
    real(kind=8) :: dsdedt, dsigy, drprim, dp, drp, dairep
    real(kind=8) :: nrj, dnrj, valres(3), devres(3)
    character(len=16) :: nomres(3)
    character(len=8) :: para_type
    aster_logical :: trac, line, puis
! --------------------------------------------------------------------------------------------------
    common/rconm2/alfafa, unsurn, sieleq
    real(kind=8) :: alfafa, unsurn, sieleq
! ----------------------------------------------------------------------
!
    p = p_arg
    trac = (compor(RELA_NAME) (1:14) .eq. 'ELAS_VMIS_TRAC')
    line = (compor(RELA_NAME) (1:14) .eq. 'ELAS_VMIS_LINE')
    puis = (compor(RELA_NAME) (1:14) .eq. 'ELAS_VMIS_PUIS')
!
!====================================================================
! -  LECTURE DE E, NU, ALPHA ET DERIVEES / TEMPERATRURE
!====================================================================
    call rcvarc(' ', 'TEMP', 'REF', fami, kpg, &
                ksp, tref, iret2)
!
    if (iret2 .ne. 0) tref = 0.d0
!
    if (fami(1:4) .eq. 'XFEM') then
        call rcvarc(' ', 'TEMP', poum, fami, kpg, &
                    1, temp, iret1)
    else
        call rcvarc(' ', 'TEMP', poum, 'RIGI', kpg, &
                    1, temp, iret1)
    end if
!
    if (iret1 .ne. 0) temp = 0.d0
!
    nomres(1) = 'E'
    nomres(2) = 'NU'
    nomres(3) = 'ALPHA'
    call rcvad2(fami, kpg, ksp, poum, imate, &
                'ELAS', 3, nomres, valres, devres, &
                icodre)
!
    if (icodre(3) .ne. 0) then
        valres(3) = 0.d0
        devres(3) = 0.d0
    end if
!
    e = valres(1)
    nu = valres(2)
    alpha = valres(3)
!
    de = devres(1)
    dnu = devres(2)
    dalpha = devres(3)
!
    demu = e/(1.d0+nu)
    demudt = ((1.d0+nu)*de-e*dnu)/(1.d0+nu)**2
!
    k = e/(1.d0-2.d0*nu)/3.d0
    dk = (de+6.d0*k*dnu)/(1.d0-2.d0*nu)/3.d0
!
    k3 = 3.d0*k
!
!====================================================================
! - LECTURE DES CARACTERISTIQUES DE NON LINEARITE DU MATERIAU
!====================================================================
!
!=================================================
! CAS NON LINEAIRE
!=================================================
    if (nonlin) then
        if (line) then
!
            nomres(1) = 'D_SIGM_EPSI'
            nomres(2) = 'SY'
!
            call rcvad2(fami, kpg, ksp, poum, imate, &
                        'ECRO_LINE', 2, nomres, valres, devres, &
                        icodre)
            if (icodre(1) .ne. 0) then
                call utmess('F', 'ALGORITH7_74')
            end if
            if (icodre(2) .ne. 0) then
                call utmess('F', 'ALGORITH7_75')
            end if
!
            dsde = valres(1)
            sigy = valres(2)
            dsdedt = devres(1)
            dsigy = devres(2)
!
            rprim = e*dsde/(e-dsde)
            drprim = (de*dsde+e*dsdedt+rprim*(dsdedt-de))/(e-dsde)
!
            p = (demu*epseq-sigy)/(rprim+1.5d0*demu)
            dp = (demudt*epseq-dsigy-p*(drprim+1.5d0*demudt))/(rprim+1.5d0*demu)
!
            rp = sigy+rprim*p
            drp = dsigy+drprim*p+rprim*dp
!
            airep = 0.5d0*(sigy+rp)*p
            dairep = 0.5d0*((dsigy+drp)*p+(sigy+rp)*dp)
!
        else if (trac) then
            sieleq = demu*epseq
            call rctype(imate, 1, 'TEMP', [temp], para_vale, &
                        para_type)
            if ((para_type .eq. 'TEMP') .and. (iret1 .eq. 1)) then
                call utmess('F', 'COMPOR5_5', sk=para_type)
            end if
            call rctrac(imate, 1, 'SIGM', para_vale, jprol, &
                        jvale, nbvale, e)
            call rcfonc('S', 1, jprol, jvale, nbvale, &
                        sigy=sigy)
            call rcfonc('E', 1, jprol, jvale, nbvale, &
                        e=e, nu=nu, p=0.d0, rp=rp, rprim=rprim, &
                        airerp=airep, sieleq=sieleq, dp=p)
            dp = 0.d0
            drp = 0.d0
            dairep = 0.d0
        else if (puis) then
            nomres(1) = 'SY'
!
            call rcvalb(fami, kpg, ksp, poum, imate, &
                        ' ', 'ECRO_PUIS', 0, ' ', [0.d0], &
                        1, nomres, valres, icodre, 1)
            sigy = valres(1)
            coco = e/alfafa/sigy
            rp = sigy*(coco*p)**unsurn+sigy
            airep = sigy*p+(1.d0/(1.d0+unsurn))*sigy*(coco**unsurn)*( &
                    p**(1+unsurn))
            dp = 0.d0
            drp = 0.d0
            dairep = 0.d0
        end if
    end if
!
!=====================================================================
!  CALCUL DE L'ENERGIE LIBRE ET DE LA DERIVEE /TEMPERATURE
!=====================================================================
!
!
    nrj = 0.5d0*k*divu*divu
!
!    POUR COMPARER AVEC CALC_K_G, enlever le terme constant
!    nrj = 0.5d0*k*divu*divu - 9/2*k*alpha*alpha*(temp- tref)*(temp- tref)
    if (iret1 .eq. 0) then
        if (iret2 .eq. 1) then
            call utmess('F', 'COMPOR5_43')
        else
            dnrj = 0.5d0*dk*divu*divu-k3*divu*(alpha+dalpha*(temp-tref))
        end if
    else
        dnrj = 0.5d0*dk*divu*divu-k3*divu*alpha
    end if
!
!    POUR COMPARER AVEC CALC_K_G, enlever le terme constant
!    dnrj = dnrj - 9*k*alpha*alpha*(temp- tref)
!
    if (nonlin) then
        if (inco) then
            ener(1) = 0.5d0*k*gonf*gonf+rp*rp/demu/3.d0+airep
        else
            ener(1) = nrj+rp*rp/demu/3.d0+airep
        end if
        ener(2) = dnrj+rp*(drp-demudt*rp/demu/2.d0)/demu/1.5d0+dairep
    else
        if (inco) then
            ener(1) = 0.5d0*k*gonf*gonf+demu*epseq*epseq/3.d0
        else
            ener(1) = nrj+demu*epseq*epseq/3.d0
        end if
        ener(2) = dnrj+demudt*epseq*epseq/3.d0
    end if
!
end subroutine
