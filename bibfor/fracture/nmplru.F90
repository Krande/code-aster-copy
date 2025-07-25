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

subroutine nmplru(fami, kpg, ksp, poum, ndim, &
                  typmod, imate, compor, ppg, eps, &
                  epsp, rp, ener)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/rcfonc.h"
#include "asterfort/rctrac.h"
#include "asterfort/rcvad2.h"
#include "asterfort/rcvarc.h"
#include "asterfort/utmess.h"
!
    integer(kind=8) :: kpg, ksp, ndim, imate
    character(len=*) :: fami, poum
    character(len=8) :: typmod(*)
    character(len=16) :: compor(*)
    real(kind=8) :: ppg, eps(6), epsp(6), ener(2)
!.......................................................................
!
!     REALISE LE CALCUL DE L'ENERGIE LIBRE ET DE LA DERIVEE DE L'ENERGIE
!             LIBRE PAR RAPPORT A LA TEMPERATURE (POUR LE CALCUL DE G)
!             EN PLASTICITE
!

! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  TYPMOD  : TYPE DE MODELISATION
! IN  IMATE   : NATURE DU MATERIAU
! IN  COMPOR  : COMPORTEMENT
! IN  PPG     : DEFORMATION PLASTIQUE CUMULEE
! IN  EPS     : DEFORMATION TOTALE
! IN  EPSP    : DEFORMATION PLASTIQUE
!
! OUT RP      :
! OUT ENER(1) : ENRGIE LIBRE
! OUT ENER(1) : DERIVEE DE L'ENERGIE LIBRE / A LA TEMPERATURE
!.......................................................................
!
    integer(kind=8) :: icodre(3)
    character(len=16) :: nomres(3)
!
    real(kind=8) :: e, nu, demu, k, k3, alpha
    real(kind=8) :: de, dnu, demudt, dk, dalpha
    real(kind=8) :: dsde, sigy, rprim, rp, airep
    real(kind=8) :: dsdedt, dsigy, drprim, drp, dairep
    real(kind=8) :: nrj, dnrj, valres(3), devres(3)
    real(kind=8) :: epsth(6), epsdv(6), epseq, kron(6)
    real(kind=8) :: ther, divu, epsmo, temp, tref
!
    integer(kind=8) :: i, jprol, jvale, nbval
!
    aster_logical :: cp, trac, line, elas
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: iret1, iret2
!-----------------------------------------------------------------------
    data kron/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/
!
    cp = typmod(1) .eq. 'C_PLAN'
    trac = compor(1) (1:14) .eq. 'VMIS_ISOT_TRAC'
    line = compor(1) (1:14) .eq. 'VMIS_ISOT_LINE'
    elas = compor(1) (1:16) .eq. 'ELAS            '
!
! -  LECTURE DE E, NU, ALPHA ET DERIVEES / TEMPERATRURE
!
    call rcvarc(' ', 'TEMP', poum, fami, kpg, &
                ksp, temp, iret1)
    call rcvarc(' ', 'TEMP', 'REF', fami, 1, &
                1, tref, iret2)
    if (iret1 .eq. 1) temp = 0.d0
    if (iret2 .eq. 1) tref = 0.d0
!
!
    nomres(1) = 'E'
    nomres(2) = 'NU'
    nomres(3) = 'ALPHA'
    call rcvad2(fami, kpg, ksp, poum, imate, &
                'ELAS', 3, nomres, valres, devres, &
                icodre)
!
    if (iret1 .eq. 0) then
        if ((iret2 .ge. 1) .or. (icodre(3) .ne. 0)) then
            call utmess('F', 'CALCULEL_15')
        else
            alpha = valres(3)
            dalpha = devres(3)
!          CALL RCVAD2 (FAMI,KPG,KSP,POUM,IMATE,'ELAS',3,
!     &             NOMRES,VALRES,DEVRES,ICODRE)
!
        end if
    else
        if (icodre(3) .eq. 0) then
            alpha = valres(3)
            dalpha = devres(3)
        else
            alpha = 0
            dalpha = 0
        end if
    end if
!
    e = valres(1)
    nu = valres(2)
!
    de = devres(1)
    dnu = devres(2)
!
    demu = e/(1.d0+nu)
    demudt = ((1.d0+nu)*de-e*dnu)/(1.d0+nu)**2
!
    k = e/(1.d0-2.d0*nu)/3.d0
    dk = (de+2.d0*k*dnu)/(1.d0-2.d0*nu)/3.d0
!
    k3 = 3.d0*k
!
! - LECTURE DES CARACTERISTIQUES DE NON LINEARITE DU MATERIAU
!
    airep = 0.d0
    dairep = 0.d0
!
    if (line) then
        nomres(1) = 'D_SIGM_EPSI'
        nomres(2) = 'SY'
        call rcvad2(fami, kpg, ksp, poum, imate, &
                    'ECRO_LINE', 2, nomres, valres, devres, &
                    icodre)
        if (icodre(1) .ne. 0) then
            call utmess('F', 'ALGORITH7_74')
        end if
        if (icodre(2) .ne. 0) then
            call utmess('F', 'ALGORITH7_75')
        end if
        dsde = valres(1)
        sigy = valres(2)
        dsdedt = devres(1)
        dsigy = devres(2)
!
        rprim = e*dsde/(e-dsde)
        drprim = (de*dsde+e*dsdedt+rprim*(dsdedt-de))/(e-dsde)
!
        rp = sigy+rprim*ppg
        drp = dsigy+drprim*ppg
!
        airep = 0.5d0*(sigy+rp)*ppg
        dairep = 0.5d0*(dsigy+drp)*ppg
!
    else if (trac) then
        call rctrac(imate, 1, 'SIGM', temp, jprol, &
                    jvale, nbval, e)
        call rcfonc('V', 1, jprol, jvale, nbval, &
                    p=ppg, rp=rp, rprim=rprim, airerp=airep)
        dairep = 0.d0
!
    else if (elas) then
        rp = 0.d0
!
    else
        ASSERT(.false.)
!
    end if
!
! - CALCUL DE EPSMO ET EPSDV
    if ((iret1+iret2) .eq. 0) then
        ther = alpha*(temp-tref)
    else
        ther = 0.d0
    end if
!
    if (cp) eps(3) = -nu/(1.d0-nu)*(eps(1)+eps(2))+(1.d0+nu)/(1.d0-nu)*ther
    divu = 0.d0
    do i = 1, 3
        epsth(i) = eps(i)-epsp(i)-ther
        epsth(i+3) = eps(i+3)-epsp(i+3)
        divu = divu+epsth(i)
    end do
    epsmo = divu/3.d0
    do i = 1, 2*ndim
        epsdv(i) = epsth(i)-epsmo*kron(i)
    end do
!
! - CALCUL DE LA CONTRAINTE ELASTIQUE EQUIVALENTE
    epseq = 0.d0
    do i = 1, 2*ndim
        epseq = epseq+epsdv(i)*epsdv(i)
    end do
    epseq = sqrt(1.5d0*epseq)
!
!  CALCUL DE L'ENERGIE LIBRE ET DE LA DERIVEE /TEMPERATURE
!
    nrj = 0.5d0*k*divu*divu+demu*epseq*epseq/3.d0
    dnrj = 0.5d0*dk*divu*divu-k3*divu*(alpha+dalpha*(temp-tref))+demudt*epseq*epseq/3.d0
!
    ener(1) = nrj+airep
    ener(2) = dnrj+dairep
!
end subroutine
