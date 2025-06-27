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

subroutine cjsela(mod, crit, materf, deps, sigd, &
                  sigf, nvi, vind, vinf, iret)
    implicit none
!       INTEGRATION ELASTIQUE NON LINEAIRE DE LA LOI CJS
!       IN  MOD    :  MODELISATION
!           CRIT   : CRITERES DE CONVERGENCE
!           MATERF :  COEFFICIENTS MATERIAU A T+DT
!           SIGD   :  CONTRAINTE  A T
!           DEPS   :  INCREMENT DE DEFORMATION
!       OUT SIGF   :  CONTRAINTE A T+DT
!           IRET   : CODE RETOUR DE  L'INTEGRATION DE LA LOI CJS
!                              IRET=0 => PAS DE PROBLEME
!                              IRET=1 => ECHEC
!       ---------------------------------------------------------------
#include "asterf_types.h"
#include "asterfort/cjsci1.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: ndt, ndi, nvi, iret
    real(kind=8) :: coef, e, nu, al, la, mu, hook(6, 6), i1
    real(kind=8) :: deps(6), dsig(6), sigd(6), sigf(6)
    real(kind=8) :: vind(*), vinf(*)
    real(kind=8) :: materf(14, 2), crit(*)
    character(len=8) :: mod
    real(kind=8) :: zero, un, d12, deux, trois, pa, qinit
    aster_logical :: tract
    integer(kind=8) :: i, j
!
    common/tdim/ndt, ndi
!
    data zero/0.d0/
    data d12/.5d0/
    data un/1.d0/
    data deux/2.d0/
    data trois/3.d0/
!
!       ---------------------------------------------------------------
    pa = materf(12, 2)
    qinit = materf(13, 2)
!
!--->   CALCUL DE I1=TR(SIG) A T+DT PAR METHODE DE LA SECANTE
!       OU EXPLICITEMENT SI NIVEAU CJS1
!
    call cjsci1(crit, materf, deps, sigd, i1, &
                tract, iret)
    if (iret .eq. 1) goto 9999
!
!--->   EN CAS D'ENTREE EN TRACTION, LES CONTRAINTES SONT
!       RAMENEES SUR L'AXE HYDROSTATIQUE A DES VALEURS FAIBLES
!       ( EGALES A PA/100.0 SOIT -1 KPA )
!
    if (tract) then
        do i = 1, ndi
            sigf(i) = -qinit/3.d0+pa/100.0d0
        end do
        do i = ndi+1, ndt
            sigf(i) = zero
        end do
        goto 9999
    end if
!
!
!                         I1+QINIT
!--->   CALCUL DU COEF  (-----------)**N ET MODULE_YOUNG A T+DT
!                        3 PA
!
!
    coef = ((i1+qinit)/trois/pa)**materf(3, 2)
    e = materf(1, 1)*coef
    nu = materf(2, 1)
    al = e*(un-nu)/(un+nu)/(un-deux*nu)
    la = nu*e/(un+nu)/(un-deux*nu)
    mu = e*d12/(un+nu)
!
!--->   OPERATEUR DE RIGIDITE
!
    hook(:, :) = zero
!
! - 3D/DP/AX
    if (mod(1:2) .eq. '3D' .or. mod(1:6) .eq. 'D_PLAN' .or. mod(1:4) .eq. 'AXIS') then
        do i = 1, ndi
            do j = 1, ndi
                if (i .eq. j) hook(i, j) = al
                if (i .ne. j) hook(i, j) = la
            end do
        end do
        do i = ndi+1, ndt
            do j = ndi+1, ndt
                if (i .eq. j) hook(i, j) = deux*mu
            end do
        end do
!
! - CP/1D
    else if (mod(1:6) .eq. 'C_PLAN' .or. mod(1:2) .eq. '1D') then
        call utmess('F', 'ALGORITH2_15')
    end if
!
!
!--->   INCREMENTATION DES CONTRAINTES  SIGF = SIGD + HOOK DEPS
!
    dsig(1:ndt) = matmul(hook(1:ndt, 1:ndt), deps(1:ndt))
    sigf(1:ndt) = sigd(1:ndt)+dsig(1:ndt)
!
!
9999 continue
!
!--->   VINF = VIND, ETAT A T+DT = ELASTIQUE = 0
!
    vinf(1:nvi-1) = vind(1:nvi-1)
    vinf(nvi) = 0.d0
!
end subroutine
