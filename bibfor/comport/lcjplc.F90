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

subroutine lcjplc(rela_comp, mod, nmat, &
                  mater, timed, timef, compor, nbcomm, &
                  cpmono, pgl, nfs, nsg, toutms, &
                  hsr, nr, nvi, epsd, deps, &
                  itmax, toler, sigf, vinf, sigd, &
                  vind, dsde, drdy, option, iret)
! aslint: disable=W1504
    implicit none
!       MATRICE SYMETRIQUE DE COMPORTEMENT TANGENT ELASTO-PLASTIQUE OU
!       VISCO-PLASTIQUE COHERENT A T+DT OU T
!       COHERENT A T+DT OU T
!       IN  LOI    :  MODELE DE COMPORTEMENT
!           MOD    :  TYPE DE MODELISATION
!           NMAT   :  DIMENSION MATER
!           MATER  :  COEFFICIENTS MATERIAU
!       OUT DSDE   :  MATRICE DE COMPORTEMENT TANGENT = DSIG/DEPS
!       ----------------------------------------------------------------
#include "asterfort/cvmjpl.h"
#include "asterfort/lcmmjp.h"
#include "asterfort/lcoptg.h"
#include "asterfort/lkijpl.h"
#include "asterfort/srijpl.h"
#include "asterfort/Behaviour_type.h"
    integer(kind=8) :: nmat, nr, nvi, itmax, iret, nfs, nsg, ndt, ndi, n2
    real(kind=8) :: dsde(6, 6), epsd(*), deps(*), toler
    real(kind=8) :: mater(nmat, 2)
    real(kind=8) :: toutms(nfs, nsg, 6), hsr(nsg, nsg)
    character(len=8) :: mod
    character(len=16) :: option
    common/tdim/ndt, ndi
!
    integer(kind=8) :: nbcomm(nmat, 3)
    real(kind=8) :: sigf(*), sigd(*), vind(*), vinf(*), timed, timef, pgl(3, 3)
    real(kind=8) :: drdy(nr, nr)
    character(len=16), intent(in) :: rela_comp
    character(len=16), intent(in) :: compor(COMPOR_SIZE)
    character(len=24) :: cpmono(5*nmat+1)
!       ----------------------------------------------------------------
    iret = 0
    if (rela_comp .eq. 'VISCOCHAB') then
        call cvmjpl(mod, nmat, mater, timed, timef, &
                    epsd, deps, sigf, vinf, sigd, &
                    vind, nvi, nr, dsde)
    else if ((rela_comp .eq. 'MONOCRISTAL')) then
        call lcmmjp(mod, nmat, mater, timed, timef, &
                    compor, nbcomm, cpmono, pgl, nfs, &
                    nsg, toutms, hsr, nr, nvi, sigd, &
                    itmax, toler, vinf, vind, dsde, &
                    drdy, option, iret)
    else if (rela_comp .eq. 'LETK') then
        call lkijpl(nmat, mater, sigf, nr, drdy, &
                    dsde)
    else if (rela_comp .eq. 'LKR') then
        call srijpl(nmat, mater, sigf, nr, drdy, dsde)
    else if (rela_comp .eq. 'HAYHURST') then
        n2 = nr-ndt
        call lcoptg(nmat, mater, nr, n2, drdy, &
                    0, dsde, iret)
    else
        n2 = nr-ndt
        call lcoptg(nmat, mater, nr, n2, drdy, &
                    1, dsde, iret)
    end if
!
end subroutine
