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

subroutine epsvmc(fami, nno, ndim, nbsig, npg, &
                  j_poids, j_vf, j_dfde, xyz, disp, &
                  time, angl_naut, nharm, option, epsi)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dmatmc.h"
#include "asterfort/eps1mc.h"
#include "asterfort/eps2mc.h"
#include "asterfort/epslmc.h"
#include "asterfort/epthmc.h"
#include "asterfort/lteatt.h"
#include "asterfort/jevech.h"
!
!
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: nno
    integer(kind=8), intent(in) :: ndim
    integer(kind=8), intent(in) :: nbsig
    integer(kind=8), intent(in) :: npg
    integer(kind=8), intent(in) :: j_poids
    integer(kind=8), intent(in) :: j_vf
    integer(kind=8), intent(in) :: j_dfde
    real(kind=8), intent(in) :: xyz(1)
    real(kind=8), intent(in) :: disp(1)
    real(kind=8), intent(in) :: time
    real(kind=8), intent(in) :: angl_naut(3)
    real(kind=8), intent(in) :: nharm
    character(len=16), intent(in) :: option
    real(kind=8), intent(out) :: epsi(1)
!
! --------------------------------------------------------------------------------------------------
!
! Compute mechanical strains or total strains (depend on option)
!
! Mechanical strains = total strains - command variables strains
!
! --------------------------------------------------------------------------------------------------
!
! In  fami         : Gauss family for integration point rule
! In  nno          : number of nodes
! In  ndim         : dimension of space
! In  nbsig        : number of stress tensor components
! In  npg          : number of Gauss points
! In  j_poids      : JEVEUX adress to weight of Gauss points
! In  j_vf         : JEVEUX adress to shape functions
! In  j_dfde       : JEVEUX adress to derivatives of shape functions
! In  xyz          : coordinates of element
! In  disp         : displacements of element
! In  time         : current time
! In  j_mater      : coded material address
! In  angl_naut    : nautical angles (for non-isotropic materials)
! In  nharm        : Fourier mode
! In  option       : name of option to compute
! Out epsi         : mechanical strains or total strains
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: epsi_varc(162), epsi_tota_g(162), epsi_tota(162)
    real(kind=8) :: d(4, 4)
    real(kind=8) :: zero, un, deux
    integer(kind=8) :: i, kpg, imate
    aster_logical :: l_modi_cp
!
! --------------------------------------------------------------------------------------------------
!
    zero = 0.d0
    un = 1.d0
    deux = 2.d0
    ASSERT(nbsig*npg .le. 162)
    epsi(1:nbsig*npg) = zero
    epsi_tota(1:nbsig*npg) = zero
    epsi_tota_g(1:nbsig*npg) = zero
    epsi_varc(1:nbsig*npg) = zero
    if (option(4:4) .eq. 'L') then
        call epslmc(nno, ndim, nbsig, &
                    npg, j_poids, j_vf, &
                    j_dfde, xyz, disp, &
                    epsi)
    else
!
! - Total strains: first order (small strains)
!
        call eps1mc(nno, ndim, nbsig, npg, j_poids, &
                    j_vf, j_dfde, xyz, disp, nharm, &
                    epsi_tota)
!
! - Total strains: second order (large strains)
!
        if (option(4:4) .eq. 'G') then
            call eps2mc(nno, ndim, nbsig, npg, j_poids, &
                        j_vf, j_dfde, xyz, disp, epsi_tota_g)
        end if
!
! - Total strains
!
        do i = 1, nbsig*npg
            epsi(i) = epsi_tota(i)+epsi_tota_g(i)
        end do
    end if
!
! - Compute variable commands strains (thermics, drying, etc.)
!
    if (option(1:4) .eq. 'EPME' .or. option(1:4) .eq. 'EPMG' .or. lteatt('C_PLAN', 'OUI')) then
        call jevech('PMATERC', 'L', imate)
        call epthmc(fami, nno, ndim, nbsig, npg, &
                    zr(j_vf), angl_naut, time, zi(imate), &
                    option, epsi_varc)
    end if
!
! - Mechanical strains
!
    if (option(1:4) .eq. 'EPME' .or. option(1:4) .eq. 'EPMG') then
        do i = 1, nbsig*npg
            epsi(i) = epsi_tota(i)+epsi_tota_g(i)-epsi_varc(i)
        end do
    end if
!
! - 2D model
!
    if (lteatt('C_PLAN', 'OUI')) then
!
! ----- Plane stress
!
        do kpg = 1, npg
!
! --------- il s'agit de calculer EPS33 : pour cela il faut donner la
! --------- condition SIG33=0 dans l'expression complete de la loi de
! --------- Hooke c'est à dire avec la loi 3D :
! --------- Eps33= -1/D33 (D13.Eps11 +D12.Eps22), ce qui donne (en
! --------- isotrope) l'expression classique :
! --------- Eps33 = -Nu / (1-Nu) * (Eps11 + Eps22).
! --------- voir issue12540
!
            l_modi_cp = .true.
!
! --------- Hooke matrix for iso-parametric elements
!
            call dmatmc(fami, zi(imate), time, '+', kpg, &
                        1, angl_naut, nbsig, d, &
                        l_modi_cp)
!
            if (option(1:4) .eq. 'EPME' .or. option(1:4) .eq. 'EPMG') then
                epsi(nbsig*(kpg-1)+3) = -un/d(3, 3)* &
                                        (d(3, 1)*epsi(nbsig*(kpg-1)+1)+ &
                                         d(3, 2)*epsi(nbsig*(kpg-1)+2)+ &
                                         d(3, 4)*epsi(nbsig*(kpg-1)+4)*deux)
            else
                epsi(nbsig*(kpg-1)+3) = -un/d(3, 3)* &
                                        (d(3, 1)*(epsi(nbsig*(kpg-1)+1)- &
                                                  epsi_varc(nbsig*(kpg-1)+1))+ &
                                         d(3, 2)*(epsi(nbsig*(kpg-1)+2)- &
                                                  epsi_varc(nbsig*(kpg-1)+2))+ &
                                         d(3, 4)*(epsi(nbsig*(kpg-1)+4)- &
                                                  epsi_varc(nbsig*(kpg-1)+4))*deux)+ &
                                        epsi_varc(nbsig*(kpg-1)+3)
            end if
        end do
    else if (lteatt('D_PLAN', 'OUI')) then
!
! ----- Plane strain: EPZZ = 0
!
        do kpg = 1, npg
            epsi(nbsig*(kpg-1)+3) = zero
        end do
    end if
!
end subroutine
