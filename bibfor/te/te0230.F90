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

subroutine te0230(option, nomte)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/r8vide.h"
#include "asterfort/bmatmc.h"
#include "asterfort/btdbmc.h"
#include "asterfort/dmatmc.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/nbsigm.h"
#include "asterfort/rcangm.h"
#include "asterfort/get_elas_id.h"
#include "asterfort/pmfmats.h"
#include "asterfort/rcvalb.h"
#include "asterfort/tecach.h"
!
    character(len=16), intent(in) :: option
    character(len=16), intent(in) :: nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 2D
! Option: RIGI_MECA_HYST
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nbres, nbpar
    parameter(nbres=2)
    parameter(nbpar=3)
!
    integer(kind=8) :: i, igau, imate, imatuu, j, iret, idrigi(2), rigi
    integer(kind=8) :: k, nbinco, nbsig, ndim, nno
    integer(kind=8) :: nnos, npg1
    integer(kind=8) :: icodre(nbres), elas_id
    integer(kind=8) :: igeom, ipoids, ivf, idfde, idim
!
    real(kind=8) :: b(486), jacgau
    real(kind=8) :: btdbi(81, 81), di(36), eta
    real(kind=8) :: angl_naut(3), instan, nharm
    real(kind=8) :: bary(3)
    real(kind=8) :: valres(nbres)
!
    character(len=4) :: fami
    character(len=8) :: nompar(nbpar), nomat
    character(len=16) :: nomres(nbres)
    character(len=32) :: phenom
!
! --------------------------------------------------------------------------------------------------
!
!
! - Finite element informations
!
    fami = 'RIGI'
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, nnos=nnos, &
                     npg=npg1, jpoids=ipoids, jvf=ivf, jdfde=idfde)
!
! - Initializations
!
    instan = r8vide()
    nbinco = ndim*nno
    nharm = 0.d0
    btdbi(:, :) = 0.d0
    bary(:) = 0.d0
!
! - Number of stress components
!
    nbsig = nbsigm()
!
! - Geometry
!
    call jevech('PGEOMER', 'L', igeom)
!
! - Material parameters
!
    call jevech('PMATERC', 'L', imate)
!
! - Dilatation coefficients
!
    call pmfmats(nomat)
!
! - Get type of elasticity (Isotropic/Orthotropic/Transverse isotropic)
!
    call get_elas_id(zi(imate), elas_id, phenom)
!
! - Orthotropic parameters
!
    do i = 1, nno
        do idim = 1, ndim
            bary(idim) = bary(idim)+zr(igeom+idim+ndim*(i-1)-1)/nno
        end do
    end do
    call rcangm(ndim, bary, angl_naut)
!
    nompar(1) = 'X'
    nompar(2) = 'Y'
    nompar(3) = 'Z'
!
! - Get RIGI_MECA real part
!
    call tecach('ONO', 'PRIGIEL', 'L', iret, nval=2, itab=idrigi)
!
! - Compute RIGI_MECA imaginary part
!
!
! - Case of a viscoelastic materials
!
    if (phenom .eq. 'ELAS_VISCO' .or. &
        phenom .eq. 'ELAS_VISCO_ISTR' .or. &
        phenom .eq. 'ELAS_VISCO_ORTH') then
        do igau = 1, npg1
            !
            ! ----- Compute matrix [B]: displacement -> strain (first order)
            !
            call bmatmc(igau, nbsig, zr(igeom), ipoids, ivf, &
                        idfde, nno, nharm, jacgau, b)
            !
            ! ---------- Compute Hooke matrix [D]
            !
            call dmatmc(fami, zi(imate), instan, '+', &
                        igau, 1, angl_naut, nbsig, &
                        di_=di)
            !
            ! --------- Compute rigidity matrix [K] = [B]Tx[D]x[B]
            !
            call btdbmc(b, di, jacgau, ndim, nno, &
                        nbsig, elas_id, btdbi)
        end do
!
! ---- Case of an elastic material
!
    else

        nomres(1) = 'AMOR_HYST'
        valres(1) = 0.d0
        call rcvalb('RIGI', 1, 1, '+', zi(imate), &
                    nomat, phenom, ndim, nompar, bary, &
                    1, nomres, valres, icodre, 0, &
                    nan='NON')
        eta = valres(1)
!
    end if
!
! - Set matrix in output field
!
    rigi = idrigi(1)
    call jevech('PMATUUC', 'E', imatuu)
    k = 0
    do i = 1, nbinco
        do j = 1, i
            k = k+1
            if (phenom .eq. 'ELAS_VISCO' .or. &
                phenom .eq. 'ELAS_VISCO_ISTR' .or. &
                phenom .eq. 'ELAS_VISCO_ORTH') then
                zc(imatuu+k-1) = dcmplx(zr(rigi+k-1), btdbi(i, j))
            else
                zc(imatuu+k-1) = dcmplx(zr(rigi+k-1), eta*zr(rigi+k-1))
            end if
        end do
    end do
!
end subroutine
