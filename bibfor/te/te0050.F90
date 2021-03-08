! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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

subroutine te0050(option, nomte)
!
implicit none
!
#include "jeveux.h"
#include "asterfort/bmatmc.h"
#include "asterfort/btdbmc.h"
#include "asterfort/dmatmc.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/nbsigm.h"
#include "asterfort/ortrep.h"
#include "asterfort/get_elas_id.h"
#include "asterfort/pmfmats.h"
#include "asterfort/rccoma.h"
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
! Elements: 3D
! Option: RIGI_MECA_HYST
!
! --------------------------------------------------------------------------------------------------
!
    integer :: nbres, nbpar
    parameter  ( nbres=2 )
    parameter  ( nbpar=3 )
!
    integer :: i, idecno, idecpg, igau, imate, imatuu, j, iret, idrigi(2), rigi
    integer :: k, nbinco, nbsig, ndim, nno
    integer :: nnos, npg1
    integer :: icodre(nbres), elas_id
    integer :: igeom, ipoids, ivf, idfde, idim, nbval
!
    real(kind=8) :: b(486), jacgau
    real(kind=8) :: btdbi(81, 81), di(36), eta
    real(kind=8) :: repere(7), xyzgau(3), instan, nharm
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
    call elrefe_info(fami=fami,ndim=ndim,nno=nno,nnos=nnos,&
                        npg=npg1,jpoids=ipoids,jvf=ivf,jdfde=idfde)
!
! - Initializations
!
    instan     = 0.d0
    nbinco     = ndim*nno
    nharm      = 0.d0
    btdbi(:,:) = 0.d0
    xyzgau(:)  = 0.d0
    bary(:)    = 0.d0
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
    call pmfmats(zi(imate), nomat)
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
    call ortrep(ndim, bary, repere)
!
    nompar(1)='X'
    nompar(2)='Y'
    nompar(3)='Z'
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
    if (phenom .eq. 'ELAS_VISCO'.or. &
        phenom .eq. 'ELAS_VISCO_ISTR' .or.&
        phenom .eq. 'ELAS_VISCO_ORTH') then

!
! - Number of stress components
!
        nbsig     = nbsigm()
!
        do igau = 1, npg1
!
            idecpg = nno* (igau-1) - 1
    !
    ! ----- Coordinates for current Gauss point
    !
            xyzgau(:) = 0.d0
            do i = 1, nno
                idecno = 3* (i-1) - 1
                xyzgau(1) = xyzgau(1) + zr(ivf+i+idecpg)*zr(igeom+1+idecno)
                xyzgau(2) = xyzgau(2) + zr(ivf+i+idecpg)*zr(igeom+2+idecno)
                xyzgau(3) = xyzgau(3) + zr(ivf+i+idecpg)*zr(igeom+3+idecno)
            end do
    !
    ! ----- Compute matrix [B]: displacement -> strain (first order)
    !
            call bmatmc(igau, nbsig, zr(igeom), ipoids, ivf,&
                        idfde, nno, nharm, jacgau, b)
    !
    ! ---------- Compute Hooke matrix [D]
    !
                call dmatmc(fami, zi(imate), instan, '+',&
                            igau, 1, repere, xyzgau, nbsig,&
                            di_=di)
    !
    ! --------- Compute rigidity matrix [K] = [B]Tx[D]x[B]
    !
                call btdbmc(b, di, jacgau, ndim, nno,&
                            nbsig, elas_id, btdbi)
        enddo
!
! ---- Case of an elastic material
!
    else

        nomres(1)='AMOR_HYST'
        valres(1) = 0.d0
        call rcvalb('RIGI', 1, 1, '+', zi(imate),&
                    nomat, phenom, ndim, nompar, bary,&
                    1, nomres, valres, icodre, 0,&
                    nan='NON')
        eta = valres(1)
!
    endif
!
! - Set matrix in output field
!
    rigi = idrigi(1)
    call jevech('PMATUUC', 'E', imatuu)
    if (phenom .eq. 'ELAS_VISCO'.or. &
        phenom .eq. 'ELAS_VISCO_ISTR' .or.&
        phenom .eq. 'ELAS_VISCO_ORTH') then
        k = 0
        do i = 1, nbinco
            do j = 1, i
                k = k + 1
                zc(imatuu+k-1) = dcmplx(zr(rigi+k-1), btdbi(i,j))
            end do
        end do
    else
        nbval = idrigi(2)
        do k = 1, nbval
            zc(imatuu+k-1) = dcmplx(zr(rigi+k-1), eta*zr(rigi+k-1))
        end do
    endif
!
end subroutine
