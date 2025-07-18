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
subroutine te0082(option, nomte)
!
    use THM_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/thmGetGene.h"
#include "asterfort/thmGetElemModel.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/pmavec.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/vecma.h"
#include "blas/ddot.h"
!
    character(len=16) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
!
!    - CALCULE DES MATRICES ELEMENTAIRES
!                          OPTION : 'MASS_MECA'
!    - CALCULE DES VECTEURS ELEMENTAIRES
!                          OPTION : 'M_GAMMA'
!    - CALCULE DES GRANDEURS ELEMENTAIRES
!                          OPTION : 'ECIN_ELEM'
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: phenom
    character(len=3) :: stopz
    integer(kind=8) :: icodre(1)
    real(kind=8) :: valres(1), poids, r, vfi, vfj
    real(kind=8) :: matp(18, 18), matv(171), masvit(18), masdep(18)
    real(kind=8) :: vect1(18), vect2(18)
    integer(kind=8) :: nno, kp, nnos, npg2, ii, jj, i, j, k, imatuu
    integer(kind=8) :: l, n1, n2, i2, j2
    integer(kind=8) :: ipoids, ivf, idfde, igeom, imate
    integer(kind=8) :: kd1, kd2, ij1, ij2, nddl, nvec, iacce, ivect
    integer(kind=8) :: idepl, ivite, ifreq, iecin
    integer(kind=8) :: mecani(5), press1(7), press2(7), tempe(5), second(5)
    integer(kind=8) :: idec, iret
    aster_logical :: l_axi
    aster_logical, parameter :: l_vf = ASTER_FALSE
    type(THM_DS) :: ds_thm
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami='MASS', nno=nno, nnos=nnos, npg=npg2, jpoids=ipoids, &
                     jvf=ivf, jdfde=idfde)
    nddl = 2*nno
    nvec = nddl*(nddl+1)/2
!
! - Get model of finite element
!
    call thmGetElemModel(ds_thm, l_axi_=l_axi)
!
! - Get generalized coordinates
!
    call thmGetGene(ds_thm, l_vf, 2, mecani, press1, &
                    press2, tempe, second)
    if (lteatt('TYPMOD2', 'THM')) then
        idec = press1(1)+press2(1)+tempe(1)+second(1)*2
    else
        idec = 0
    end if
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
!
    call rccoma(zi(imate), 'ELAS', 1, phenom, icodre(1))
    call rcvalb('FPG1', 1, 1, '+', zi(imate), &
                ' ', phenom, 0, ' ', [0.d0], &
                1, 'RHO', valres, icodre(1), 1)
!
    matv(:) = 0.d0
!
    do kp = 1, npg2
        k = (kp-1)*nno
        call dfdm2d(nno, kp, ipoids, idfde, zr(igeom), &
                    poids)
        if (lteatt('AXIS', 'OUI')) then
            r = 0.0d0
            do i = 1, nno
                r = r+zr(igeom+2*(i-1))*zr(ivf+k+i-1)
            end do
            poids = poids*r
        end if
        poids = poids*valres(1)
!
        kd1 = 2
        kd2 = 1
        do i = 1, 2*nno, 2
            kd1 = kd1+2*i-3
            kd2 = kd2+2*i-1
            ii = (i+1)/2
            do j = 1, i, 2
                jj = (j+1)/2
                ij1 = kd1+j-1
                ij2 = kd2+j-1
                vfi = zr(ivf+k+ii-1)
                vfj = zr(ivf+k+jj-1)
                matv(ij1) = matv(ij1)+poids*vfi*vfj
                matv(ij2+1) = matv(ij2+1)+poids*vfi*vfj
            end do
        end do
    end do
!
    if (option .eq. 'MASS_MECA') then
        call jevech('PMATUUR', 'E', imatuu)
        if (idec .eq. 0) then
            do i = 1, nvec
                zr(imatuu+i-1) = matv(i)
            end do
        else
            do k = 1, nno
                do n1 = 1, 2
                    i = 2*k+n1-2
                    if (k .le. nnos) then
                        i2 = i+idec*(k-1)
                    else
                        i2 = i+idec*nnos
                    end if
                    do l = 1, nno
                        do n2 = 1, 2
                            j = 2*l+n2-2
                            if (j .gt. i) goto 105
                            if (l .le. nnos) then
                                j2 = j+idec*(l-1)
                            else
                                j2 = j+idec*nnos
                            end if
                            zr(imatuu+i2*(i2-1)/2+j2-1) = matv(i*(i-1)/2+j)
                        end do
                    end do
105                 continue
                end do
            end do
        end if
!
    else if (option .eq. 'M_GAMMA') then
        call jevech('PACCELR', 'L', iacce)
        call jevech('PVECTUR', 'E', ivect)
        call vecma(matv, nvec, matp, nddl)
        if (idec .eq. 0) then
            call pmavec('ZERO', nddl, matp, zr(iacce), zr(ivect))
        else
            do k = 1, nddl
                vect1(k) = 0.0d0
                vect2(k) = 0.0d0
            end do
            do k = 1, nno
                do n1 = 1, 2
                    i = 2*k+n1-2
                    if (k .le. nnos) then
                        i2 = i+idec*(k-1)
                    else
                        i2 = i+idec*nnos
                    end if
                    vect1(i) = zr(iacce+i2-1)
                end do
            end do
            call pmavec('ZERO', nddl, matp, vect1, vect2)
            do k = 1, nno
                do n1 = 1, 2
                    i = 2*k+n1-2
                    if (k .le. nnos) then
                        i2 = i+idec*(k-1)
                    else
                        i2 = i+idec*nnos
                    end if
                    zr(ivect+i2-1) = vect2(i)
                end do
            end do
        end if
!
! OPTION ECIN_ELEM : CALCUL DE L'ENERGIE CINETIQUE
!
    else if (option .eq. 'ECIN_ELEM') then
        stopz = 'ONO'
        call tecach(stopz, 'PVITESR', 'L', iret, iad=ivite)
! IRET NE PEUT VALOIR QUE 0 (TOUT EST OK) OU 2 (CHAMP NON FOURNI)
        if (iret .eq. 0) then
            call jevech('PENERCR', 'E', iecin)
            call vecma(matv, nvec, matp, nddl)
            call pmavec('ZERO', nddl, matp, zr(ivite), masvit)
            b_n = to_blas_int(nddl)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            zr(iecin) = .5d0*ddot(b_n, zr(ivite), b_incx, masvit, b_incy)
        else
            call tecach(stopz, 'PDEPLAR', 'L', iret, iad=idepl)
            if (iret .eq. 0) then
                call jevech('PENERCR', 'E', iecin)
                call jevech('POMEGA2', 'L', ifreq)
                call vecma(matv, nvec, matp, nddl)
                call pmavec('ZERO', nddl, matp, zr(idepl), masdep)
                b_n = to_blas_int(nddl)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                zr(iecin) = .5d0*ddot(b_n, zr(idepl), b_incx, masdep, b_incy)*zr(ifreq)
            else
                call utmess('F', 'ELEMENTS2_1', sk=option)
            end if
        end if
!
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
