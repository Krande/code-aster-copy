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
subroutine te0488(option, nomte)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/dxqpgl.h"
#include "asterfort/dxtpgl.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/subaco.h"
#include "asterfort/sumetr.h"
#include "asterfort/tecach.h"
#include "asterfort/utpvlg.h"
#include "asterfort/vectan.h"
#include "asterfort/vectgt.h"
#include "asterfort/lteatt.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: all
!
! Options: COOR_ELGA
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: jv_poids, jv_vf, jv_geom, jv_coopg, jv_dfde
    integer(kind=8) :: jv_nbsp_i, jv_cacoqu
    integer(kind=8) :: nno, kp, npg, ino, ndim
    real(kind=8) :: xx, yy, zz, poids, cova(3, 3), metr(2, 2), jac
    integer(kind=8) ::  nbc, decpo, ic, ispc, jtab(7), nbsp, iret
    real(kind=8) :: epais, excen, bas, epc, pgl(3, 3), gm2(3)
    integer(kind=8) :: lzi, lzr, nb1, nb2
    real(kind=8) :: vectg(2, 3), vectt(3, 3)
    real(kind=8) :: vecta(9, 2, 3), vectn(9, 3), vectpt(9, 2, 3), hh
    real(kind=8), parameter :: zero = 0.d0
    aster_logical :: l_coq3d, l_grille, l_solid_shell
    real(kind=8), parameter :: gm1(3) = (/0.d0, 0.d0, 1.d0/)
    real(kind=8), parameter :: poidc(3) = (/0.16666666666666666d0, 0.66666666666666663d0, &
                                            0.16666666666666666d0/)
!
! --------------------------------------------------------------------------------------------------
!
    l_coq3d = lteatt('MODELI', 'CQ3')
    l_grille = lteatt('MODELI', 'GRC')
    l_solid_shell = lteatt('MODELI', 'SSH')
    if (l_coq3d) then
        call elrefe_info(fami='MASS', ndim=ndim, nno=nno, npg=npg, &
                         jpoids=jv_poids, jvf=jv_vf, jdfde=jv_dfde)
    else
        call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, npg=npg, &
                         jpoids=jv_poids, jvf=jv_vf, jdfde=jv_dfde)
    end if
!
! - Access to input fields
!
    call jevech('PGEOMER', 'L', jv_geom)
!
! - Access to output fields
!
    call tecach('OOO', 'PCOORPG', 'E', iret, nval=7, itab=jtab)
    jv_coopg = jtab(1)
    nbsp = jtab(7)
!
! - Get parameters for structural elements
!
    if (l_grille) then
        nbsp = 1
        call jevech('PCACOQU', 'L', jv_cacoqu)
        excen = zr(jv_cacoqu+3)
        if (nno .eq. 3) then
            call dxtpgl(zr(jv_geom), pgl)
        else if (nno .eq. 4) then
            call dxqpgl(zr(jv_geom), pgl)
        end if
        call utpvlg(1, 3, pgl, gm1, gm2)
    end if
    if (nbsp .ne. 1) then
        call jevech('PNBSP_I', 'L', jv_nbsp_i)
        nbc = zi(jv_nbsp_i)
        call jevech('PCACOQU', 'L', jv_cacoqu)
        epais = zr(jv_cacoqu)
        if (l_coq3d) then
            excen = zr(jv_cacoqu+5)
            call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
            nb1 = zi(lzi-1+1)
            nb2 = zi(lzi-1+2)
            npg = zi(lzi-1+4)
            call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)
            call vectan(nb1, nb2, zr(jv_geom), zr(lzr), vecta, vectn, vectpt)
        else
            excen = zr(jv_cacoqu+4)
            if (nno .eq. 3) then
                call dxtpgl(zr(jv_geom), pgl)
            else if (nno .eq. 4) then
                call dxqpgl(zr(jv_geom), pgl)
            end if
            call utpvlg(1, 3, pgl, gm1, gm2)
        end if
        bas = -epais/2.d0+excen
        epc = epais/nbc
    end if
!
    do kp = 1, npg
        xx = zero
        yy = zero
        zz = zero
        if (l_solid_shell) then
            do ino = 1, nno-1
                xx = xx+zr(jv_geom+3*(ino-1)+0)*zr(jv_vf+(kp-1)*nno+ino-1)
                yy = yy+zr(jv_geom+3*(ino-1)+1)*zr(jv_vf+(kp-1)*nno+ino-1)
                zz = zz+zr(jv_geom+3*(ino-1)+2)*zr(jv_vf+(kp-1)*nno+ino-1)
            end do
        else
            do ino = 1, nno
                xx = xx+zr(jv_geom+3*(ino-1)+0)*zr(jv_vf+(kp-1)*nno+ino-1)
                yy = yy+zr(jv_geom+3*(ino-1)+1)*zr(jv_vf+(kp-1)*nno+ino-1)
                zz = zz+zr(jv_geom+3*(ino-1)+2)*zr(jv_vf+(kp-1)*nno+ino-1)
            end do
        end if
        if (ndim .eq. 3) then
            call dfdm3d(nno, kp, jv_poids, jv_dfde, zr(jv_geom), poids)
        else if (ndim .eq. 2) then
            call subaco(nno, zr(jv_dfde+(kp-1)*ndim*nno), zr(jv_geom), cova)
            call sumetr(cova, metr, jac)
            poids = jac*zr(jv_poids-1+kp)
        else
            ASSERT(ASTER_FALSE)
        end if
!
        if (nbsp .ne. 1) then
            decpo = 4*3*nbc*(kp-1)
            if (l_coq3d) then
                call vectgt(1, nb1, zr(jv_geom), zero, kp, &
                            zr(lzr), epais, vectn, vectg, vectt)
                gm2(1) = vectt(3, 1)
                gm2(2) = vectt(3, 2)
                gm2(3) = vectt(3, 3)
            end if
            do ic = 1, nbc
                do ispc = 1, 3
                    hh = bas+dble(ic-1)*epc+dble(ispc-1)*epc/2.d0
                    zr(jv_coopg+decpo+(ic-1)*12+(ispc-1)*4+0) = xx+hh*gm2(1)
                    zr(jv_coopg+decpo+(ic-1)*12+(ispc-1)*4+1) = yy+hh*gm2(2)
                    zr(jv_coopg+decpo+(ic-1)*12+(ispc-1)*4+2) = zz+hh*gm2(3)
                    zr(jv_coopg+decpo+(ic-1)*12+(ispc-1)*4+3) = poids*epc*poidc(ispc)
                end do
            end do
        else if (l_grille) then
            zr(jv_coopg+4*(kp-1)+0) = xx+excen*gm2(1)
            zr(jv_coopg+4*(kp-1)+1) = yy+excen*gm2(2)
            zr(jv_coopg+4*(kp-1)+2) = zz+excen*gm2(3)
            zr(jv_coopg+4*(kp-1)+3) = poids
        else
            zr(jv_coopg+4*(kp-1)+0) = xx
            zr(jv_coopg+4*(kp-1)+1) = yy
            zr(jv_coopg+4*(kp-1)+2) = zz
            zr(jv_coopg+4*(kp-1)+3) = poids
        end if
    end do
!
end subroutine
