! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
subroutine pipeei(ndim, axi, nno1, nno2, npg, &
                  wref, vff1, vff2, dffr2, geom, &
                  ang, mat, compor, lgpg, ddlm, &
                  ddld, ddl0, ddl1, dtau, vim, &
                  iu, im, copilo)
!
!
!
! aslint: disable=W1504
    implicit none
#include "asterf_types.h"
#include "asterc/r8vide.h"
#include "asterfort/eicine.h"
#include "asterfort/pipeex.h"
#include "asterfort/pipeou.h"
#include "asterfort/pipetc.h"
#include "asterfort/r8inir.h"
#include "asterfort/utmess.h"
    aster_logical :: axi
    integer :: ndim, nno1, nno2, npg, mat, lgpg, iu(3, 18), im(3, 9)
    real(kind=8) :: vff1(nno1, npg), vff2(nno2, npg), geom(ndim, nno2)
    real(kind=8) :: wref(npg)
    real(kind=8) :: ddlm(2*nno1*ndim+nno2*ndim), ddld(2*nno1*ndim+nno2*ndim)
    real(kind=8) :: ddl0(2*nno1*ndim+nno2*ndim), ddl1(2*nno1*ndim+nno2*ndim)
    real(kind=8) :: vim(lgpg, npg), dffr2(ndim-1, nno2, npg), ang(*)
    real(kind=8) :: dtau, copilo(5, npg)
    character(len=16) :: compor
!
!-----------------------------------------------------------------------
!
!  PILOTAGE PRED_ELAS POUR LES ELEMENTS D'INTERFACE
!
!-----------------------------------------------------------------------
    integer :: g, n, i, j, kk
    real(kind=8) :: mup(3), sup(3), mud(3), sud(3), wg, b(3, 3, 18)
!-----------------------------------------------------------------------
!
!
    call r8inir(3, 0.d0, sup, 1)
    call r8inir(3, 0.d0, sud, 1)
    call r8inir(3, 0.d0, mup, 1)
    call r8inir(3, 0.d0, mud, 1)
    call r8inir(5*npg, 0.d0, copilo, 1)
!
    do g = 1, npg
!
! -- INITIALISATION DES ELEMENTS CINEMATIQUES
!
        call eicine(ndim, axi, nno1, nno2, vff1(1, g), &
                    vff2(1, g), wref(g), dffr2(1, 1, g), geom, ang, &
                    wg, b)
!
        do i = 1, ndim
            sup(i) = 0.d0
            sud(i) = 0.d0
            do j = 1, ndim
                do n = 1, 2*nno1
                    kk = iu(j, n)
                    sup(i) = sup(i)+b(i, j, n)*(ddlm(kk)+ddld(kk)+ddl0(kk))
                    sud(i) = sud(i)+b(i, j, n)*ddl1(kk)
                end do
            end do
        end do
!
        do i = 1, ndim
            mup(i) = 0.d0
            mud(i) = 0.d0
            do n = 1, nno2
                kk = im(i, n)
                mup(i) = mup(i)+vff2(n, g)*(ddlm(kk)+ddld(kk)+ddl0(kk))
                mud(i) = mud(i)+vff2(n, g)*ddl1(kk)
            end do
        end do
!
!
! -- APPEL DU PILOTAGE PRED_ELAS SPECIFIQUE A LA LOI DE COMPORTEMENT
!
        copilo(5, g) = r8vide()
        if (compor .eq. 'CZM_TAC_MIX') then
            call pipetc(mat, sup, sud, mup, mud, &
                        vim(1, g), dtau, copilo(1, g))
        else if (compor .eq. 'CZM_OUV_MIX') then
            call pipeou(mat, sup, sud, mup, mud, &
                        vim(1, g), dtau, copilo(1, g))
        else if (compor .eq. 'CZM_EXP_MIX') then
            call pipeex(mat, sup, sud, mup, mud, &
                        vim(1, g), dtau, copilo(1, g))
        else
            call utmess('F', 'MECANONLINE_59')
        end if
!
    end do
!
end subroutine
