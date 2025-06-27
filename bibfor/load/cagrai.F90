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
subroutine cagrai(load, mesh, model, valeType)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/assert.h"
#include "asterfort/alcart.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nocart.h"
#include "asterfort/reliem.h"
!
    character(len=8), intent(in) :: load, mesh, model
    character(len=4), intent(in) :: valeType
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Treatment of load PRE_GRAD_TEMP
!
! --------------------------------------------------------------------------------------------------
!
! In  load             : load
! In  mesh             : mesh
! In  model            : model
! In  valeType         : affected value type (real, complex or function)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: keywordFact = 'PRE_GRAD_TEMP'
    integer(kind=8) :: nchgi, jvalv, nx, ny, nz, i, iocc, nbtou, nbma, jma
    real(kind=8) :: grx, gry, grz
    character(len=8) :: k8b, typmcl(2), grxf, gryf, grzf
    character(len=16) :: motcle(2)
    character(len=19) :: carte
    character(len=24) :: mesmai
    character(len=8), pointer :: ncmp(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    call getfac(keywordFact, nchgi)
!
    carte = load//'.CHTH.'//'GRAIN'
!
    if (valeType .eq. 'REEL') then
        call alcart('G', carte, mesh, 'FLUX_R')
    else if (valeType .eq. 'FONC') then
        call alcart('G', carte, mesh, 'FLUX_F')
    else
        ASSERT(ASTER_FALSE)
    end if
!
    call jeveuo(carte//'.NCMP', 'E', vk8=ncmp)
    call jeveuo(carte//'.VALV', 'E', jvalv)
!
    ncmp(1) = 'FLUX'
    ncmp(2) = 'FLUY'
    ncmp(3) = 'FLUZ'
    if (valeType .eq. 'REEL') then
        do i = 1, 3
            zr(jvalv-1+i) = 0.d0
        end do
    else if (valeType .eq. 'FONC') then
        do i = 1, 3
            zk8(jvalv-1+i) = '&FOZERO'
        end do
    end if
    call nocart(carte, 1, 3)
!
    mesmai = '&&CAGRAI.MES_MAILLES'
    motcle(1) = 'GROUP_MA'
    motcle(2) = 'MAILLE'
    typmcl(1) = 'GROUP_MA'
    typmcl(2) = 'MAILLE'
!
    do iocc = 1, nchgi
!
        if (valeType .eq. 'REEL') then
            call getvr8(keywordFact, 'FLUX_X', iocc=iocc, scal=grx, nbret=nx)
            call getvr8(keywordFact, 'FLUX_Y', iocc=iocc, scal=gry, nbret=ny)
            call getvr8(keywordFact, 'FLUX_Z', iocc=iocc, scal=grz, nbret=nz)
            do i = 1, 3
                zr(jvalv-1+i) = 0.d0
            end do
            if (nx .ne. 0) zr(jvalv-1+1) = grx
            if (ny .ne. 0) zr(jvalv-1+2) = gry
            if (nz .ne. 0) zr(jvalv-1+3) = grz
        else if (valeType .eq. 'FONC') then
            call getvid(keywordFact, 'FLUX_X', iocc=iocc, scal=grxf, nbret=nx)
            call getvid(keywordFact, 'FLUX_Y', iocc=iocc, scal=gryf, nbret=ny)
            call getvid(keywordFact, 'FLUX_Z', iocc=iocc, scal=grzf, nbret=nz)
            do i = 1, 3
                zk8(jvalv-1+i) = '&FOZERO'
            end do
            if (nx .ne. 0) zk8(jvalv-1+1) = grxf
            if (ny .ne. 0) zk8(jvalv-1+2) = gryf
            if (nz .ne. 0) zk8(jvalv-1+3) = grzf
        end if
!
        call getvtx(keywordFact, 'TOUT', iocc=iocc, scal=k8b, nbret=nbtou)
        if (nbtou .ne. 0) then
            call nocart(carte, 1, 3)
!
        else
            call reliem(model, mesh, 'NU_MAILLE', keywordFact, iocc, &
                        2, motcle, typmcl, mesmai, nbma)
            if (nbma .eq. 0) cycle
            call jeveuo(mesmai, 'L', jma)
            call nocart(carte, 3, 3, mode='NUM', nma=nbma, limanu=zi(jma))
            call jedetr(mesmai)
        end if
!
    end do
!
    call jedema()
end subroutine
