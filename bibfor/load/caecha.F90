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

subroutine caecha(load, model, mesh, geomDime, valeType)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/alcart.h"
#include "asterfort/assert.h"
#include "asterfort/char_affe_neum.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nocart.h"
#include "asterfort/tecart.h"
!
    character(len=8), intent(in) :: load, mesh, model
    integer(kind=8), intent(in) :: geomDime
    character(len=4), intent(in) :: valeType
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Treatment of loads ECHANGE
!
! --------------------------------------------------------------------------------------------------
!
! In  load             : load
! In  model            : model
! In  mesh             : mesh
! In  geomDime         : dimension of space
! In  valeType         : affected value type (real, complex or function)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: keywordFact = 'ECHANGE'
    integer(kind=8) :: necha, ncmp, jvalv1, jvalv2, n, ncmp1, ncmp2, ncmps(2)
    integer(kind=8) :: iocc
    real(kind=8) :: r8b
    character(len=8) :: k8b
    character(len=19) :: carte1, carte2, cartes(2)
    character(len=8), pointer :: vncmp1(:) => null()
    character(len=8), pointer :: vncmp2(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    call getfac(keywordFact, necha)
!
    carte1 = load//'.CHTH.COEFH'
    carte2 = load//'.CHTH.T_EXT'
!
    if (valeType .eq. 'REEL') then
        call alcart('G', carte1, mesh, 'COEH_R')
        call alcart('G', carte2, mesh, 'TEMP_R')
    else if (valeType .eq. 'FONC') then
        call alcart('G', carte1, mesh, 'COEH_F')
        call alcart('G', carte2, mesh, 'TEMP_F')
    else
        ASSERT(ASTER_FALSE)
    end if
!
    call jeveuo(carte1//'.NCMP', 'E', vk8=vncmp1)
    call jeveuo(carte1//'.VALV', 'E', jvalv1)
    call jeveuo(carte2//'.NCMP', 'E', vk8=vncmp2)
    call jeveuo(carte2//'.VALV', 'E', jvalv2)
!
! --- STOCKAGE DE FLUX NULS SUR TOUT LE MAILLAGE
!
    ncmp = 3
    vncmp1(1) = 'H'
    vncmp1(2) = 'H_INF'
    vncmp1(3) = 'H_SUP'
    vncmp2(1) = 'TEMP'
    vncmp2(2) = 'TEMP_INF'
    vncmp2(3) = 'TEMP_SUP'
    if (valeType .eq. 'REEL') then
        zr(jvalv1-1+1) = 0.d0
        zr(jvalv1-1+2) = 0.d0
        zr(jvalv1-1+3) = 0.d0
        zr(jvalv2-1+1) = 0.d0
        zr(jvalv2-1+2) = 0.d0
        zr(jvalv2-1+3) = 0.d0
    else
        zk8(jvalv1-1+1) = '&FOZERO'
        zk8(jvalv1-1+2) = '&FOZERO'
        zk8(jvalv1-1+3) = '&FOZERO'
        zk8(jvalv2-1+1) = '&FOZERO'
        zk8(jvalv2-1+2) = '&FOZERO'
        zk8(jvalv2-1+3) = '&FOZERO'
    end if
    call nocart(carte1, 1, ncmp)
    call nocart(carte2, 1, ncmp)
!
! --- STOCKAGE DANS LA CARTE
!
    do iocc = 1, necha
        ncmp1 = 0
        ncmp2 = 0
        if (valeType .eq. 'REEL') then
            call getvr8(keywordFact, 'COEF_H', iocc=iocc, scal=r8b, nbret=n)
            if (n .eq. 1) then
                ncmp1 = ncmp1+1
                vncmp1(ncmp1) = 'H'
                zr(jvalv1-1+ncmp1) = r8b
            end if
            call getvr8(keywordFact, 'COEF_H_INF', iocc=iocc, scal=r8b, nbret=n)
            if (n .eq. 1) then
                ncmp1 = ncmp1+1
                vncmp1(ncmp1) = 'H_INF'
                zr(jvalv1-1+ncmp1) = r8b
            end if
            call getvr8(keywordFact, 'COEF_H_SUP', iocc=iocc, scal=r8b, nbret=n)
            if (n .eq. 1) then
                ncmp1 = ncmp1+1
                vncmp1(ncmp1) = 'H_SUP'
                zr(jvalv1-1+ncmp1) = r8b
            end if
            call getvr8(keywordFact, 'TEMP_EXT', iocc=iocc, scal=r8b, nbret=n)
            if (n .eq. 1) then
                ncmp2 = ncmp2+1
                vncmp2(ncmp2) = 'TEMP'
                zr(jvalv2-1+ncmp2) = r8b
            end if
            call getvr8(keywordFact, 'TEMP_EXT_INF', iocc=iocc, scal=r8b, nbret=n)
            if (n .eq. 1) then
                ncmp2 = ncmp2+1
                vncmp2(ncmp2) = 'TEMP_INF'
                zr(jvalv2-1+ncmp2) = r8b
            end if
            call getvr8(keywordFact, 'TEMP_EXT_SUP', iocc=iocc, scal=r8b, nbret=n)
            if (n .eq. 1) then
                ncmp2 = ncmp2+1
                vncmp2(ncmp2) = 'TEMP_SUP'
                zr(jvalv2-1+ncmp2) = r8b
            end if
        else
            call getvid(keywordFact, 'COEF_H', iocc=iocc, scal=k8b, nbret=n)
            if (n .eq. 1) then
                ncmp1 = ncmp1+1
                vncmp1(ncmp1) = 'H'
                zk8(jvalv1-1+ncmp1) = k8b
            end if
            call getvid(keywordFact, 'COEF_H_INF', iocc=iocc, scal=k8b, nbret=n)
            if (n .eq. 1) then
                ncmp1 = ncmp1+1
                vncmp1(ncmp1) = 'H_INF'
                zk8(jvalv1-1+ncmp1) = k8b
            end if
            call getvid(keywordFact, 'COEF_H_SUP', iocc=iocc, scal=k8b, nbret=n)
            if (n .eq. 1) then
                ncmp1 = ncmp1+1
                vncmp1(ncmp1) = 'H_SUP'
                zk8(jvalv1-1+ncmp1) = k8b
            end if
            call getvid(keywordFact, 'TEMP_EXT', iocc=iocc, scal=k8b, nbret=n)
            if (n .eq. 1) then
                ncmp2 = ncmp2+1
                vncmp2(ncmp2) = 'TEMP'
                zk8(jvalv2-1+ncmp2) = k8b
            end if
            call getvid(keywordFact, 'TEMP_EXT_INF', iocc=iocc, scal=k8b, nbret=n)
            if (n .eq. 1) then
                ncmp2 = ncmp2+1
                vncmp2(ncmp2) = 'TEMP_INF'
                zk8(jvalv2-1+ncmp2) = k8b
            end if
            call getvid(keywordFact, 'TEMP_EXT_SUP', iocc=iocc, scal=k8b, nbret=n)
            if (n .eq. 1) then
                ncmp2 = ncmp2+1
                vncmp2(ncmp2) = 'TEMP_SUP'
                zk8(jvalv2-1+ncmp2) = k8b
            end if
        end if
!
        cartes(1) = carte1
        cartes(2) = carte2
        ncmps(1) = ncmp1
        ncmps(2) = ncmp2
        call char_affe_neum(model, mesh, geomDime, keywordFact, iocc, 2, &
                            cartes, ncmps)
!
    end do
    call tecart(carte1)
    call tecart(carte2)
!
    call jedema()
end subroutine
