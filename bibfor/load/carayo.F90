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
subroutine carayo(load, mesh, model, valeType)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/alcart.h"
#include "asterfort/assert.h"
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
! Treatment of load RAYONNEMENT
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : mesh
! In  load             : load
! In  model            : model
! In  valeType         : affected value type (real, complex or function)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: keywordFact = 'RAYONNEMENT'
    integer(kind=8) :: nrayo, jvalv, ncmp, n, iocc, nbtou, nbma, jma
    character(len=8) :: k8b, typmcl(2)
    character(len=16) :: motcle(2)
    character(len=19) :: carte
    character(len=24) :: mesmai
    character(len=8), pointer :: vncmp(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    call getfac(keywordFact, nrayo)
!
    carte = load//'.CHTH.RAYO'
!
    if (valeType .eq. 'REEL') then
        call alcart('G', carte, mesh, 'RAYO_R')
    else if (valeType .eq. 'FONC') then
        call alcart('G', carte, mesh, 'RAYO_F')
    else
        ASSERT(ASTER_FALSE)
    end if
!
    call jeveuo(carte//'.NCMP', 'E', vk8=vncmp)
    call jeveuo(carte//'.VALV', 'E', jvalv)
!
! --- STOCKAGE DE FLUX NULS SUR TOUT LE MAILLAGE
!
    ncmp = 3
    vncmp(1) = 'SIGMA'
    vncmp(2) = 'EPSIL'
    vncmp(3) = 'TPINF'
    if (valeType .eq. 'REEL') then
        zr(jvalv-1+1) = 0.d0
        zr(jvalv-1+2) = 0.d0
        zr(jvalv-1+3) = 0.d0
    else
        zk8(jvalv-1+1) = '&FOZERO'
        zk8(jvalv-1+2) = '&FOZERO'
        zk8(jvalv-1+3) = '&FOZERO'
    end if
    call nocart(carte, 1, ncmp)
!
    mesmai = '&&CARAYO.MES_MAILLES'
    motcle(1) = 'GROUP_MA'
    motcle(2) = 'MAILLE'
    typmcl(1) = 'GROUP_MA'
    typmcl(2) = 'MAILLE'
!
! --- STOCKAGE DANS LA CARTE
!
    do iocc = 1, nrayo
        if (valeType .eq. 'REEL') then
            call getvr8(keywordFact, 'SIGMA', iocc=iocc, scal=zr(jvalv), nbret=n)
            call getvr8(keywordFact, 'EPSILON', iocc=iocc, scal=zr(jvalv+1), nbret=n)
            call getvr8(keywordFact, 'TEMP_EXT', iocc=iocc, scal=zr(jvalv+2), nbret=n)
        else
            call getvid(keywordFact, 'SIGMA', iocc=iocc, scal=zk8(jvalv), nbret=n)
            call getvid(keywordFact, 'EPSILON', iocc=iocc, scal=zk8(jvalv+1), nbret=n)
            call getvid(keywordFact, 'TEMP_EXT', iocc=iocc, scal=zk8(jvalv+2), nbret=n)
        end if
!
        call getvtx(keywordFact, 'TOUT', iocc=iocc, scal=k8b, nbret=nbtou)
        if (nbtou .ne. 0) then
            call nocart(carte, 1, ncmp)
!
        else
            call reliem(model, mesh, 'NU_MAILLE', keywordFact, iocc, &
                        2, motcle, typmcl, mesmai, nbma)
            if (nbma .eq. 0) cycle
            call jeveuo(mesmai, 'L', jma)
            call nocart(carte, 3, ncmp, mode='NUM', nma=nbma, &
                        limanu=zi(jma))
            call jedetr(mesmai)
        end if
!
    end do
!
    call jedema()
!
end subroutine
