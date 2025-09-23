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
subroutine caflnl(load, mesh, model)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/alcart.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nocart.h"
#include "asterfort/reliem.h"
#include "asterfort/utmess.h"
!
    character(len=8), intent(in) :: load, mesh, model
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Treatment of loads 'FLUX_NL'
!
! --------------------------------------------------------------------------------------------------
!
! In  load             : load
! In  model            : model
! In  mesh             : mesh
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: keywordFact = 'FLUX_NL'
    integer(kind=8) :: nflux, jvalv, nf, iocc, nbtou, nbma, jma, ncmp, lprol
    character(len=8) :: k8b, typmcl(2)
    character(len=16) :: motcle(2)
    character(len=19) :: carte
    character(len=24) :: mesmai, prol, nompar
    character(len=8), pointer :: vncmp(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    call getfac(keywordFact, nflux)
!
    carte = load//'.CHTH.FLUNL'
    call alcart('G', carte, mesh, 'FLUN_F')
    call jeveuo(carte//'.NCMP', 'E', vk8=vncmp)
    call jeveuo(carte//'.VALV', 'E', jvalv)
!
! --- STOCKAGE DE FLUX NULS SUR TOUT LE MAILLAGE
!
    ncmp = 3
    vncmp(1) = 'FLUN'
    vncmp(2) = 'FLUN_INF'
    vncmp(3) = 'FLUN_SUP'
    zk8(jvalv-1+1) = '&FOZERO'
    zk8(jvalv-1+2) = '&FOZERO'
    zk8(jvalv-1+3) = '&FOZERO'
    call nocart(carte, 1, ncmp)
!
    mesmai = '&&CAFLNL.MES_MAILLES'
    motcle(1) = 'GROUP_MA'
    motcle(2) = 'MAILLE'
    typmcl(1) = 'GROUP_MA'
    typmcl(2) = 'MAILLE'
!
! --- STOCKAGE DANS LES CARTES
!
    do iocc = 1, nflux
!
        call getvid(keywordFact, 'FLUN', iocc=iocc, scal=zk8(jvalv), nbret=nf)

        prol = zk8(jvalv)//'           .PROL'
        call jeveuo(prol, 'L', lprol)
        nompar = zk24(lprol+2)
        if (nompar(1:4) .ne. 'TEMP' .and. nompar(1:4) .ne. 'SECH') then
            call utmess('F', 'CHARGES2_3', sk=zk8(jvalv))
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
