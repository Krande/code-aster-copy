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
! person_in_charge: sylvie.granet at edf.fr
!
subroutine cafthm(load, mesh, model, valeType)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/alcart.h"
#include "asterfort/assert.h"
#include "asterfort/exixfe.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
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
    character(len=4), intent(in) :: valeType
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Treatment of load FLUX_THM_REP
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
    character(len=16), parameter :: keywordFact = 'FLUX_THM_REP'
    integer(kind=8) :: n1, n2, n3, n4, nflux, jvalv, iocc
    integer(kind=8) :: nbtou, nbma, jma, ncmp
    integer(kind=8) :: iret, nfiss, jnfis
    character(len=8) :: k8b, typmcl(2)
    character(len=16) :: motcle(2)
    character(len=19) :: carte
    character(len=24) :: mesmai, lismai
    character(len=8), pointer :: vncmp(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    call getfac(keywordFact, nflux)
    if (nflux .eq. 0) goto 99
!
    carte = load//'.CHME.FLUX '
!
    call exixfe(model, iret)
!
    if (valeType .eq. 'REEL') then
        call alcart('G', carte, mesh, 'FTHM_R')
    else if (valeType .eq. 'FONC') then
        call alcart('G', carte, mesh, 'FTHM_F')
    else
        ASSERT(ASTER_FALSE)
    end if
!
    call jeveuo(carte//'.NCMP', 'E', vk8=vncmp)
    call jeveuo(carte//'.VALV', 'E', jvalv)
!
! --- STOCKAGE DE FORCES NULLES SUR TOUT LE MAILLAGE
!
    ncmp = 3
    vncmp(1) = 'PFLU1'
    vncmp(2) = 'PFLU2'
    vncmp(3) = 'PTHER'
!
    if (valeType .eq. 'FONC') then
        ncmp = 4
        vncmp(4) = 'PFLUF'
    end if
!
    if (valeType .eq. 'REEL') then
        zr(jvalv) = 0.d0
        zr(jvalv+1) = 0.d0
        zr(jvalv+2) = 0.d0
    else
        zk8(jvalv) = '&FOZERO'
        zk8(jvalv+1) = '&FOZERO'
        zk8(jvalv+2) = '&FOZERO'
        zk8(jvalv+3) = '&FOZERO'
    end if
    call nocart(carte, 1, ncmp)
!
    mesmai = '&&CAFTHM.MAILLES_INTE'
    lismai = '&&CAFTHM.NUM_MAILLES'
    motcle(1) = 'GROUP_MA'
    motcle(2) = 'MAILLE'
    typmcl(1) = 'GROUP_MA'
    typmcl(2) = 'MAILLE'
!
! --- STOCKAGE DANS LA CARTE
!
    do iocc = 1, nflux
        if (valeType .eq. 'REEL') then
            call getvr8(keywordFact, 'FLUN_HYDR1', iocc=iocc, scal=zr(jvalv), nbret=n1)
            call getvr8(keywordFact, 'FLUN_HYDR2', iocc=iocc, scal=zr(jvalv+1), nbret=n2)
            call getvr8(keywordFact, 'FLUN', iocc=iocc, scal=zr(jvalv+2), nbret=n3)
        else
            call getvid(keywordFact, 'FLUN_HYDR1', iocc=iocc, scal=zk8(jvalv), nbret=n1)
            call getvid(keywordFact, 'FLUN_HYDR2', iocc=iocc, scal=zk8(jvalv+1), nbret=n2)
            call getvid(keywordFact, 'FLUN', iocc=iocc, scal=zk8(jvalv+2), nbret=n3)
            call getvid(keywordFact, 'FLUN_FRAC', iocc=iocc, scal=zk8(jvalv+3), nbret=n4)
        end if
!
! --- TEST SUR LES CAL
!
!
        call getvtx(keywordFact, 'TOUT', iocc=iocc, scal=k8b, nbret=nbtou)
!
        nfiss = 0
        if (iret .ne. 0) then
            call jeveuo(model//'.NFIS', 'L', jnfis)
            nfiss = zi(jnfis)
        end if
!
        if (nbtou .ne. 0) then
!
            call nocart(carte, 1, ncmp)
        else
            if (nfiss .ne. 0) then
!           LES FLUX POUR LA DEUXIEME PRESSION PRE2 ET POUR LA THERMIQUE
!           NE SONT PAS AUTORISES EN HM-XFEM
                if ((n2 .ne. 0) .and. (n3 .ne. 0)) call utmess('F', 'XFEM_48')
            end if
!           LES FLUX DE FLUIDE DANS LES FRACTURES NE SONT AUTORISES QU'EN HM_XFEM
            if ((valeType .eq. 'FONC') .and. (nfiss .eq. 0) .and. (n4 .ne. 0)) then
                call utmess('F', 'XFEM_28')
            end if
            call reliem(model, mesh, 'NU_MAILLE', keywordFact, iocc, &
                        2, motcle, typmcl, mesmai, nbma)
            if (nbma .ne. 0) then
                call jeveuo(mesmai, 'L', jma)
                call nocart(carte, 3, ncmp, mode='NUM', nma=nbma, &
                            limanu=zi(jma))
                call jedetr(mesmai)
            end if
        end if
!
    end do
99  continue
!
!
    call jedema()
end subroutine
