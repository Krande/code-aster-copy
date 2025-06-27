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
subroutine casour(load, mesh, model, geomDime, valeType)
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
!
    character(len=8), intent(in) :: load, mesh, model
    character(len=4), intent(in) :: valeType
    integer(kind=8), intent(in) :: geomDime
!
!
! BUT : STOCKAGE DES SOURCES DANS UNE CARTE ALLOUEE SUR LE
!       LIGREL DU MODELE
!
!-----------------------------------------------------------------------
    integer(kind=8) :: nsour, jvalv, n1, ncmp, iocc
    character(len=16) :: motclf
    character(len=19) :: carte
    character(len=19) :: cartes(1)
    integer(kind=8) :: ncmps(1)
    character(len=8), pointer :: vncmp(:) => null()
!     ------------------------------------------------------------------
    call jemarq()
!
    motclf = 'SOURCE'
    call getfac(motclf, nsour)
!
    carte = load//'.CHTH.SOURE'
!
    if (valeType .eq. 'REEL') then
        call alcart('G', carte, mesh, 'SOUR_R')
    else if (valeType .eq. 'FONC') then
        call alcart('G', carte, mesh, 'SOUR_F')
    else
        ASSERT(.false.)
    end if
!
    call jeveuo(carte//'.NCMP', 'E', vk8=vncmp)
    call jeveuo(carte//'.VALV', 'E', jvalv)
!
! --- STOCKAGE DE SOURCES NULLES SUR TOUT LE MAILLAGE
!
    ncmp = 1
    vncmp(1) = 'SOUR'
    if (valeType .eq. 'REEL') then
        zr(jvalv) = 0.d0
    else
        zk8(jvalv) = '&FOZERO'
    end if
    call nocart(carte, 1, ncmp)
!
! --- STOCKAGE DANS LA CARTE
!
    do iocc = 1, nsour
!
        if (valeType .eq. 'REEL') then
            call getvr8(motclf, 'SOUR', iocc=iocc, scal=zr(jvalv), nbret=n1)
        else
            call getvid(motclf, 'SOUR', iocc=iocc, scal=zk8(jvalv), nbret=n1)
        end if
!
        cartes(1) = carte
        ncmps(1) = ncmp
        call char_affe_neum(model, mesh, geomDime, motclf, iocc, 1, &
                            cartes, ncmps)
!
    end do
!
    call jedema()
end subroutine
