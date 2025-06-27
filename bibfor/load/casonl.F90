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
subroutine casonl(load, mesh, model, geomDime)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/alcart.h"
#include "asterfort/char_affe_neum.h"
#include "asterfort/getvid.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nocart.h"
!
    character(len=8), intent(in) :: load, mesh, model
    integer(kind=8), intent(in) :: geomDime
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Treatment of load SOUR_NL
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : mesh
! In  load             : load
! In  model            : model
! In  geomDime         : space dimension
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: keywordFact = 'SOUR_NL'
    integer(kind=8) :: nsour, jvalv, n1, ncmp, iocc
    character(len=19) :: carte
    character(len=19) :: cartes(1)
    integer(kind=8) :: ncmps(1)
    character(len=8), pointer :: vncmp(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    ncmp = 1
    call getfac(keywordFact, nsour)
!
    carte = load//'.CHTH.SOUNL'
    call alcart('G', carte, mesh, 'SOUR_F')
    call jeveuo(carte//'.NCMP', 'E', vk8=vncmp)
    call jeveuo(carte//'.VALV', 'E', jvalv)
    vncmp(1) = 'SOUR'
!
!
! --- DEFAUT: STOCKAGE D'UNE SOURCE NULLE SUR TOUT LE MAILLAGE
!
    zk8(jvalv) = '&FOZERO'
    call nocart(carte, 1, ncmp)
!
!
! --- STOCKAGE DES FONCTIONS SOURCES DANS LA CARTE
!
    do iocc = 1, nsour
        call getvid(keywordFact, 'SOUR', iocc=iocc, scal=zk8(jvalv), nbret=n1)
!
        cartes(1) = carte
        ncmps(1) = ncmp
        call char_affe_neum(model, mesh, geomDime, keywordFact, iocc, 1, &
                            cartes, ncmps)
    end do
!
    call jedema()
end subroutine
