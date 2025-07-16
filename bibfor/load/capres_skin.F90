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
subroutine capres_skin(load, mesh, model, geomDime, valeType, nbOccPresRep)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/alcart.h"
#include "asterfort/assert.h"
#include "asterfort/char_affe_neum.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nocart.h"
#include "asterfort/utmess.h"
#include "asterfort/xtmafi.h"
#include "asterfort/xvelfm.h"
!
    character(len=8), intent(in) :: load, mesh, model
    integer(kind=8), intent(in) :: geomDime
    character(len=4), intent(in) :: valeType
    integer(kind=8), intent(in) :: nbOccPresRep
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Treatment of loads PRES_REP on skin elements
!
! --------------------------------------------------------------------------------------------------
!
! In  load             : load
! In  model            : model
! In  mesh             : mesh
! In  geomDime         : dimension of space
! In  valeType         : affected value type (real, complex or function)
! In  nbOccPresRep     : number of occurrences for PRES_REP
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: keywordFact = 'PRES_REP'
    integer(kind=8), parameter :: nbCmp = 2
    character(len=8), parameter :: cmpName(nbCmp) = (/'PRES', 'CISA'/)
    character(len=8), pointer :: mapCmpName(:) => null()
    real(kind=8), pointer :: mapCmpValeR(:) => null()
    character(len=8), pointer :: mapCmpValeK(:) => null()
    integer(kind=8) :: iocc, nbCisa, nbPres, nbma, jma, nbCrack
    integer(kind=8), parameter :: nfismx = 100
    character(len=8) :: fiss(nfismx)
    character(len=19) :: carte
    character(len=24), parameter :: mesmai = '&&CAPRES.MES_MAILLES'
    character(len=24), parameter :: lismai = '&&CAPRES.NUM_MAILLES'
    integer(kind=8), parameter :: mapListNb = 1
    character(len=19) :: mapListName(mapListNb)
    integer(kind=8) :: mapListNbCmp(mapListNb)
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Name of CARTE
    carte = load//'.CHME.PRESS'

! - Pre-allocation of CARTE
    if (valeType .eq. 'REEL') then
        call alcart('G', carte, mesh, 'PRES_R')
    else if (valeType .eq. 'FONC') then
        call alcart('G', carte, mesh, 'PRES_F')
    else
        ASSERT(ASTER_FALSE)
    end if

! - Set name of components
    call jeveuo(carte//'.NCMP', 'E', vk8=mapCmpName)
    mapCmpName(1) = cmpName(1)
    mapCmpName(2) = cmpName(2)

! - STOCKAGE DE FORCES NULLES SUR TOUT LE MAILLAGE
    if (valeType .eq. 'REEL') then
        call jeveuo(carte//'.VALV', 'E', vr=mapCmpValeR)
        mapCmpValeR(1) = 0.d0
        mapCmpValeR(2) = 0.d0
    else
        call jeveuo(carte//'.VALV', 'E', vk8=mapCmpValeK)
        mapCmpValeK(1) = '&FOZERO'
        mapCmpValeK(2) = '&FOZERO'
    end if
    call nocart(carte, 1, nbCmp)

! - Set values in CARTE
    do iocc = 1, nbOccPresRep
        if (valeType .eq. 'REEL') then
            call getvr8(keywordFact, 'PRES', iocc=iocc, scal=mapCmpValeR(1), nbret=nbPres)
            call getvr8(keywordFact, 'CISA_2D', iocc=iocc, scal=mapCmpValeR(2), nbret=nbCisa)
        else
            call getvid(keywordFact, 'PRES', iocc=iocc, scal=mapCmpValeK(1), nbret=nbPres)
            call getvid(keywordFact, 'CISA_2D', iocc=iocc, scal=mapCmpValeK(2), nbret=nbCisa)
        end if
        if (nbCisa .ne. 0 .and. geomDime .eq. 3) then
            call utmess('F', 'CHARGES6_94')
        end if
        call getvid(keywordFact, 'FISSURE', iocc=iocc, nbval=0, nbret=nbCrack)
        if (nbCrack .ne. 0) then
            nbCrack = -nbCrack
            call getvid(keywordFact, 'FISSURE', iocc=iocc, nbval=nbCrack, vect=fiss)
!           VERIFICATION DE LA COHERENCE ENTRE LES FISSURES ET LE MODELE
            call xvelfm(nbCrack, fiss, model)
!           RECUPERATION DES MAILLES PRINCIPALES X-FEM FISSUREES
            call xtmafi(geomDime, fiss, nbCrack, lismai, &
                        mesmai, nbma, model=model)
            call jeveuo(mesmai, 'L', jma)
            call nocart(carte, 3, nbCmp, mode='NOM', nma=nbma, &
                        limano=zk8(jma))
            call jedetr(mesmai)
            call jedetr(lismai)
        else
            mapListName(1) = carte
            mapListNbCmp(1) = nbCmp
            call char_affe_neum(model, mesh, geomDime, &
                                keywordFact, iocc, &
                                mapListNb, mapListName, mapListNbCmp)
        end if
    end do
!
    call jedema()
end subroutine
