! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
subroutine capres(load, ligrmo, mesh, ndim, valeType, nbOcc)
!
implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/char_affe_neum.h"
#include "asterfort/char_crea_cart.h"
#include "asterfort/dismoi.h"
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
character(len=8), intent(in) :: load, mesh
character(len=19), intent(in) :: ligrmo
integer, intent(in) :: ndim
character(len=4), intent(in) :: valeType
integer, intent(in) :: nbOcc
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Treatment of load PRES_REP
!
! --------------------------------------------------------------------------------------------------
!
! In  load             : load
! In  mesh             : mesh
! In  valeType         : affected value type (real, complex or function)
! In  ligrmo           : list of elements in model
! In  ndim             : dimension of space (2 or 3)
! In  nbOcc            : number of factor keywords
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: keywordfact = 'PRES_REP'
    integer, parameter :: nfismx = 100
    integer :: nbCell
    integer :: nbCisa, nfiss, nbret
    integer :: jvCell, jvalv, iocc
    character(len=8) :: fiss(nfismx), model
    character(len=24), parameter :: mesmai = '&&CAPRES.MES_MAILLES'
    character(len=24), parameter :: lismai = '&&CAPRES.NUM_MAILLES'
    character(len=19) :: map(LOAD_MAP_NBMAX)
    integer :: nbMap, nbCmp(LOAD_MAP_NBMAX)
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Initializations
!
    call dismoi('NOM_MODELE', ligrmo, 'LIGREL', repk=model)
    ASSERT(valeType .eq. 'REEL' .or. valeType .eq. 'FONC')
!
! - Creation and initialization to zero of <CARTE>
!
    call char_crea_cart('MECANIQUE', keywordfact, load, mesh, valeType,&
                        nbMap, map, nbCmp)
    ASSERT(nbMap .eq. 1)
    call jeveuo(map(1)//'.VALV', 'E', jvalv)
!
! - Loop on factor keyword
!
    do iocc = 1, nbOcc


! ----- Get values of load
        if (valeType .eq. 'REEL') then
            call getvr8(keywordfact, 'PRES', iocc=iocc, scal=zr(jvalv), nbret=nbret)
            call getvr8(keywordfact, 'CISA_2D', iocc=iocc, scal=zr(jvalv+1), nbret=nbCisa)
        else
            call getvid(keywordfact, 'PRES', iocc=iocc, scal=zk8(jvalv), nbret=nbret)
            call getvid(keywordfact, 'CISA_2D', iocc=iocc, scal=zk8(jvalv+1), nbret=nbCisa)
        endif
        if (nbCisa .ne. 0 .and. ndim .eq. 3) then
            call utmess('F', 'CHARGES6_94')
        endif
        call getvid(keywordfact, 'FISSURE', iocc=iocc, nbval=0, nbret=nfiss)
!
        if (nfiss .ne. 0) then
            nfiss = -nfiss
            call getvid(keywordfact, 'FISSURE', iocc=iocc, nbval=nfiss, vect=fiss, nbret=nbret)
!           VERIFICATION DE LA COHERENCE ENTRE LES FISSURES ET LE MODELE
            call xvelfm(nfiss, fiss, model)
!           RECUPERATION DES MAILLES PRINCIPALES X-FEM FISSUREES
            call xtmafi(ndim, fiss, nfiss, lismai,&
                        mesmai, nbCell, model = model)
            call jeveuo(mesmai, 'L', jvCell)
            call nocart(map(1), 3, nbCmp(1), mode='NOM', nma=nbCell,&
                        limano=zk8(jvCell))
            call jedetr(mesmai)
            call jedetr(lismai)
        else
            call char_affe_neum(model      , mesh, ndim,&
                                keywordfact, iocc,&
                                nbMap      , map , nbCmp)
        endif
    end do
!
    call jedema()
end subroutine
