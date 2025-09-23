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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine charth(load, valeType)
!
    implicit none
!
#include "asterfort/adalig.h"
#include "asterfort/assert.h"
#include "asterfort/caddli.h"
#include "asterfort/caechp.h"
#include "asterfort/cagene.h"
#include "asterfort/cagrou.h"
#include "asterfort/caliag.h"
#include "asterfort/caliai.h"
#include "asterfort/calich.h"
#include "asterfort/calirc.h"
#include "asterfort/cbconv.h"
#include "asterfort/cbecha.h"
#include "asterfort/cbflnl.h"
#include "asterfort/cbflux.h"
#include "asterfort/cbgrai.h"
#include "asterfort/cbprca.h"
#include "asterfort/cbrayo.h"
#include "asterfort/cbsonl.h"
#include "asterfort/cbsour.h"
#include "asterfort/cormgi.h"
#include "asterfort/dismoi.h"
#include "asterfort/initel.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lisnnl.h"
#include "asterfort/utmess.h"
#include "asterfort/verif_affe.h"
!
    character(len=8), intent(in) :: load
    character(len=4), intent(in) :: valeType
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Treatment of loads for AFFE_CHAR_THER_*
!
! --------------------------------------------------------------------------------------------------
!
! In  load             : load
! In  valeType         : affected value type (real, complex or function)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: phenom = 'THERMIQUE'
    character(len=4), parameter :: phenomS = 'THER'
    character(len=16), parameter :: command = 'AFFE_CHAR_THER'
    character(len=16), parameter :: keywFactEnforceDOF = 'TEMP_IMPO'
    character(len=16), parameter :: keywFactEnforceSECH = 'SECH_IMPO'
    integer(kind=8) :: geomDime, iret
    character(len=8) :: mesh, model
    character(len=13) :: loadDescBase
    character(len=19) :: loadLigrel
    character(len=8), pointer :: loadLigrelLgrf(:) => null()
!
! --------------------------------------------------------------------------------------------------
!

! - Mesh, Ligrel for model, dimension of model
    call cagene(load, command, model, mesh, geomDime)
    if (geomDime .gt. 3) then
        call utmess('A', 'CHARGES2_4')
    end if

! - Get Ligrel for load
    call lisnnl(phenom, load, loadDescBase)
    loadLigrel = loadDescBase//'.LIGRE'

    if (valeType .eq. 'REEL') then

! ----- SOURCE
        call cbsour(load, mesh, model, geomDime, valeType)

! ----- CONVECTION
        call cbconv(load)

! ----- FLUX_REP
        call cbflux(load, mesh, model, geomDime, valeType)

! ----- RAYONNEMENT
        call cbrayo(load, mesh, model, valeType)

! ----- ECHANGE
        call cbecha(load, mesh, model, geomDime, valeType)

! ----- ECHANGE_PAROI
        call caechp(load, loadLigrel, mesh, model, geomDime, valeType)

! ----- EVOL_CHAR
        call cbprca(phenom, load)

! ----- GRADIENT INITIAL
        call cbgrai(load, mesh, model, valeType)

! ----- TEMP_IMPO
        call caddli(keywFactEnforceDOF, load, mesh, model, valeType)

! ----- SECH_IMPO
        call caddli(keywFactEnforceSECH, load, mesh, model, valeType)

! ----- LIAISON_DDL
        call caliai(valeType, load, phenomS)

! ----- LIAISON_GROUP
        call caliag(valeType, load, phenomS)

! ----- LIAISON_UNIF
        call cagrou(load, mesh, valeType, phenomS)

! ----- LIAISON_CHAMNO
        call calich(load, phenomS)

! ----- LIAISON_MAIL
        call calirc(phenom, load, model)

    else if (valeType .eq. 'FONC') then

! ----- SOURCE
        call cbsour(load, mesh, model, geomDime, valeType)

! ----- SOUR_NL
        call cbsonl(load, mesh, model, geomDime)

! ----- CONVECTION
        call cbconv(load)

! ----- FLUX_REP
        call cbflux(load, mesh, model, geomDime, valeType)

! ----- FLUX_NL
        call cbflnl(load, mesh, model)

! ----- RAYONNEMENT
        call cbrayo(load, mesh, model, valeType)

! ----- ECHANGE
        call cbecha(load, mesh, model, geomDime, valeType)

! ----- ECHANGE_PAROI
        call caechp(load, loadLigrel, mesh, model, geomDime, valeType)

! ----- GRADIENT INITIAL
        call cbgrai(load, mesh, model, valeType)

! ----- TEMP_IMPO
        call caddli(keywFactEnforceDOF, load, mesh, model, valeType)

! ----- SECH_IMPO
        call caddli(keywFactEnforceSECH, load, mesh, model, valeType)

! ----- LIAISON_DDL
        call caliai(valeType, load, phenomS)

! ----- LIAISON_GROUP
        call caliag(valeType, load, phenomS)

! ----- LIAISON_UNIF
        call cagrou(load, mesh, valeType, phenomS)
    else
        ASSERT(ASTER_FALSE)

    end if

! - Update loads <LIGREL>
    call jeexin(loadLigrel//'.LGRF', iret)
    if (iret .ne. 0) then
        call adalig(loadLigrel)
        call cormgi('G', loadLigrel)
        call jeecra(loadLigrel//'.LGRF', 'DOCU', cval=phenomS)
        call initel(loadLigrel)
        call jeveuo(loadLigrel//'.LGRF', 'E', vk8=loadLigrelLgrf)
        call dismoi('PARTITION', model, "MODELE", repk=loadLigrelLgrf(2))
    end if

! - Audit assignments
    call verif_affe(modele=model, sd=load)
!
end subroutine
