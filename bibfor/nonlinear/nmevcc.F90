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
subroutine nmevcc(sddisc, nume_inst, ds_contact, i_echec, i_echec_acti)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cfdisi.h"
#include "asterfort/cfmmvd.h"
#include "asterfort/diinst.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utdidt.h"
!
    character(len=19), intent(in) :: sddisc
    integer(kind=8), intent(in) :: nume_inst
    type(NL_DS_Contact), intent(in) :: ds_contact
    integer(kind=8), intent(in) :: i_echec
    integer(kind=8), intent(out) :: i_echec_acti
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME - EVENEMENTS)
!
! DETECTION DE L'EVENEMENT COLLISION - CAS CONTINU
!
! --------------------------------------------------------------------------------------------------
!
! In  sddisc           : datastructure for time discretization TEMPORELLE
! In  ds_contact       : datastructure for contact management
! IN  NUMINS : NUMERO D'INSTANT
! IN  IECHEC : OCCURRENCE DE L'ECHEC
! OUT IEVDAC : VAUT IECHEC SI EVENEMENT DECLENCHE
!                   0 SINON
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    character(len=24) :: ctevco
    integer(kind=8) :: jctevc
    real(kind=8) :: etacin, etacfi, etacol
    real(kind=8) :: fincol, subdur
    integer(kind=8) :: ntpc, iptc
    real(kind=8) :: instam, instap
    aster_logical :: levent
    integer(kind=8) :: zeven
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<MECANONLINE> ... COLLISION'
    end if
!
! - Initializations
!
    i_echec_acti = 0
    levent = .false.
!
! - Get contact parameters
!
    ntpc = cfdisi(ds_contact%sdcont_defi, 'NTPC')
!
! - Get time dicretization parameters
!
    ASSERT(nume_inst .gt. 0)
    instap = diinst(sddisc, nume_inst)
    instam = diinst(sddisc, nume_inst-1)
!
! - Get event parameters
!
    call utdidt('L', sddisc, 'ECHE', 'SUBD_DUREE', index_=i_echec, &
                valr_=subdur)
!
! - Access to contact datastructures
!
    ctevco = ds_contact%sdcont_solv(1:14)//'.EVENCO'
    call jeveuo(ctevco, 'E', jctevc)
    zeven = cfmmvd('ZEVEN')
!
! --- STATUT DE LA COLLISION
!
    do iptc = 1, ntpc
        etacin = zr(jctevc+zeven*(iptc-1)+1-1)
        etacfi = zr(jctevc+zeven*(iptc-1)+2-1)
        etacol = zr(jctevc+zeven*(iptc-1)+3-1)
        fincol = zr(jctevc+zeven*(iptc-1)+4-1)
!
! ----- COLLISION EN COURS CE N'EST PAS UN EVENEMENT
!
        if (etacol .eq. 1.d0) then
            if (instap .gt. fincol) then
                etacol = 0.d0
                fincol = 0.d0
            end if
        else if (etacol .eq. 0.d0) then
!
! ------- COLLISION QUI S'ACTIVE: C'EST UN EVENEMENT
!
            if ((etacin .eq. 1.d0) .or. (etacfi .eq. 1.d0)) then
                etacol = 1.d0
                fincol = instam+subdur
                levent = .true.
            end if
        end if
        zr(jctevc+zeven*(iptc-1)+3-1) = etacol
        zr(jctevc+zeven*(iptc-1)+4-1) = fincol
    end do
!
! --- ACTIVATION EVENEMENT
!
    if (levent) then
        i_echec_acti = i_echec
    end if
!
    call jedema()
end subroutine
