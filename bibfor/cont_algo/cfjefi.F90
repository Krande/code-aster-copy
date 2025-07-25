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

subroutine cfjefi(mesh, disp_iter, ds_contact)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/caladu.h"
#include "asterfort/cfdisd.h"
#include "asterfort/cfdisl.h"
#include "asterfort/cfimp1.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=8), intent(in) :: mesh
    character(len=19), intent(in) :: disp_iter
    type(NL_DS_Contact), intent(in) :: ds_contact
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Solve
!
! Discrete methods - Compute final gaps
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  disp_iter        : displacement iteration
! In  ds_contact       : datastructure for contact management
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: iliai, jdecal, nbddl
    real(kind=8) :: jeuini, jeuold, jeuinc
    aster_logical :: l_pena_cont, l_frot
    character(len=24) :: sdcont_apcoef, sdcont_apddl, sdcont_appoin
    integer(kind=8) :: japcoe, japddl, japptr
    character(len=24) :: sdcont_apcofr
    integer(kind=8) :: japcof
    character(len=24) :: sdcont_jeuite, sdcont_jeux
    integer(kind=8) :: jjeuit, jjeux
    integer(kind=8) :: nbliai, nb_equa, model_ndim
    real(kind=8), pointer :: vale(:) => null()
!
! ----------------------------------------------------------------------
!
    call jemarq()
    call infdbg('CONTACT', ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<CONTACT> ...... CALCUL DES JEUX FINAUX'
    end if
!
! - Get contact parameters
!
    nbliai = cfdisd(ds_contact%sdcont_solv, 'NBLIAI')
    nb_equa = cfdisd(ds_contact%sdcont_solv, 'NEQ')
    model_ndim = cfdisd(ds_contact%sdcont_solv, 'NDIM')
    l_pena_cont = cfdisl(ds_contact%sdcont_defi, 'CONT_PENA')
    l_frot = cfdisl(ds_contact%sdcont_defi, 'FROT_DISCRET')
!
! - Access to contact datastructures
!
    sdcont_appoin = ds_contact%sdcont_solv(1:14)//'.APPOIN'
    sdcont_apddl = ds_contact%sdcont_solv(1:14)//'.APDDL'
    sdcont_apcoef = ds_contact%sdcont_solv(1:14)//'.APCOEF'
    call jeveuo(sdcont_appoin, 'L', japptr)
    call jeveuo(sdcont_apddl, 'L', japddl)
    call jeveuo(sdcont_apcoef, 'L', japcoe)
    if (l_frot) then
        sdcont_apcofr = ds_contact%sdcont_solv(1:14)//'.APCOFR'
        call jeveuo(sdcont_apcofr, 'L', japcof)
    end if
    sdcont_jeuite = ds_contact%sdcont_solv(1:14)//'.JEUITE'
    sdcont_jeux = ds_contact%sdcont_solv(1:14)//'.JEUX'
    call jeveuo(sdcont_jeux, 'L', jjeux)
    call jeveuo(sdcont_jeuite, 'E', jjeuit)
!
! - Access to displacements
!
    call jeveuo(disp_iter(1:19)//'.VALE', 'L', vr=vale)
!
! - Gap update
!
    do iliai = 1, nbliai
        jeuini = zr(jjeux+3*(iliai-1)+1-1)
        if (l_pena_cont) then
            zr(jjeuit+3*(iliai-1)+1-1) = jeuini
        else
            jdecal = zi(japptr+iliai-1)
            nbddl = zi(japptr+iliai)-zi(japptr+iliai-1)
            call caladu(nb_equa, nbddl, zr(japcoe+jdecal), zi(japddl+jdecal), vale, &
                        jeuinc)
            jeuold = zr(jjeuit+3*(iliai-1)+1-1)
            zr(jjeuit+3*(iliai-1)+1-1) = jeuold-jeuinc
        end if
    end do
!
! - Print
!
    if (niv .ge. 2) then
        call cfimp1('FIN', mesh, ds_contact%sdcont_defi, ds_contact%sdcont_solv, ifm)
    end if
!
    call jedema()
!
end subroutine
