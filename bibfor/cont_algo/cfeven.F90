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

subroutine cfeven(phase, ds_contact)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cfdisd.h"
#include "asterfort/cfdisl.h"
#include "asterfort/cfmmvd.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=3), intent(in) :: phase
    type(NL_DS_Contact), intent(in) :: ds_contact
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (METHODES DISCRETES - EVENT-DRIVEN)
!
! DETECTION D'UNE COLLISION
!
! ----------------------------------------------------------------------
!
! IN  PHASE  : PHASE DE DETECTION
!              'INI' - AU DEBUT DU PAS DE TEMPS
!              'FIN' - A LA FIN DU PAS DE TEMPS
! In  ds_contact       : datastructure for contact management
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: nbliai, nbliac, ip
    integer(kind=8) :: iliai, iliac
    character(len=24) :: numlia, ctevco
    integer(kind=8) :: jnumli, jctevc
    character(len=19) :: liac
    integer(kind=8) :: jliac
    integer(kind=8) :: zeven
    aster_logical :: lactif
    real(kind=8) :: etacin, etacfi
    aster_logical :: lexiv
!
! ----------------------------------------------------------------------
!
    call jemarq()
    call infdbg('CONTACT', ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<CONTACT> ...... GESTION DES JEUX POUR EVENT-DRIVEN'
    end if
!
! --- PARAMETRES
!
    nbliai = cfdisd(ds_contact%sdcont_solv, 'NBLIAI')
    nbliac = cfdisd(ds_contact%sdcont_solv, 'NBLIAC')
    zeven = cfmmvd('ZEVEN')
!
! --- UNE ZONE EN MODE SANS CALCUL: ON NE PEUT RIEN FAIRE
!
    lexiv = cfdisl(ds_contact%sdcont_defi, 'EXIS_VERIF')
    if (lexiv) goto 999
!
! --- ACCES OBJETS DU CONTACT
!
    liac = ds_contact%sdcont_solv(1:14)//'.LIAC'
    numlia = ds_contact%sdcont_solv(1:14)//'.NUMLIA'
    ctevco = ds_contact%sdcont_solv(1:14)//'.EVENCO'
    call jeveuo(liac, 'L', jliac)
    call jeveuo(numlia, 'L', jnumli)
    call jeveuo(ctevco, 'E', jctevc)
!
! --- DETECTION
!
    do iliai = 1, nbliai
        ip = zi(jnumli+4*(iliai-1)+1-1)
        etacin = zr(jctevc+zeven*(ip-1)+1-1)
        etacfi = zr(jctevc+zeven*(ip-1)+2-1)
        lactif = .false.
!
! ----- LA LIAISON EST-ELLE ACTIVE ?
!
        do iliac = 1, nbliac
            if (zi(jliac-1+iliac) .eq. iliai) lactif = .true.
        end do
!
! ----- CHANGEMENT STATUT
!
        if (lactif) then
            if (phase .eq. 'INI') then
                etacin = 1.d0
            else if (phase .eq. 'FIN') then
                etacfi = 1.d0
            else
                ASSERT(.false.)
            end if
        else
            if (phase .eq. 'INI') then
                etacin = 0.d0
            else if (phase .eq. 'FIN') then
                etacfi = 0.d0
            else
                ASSERT(.false.)
            end if
        end if
        zr(jctevc+zeven*(ip-1)+1-1) = etacin
        zr(jctevc+zeven*(ip-1)+2-1) = etacfi
    end do
!
999 continue
!
    call jedema()
end subroutine
