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

subroutine nmadcp(sddisc, ds_contact, i_event_acti, retpen)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cfdisd.h"
#include "asterfort/cfdisi.h"
#include "asterfort/cfmmco.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utdidt.h"
#include "asterfort/utmess.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    integer(kind=8), intent(in) :: i_event_acti
    type(NL_DS_Contact), intent(in) :: ds_contact
    character(len=19), intent(in) :: sddisc
    integer(kind=8), intent(out) :: retpen
!
! ----------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME)
!
! GESTION DE L'ACTION ADAPTATION COEF. PENALISATION
!
! ----------------------------------------------------------------------
!
! In  ds_contact       : datastructure for contact management
! In  sddisc           : datastructure for time discretization
! IN  i_event_acti     : INDICE DE L'EVENEMENT ACTIF
! OUT RETPEN : CODE RETOUR ADAPTATION PENALISATION
!               0 ON N'A PAS ADAPTE
!               1 ON A ADAPTE
!
! ----------------------------------------------------------------------
!
    real(kind=8) :: pene_maxi, coefpn, newcoe
    real(kind=8) :: coef_maxi
    real(kind=8) :: jeumin, jeumax, jeufin
    integer(kind=8) :: nbliai, nb_cont_zone
    integer(kind=8) :: iliai, i_zone
    character(len=24) :: jeuite, numlia
    integer(kind=8) :: jjeuit, jnumli
    character(len=24) :: ctevpe
    integer(kind=8) :: jctevp
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATIONS
!
    retpen = 1
    call utdidt('L', sddisc, 'ECHE', 'PENE_MAXI', index_=i_event_acti, &
                valr_=pene_maxi)
    call utdidt('L', sddisc, 'ECHE', 'COEF_MAXI', index_=i_event_acti, &
                valr_=coef_maxi)
!
! - Get contact parameters
!
    nbliai = cfdisd(ds_contact%sdcont_solv, 'NBLIAI')
    nb_cont_zone = cfdisi(ds_contact%sdcont_defi, 'NZOCO')
!
! --- ACCES OBJETS DU CONTACT
!
    jeuite = ds_contact%sdcont_solv(1:14)//'.JEUITE'
    numlia = ds_contact%sdcont_solv(1:14)//'.NUMLIA'
    call jeveuo(jeuite, 'L', jjeuit)
    call jeveuo(numlia, 'L', jnumli)
    ctevpe = ds_contact%sdcont_solv(1:14)//'.EVENPE'
    call jeveuo(ctevpe, 'E', jctevp)
!
! --- DETECTION PENETRATION MAXIMUM/MINIMUM
!
    do iliai = 1, nbliai
        jeufin = zr(jjeuit+3*(iliai-1)+1-1)
        i_zone = zi(jnumli+4*(iliai-1)+4-1)
        jeumin = zr(jctevp+3*(i_zone-1)+1-1)
        jeumax = zr(jctevp+3*(i_zone-1)+2-1)
        if (jeufin .le. 0.d0) then
            jeufin = abs(jeufin)
            jeumax = max(jeumax, jeufin)
        else
            jeumin = max(jeumin, jeufin)
        end if
        zr(jctevp+3*(i_zone-1)+1-1) = jeumin
        zr(jctevp+3*(i_zone-1)+2-1) = jeumax
        zr(jctevp+3*(i_zone-1)+3-1) = jeufin
    end do
!
! --- DETECTION PENETRATION MAXIMUM
!
    do i_zone = 1, nb_cont_zone
        call cfmmco(ds_contact, i_zone, 'E_N', 'L', coefpn)
        if (jeumax .gt. pene_maxi) then
            newcoe = coefpn*2.d0
            if (newcoe .gt. coef_maxi) then
                newcoe = coef_maxi
                retpen = 0
            end if
            call cfmmco(ds_contact, i_zone, 'E_N', 'E', newcoe)
        end if
        if (retpen .eq. 1) then
            call utmess('I', 'MECANONLINE10_46', si=i_zone, sr=newcoe)
        end if
    end do
!
! --- AFFICHAGE
!
    if (retpen .eq. 0) then
        call utmess('I', 'MECANONLINE10_44')
    else if (retpen .eq. 1) then
        call utmess('I', 'MECANONLINE10_45')
    else
        ASSERT(.false.)
    end if
!
    call jedema()
end subroutine
