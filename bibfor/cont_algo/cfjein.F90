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

subroutine cfjein(mesh, ds_contact, disp_cumu_inst)
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
#include "asterfort/cfdisr.h"
#include "asterfort/cfimp2.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=8), intent(in) :: mesh
    type(NL_DS_Contact), intent(in) :: ds_contact
    character(len=19), intent(in) :: disp_cumu_inst
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Solve
!
! Discrete methods - Compute initial gaps
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  ds_contact       : datastructure for contact management
! In  disp_cumu_inst   : displacement increment from beginning of current time
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: nbddl, jdecal
    integer(kind=8) :: iliai, ialarm
    character(len=19) :: ddepl0
    character(len=24) :: jeuite, jeux
    integer(kind=8) :: jjeuit, jjeux
    character(len=24) :: apddl, apcoef, appoin
    integer(kind=8) :: japddl, japcoe, japptr
    character(len=24) :: apcofr
    integer(kind=8) :: japcof
    integer(kind=8) :: nbliai, neq, ndimg, nesmax
    aster_logical :: lgliss, lctfd, l_first_geom
    real(kind=8) :: aljeu
    real(kind=8) :: jeuold, jeuini, jexini, jeyini
    real(kind=8) :: val1, val2, val
    real(kind=8), pointer :: ddep0(:) => null()
    real(kind=8), pointer :: depde(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infdbg('CONTACT', ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<CONTACT> ......... CALCUL DES JEUX INITIAUX'
    end if
!
! - Get contact parameters
!
    ndimg = cfdisd(ds_contact%sdcont_solv, 'NDIM')
    nbliai = cfdisd(ds_contact%sdcont_solv, 'NBLIAI')
    neq = cfdisd(ds_contact%sdcont_solv, 'NEQ')
    nesmax = cfdisd(ds_contact%sdcont_solv, 'NESMAX')
    lctfd = cfdisl(ds_contact%sdcont_defi, 'FROT_DISCRET')
    lgliss = cfdisl(ds_contact%sdcont_defi, 'CONT_DISC_GLIS')
    aljeu = cfdisr(ds_contact%sdcont_defi, 'ALARME_JEU')
!
! --- INITIALISATIONS
!
    ialarm = 0
!
! --- LECTURE DES STRUCTURES DE DONNEES DE CONTACT
!
    appoin = ds_contact%sdcont_solv(1:14)//'.APPOIN'
    apddl = ds_contact%sdcont_solv(1:14)//'.APDDL'
    apcoef = ds_contact%sdcont_solv(1:14)//'.APCOEF'
    call jeveuo(appoin, 'L', japptr)
    call jeveuo(apddl, 'L', japddl)
    call jeveuo(apcoef, 'L', japcoe)
!
    jeuite = ds_contact%sdcont_solv(1:14)//'.JEUITE'
    jeux = ds_contact%sdcont_solv(1:14)//'.JEUX'
    call jeveuo(jeuite, 'E', jjeuit)
    call jeveuo(jeux, 'E', jjeux)
!
    if (lctfd) then
        apcofr = ds_contact%sdcont_solv(1:14)//'.APCOFR'
        call jeveuo(apcofr, 'L', japcof)
    end if
!
! - Get geometric loop state
!
    l_first_geom = ds_contact%l_first_geom
!
! --- ACCES AUX CHAMPS DE TRAVAIL
! --- DDEPL0: INCREMENT DE SOLUTION SANS CORRECTION DU CONTACT
!
    ddepl0 = ds_contact%sdcont_solv(1:14)//'.DEL0'
    call jeveuo(ddepl0(1:19)//'.VALE', 'L', vr=ddep0)
!
! --- INCREMENT DE DEPLACEMENT DEPUIS LE DEBUT DU PAS DE TEMPS
!
    call jeveuo(disp_cumu_inst(1:19)//'.VALE', 'L', vr=depde)
!
! --- ON CALCULE LE NOUVEAU JEU : AJEU+ = AJEU/I/N - A.DDEPLA
!
    do iliai = 1, nbliai
!
! ----- ACCES TABLEAU LIAISONS
!
        jdecal = zi(japptr+iliai-1)
        nbddl = zi(japptr+iliai)-zi(japptr+iliai-1)
!
! ----- INCR. DE JEU SANS CORRECTION [A].{DDEPL0}
!
        call caladu(neq, nbddl, zr(japcoe+jdecal), zi(japddl+jdecal), ddep0, &
                    val)
!
! ----- JEU AVANT L'ITERATION DE NEWTON {JEU(DEPTOT)}
!
        jeuold = zr(jjeuit+3*(iliai-1)+1-1)
!
! ----- JEU SANS CORRECTION DU CONTACT: {JEU(DEPTOT)} - [A].{DDEPL0}
!
        jeuini = jeuold-val
!
! ----- SAUVEGARDES
!
        zr(jjeux+3*(iliai-1)+1-1) = jeuini
!
! ----- JEUX TANGENTS
!
        if (lctfd) then
!
! ------- INCR. JEUX TANGENTS SANS CORRECTION
!
            call caladu(neq, nbddl, zr(japcof+jdecal), zi(japddl+jdecal), ddep0, &
                        val1)
!
! ------- INCR. JEUX TANGENTS DEPUIS LE DEBUT DU PAS DE TEMPS
!
            call caladu(neq, nbddl, zr(japcof+jdecal), zi(japddl+jdecal), depde, &
                        val2)
!
! ------- INCR. DE JEU DEPUIS LE DEBUT DU PAS DE TEMPS SANS CORR.
!
            jexini = val1+val2
!
            zr(jjeux+3*(iliai-1)+2-1) = jexini
!
            if (ndimg .eq. 3) then
!
! --------- INCR. DE JEU SANS CORRECTION
!
                call caladu(neq, nbddl, zr(japcof+jdecal+30*nesmax), zi(japddl+jdecal), ddep0, &
                            val1)
!
! --------- INCR. DE JEU DEPUIS LE DEBUT DU PAS DE TEMPS
!
                call caladu(neq, nbddl, zr(japcof+jdecal+30*nesmax), zi(japddl+jdecal), depde, &
                            val2)
!
! --------- INCR. DE JEU DEPUIS LE DEBUT DU PAS DE TEMPS SANS CORR.
!
                jeyini = val1+val2
                zr(jjeux+3*(iliai-1)+3-1) = jeyini
            end if
        end if
!
    end do
!
! --- ALARME SI DECOLLEMENT ALORS QUE GLISSIERE
!
    do iliai = 1, nbliai
        if (lgliss) then
            jeuini = zr(jjeux+3*(iliai-1)+1-1)
            if (jeuini .gt. aljeu) then
                ialarm = ialarm+1
                if (ialarm .eq. 1) then
                    call utmess('A', 'CONTACT_9')
                end if
                call cfimp2(ds_contact%sdcont_defi, ds_contact%sdcont_solv, mesh, iliai, &
                            'ALJ')
            end if
        end if
    end do
!
    call jedema()
!
end subroutine
