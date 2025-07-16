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
subroutine vrcomp(comporCurrZ, variZ, ligrelCurrZ, iret, &
                  comporPrevZ_, typeStop_, verbose_, &
                  lModiVari_)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jemarq.h"
#include "asterfort/jedema.h"
#include "asterfort/utmess.h"
#include "asterfort/vrcomp_chck_cmp.h"
#include "asterfort/vrcomp_chck_rela.h"
#include "asterfort/vrcomp_prep.h"
!
    character(len=*), intent(in) :: comporCurrZ, variZ, ligrelCurrZ
    integer(kind=8), intent(out) :: iret
    character(len=*), optional, intent(in) :: comporPrevZ_
    character(len=1), optional, intent(in) :: typeStop_
    aster_logical, intent(in), optional :: verbose_
    aster_logical, intent(out), optional :: lModiVari_
!
! --------------------------------------------------------------------------------------------------
!
! Check compatibility of comportments
!
! Is internal variable field is correct with comportment definition ?
!
! --------------------------------------------------------------------------------------------------
!
! In  vari          : internal variable field
! In  comporCurr   : current comportment
! In  ligrelCurr   : current LIGREL
! In  comporPrev   : previous comportment
!
! --------------------------------------------------------------------------------------------------
!
!
! ------------------------------------------------------------------
! BUT: VERIFIER LA COHERENCE DU CHAMP DE VARIABLES INTERNES "-" AVEC
!      LE COMPORTEMENT CHOISI.
!      ATTENTION : ON MODIFIE PARFOIS VARMOI POUR TENIR COMPTE DES
!      POSSIBLES CHANGEMENTS DE MODELE ET/OU DE COMPORTEMENT
!
!
!      VARMOI EST LE CHAMP DE VARIABLES INTERNES A L'INSTANT "-"
!      COMPOM EST LA CARTE DE COMPORTEMENT DE VARMOI (OU ' ')
!      COMPOP EST LA CARTE DE COMPOPTEMENT A L'INSTANT "+"
!      LIGREP EST LE LIGREL DU MODELE DE L'INSTANT "+"
!
!      - SI COMPOM = ' ' :
!           ON SE CONTENTE DE COMPARER LE NOMBRE DE V.I. DE VARMOI
!           AVEC LE NOMBRE ATTENDU DANS COMPOP.
!      - SI COMPOM /= ' ':
!           ON PEUT ALORS COMPARER LE NOM DES COMPORTEMENTS DES
!           INSTANTS "+" ET "-"
!           ON EXIGE QUE CES NOMS SOIENT IDENTIQUES OU BIEN QUE :
!             "-" : 'ELAS' OU 'SANS'  -> "+" : N'IMPORTE QUOI
!             "+" : 'ELAS' OU 'SANS'  -> "-" : N'IMPORTE QUOI
!
!
! ------------------------------------------------------------------
!     ARGUMENTS:
! COMPOM   IN/JXIN  K19 : CARTE DE COMPOPTEMENT "-"
! COMPOP   IN/JXIN  K19 : CARTE DE COMPOPTEMENT "+"
! COMPOP   EST AUSSI LE NOM DU CHAM_ELEM_S DE DCEL_I PERMETTANT DE
!          DE CONNAITRE LE NOMBRE DE SOUS-POINTS ET LE NOMBRE DE VARI
! VARMOI   IN/JXVAR K19 : SD CHAM_ELEM   (VARI_R) "-"
! TYPE_STOP IN   : COMPORTEMENT SI PROBLEME DETECTE
!                  'A' ALARME
!                  'E' ERREUR
! IRET      OUT  : CODE RETOUR
!                   0 SI PAS DE PROBLEME
!
! REMARQUES :
!  - VARMOI EST PARFOIS MODIFIE POUR ETRE COHERENT AVEC COMPOP
!           ON LE RECREE ALORS SUR LA BASE VOLATILE
!  - ON VERIFIE EGALEMENT LE NOMBRE DES SOUS-POINTS
!
!-----------------------------------------------------------------------
!
    aster_logical :: lModiVari
    character(len=1) :: typeStop
    character(len=19) :: variRedu, vari
    character(len=19) :: comporCurr, comporPrev
    character(len=19) :: comporCurrRedu, comporPrevRedu
! - Beahviours can been mixed with each other
    character(len=48), parameter :: comp_comb_1 = &
                                    'LEMAITRE        VMIS_ISOT_LINE  VMIS_ISOT_TRAC'
! - Beahviours can been mixed with all other ones
    character(len=48), parameter :: comp_comb_2 = &
                                    'ELAS            SANS            KIT_CG'
    aster_logical :: newBehaviourOnCell, nbSpgDifferent
    aster_logical :: inconsistentBehaviour, nbVariDifferent
    character(len=19) :: ligrelCurr, ligrelPrev
    aster_logical :: verbose
    character(len=8) :: meshCompor, meshField, mesh
    integer(kind=8) :: nbCell
    character(len=8), pointer :: cesk(:) => null()
    integer(kind=8), pointer :: cesd(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
!     -- MODIF : .TRUE. => IL FAUT MODIFIER VARMOI CAR CERTAINES
!        MAILLES ONT DISPARU OU SONT NOUVELLES OU ONT CHANGE DE
!        COMPORTEMENT
    lModiVari = ASTER_FALSE
    nbVariDifferent = ASTER_FALSE
    nbSpgDifferent = ASTER_FALSE
    newBehaviourOnCell = ASTER_FALSE
    inconsistentBehaviour = ASTER_FALSE
    verbose = ASTER_FALSE
    if (present(verbose_)) then
        verbose = verbose_
    end if
    comporPrev = ' '
    if (present(comporPrevZ_)) then
        comporPrev = comporPrevZ_
    end if
    comporCurr = comporCurrZ
    ligrelCurr = ligrelCurrZ
    vari = variZ

! - Error management
    iret = 0
    if (present(typeStop_)) then
        typeStop = typeStop_
    else
        typeStop = 'E'
    end if

! - Acces to reduced CARTE DCEL_I (see CESVAR) on current comportement
    call jeveuo(comporCurr//'.CESD', 'L', vi=cesd)
    call jeveuo(comporCurr//'.CESK', 'L', vk8=cesk)
    meshCompor = cesk(1)
    nbCell = cesd(1)

! - Get LIGREL from VARI_ELGA field
    call dismoi('NOM_LIGREL', vari, 'CHAM_ELEM', repk=ligrelPrev)

! - Check meshes
    call dismoi('NOM_MAILLA', vari, 'CHAMP', repk=meshField)
    if (meshCompor .ne. meshField) then
        call utmess('F', 'COMPOR6_1')
    end if
    mesh = meshCompor

! - Prepare fields
    if (present(comporPrevZ_)) then
        call vrcomp_prep(vari, variRedu, &
                         comporCurr, comporCurrRedu, &
                         comporPrev, comporPrevRedu)
    else
        call vrcomp_prep(vari, variRedu, &
                         comporCurr, comporCurrRedu)
        comporPrevRedu = ' '
    end if

! - Check if comportments are the same (or compatible)
    if (present(comporPrevZ_)) then
        call vrcomp_chck_rela(mesh, nbCell, &
                              comporCurrRedu, comporPrevRedu, &
                              ligrelCurr, ligrelPrev, &
                              comp_comb_1, comp_comb_2, verbose, &
                              newBehaviourOnCell, inconsistentBehaviour, &
                              lModiVari)
    end if
    if (newBehaviourOnCell) then
        iret = 1
        call utmess(typeStop, 'COMPOR6_2')
    end if
    if (inconsistentBehaviour) then
        iret = 1
        call utmess(typeStop, 'COMPOR6_3')
    end if

! - Check if elements have the same number of internal variables and Gauss-subpoints
    call vrcomp_chck_cmp(mesh, nbCell, &
                         comporCurr, comporCurrRedu, comporPrevRedu, &
                         variRedu, comp_comb_2, &
                         ligrelCurr, ligrelPrev, &
                         verbose, &
                         nbSpgDifferent, nbVariDifferent, lModiVari)
    if (nbSpgDifferent) then
        iret = 1
        call utmess(typeStop, 'COMPOR6_4')
    end if
    if (nbVariDifferent) then
        iret = 1
        call utmess(typeStop, 'COMPOR6_5')
    end if

! - Clean
    call detrsd('CHAM_ELEM_S', variRedu)
    if (present(comporPrevZ_)) then
        call detrsd('CHAM_ELEM_S', comporPrevRedu)
    end if
    call detrsd('CHAM_ELEM_S', comporCurrRedu)

! - Output
    if (present(lModiVari_)) then
        lModiVari_ = lModiVari
    end if
!
    call jedema()
!
end subroutine
