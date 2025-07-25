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
subroutine nmarc0(result, modele, ds_material, carele, fonact, &
                  sdcrit, sddyna, ds_errorindic, &
                  sdpilo, listLoadResu, numarc, time_curr)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/isfonc.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/ndaram.h"
#include "asterfort/ndynlo.h"
#include "asterfort/ndynre.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rssepa.h"
#include "asterfort/Behaviour_type.h"
!
    character(len=8) :: result
    integer(kind=8) :: numarc
    integer(kind=8) :: fonact(*)
    real(kind=8) :: time_curr
    character(len=19) :: sddyna, sdpilo
    character(len=19) :: sdcrit
    character(len=24) :: modele, carele, listLoadResu
    type(NL_DS_ErrorIndic), intent(in) :: ds_errorindic
    type(NL_DS_Material), intent(in) :: ds_material
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME - ARCHIVAGE)
!
! ARCHIVAGE DES PARAMETRES
!
! --------------------------------------------------------------------------------------------------
!
! IN  RESULT : NOM UTILISATEUR DU CONCEPT RESULTAT
! IN  MODELE : NOM DU MODELE
! In  ds_material      : datastructure for material parameters
! IN  SDCRIT : VALEUR DES CRITERES DE CONVERGENCE
! IN  SDPILO : SD PILOTAGE
! In  ds_errorindic    : datastructure for error indicator
! IN  COMPOR : CARTE DECRIVANT LE TYPE DE COMPORTEMENT
! IN  CARELE : CARACTERISTIQUES DES ELEMENTS DE STRUCTURE
! IN  FONACT : FONCTIONNALITES ACTIVEES (VOIR NMFONC)
! IN  SDDYNA : SD DEDIEE A LA DYNAMIQUE
! In  list_load_resu : name of list of loads saved in results datastructure
! IN  NUMARC : NUMERO D'ARCHIVAGE
! IN  INSTAN : VALEUR DE L'INSTANT
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: valk
    aster_logical :: lerrt, lthm, lpilo, ldyna
    aster_logical :: lexge
    character(len=24) :: typsel, typpil
    real(kind=8) :: valr, coef, time_prev
    integer(kind=8) :: jv_para
    real(kind=8), pointer :: plir(:) => null()
    real(kind=8), pointer :: crtr(:) => null()
    character(len=24), pointer :: pltk(:) => null()
    character(len=16), pointer :: crde(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! --- FONCTIONNALITES ACTIVEES
!
    ldyna = ndynlo(sddyna, 'DYNAMIQUE')
    lexge = ndynlo(sddyna, 'EXPL_GENE')
    lerrt = isfonc(fonact, 'ERRE_TEMPS_THM')
    lthm = isfonc(fonact, 'THM')
    lpilo = isfonc(fonact, 'PILOTAGE')
!
! --- ARCHIVAGE DE THETA EN THM
!
    if (lthm) then
        call rsadpa(result, 'E', 1, 'PARM_THETA', numarc, 0, sjv=jv_para)
        zr(jv_para) = ds_errorindic%parm_theta
    end if
!
! --- ARCHIVAGE DE L'INSTANT
!
    call rsadpa(result, 'E', 1, 'INST', numarc, 0, sjv=jv_para)
    zr(jv_para) = time_curr
!
! --- ARCHIVAGE DE L'INSTANT PRECEDENT
!
    if (ldyna) then
        time_prev = ndynre(sddyna, 'INST_PREC')
    else
        time_prev = time_curr
    end if

    call rsadpa(result, 'E', 1, 'INST_PREC', numarc, 0, sjv=jv_para)
    zr(jv_para) = time_prev
!
! --- ARCHIVAGE DU MODELE, MATERIAU, CARA_ELEM ET DE LA SD CHARGE
!
    call rssepa(result, numarc, modele(1:8), ds_material%mater(1:8), carele(1:8), &
                listLoadResu)
!
! --- ARCHIVAGE DES CRITERES DE CONVERGENCE
!
    call jeveuo(sdcrit//'.CRTR', 'L', vr=crtr)
    call jeveuo(sdcrit//'.CRDE', 'L', vk16=crde)
    valr = crtr(1)
    valk = crde(1)
    call rsadpa(result, 'E', 1, valk, numarc, 0, sjv=jv_para)
    zi(jv_para) = nint(valr)
    valr = crtr(5)
    valk = crde(5)
    call rsadpa(result, 'E', 1, valk, numarc, 0, sjv=jv_para)
    zr(jv_para) = valr
    valr = crtr(6)
    valk = crde(6)
    call rsadpa(result, 'E', 1, valk, numarc, 0, sjv=jv_para)
    zr(jv_para) = valr
!
! --- ARCHIVAGE DES INDICATEURS D'ERREUR EN TEMPS EN THM UNIQUEMENT
!
    if (lerrt) then
        call rsadpa(result, 'E', 1, 'ERRE_TPS_LOC', numarc, 0, sjv=jv_para)
        zr(jv_para) = ds_errorindic%erre_thm_loca
        call rsadpa(result, 'E', 1, 'ERRE_TPS_GLOB', numarc, 0, sjv=jv_para)
        zr(jv_para) = ds_errorindic%erre_thm_glob
    end if
!
! --- ARCHIVAGE DE COEF_MULT SI PILOTAGE
!
    if (lpilo) then
        call jeveuo(sdpilo(1:19)//'.PLTK', 'L', vk24=pltk)
        typpil = pltk(1)
        typsel = pltk(6)
        if ((typpil .eq. 'LONG_ARC') .and. (typsel .eq. 'ANGL_INCR_DEPL')) then
            call jeveuo(sdpilo(1:19)//'.PLIR', 'L', vr=plir)
            coef = plir(1)
            call rsadpa(result, 'E', 1, 'COEF_MULT', numarc, 0, sjv=jv_para)
            zr(jv_para) = coef
        end if
    end if
!
! --- ARCHIVAGE DEPL/VITE/ACCE GENERALISES
!
    if (lexge) then
        call ndaram(result, sddyna, numarc)
    end if
!
    call jedema()
end subroutine
