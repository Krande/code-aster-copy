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
subroutine nmcrch(numeDof, listFuncActi, sddyna, nlDynaDamping, &
                  hval_incr, hval_algo, hval_veasse)
!
    use NonLin_Datastructure_type
    use NonLinearDyna_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/isfonc.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/ndynkk.h"
#include "asterfort/ndynlo.h"
#include "asterfort/nmchex.h"
#include "asterfort/vtcreb.h"
!
    character(len=24), intent(in) :: numeDof
    integer(kind=8), intent(in) :: listFuncActi(*)
    character(len=19), intent(in) :: sddyna
    type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
    character(len=19), intent(in) :: hval_incr(*), hval_algo(*), hval_veasse(*)
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Initializations
!
! Create vectors
!
! --------------------------------------------------------------------------------------------------
!
! In  sddyna           : name of datastructure for dynamic parameters
! In  nlDynaDamping    : damping parameters
! IN  NUMEDD : NUME_DDL
! IN  VALINC : VARIABLE CHAPEAU POUR INCREMENTS VARIABLES
! IN  SOLALG : VARIABLE CHAPEAU POUR INCREMENTS SOLUTIONS
! IN  VEASSE : VARIABLE CHAPEAU POUR NOM DES VECT_ASSE
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: ldyna, lDampModal, lmpas, lrefe, lSuperElement, lmuap, lviss
    aster_logical :: lsstf
    aster_logical :: ldidi, lpilo, lener
    character(len=19) :: depplu, vitplu, accplu
    character(len=19) :: depmoi, vitmoi, accmoi
    character(len=19) :: fexmoi, fammoi, flimoi, fnomoi
    character(len=19) :: fexplu, famplu, fliplu, fnoplu
    character(len=19) :: depso1, depso2
    character(len=19) :: depdel, depold, ddepla, deppr1, deppr2
    character(len=19) :: vitdel, vitold, dvitla, vitpr1, vitpr2
    character(len=19) :: accdel, accold, daccla, accpr1, accpr2
    character(len=19) :: depkm1, vitkm1, acckm1, romkm1, romk
    character(len=19) :: depent, vitent, accent
    character(len=19) :: depabs, vitabs, accabs
    character(len=19) :: cndyna, cnamod
    character(len=19) :: cnfext
    character(len=19) :: cnfedo, cnfsdo, cndidi
    character(len=19) :: cndido, cncine, cndiri
    character(len=19) :: cnondp, cnviss, cnhyst, cnfint
    character(len=19) :: cnsstf, cnsstr
    character(len=19) :: cnfepi, cndipi, cnrefe, cneltc, cneltf
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Active functionnalities
    lSuperElement = isfonc(listFuncActi, 'MACR_ELEM_STAT')
    ldyna = ndynlo(sddyna, 'DYNAMIQUE')
    lDampModal = nlDynaDamping%lDampModal
    lmpas = ndynlo(sddyna, 'MULTI_PAS')
    ldidi = isfonc(listFuncActi, 'DIDI')
    lpilo = isfonc(listFuncActi, 'PILOTAGE')
    lsstf = isfonc(listFuncActi, 'SOUS_STRUC')
    lrefe = isfonc(listFuncActi, 'RESI_REFE')
    lviss = ndynlo(sddyna, 'VECT_ISS')
    lener = isfonc(listFuncActi, 'ENERGIE')
    lmuap = ndynlo(sddyna, 'MULTI_APPUI')

! - For energy
    if (lener) then
        call nmchex(hval_incr, 'VALINC', 'FEXMOI', fexmoi)
        call nmchex(hval_incr, 'VALINC', 'FAMMOI', fammoi)
        call nmchex(hval_incr, 'VALINC', 'FLIMOI', flimoi)
        call nmchex(hval_incr, 'VALINC', 'FNOMOI', fnomoi)
        call nmchex(hval_incr, 'VALINC', 'FEXPLU', fexplu)
        call nmchex(hval_incr, 'VALINC', 'FAMPLU', famplu)
        call nmchex(hval_incr, 'VALINC', 'FLIPLU', fliplu)
        call nmchex(hval_incr, 'VALINC', 'FNOPLU', fnoplu)
        call vtcreb(fexmoi, 'V', 'R', nume_ddlz=numeDof)
        call vtcreb(fammoi, 'V', 'R', nume_ddlz=numeDof)
        call vtcreb(flimoi, 'V', 'R', nume_ddlz=numeDof)
        call vtcreb(fnomoi, 'V', 'R', nume_ddlz=numeDof)
        call vtcreb(fexplu, 'V', 'R', nume_ddlz=numeDof)
        call vtcreb(famplu, 'V', 'R', nume_ddlz=numeDof)
        call vtcreb(fliplu, 'V', 'R', nume_ddlz=numeDof)
        call vtcreb(fnoplu, 'V', 'R', nume_ddlz=numeDof)
    end if

! --- CREATION DES CHAMPS DE BASE - ETAT EN T-
    call nmchex(hval_incr, 'VALINC', 'DEPMOI', depmoi)
    call nmchex(hval_incr, 'VALINC', 'VITMOI', vitmoi)
    call nmchex(hval_incr, 'VALINC', 'ACCMOI', accmoi)
    call vtcreb(depmoi, 'V', 'R', nume_ddlz=numeDof)
    if (ldyna) then
        call vtcreb(vitmoi, 'V', 'R', nume_ddlz=numeDof)
        call vtcreb(accmoi, 'V', 'R', nume_ddlz=numeDof)
    end if
!
! --- CREATION DES CHAMPS DE BASE - ETAT EN T+
!
    call nmchex(hval_incr, 'VALINC', 'DEPPLU', depplu)
    call nmchex(hval_incr, 'VALINC', 'VITPLU', vitplu)
    call nmchex(hval_incr, 'VALINC', 'ACCPLU', accplu)
    call vtcreb(depplu, 'V', 'R', nume_ddlz=numeDof)
    if (ldyna) then
        call vtcreb(vitplu, 'V', 'R', nume_ddlz=numeDof)
        call vtcreb(accplu, 'V', 'R', nume_ddlz=numeDof)
    end if
!
! --- CREATION DES CHAMPS DE BASE - POUTRES EN GRANDES ROTATIONS
!
    call nmchex(hval_incr, 'VALINC', 'DEPKM1', depkm1)
    call nmchex(hval_incr, 'VALINC', 'VITKM1', vitkm1)
    call nmchex(hval_incr, 'VALINC', 'ACCKM1', acckm1)
    call nmchex(hval_incr, 'VALINC', 'ROMKM1', romkm1)
    call nmchex(hval_incr, 'VALINC', 'ROMK  ', romk)
    call vtcreb(depkm1, 'V', 'R', nume_ddlz=numeDof)
    call vtcreb(vitkm1, 'V', 'R', nume_ddlz=numeDof)
    call vtcreb(acckm1, 'V', 'R', nume_ddlz=numeDof)
    call vtcreb(romkm1, 'V', 'R', nume_ddlz=numeDof)
    call vtcreb(romk, 'V', 'R', nume_ddlz=numeDof)
!
! --- CREATION DES CHAMPS DE BASE - INCREMENTS SOLUTIONS
!
    call nmchex(hval_algo, 'SOLALG', 'DEPDEL', depdel)
    call nmchex(hval_algo, 'SOLALG', 'DDEPLA', ddepla)
    call nmchex(hval_algo, 'SOLALG', 'DEPPR1', deppr1)
    call nmchex(hval_algo, 'SOLALG', 'DEPPR2', deppr2)
    call nmchex(hval_algo, 'SOLALG', 'DEPOLD', depold)
    call vtcreb(depdel, 'V', 'R', nume_ddlz=numeDof)
    call vtcreb(ddepla, 'V', 'R', nume_ddlz=numeDof)
    call vtcreb(depold, 'V', 'R', nume_ddlz=numeDof)
    call vtcreb(deppr1, 'V', 'R', nume_ddlz=numeDof)
    call vtcreb(deppr2, 'V', 'R', nume_ddlz=numeDof)
    if (ldyna) then
        call nmchex(hval_algo, 'SOLALG', 'VITDEL', vitdel)
        call nmchex(hval_algo, 'SOLALG', 'DVITLA', dvitla)
        call nmchex(hval_algo, 'SOLALG', 'VITPR1', vitpr1)
        call nmchex(hval_algo, 'SOLALG', 'VITPR2', vitpr2)
        call nmchex(hval_algo, 'SOLALG', 'VITOLD', vitold)
        call vtcreb(vitdel, 'V', 'R', nume_ddlz=numeDof)
        call vtcreb(dvitla, 'V', 'R', nume_ddlz=numeDof)
        call vtcreb(vitold, 'V', 'R', nume_ddlz=numeDof)
        call vtcreb(vitpr1, 'V', 'R', nume_ddlz=numeDof)
        call vtcreb(vitpr2, 'V', 'R', nume_ddlz=numeDof)
        call nmchex(hval_algo, 'SOLALG', 'ACCDEL', accdel)
        call nmchex(hval_algo, 'SOLALG', 'DACCLA', daccla)
        call nmchex(hval_algo, 'SOLALG', 'ACCPR1', accpr1)
        call nmchex(hval_algo, 'SOLALG', 'ACCPR2', accpr2)
        call nmchex(hval_algo, 'SOLALG', 'ACCOLD', accold)
        call vtcreb(accdel, 'V', 'R', nume_ddlz=numeDof)
        call vtcreb(daccla, 'V', 'R', nume_ddlz=numeDof)
        call vtcreb(accold, 'V', 'R', nume_ddlz=numeDof)
        call vtcreb(accpr1, 'V', 'R', nume_ddlz=numeDof)
        call vtcreb(accpr2, 'V', 'R', nume_ddlz=numeDof)
    end if
!
! --- REACTIONS D'APPUI BT.LAMBDA
!
    call nmchex(hval_veasse, 'VEASSE', 'CNDIRI', cndiri)
    call vtcreb(cndiri, 'V', 'R', nume_ddlz=numeDof)
!
! --- VECTEURS SOLUTION
!
    call nmchex(hval_algo, 'SOLALG', 'DEPSO1', depso1)
    call nmchex(hval_algo, 'SOLALG', 'DEPSO2', depso2)
    call vtcreb(depso1, 'V', 'R', nume_ddlz=numeDof)
    call vtcreb(depso2, 'V', 'R', nume_ddlz=numeDof)
!
! --- CREATION DES CHAMPS DE BASE - DEPL/VITE/ACCE D'ENTRAINEMENT
!
    if (ldyna) then
        call ndynkk(sddyna, 'DEPENT', depent)
        call ndynkk(sddyna, 'VITENT', vitent)
        call ndynkk(sddyna, 'ACCENT', accent)
        call vtcreb(depent, 'V', 'R', nume_ddlz=numeDof)
        call vtcreb(vitent, 'V', 'R', nume_ddlz=numeDof)
        call vtcreb(accent, 'V', 'R', nume_ddlz=numeDof)
    end if
!
! --- CREATION DES CHAMPS DEPL/VITE/ACCE ABSOLUS POUR LE MULTI-APPUIS
!
    if (lmuap) then
        call ndynkk(sddyna, 'DEPABS', depabs)
        call ndynkk(sddyna, 'VITABS', vitabs)
        call ndynkk(sddyna, 'ACCABS', accabs)
        call vtcreb(depabs, 'V', 'R', nume_ddlz=numeDof)
        call vtcreb(vitabs, 'V', 'R', nume_ddlz=numeDof)
        call vtcreb(accabs, 'V', 'R', nume_ddlz=numeDof)
    end if
!
! --- FORCES DE SOL
!
    if (lviss) then
        call nmchex(hval_veasse, 'VEASSE', 'CNVISS', cnviss)
        call vtcreb(cnviss, 'V', 'R', nume_ddlz=numeDof)
    end if

! --- SECOND MEMBRE
    call nmchex(hval_veasse, 'VEASSE', 'CNFEDO', cnfedo)
    call vtcreb(cnfedo, 'V', 'R', nume_ddlz=numeDof)
    call nmchex(hval_veasse, 'VEASSE', 'CNFSDO', cnfsdo)
    call vtcreb(cnfsdo, 'V', 'R', nume_ddlz=numeDof)
    call nmchex(hval_veasse, 'VEASSE', 'CNDIDO', cndido)
    if (ldidi) then
        call nmchex(hval_veasse, 'VEASSE', 'CNDIDI', cndidi)
        call vtcreb(cndidi, 'V', 'R', nume_ddlz=numeDof)
    end if
    if (lpilo) then
        call nmchex(hval_veasse, 'VEASSE', 'CNFEPI', cnfepi)
        call vtcreb(cnfepi, 'V', 'R', nume_ddlz=numeDof)
        call nmchex(hval_veasse, 'VEASSE', 'CNDIPI', cndipi)
        call vtcreb(cndipi, 'V', 'R', nume_ddlz=numeDof)
    end if
!
! --- PAS VRAIMENT DES VECT_ELEM MAIS DES CHAM_NO A CREER
!
    call nmchex(hval_veasse, 'VEASSE', 'CNFEXT', cnfext)
    call vtcreb(cnfext, 'V', 'R', nume_ddlz=numeDof)
    if (ldyna) then
        call nmchex(hval_veasse, 'VEASSE', 'CNDYNA', cndyna)
        call vtcreb(cndyna, 'V', 'R', nume_ddlz=numeDof)
        call nmchex(hval_veasse, 'VEASSE', 'CNHYST', cnhyst)
        call vtcreb(cnhyst, 'V', 'R', nume_ddlz=numeDof)
        if (lmpas) then
            call ndynkk(sddyna, 'OLDP_CNFEDO', cnfedo)
            call vtcreb(cnfedo, 'V', 'R', nume_ddlz=numeDof)
            call ndynkk(sddyna, 'OLDP_CNFSDO', cnfsdo)
            call vtcreb(cnfsdo, 'V', 'R', nume_ddlz=numeDof)
            call ndynkk(sddyna, 'OLDP_CNDIDO', cndido)
            call vtcreb(cndido, 'V', 'R', nume_ddlz=numeDof)
            call ndynkk(sddyna, 'OLDP_CNDIDI', cndidi)
            call vtcreb(cndidi, 'V', 'R', nume_ddlz=numeDof)
            call ndynkk(sddyna, 'OLDP_CNFINT', cnfint)
            call vtcreb(cnfint, 'V', 'R', nume_ddlz=numeDof)
            call ndynkk(sddyna, 'OLDP_CNONDP', cnondp)
            call vtcreb(cnondp, 'V', 'R', nume_ddlz=numeDof)
            call ndynkk(sddyna, 'OLDP_CNSSTF', cnsstf)
            call vtcreb(cnsstf, 'V', 'R', nume_ddlz=numeDof)
            call ndynkk(sddyna, 'OLDP_CNCINE', cncine)
            call vtcreb(cncine, 'V', 'R', nume_ddlz=numeDof)
            call ndynkk(sddyna, 'OLDP_CNVISS', cnviss)
            call vtcreb(cnviss, 'V', 'R', nume_ddlz=numeDof)
            call ndynkk(sddyna, 'OLDP_CNHYST', cnhyst)
            call vtcreb(cnhyst, 'V', 'R', nume_ddlz=numeDof)
            call ndynkk(sddyna, 'OLDP_CNSSTR', cnsstr)
            call vtcreb(cnsstr, 'V', 'R', nume_ddlz=numeDof)
            call ndynkk(sddyna, 'OLDP_CNELTC', cneltc)
            call vtcreb(cneltc, 'V', 'R', nume_ddlz=numeDof)
            call ndynkk(sddyna, 'OLDP_CNELTF', cneltf)
            call vtcreb(cneltf, 'V', 'R', nume_ddlz=numeDof)
        end if
    end if

! --- FORCES ISSUES DES MACRO-ELEMENTS STATIQUES
    if (lSuperElement) then
        call nmchex(hval_veasse, 'VEASSE', 'CNSSTR', cnsstr)
        call vtcreb(cnsstr, 'V', 'R', nume_ddlz=numeDof)
    end if

! --- CALCUL PAR SOUS-STRUCTURATION
    if (lsstf) then
        call nmchex(hval_veasse, 'VEASSE', 'CNSSTF', cnsstf)
        call vtcreb(cnsstf, 'V', 'R', nume_ddlz=numeDof)
    end if

! - Modal damping
    if (lDampModal) then
        call nmchex(hval_veasse, 'VEASSE', 'CNAMOD', cnamod)
        call vtcreb(cnamod, 'V', 'R', nume_ddlz=numeDof)
    end if
!
! --- RESIDU DE REFERENCE
!
    if (lrefe) then
        call nmchex(hval_veasse, 'VEASSE', 'CNREFE', cnrefe)
        call vtcreb(cnrefe, 'V', 'R', nume_ddlz=numeDof)
    end if
!
! --- CREATION DE CHAMPS NODAUX PARTAGES (PASSES EN SOUTERRAIN)
!      OBJECTIFS :
!         NE PAS FRAGMENTER LA MEMOIRE
!      REGLES :
!         CNZERO : LECTURE SEULE -> IL VAUT TJRS 0
!         CNTMPX : NE TRANSITENT PAS D'UNE ROUTINE A L'AUTRE
!
    call vtcreb('&&CNPART.ZERO', 'V', 'R', nume_ddlz=numeDof)
    call vtcreb('&&CNPART.CHP1', 'V', 'R', nume_ddlz=numeDof)
    call vtcreb('&&CNPART.CHP2', 'V', 'R', nume_ddlz=numeDof)
    call vtcreb('&&CNPART.CHP3', 'V', 'R', nume_ddlz=numeDof)
    call vtcreb('&&CNREPL.CHP1', 'V', 'R', nume_ddlz=numeDof)
    call vtcreb('&&CNREPL.CHP2', 'V', 'R', nume_ddlz=numeDof)
    call vtcreb('&&CNREPL.CHP3', 'V', 'R', nume_ddlz=numeDof)
    call vtcreb('&&CNREPL.CHP4', 'V', 'R', nume_ddlz=numeDof)
    call vtcreb('&&CNCETA.CHP0', 'V', 'R', nume_ddlz=numeDof)
    call vtcreb('&&CNCETA.CHP1', 'V', 'R', nume_ddlz=numeDof)
    call vtcreb('&&CNCETA.CHP2', 'V', 'R', nume_ddlz=numeDof)
    call vtcreb('&&CNCHAR.FFDO', 'V', 'R', nume_ddlz=numeDof)
    call vtcreb('&&CNCHAR.FFPI', 'V', 'R', nume_ddlz=numeDof)
    call vtcreb('&&CNCHAR.DFDO', 'V', 'R', nume_ddlz=numeDof)
    call vtcreb('&&CNCHAR.DFPI', 'V', 'R', nume_ddlz=numeDof)
    call vtcreb('&&CNCHAR.FVDO', 'V', 'R', nume_ddlz=numeDof)
    call vtcreb('&&CNCHAR.FVDY', 'V', 'R', nume_ddlz=numeDof)
    call vtcreb('&&CNCHAR.DUMM', 'V', 'R', nume_ddlz=numeDof)
    call vtcreb('&&CNCHAR.CINE', 'V', 'R', nume_ddlz=numeDof)
    call vtcreb('&&CNCHAR.DONN', 'V', 'R', nume_ddlz=numeDof)
    call vtcreb('&&CNCHAR.PILO', 'V', 'R', nume_ddlz=numeDof)
    call vtcreb('&&CNCHAR.DIDI', 'V', 'R', nume_ddlz=numeDof)
!
    call jedema()
end subroutine
