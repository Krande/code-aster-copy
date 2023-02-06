! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
subroutine nmmatr(phaseType, listFuncActi, listLoad, numeDof, &
                  sddyna, nlDynaDamping, &
                  numeTime, ds_contact, meelem, measse, &
                  matrAsse)
!
    use NonLin_Datastructure_type
    use NonLinearDyna_type
    use NonLinearDyna_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/ascoma.h"
#include "asterfort/cfdisl.h"
#include "asterfort/detrsd.h"
#include "asterfort/infdbg.h"
#include "asterfort/isfonc.h"
#include "asterfort/lccmst.h"
#include "asterfort/mtcmbl.h"
#include "asterfort/mtdefs.h"
#include "asterfort/ndynlo.h"
#include "asterfort/ndynre.h"
#include "asterfort/nmasfr.h"
#include "asterfort/nmasun.h"
#include "asterfort/nmchex.h"
#include "asterfort/NonLinear_type.h"
#include "asterfort/utmess.h"
!
    integer, intent(in) :: phaseType, listFuncActi(*)
    character(len=19), intent(in) :: listLoad
    character(len=24), intent(in) :: numeDof
    character(len=19), intent(in) :: sddyna
    type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
    integer, intent(in) :: numeTime
    type(NL_DS_Contact), intent(in) :: ds_contact
    character(len=19), intent(in) :: meelem(*), measse(*)
    character(len=19), intent(inout) :: matrAsse
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME - CALCUL)
!
! ASSEMBLAGE DE LA MATRICE GLOBALE
!
! --------------------------------------------------------------------------------------------------
!
! In  phaseType        : name of current phase of algorithm
! In  ds_contact       : datastructure for contact management
! In  sddyna           : name of datastructure for dynamic parameters
! In  nlDynaDamping    : damping parameters
! IN  numeTime : NUMERO D'INSTANT
! IN  numeDof : NOM DE LA NUMEROTATION MECANIQUE
! IN  listLoad : SD LISTE DES CHARGES
! IN  MEASSE : VARIABLE CHAPEAU POUR NOM DES MATR_ASSE
! IN  MEELEM : VARIABLE CHAPEAU POUR NOM DES MATR_ELEM
! IO  MATASS : MATRICE ASSEMBLEE RESULTANTE
!
! --------------------------------------------------------------------------------------------------
!
    integer :: ifm, niv
    aster_logical :: lDyna, lctcd, lexpl, lDampMatrix, l_neum_undead, lShiftMass, lprem, l_cont_lac
    real(kind=8) :: coerig, coeamo, coemas, coeshi
    character(len=8) :: nomddl
    real(kind=8) :: coemat(3)
    character(len=24) :: limat(3)
    character(len=4) :: typcst(3)
    real(kind=8) :: coemam(3)
    character(len=24) :: limam(3)
    character(len=4) :: typcsm(3)
    integer :: nbmat
    character(len=19) :: rigiAsse, rigiMass, rigiDamp
    aster_logical :: lunil, l_unil_pena
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE13_41')
    end if

! - Initializations
    nomddl = ' '

! - Get name of matrices
    call nmchex(measse, 'MEASSE', 'MERIGI', rigiAsse)
    call nmchex(measse, 'MEASSE', 'MEMASS', rigiMass)
    call nmchex(measse, 'MEASSE', 'MEAMOR', rigiDamp)

! - Active functionnalites
    lctcd = isfonc(listFuncActi, 'CONT_DISCRET')
    lunil = isfonc(listFuncActi, 'LIAISON_UNILATER')
    l_neum_undead = isfonc(listFuncActi, 'NEUM_UNDEAD')
    l_cont_lac = isfonc(listFuncActi, 'CONT_LAC')
    lDyna = isfonc(listFuncActi, 'DYNAMIQUE')
    lDampMatrix = nlDynaDamping%hasMatrDamp
    lexpl = ndynlo(sddyna, 'EXPLICITE')
    lShiftMass = ndynlo(sddyna, 'COEF_MASS_SHIFT')
!
! --- PREMIER PAS DE TEMPS ?
!
    lprem = numeTime .le. 1
!
! --- SUPPRESSION ANCIENNE MATRICE ASSEMBLEE
!
    if (lDyna) then
        call detrsd('MATR_ASSE', matrAsse)
    end if
!
! --- COEFFICIENTS POUR MATRICES
!
    if (lDyna) then
        coerig = ndynre(sddyna, 'COEF_MATR_RIGI')
        coeamo = ndynre(sddyna, 'COEF_MATR_AMOR')
        coemas = ndynre(sddyna, 'COEF_MATR_MASS')
        coeshi = ndynre(sddyna, 'COEF_MASS_SHIFT')
    else
        coerig = 1.d0
    end if
    typcst(1) = 'R'
    typcst(2) = 'R'
    typcst(3) = 'R'
!
! --- DECALAGE DE LA MATRICE MASSE (COEF_MASS_SHIFT)
!
    if (lShiftMass .and. lprem .and. (phaseType .eq. PRED_EULER)) then
        typcsm(1) = 'R'
        typcsm(2) = 'R'
        coemam(1) = 1.d0
        coemam(2) = coeshi
        limam(1) = rigiMass
        limam(2) = rigiAsse
        if (lexpl) then
            call mtcmbl(2, typcsm, coemam, limam, rigiMass, &
                        ' ', ' ', 'ELIM=')
        else
            call mtcmbl(2, typcsm, coemam, limam, rigiMass, &
                        'LAGR', ' ', 'ELIM=')
        end if
    end if
!
! --- MATRICES ET COEFFICIENTS
!
    if (ldyna) then
        if (phaseType .eq. ACCEL_INIT) then
            limat(1) = rigiMass
            nbmat = 1
            coemat(1) = 1.d0
        else
            if (lexpl) then
                limat(1) = rigiMass
                nbmat = 1
                coemat(1) = coemas
            else
                coemat(1) = coerig
                coemat(2) = coemas
                coemat(3) = coeamo
                limat(1) = rigiAsse
                limat(2) = rigiMass
                limat(3) = rigiDamp
                if (lDampMatrix) then
                    nbmat = 3
                else
                    nbmat = 2
                end if
            end if
        end if
    end if
!
! --- DEFINITION DE LA STRUCTURE DE LA MATRICE
!
    if (ldyna) then
        if (phaseType .eq. ACCEL_INIT) then
            call mtdefs(matrAsse, rigiMass, 'V', 'R')
        else
            if (lexpl) then
                call mtdefs(matrAsse, rigiMass, 'V', 'R')
            else
                call mtdefs(matrAsse, rigiAsse, 'V', 'R')
            end if
        end if
    end if
!
! --- ASSEMBLAGE
!
    if (lDyna) then
        call mtcmbl(nbmat, typcst, coemat, limat, matrAsse, &
                    nomddl, ' ', 'ELIM=')
    else
        matrAsse = rigiAsse
    end if
    if (phaseType .eq. ACCEL_INIT) then
        goto 999
    end if
!
! --- PRISE EN COMPTE DE LA MATRICE TANGENTE DES FORCES SUIVEUSES
!
    if (l_neum_undead) then
        call ascoma(meelem, numeDof, listLoad, matrAsse)
    end if
!
! --- PRISE EN COMPTE DE LA MATRICE TANGENTE DU FROTTEMENT
!
    if (lctcd .and. (phaseType .eq. CORR_NEWTON)) then
        call nmasfr(ds_contact, matrAsse)
    end if
!
! - Special post-treatment for LAC contact method
!
    if (l_cont_lac) then
        call lccmst(ds_contact, matrAsse)
    end if
!
! --- PRISE EN COMPTE DE LA MATRICE TANGENTE DE PENALISATION
! --- AVEC LES LIAISONS UNILATERALES
!
    if (lunil .and. (phaseType .eq. CORR_NEWTON)) then
        l_unil_pena = cfdisl(ds_contact%sdcont_defi, 'UNIL_PENA')
        if (l_unil_pena) then
            call nmasun(ds_contact, matrAsse)
        end if
    end if
!
999 continue
!
end subroutine
