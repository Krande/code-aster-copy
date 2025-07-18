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
function isfonc(list_func_acti, func_name_z)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
!
    aster_logical :: isfonc
    integer(kind=8), intent(in) :: list_func_acti(*)
    character(len=*), intent(in) :: func_name_z
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Utility
!
! Is functionnality active ?
!
! --------------------------------------------------------------------------------------------------
!
!
! In  list_func_acti : list of active functionnalities
! In  func_name      : name of functionnality asked
!       RECH_LINE          :  line search
!       PILOTAGE           :  PILOTAGE
!       NEWTON_KRYLOV      :  Newton-Krylov method
!       DEBORST            :  De Borst plane stress method
!       SOUS_STRUC         :  CALCUL PAR SOUS-STRUCTURATION
!       PROJ_MODAL         :  PROJECTION MODALE
!       DYNAMIQUE          :  DYNAMIQUE
!       IMPLEX             :  METHODE IMPLEX
!       EXPLICITE          :  METHODE EXPLICITE
!       CONTACT            :  CONTACT DISCRET OU CONTINU OU XFEM OU LAC
!       CONT_DISCRET       :  CONTACT DISCRET
!       CONT_CONTINU       :  CONTACT CONTINU
!       CONT_LAC           :  CONTACT LAC
!       CONTACT_INIT       :  CONTACT INITIAL
!       FROT_DISCRET       :  FROTTEMENT DISCRET
!       FROT_CONTINU       :  FROTTEMENT CONTINU
!       BOUCLE_EXTERNE     :  PRESENCE D'UNE BOUCLE EXTERNE
!       BOUCLE_EXT_GEOM    :  PRESENCE D'UNE BOUCLE EXTERNE POUR
!                             LA GEOMETRIE
!       BOUCLE_EXT_FROT    :  PRESENCE D'UNE BOUCLE EXTERNE POUR
!                             LE FROTTEMENT
!       BOUCLE_EXT_CONT    :  PRESENCE D'UNE BOUCLE EXTERNE POUR
!                             LE CONTACT
!       FROT_NEWTON        :  NEWTON GENERALISE POUR LE CONTACT CONTINU
!                             FROTTEMENT
!       CONT_NEWTON        :  NEWTON GENERALISE POUR LE CONTACT CONTINU
!                             CONTACT
!       GEOM_NEWTON        :  NEWTON GENERALISE POUR LE CONTACT CONTINU
!                             GEOMETRIE
!       EXIS_PENA          :  AT LEAST ONE CONTACT ZONE IS ACTIVATED ON PENALISATION
!       CONT_ALL_VERIF     :  CONTACT SANS CALCUL SUR TOUTES LES ZONES
!       LIAISON_UNILATER   :  LIAISON UNILATERALE
!       ELT_CONTACT        :  ELEMENTS DE CONTACT (CONTINU/XFEM)
!       DIS_CHOC           :  ELEMENTS DIS_CHOC
!       GD_ROTA            :  ELEMENTS DE STRUCTURE EN GRANDES ROTATIONS
!       XFEM               :  ELEMENTS XFEM
!       EXI_STRX           :  ELEMENTS UTILISANT STRX (PMF)
!       RESI_REFE          :  CONVERGENCE PAR RESIDU DE REFERENCE
!       RESI_COMP          :  CONVERGENCE NORME PAR FORCE NODALE CMP
!       RESI_RELA          :  CONVERGENCE PAR RESIDU DE RELATIF GLOBAL
!       RESI_MAXI          :  CONVERGENCE PAR RESIDU DE MAXIMAL GLOBAL
!       DIRI_CINE          :  PRESENCE DE CHARGEMENTS DE DIRICHLET
!                             DE TYPE ELIMINATION (AFFE_CHAR_CINE)
!       NEUM_UNDEAD        :  undead Neumann load
!       DIDI               :  FORCE DE TYPE DIFF. DIRICHLET
!       THM                :  MODELISATION THM
!       HHO                :  MODELISATION HHO
!       MACR_ELEM_STAT     :  MACRO-ELEMENTS STATIQUES
!       ENDO_NO            :  MODELISATION ENDO AUX NOEUDS *_GVNO
!
!       ENERGIE            :  CALCUL DES ENERGIES
!       CRIT_STAB          :  CALCUL DE STABILITE
!       DDL_STAB           :  DDLS IRREVERSIBLES DANS CRIT_STAB
!       MODE_VIBR          :  CALCUL DE MODES VIBRATOIRES
!       ERRE_TEMPS_THM     :  CALCUL DE L'ERREUR EN TEMPS POUR LA THM
!       EXI_VARC           :  PRESENCE DE VARIABLES DE COMMANDES
!       REUSE              :  CONCEPT RE-ENTRANT
!       LDLT               :  SOLVEUR LDLT
!       MULT_FRONT         :  SOLVEUR MULT_FRONT
!       GCPC               :  SOLVEUR GCPC
!       MUMPS              :  SOLVEUR MUMPS
!       PETSC              :  SOLVEUR PETSC
!       LDLT_SP            :  PRECONDITIONNEUR LDLT_SP
!       MATR_DISTRIBUEE    :  MATRICES DISTRIBUEES
!       ELAS_FO            :  elastic properties are functions
!       ANNEALING          :  post-treatment for annealing
!       ETAT_INIT          :  initial state
!       ROM                :  reduced order model
!       HROM               :  hyper-reduced order model
!       HROM_CORR_EF       :  hyper-reduced order model with EF correction
!
! DERNIER NUMERO UTILISE: 71
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24) :: func_name
!
! --------------------------------------------------------------------------------------------------
!
    func_name = func_name_z
!
    if (func_name .eq. 'RECH_LINE') then
        isfonc = list_func_acti(1) .eq. 1
    else if (func_name .eq. 'PILOTAGE') then
        isfonc = list_func_acti(2) .eq. 1
    else if (func_name .eq. 'CONTACT') then
        isfonc = list_func_acti(64) .eq. 1
    else if (func_name .eq. 'LIAISON_UNILATER') then
        isfonc = list_func_acti(12) .eq. 1
    else if (func_name .eq. 'ELT_CONTACT') then
        isfonc = list_func_acti(26) .eq. 1
    else if (func_name .eq. 'CONT_DISCRET') then
        isfonc = list_func_acti(4) .eq. 1
    else if (func_name .eq. 'CONT_CONTINU') then
        isfonc = list_func_acti(5) .eq. 1
    else if (func_name .eq. 'CONT_LAC') then
        isfonc = list_func_acti(63) .eq. 1
    else if (func_name .eq. 'CONTACT_INIT') then
        isfonc = list_func_acti(17) .eq. 1
    else if (func_name .eq. 'DIS_CHOC') then
        isfonc = list_func_acti(29) .eq. 1
    else if (func_name .eq. 'FROT_DISCRET') then
        isfonc = list_func_acti(3) .eq. 1
    else if (func_name .eq. 'FROT_CONTINU') then
        isfonc = list_func_acti(10) .eq. 1
    else if (func_name .eq. 'BOUCLE_EXT_GEOM') then
        isfonc = list_func_acti(31) .eq. 1
    else if (func_name .eq. 'BOUCLE_EXT_FROT') then
        isfonc = list_func_acti(32) .eq. 1
    else if (func_name .eq. 'BOUCLE_EXT_CONT') then
        isfonc = list_func_acti(33) .eq. 1
    else if (func_name .eq. 'BOUCLE_EXTERNE') then
        isfonc = list_func_acti(34) .eq. 1
    else if (func_name .eq. 'GEOM_NEWTON') then
        isfonc = list_func_acti(55) .eq. 1
    else if (func_name .eq. 'EXIS_PENA') then
        isfonc = list_func_acti(66) .eq. 1
    else if (func_name .eq. 'FROT_NEWTON') then
        isfonc = list_func_acti(47) .eq. 1
    else if (func_name .eq. 'CONT_NEWTON') then
        isfonc = list_func_acti(53) .eq. 1
    else if (func_name .eq. 'CONT_ALL_VERIF') then
        isfonc = list_func_acti(38) .eq. 1
!
    else if (func_name .eq. 'XFEM') then
        isfonc = list_func_acti(6) .eq. 1
    else if (func_name .eq. 'DEBORST') then
        isfonc = list_func_acti(7) .eq. 1
!
    else if (func_name .eq. 'RESI_REFE') then
        isfonc = list_func_acti(8) .eq. 1
    else if (func_name .eq. 'RESI_COMP') then
        isfonc = list_func_acti(35) .eq. 1
    else if (func_name .eq. 'RESI_RELA') then
        isfonc = list_func_acti(70) .eq. 1
    else if (func_name .eq. 'RESI_MAXI') then
        isfonc = list_func_acti(71) .eq. 1
    else if (func_name .eq. 'DIRI_CINE') then
        isfonc = list_func_acti(36) .eq. 1
    else if (func_name .eq. 'NEUM_UNDEAD') then
        isfonc = list_func_acti(13) .eq. 1
    else if (func_name .eq. 'DIRI_UNDEAD') then
        isfonc = list_func_acti(60) .eq. 1
    else if (func_name .eq. 'DIDI') then
        isfonc = list_func_acti(22) .eq. 1
!
    else if (func_name .eq. 'MACR_ELEM_STAT') then
        isfonc = list_func_acti(14) .eq. 1
    else if (func_name .eq. 'GD_ROTA') then
        isfonc = list_func_acti(15) .eq. 1
    else if (func_name .eq. 'ENDO_NO') then
        isfonc = list_func_acti(40) .eq. 1
    else if (func_name .eq. 'CRIT_STAB') then
        isfonc = list_func_acti(18) .eq. 1
    else if (func_name .eq. 'DDL_STAB') then
        isfonc = list_func_acti(49) .eq. 1
    else if (func_name .eq. 'MODE_VIBR') then
        isfonc = list_func_acti(19) .eq. 1
    else if (func_name .eq. 'ERRE_TEMPS_THM') then
        isfonc = list_func_acti(21) .eq. 1
    else if (func_name .eq. 'SOUS_STRUC') then
        isfonc = list_func_acti(24) .eq. 1
    else if (func_name .eq. 'IMPLEX') then
        isfonc = list_func_acti(28) .eq. 1
    else if (func_name .eq. 'EXI_VARC') then
        isfonc = list_func_acti(30) .eq. 1
    else if (func_name .eq. 'THM') then
        isfonc = list_func_acti(37) .eq. 1
    else if (func_name .eq. 'HHO') then
        isfonc = list_func_acti(68) .eq. 1
    else if (func_name .eq. 'REUSE') then
        isfonc = list_func_acti(39) .eq. 1
!
    else if (func_name .eq. 'LDLT') then
        isfonc = list_func_acti(41) .eq. 1
    else if (func_name .eq. 'MULT_FRONT') then
        isfonc = list_func_acti(42) .eq. 1
    else if (func_name .eq. 'GCPC') then
        isfonc = list_func_acti(43) .eq. 1
    else if (func_name .eq. 'MUMPS') then
        isfonc = list_func_acti(44) .eq. 1
    else if (func_name .eq. 'PETSC') then
        isfonc = list_func_acti(45) .eq. 1
    else if (func_name .eq. 'LDLT_SP') then
        isfonc = list_func_acti(46) .eq. 1
    else if (func_name .eq. 'NEWTON_KRYLOV') then
        isfonc = list_func_acti(48) .eq. 1
    else if (func_name .eq. 'ROM') then
        isfonc = list_func_acti(61) .eq. 1
    else if (func_name .eq. 'HROM') then
        isfonc = list_func_acti(62) .eq. 1
    else if (func_name .eq. 'HROM_CORR_EF') then
        isfonc = list_func_acti(67) .eq. 1
    else if (func_name .eq. 'ENERGIE') then
        isfonc = list_func_acti(50) .eq. 1
    else if (func_name .eq. 'PROJ_MODAL') then
        isfonc = list_func_acti(51) .eq. 1
    else if (func_name .eq. 'MATR_DISTRIBUEE') then
        isfonc = list_func_acti(52) .eq. 1
    else if (func_name .eq. 'EXPLICITE') then
        isfonc = list_func_acti(54) .eq. 1
    else if (func_name .eq. 'DYNAMIQUE') then
        isfonc = list_func_acti(69) .eq. 1
    else if (func_name .eq. 'EXI_STRX') then
        isfonc = list_func_acti(56) .eq. 1
    else if (func_name .eq. 'ELAS_FO') then
        isfonc = list_func_acti(57) .eq. 1
    else if (func_name .eq. 'ANNEALING') then
        isfonc = list_func_acti(58) .eq. 1
    else if (func_name .eq. 'ETAT_INIT') then
        isfonc = list_func_acti(59) .eq. 1
!
    else
        ASSERT(.false.)
    end if
!
end function
