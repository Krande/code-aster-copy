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
subroutine nmflal(optionSpec, ds_posttimestep, &
                  mod45, l_hpp, &
                  nbFreq, coefDimSpace, matrType, optionModal, bande, &
                  nbDofExcl, nbDofStab, lModiRigi, calcLevel)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
!
    character(len=16), intent(in) :: optionSpec
    type(NL_DS_PostTimeStep), intent(in) :: ds_posttimestep
    character(len=16), intent(out) :: optionModal
    character(len=4), intent(out) :: mod45
    integer(kind=8), intent(out) :: nbFreq, nbDofExcl, nbDofStab, coefDimSpace
    character(len=16), intent(out) :: matrType
    aster_logical, intent(out) :: lModiRigi
    real(kind=8), intent(out) :: bande(2)
    aster_logical, intent(out) :: l_hpp
    character(len=16), intent(out) :: calcLevel
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME - CALCUL DE MODES)
!
! LECTURE DES OPTIONS DANS MECA_NON_LINE
!
! --------------------------------------------------------------------------------------------------
!
! In  optionSpec       : option to compute (FLAMBSTA/FLAMBDYN/VIBRDYNA)
! In  ds_posttimestep  : datastructure for post-treatment at each time step
! OUT MOD45  : TYPE DE CALCUL DE MODES PROPRES
!              'VIBR'     MODES VIBRATOIRES
!              'FLAM'     MODES DE FLAMBEMENT
!              'STAB'     MODE DE STABILITE
! OUT NFREQ  : NOMBRE DE FREQUENCES A CALCULER
! OUT CDSP   : COEFFICIENT MULTIPLICATEUR DE NFREQ -> DIM_SPACE
! Out matrType         : type of matrix
! OUT OPTMOD : OPTION DE RECHERCHE DE MODES
!               'PLUS_PETITE' LA PLUS PETITE FREQUENCE
!               'BANDE'       DANS UNE BANDE DE FREQUENCE DONNEE
! OUT l_hpp   : TYPE DE DEFORMATIONS
!                0            PETITES DEFORMATIONS (MATR. GEOM.)
!                1            GRANDES DEFORMATIONS (PAS DE MATR. GEOM.)
! OUT BANDE  : BANDE DE FREQUENCE SI OPTMOD='BANDE'
! OUT NDDLE  : NOMBRE DE DDL EXCLUS
! OUT DDLEXC : NOM DE L'OBJET JEVEUX CONTENANT LE NOM DES DDLS EXCLUS
! OUT NSTA   : NOMBRE DE DDL STAB
! OUT DDLSTA : NOM DE L'OBJET JEVEUX CONTENANT LE NOM DES DDLS STAB
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: optrig
!
! --------------------------------------------------------------------------------------------------
!
    bande(1) = 1.d-5
    bande(2) = 1.d5
    nbFreq = 0
    coefDimSpace = 0
    nbDofExcl = 0
    mod45 = ' '
    optionModal = ' '
    optrig = ' '
    matrType = ' '
    lModiRigi = ASTER_FALSE
    nbDofStab = 0
    calcLevel = 'TOUT'
    l_hpp = ASTER_TRUE

! - Get parameters
    if (optionSpec .eq. 'VIBRDYNA') then
        nbFreq = ds_posttimestep%mode_vibr%nb_eigen
        coefDimSpace = ds_posttimestep%mode_vibr%coef_dim_espace
        matrType = ds_posttimestep%mode_vibr%type_matr_rigi
        if (ds_posttimestep%mode_vibr%l_small) then
            optionModal = 'PLUS_PETITE'
        else
            optionModal = 'BANDE'
        end if
        bande = ds_posttimestep%mode_vibr%strip_bounds
        if (ds_posttimestep%mode_vibr%level .eq. 'CALIBRATION') then
            calcLevel = 'CALIBRATION'
        end if
        mod45 = 'VIBR'

    else if (optionSpec(1:5) .eq. 'FLAMB') then
        nbFreq = ds_posttimestep%crit_stab%nb_eigen
        coefDimSpace = ds_posttimestep%crit_stab%coef_dim_espace
        matrType = ds_posttimestep%crit_stab%type_matr_rigi
        if (ds_posttimestep%crit_stab%l_small) then
            optionModal = 'PLUS_PETITE'
        else
            optionModal = 'BANDE'
        end if
        bande = ds_posttimestep%crit_stab%strip_bounds
        if (ds_posttimestep%crit_stab%level .eq. 'CALIBRATION') then
            calcLevel = 'CALIBRATION'
        end if
        mod45 = 'FLAM'
        nbDofExcl = ds_posttimestep%stab_para%nb_dof_excl
        nbDofStab = ds_posttimestep%stab_para%nb_dof_stab
        lModiRigi = ds_posttimestep%stab_para%l_modi_rigi
        l_hpp = ds_posttimestep%l_hpp

    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
