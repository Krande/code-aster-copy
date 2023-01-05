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
! person_in_charge: jean-luc.flejou at edf.fr
!
module te0047_type
!
    implicit none
!
#include "asterf_types.h"
!
    type :: te0047_dscr
        ! Ce qu'il faut sauvegarder en fonction de l'option :
        !       matrice Tangente
        !       vecteur Force
        !       vecteur Sigma
        !       vecteur Variables internes
        !   lMatr       :   FULL_MECA*             RIGI_MECA*
        !   lVect       :   FULL_MECA*  RAPH_MECA  RIGI_MECA_TANG
        !   lSigm       :   FULL_MECA*  RAPH_MECA  RIGI_MECA_TANG
        !   lVari       :   FULL_MECA*  RAPH_MECA
        !   lPred       :                          RIGI_MECA_TANG
        ! Si on fait uniquement une prédiction
        !   lMatrPred   :                          RIGI*
        !
        !   nomte   : nom terme élémentaire
        !   ndim    : dimension de l'espace
        !   nbt     : nombre de terme dans la matrice
        !   nno     : nombre de noeuds de l'élément
        !   nc      : nombre de composante par noeud
        !
        !   ugm     : déplacement moins dans le repère global
        !   dug     : incrément de déplacement dans le repère global
        !
        !   ulm     : déplacement moins dans le repère local
        !   dul     : incrément de déplacement dans le repère local
        !
        !   pgl     : matrice de passage global vers local
        !
        aster_logical       :: lVect = ASTER_FALSE
        aster_logical       :: lMatr = ASTER_FALSE
        aster_logical       :: lVari = ASTER_FALSE
        aster_logical       :: lSigm = ASTER_FALSE
        aster_logical       :: lPred = ASTER_FALSE
        aster_logical       :: lMatrPred = ASTER_FALSE
        aster_logical       :: lTraceDbg = ASTER_FALSE
        !
        character(len=16)   :: option = ''
        character(len=16)   :: nomte = ''
        character(len=16)   :: rela_comp = ''
        character(len=16)   :: defo_comp = ''
        character(len=16)   :: type_comp = ''
        !
        integer             :: ndim
        integer             :: nbt
        integer             :: nno
        integer             :: nc
        !
        real(kind=8)        :: ulm(12), ugm(12)
        real(kind=8)        :: dul(12), dug(12)
        real(kind=8)        :: pgl(3, 3)
        real(kind=8)        :: TempsPlus, TempsMoins
        !
    end type te0047_dscr
!
contains

    subroutine te0047_dscr_write(D)
        type(te0047_dscr), intent(in) :: D
        !
    write (*, '(5(A16))') [character(16) ::'option', 'nomte', 'rela_comp', 'type_comp', 'defo_comp']
        write (*, '(5(A16))') D%option, D%nomte, D%rela_comp, D%type_comp, D%defo_comp
        !
        write (*, '(A)') 'Vect Matr Vari Sigm Pred MPred '
        write (*, '(6(L1,4X))') D%lVect, D%lMatr, D%lVari, D%lSigm, D%lPred, D%lMatrPred
    end

end module te0047_type
