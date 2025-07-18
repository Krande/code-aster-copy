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

subroutine xgrdhm(nomte, ndim, mecani, press1, press2, &
                  tempe, enrmec, dimdef, dimcon, nmec, &
                  np1, np2, nenr, dimenr, enrhyd, nfh)
    implicit none
!
#include "asterfort/teattr.h"
#include "asterfort/assert.h"
    integer(kind=8) :: mecani(5), press1(7), press2(7), tempe(5)
    integer(kind=8) :: dimdef, dimcon, ier, nfh
    integer(kind=8) :: ndim, nmec, np1, np2, i
    character(len=8) :: enr
    character(len=16) :: nomte
!
! DECLARATIONS POUR XFEM
    integer(kind=8) :: enrmec(3), enrhyd(3), dimenr
    integer(kind=8) :: nenr
!
! person_in_charge: daniele.colombo at ifpen.fr
!
!   TABLEAU MECANI :
!   MECANI(1) = 1 : IL Y A UNE EQUATION MECANIQUE
!               0 : SINON
!   MECANI(2) = ADRESSE DANS LES TABLEAUX DES DEFORMATIONS
!               GENERALISEES AU POINT DE GAUSS DES
!               DEFORMATIONS CORRESPONDANT A LA MECANIQUE
!   MECANI(3) = ADRESSE DANS LES TABLEAUX DES CONTRAINTES
!               GENERALISEES AU POINT DE GAUSS DES
!               CONTRAINTES CORRESPONDANT A LA MECANIQUE
!   MECANI(4) = NOMBRE DE DEFORMATIONS MECANIQUES
!   MECANI(5) = NOMBRE DE CONTRAINTES MECANIQUES
!
!   TABLEAU PRESS1 :
!   PRESS1(1) = 1 : IL Y A UNE EQUATION SUR LA PREMIERE PRESSION
!               0 : SINON
!   PRESS1(2) = NOMBRE DE PHASES POUR LE CONSTITUANT 1
!   PRESS1(3) = ADRESSE DANS LES TABLEAUX DES DEFORMATIONS
!               GENERALISEES AU POINT DE GAUSS
!   PRESS1(4) = ADRESSE DANS LES TABLEAUX DES CONTRAINTES
!               GENERALISEES AU POINT DE GAUSS DES CONTRAINTES
!               CORRESPONDANT A LA PREMIERE PHASE DU 1ER CONSTITUANT
!   PRESS1(5) = ADRESSE DANS LES TABLEAUX DES CONTRAINTES
!               GENERALISEES AU POINT DE GAUSS DES CONTRAINTES
!               CORRESPONDANT A LA DEUXIEME PHASE DU 1ER CONSTITUANT
!   PRESS1(6) = NOMBRE DE DEFORMATIONS PRESSION
!   PRESS1(7) = NOMBRE DE CONTRAINTES POUR CHAQUE PHASE DU CONSTITUANT 1
!
!   TABLEAU ENRMEC : (POUR LA MECANIQUE)
!   ENRMEC(1) = 1 : IL Y A ENRICHISSEMENT PAS FONCTION HEAVISIDE
!   ENRMEC(2) = ADRESSE DANS LES TABLEAUX DES DEFORMATIONS
!               GENERALISEES AU POINT DE GAUSS
!   ENRMEC(3) = NOMBRE DE DEFORMATIONS POUR L'ENRICHISSEMENT
!               MECANIQUE (TYPE HEAVISIDE)
!
!   TABLEAU ENRHYD : (POUR L'HYDRAULIQUE)
!   ENRHYD(1) = 1 : IL Y A ENRICHISSEMENT PAR FONCTION HEAVISIDE
!   ENRHYD(2) = ADRESSE DANS LES TABLEAUX DES DEFORMATIONS
!               GENERALISEES AU POINT DE GAUSS
!   ENRHYD(3) = NOMBRE DE DEFORMATIONS POUR L'ENRICHISSEMENT
!               HYDRAULIQUE (TYPE HEAVISIDE)
!
    integer(kind=8) :: iaux
!
!====
! 1. REPERAGE DES CALCULS A FAIRE : MECANIQUE, HYDRAULIQUE, ETC.
!====
!
! =====================================================================
! --- INITIALISATION DES GRANDEURS PRESENTES SELON LA MODELISATION ----
! --- EN HM-XFEM ----------------------------------------------------------
!======================================================================
    do i = 1, 5
        mecani(i) = 0
        tempe(i) = 0
    end do
    do i = 1, 7
        press1(i) = 0
        press2(i) = 0
    end do
    do i = 1, 3
        enrmec(i) = 0
        enrhyd(i) = 0
    end do
    np2 = 0
! =====================================================================
! --- SI MODELISATION = HM --------------------------------------------
! =====================================================================
    call teattr('C', 'MECA', enr, ier, typel=nomte)
    ASSERT(enr .eq. 'OUI')
    call teattr('C', 'THER', enr, ier, typel=nomte)
    ASSERT(enr .eq. 'NON')
    call teattr('C', 'HYDR1', enr, ier, typel=nomte)
    ASSERT(enr .eq. '1')
    call teattr('C', 'HYDR2', enr, ier, typel=nomte)
    ASSERT(enr .eq. '0')

    mecani(1) = 1
    press1(1) = 1
    press1(2) = 1

! =====================================================================
! --- ON VERIFIE LA NATURE DE L'ELEMENT HM-XFEM -----------------------
! =====================================================================
    call teattr('S', 'XFEM', enr, ier, typel=nomte)
    if (enr(1:2) .eq. 'XH') then
        enrmec(1) = 1
        enrhyd(1) = 1
    end if
!
! 2. CALCUL PREALABLE DES ADRESSES LOCALES DES VARIABLES
! =====================================================================
! 2.1. LES AUTRES VALEURS DES TABLEAUX MECANI,PRESS1,PRESS2,TEMPE -----
! --- SE DEFINISSENT AUTOMATIQUEMENT : --------------------------------
! --- NOMBRE DE DEFORMATIONS ET DE CONTRAINTES DE CHAQUE PROBLEME -----
! =====================================================================
!
!  ATTENTION LE NOMBRE DE PARAMETRES MECANIQUES AUGMENTE
!  LE SCALAIRE SIP DEVIENT UN TENSEUR A 6 COMPOSANTES
!  DANS LE CAS ISOTROPE SIP EST UN TENSEUR ISOTROPE
!
    if (mecani(1) .eq. 1) then
        mecani(4) = ndim+6
        mecani(5) = 6+6
        nmec = ndim
    end if
!
    iaux = 1
!
    if (press1(1) .eq. 1) then
        press1(6) = 1+ndim
        press1(7) = iaux+ndim
        np1 = 1
    end if
!
    if (enrmec(1) .eq. 1) then
        enrmec(3) = nfh*ndim
        nenr = ndim*nfh
    end if
!
    if (enrhyd(1) .eq. 1) then
        enrhyd(3) = nfh
    end if
!
! =====================================================================
! 2.2. ADRESSE DES SOUS-TABLEAUX DANS LES DEFORMATIONS PHYSIQUES, LES -
!      DEFORMATIONS GENERALISEES ET LES CONTRAINTES GENERALISEES ------
! =====================================================================
!
! 2.2.1. ==> DEFORMATIONS ET CONTRAINTES EN MECANIQUE
!
    if (mecani(1) .eq. 1) then
        mecani(2) = 1
        mecani(3) = 1
    end if
!
! 2.2.2. ==> DEFORMATIONS ET CONTRAINTES POUR LA PREMIERE PRESSION
!
    if (press1(1) .eq. 1) then
        press1(3) = mecani(4)+1
        press1(4) = mecani(5)+1
    end if
!
! 2.2.3. ==> DEFORMATIONS POUR L'ENRICHISSEMENT HEAVISIDE (MECA)
!
    if (enrmec(1) .eq. 1) then
        enrmec(2) = mecani(4)+press1(6)+1
    end if
!
! 2.2.4. ==> DEFORMATIONS POUR L'ENRICHISSEMENT HEAVISIDE (HYDRO)
!
    if (enrhyd(1) .eq. 1) then
        enrhyd(2) = mecani(4)+press1(6)+nmec+1
    end if
!
! =====================================================================
! 2.3. DIMENSION DES DEFORMATIONS ET CONTRAINTES ----------------------
! =====================================================================
    dimdef = mecani(4)+press1(6)
    dimcon = mecani(5)+press1(2)*press1(7)
! DIMENSION INTERMEDIAIRE UTILISEE POUR L'ASSEMBLAGE EN XFEM
    dimenr = mecani(4)+press1(6)+enrmec(3)+enrhyd(3)
!=====================================================================
!
! =====================================================================
!
end subroutine
