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
! aslint: disable=W1504
!
subroutine nmdesc(modele, numedd, &
                  numfix, ds_material, carele, &
                  ds_constitutive, lischa, ds_contact, &
                  ds_algopara, ds_system, solveu, &
                  fonact, numins, iterat, &
                  sddisc, ds_print, ds_measure, &
                  ds_algorom, sddyna, nlDynaDamping, sdnume, &
                  sderro, matass, maprec, &
                  valinc, solalg, meelem, &
                  measse, veasse, lerrit)
!
    use NonLin_Datastructure_type
    use NonLinearDyna_type
    use ROM_Datastructure_type
    use HHO_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/copisd.h"
#include "asterfort/utmess.h"
#include "asterfort/infdbg.h"
#include "asterfort/nonlinDSPrintSepLine.h"
#include "asterfort/nmacin.h"
#include "asterfort/nmassc.h"
#include "asterfort/nmchex.h"
#include "asterfort/nmcoma.h"
#include "asterfort/nmcret.h"
#include "asterfort/nmltev.h"
#include "asterfort/nmresd.h"
#include "asterfort/vtzero.h"
!
    integer(kind=8) :: numins, iterat
    type(NL_DS_AlgoPara), intent(in) :: ds_algopara
    character(len=19) :: matass, maprec
    type(NL_DS_Measure), intent(inout) :: ds_measure
    character(len=19) :: solveu, sddisc, sdnume
    character(len=19), intent(in) :: lischa, sddyna
    type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
    character(len=24) :: numedd, numfix
    character(len=24) :: modele, carele
    type(NL_DS_Material), intent(in) :: ds_material
    type(NL_DS_Constitutive), intent(in) :: ds_constitutive
    type(NL_DS_Contact), intent(inout) :: ds_contact
    type(ROM_DS_AlgoPara), intent(in) :: ds_algorom
    type(NL_DS_System), intent(in) :: ds_system
    character(len=24) :: sderro
    integer(kind=8) :: fonact(*)
    character(len=19) :: meelem(*)
    character(len=19) :: solalg(*), valinc(*)
    character(len=19) :: measse(*), veasse(*)
    type(NL_DS_Print), intent(inout) :: ds_print
    aster_logical :: lerrit
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME)
!
! CALCUL DE LA DIRECTION DE DESCENTE
!
! --------------------------------------------------------------------------------------------------
!
! IN  MODELE : MODELE
! IN  NUMEDD : NUME_DDL (VARIABLE AU COURS DU CALCUL)
! IN  NUMFIX : NUME_DDL (FIXE AU COURS DU CALCUL)
! In  ds_material      : datastructure for material parameters
! IN  CARELE : CARACTERISTIQUES DES ELEMENTS DE STRUCTURE
! In  ds_constitutive  : datastructure for constitutive laws management
! IN  LISCHA : L_CHARGES
! IN  SDDISC : SD DISCRETISATION TEMPORELLE
! IO  ds_print         : datastructure for printing parameters
! IO  ds_measure       : datastructure for measure and statistics management
! In  ds_system        : datastructure for non-linear system management
! IN  SDERRO : SD GESTION DES ERREURS
! IN  SDNUME : SD NUMEROTATION
! IN  ITERAT : NUMERO D'ITERATION DE NEWTON
! IN  NUMINS : NUMERO D'INSTANT
! In  sddyna           : name of datastructure for dynamic parameters
! In  nlDynaDamping    : damping parameters
! IO  ds_contact       : datastructure for contact management
! In  ds_algopara      : datastructure for algorithm parameters
! In  ds_algorom       : datastructure for ROM parameters
! IN  VALINC : VARIABLE CHAPEAU POUR INCREMENTS VARIABLE
! IN  SOLALG : VARIABLE CHAPEAU POUR INCREMENTS SOLUTIONS
! IN  VEASSE : VARIABLE CHAPEAU POUR NOM DES VECT_ASSE
! IN  MEELEM : VARIABLE CHAPEAU POUR NOM DES MATR_ELEM
! IN  MEASSE : VARIABLE CHAPEAU POUR NOM DES MATR_ASSE
! OUT LERRIT : .TRUE. SI ERREUR PENDANT L'ITERATION
!
! --------------------------------------------------------------------------------------------------
!
    character(len=19) :: cncine, depdel, cndonn, cnpilo, cncind
    integer(kind=8) :: faccvg, rescvg, ldccvg
    real(kind=8) :: r8bid
    integer(kind=8) :: ifm, niv
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call nonlinDSPrintSepLine()
        call utmess('I', 'MECANONLINE13_67')
    end if
!
! --- INITIALISATIONS
!
    cndonn = '&&CNCHAR.DONN'
    cnpilo = '&&CNCHAR.PILO'
    cncind = '&&CNCHAR.CINE'
    call vtzero(cndonn)
    call vtzero(cnpilo)
    call vtzero(cncind)
!
! --- INITIALISATIONS CODES RETOURS
!
    ldccvg = -1
    faccvg = -1
    rescvg = -1
!
! --- DECOMPACTION DES VARIABLES CHAPEAUX
!
    call nmchex(veasse, 'VEASSE', 'CNCINE', cncine)
    call nmchex(solalg, 'SOLALG', 'DEPDEL', depdel)
!
! --- CALCUL DE LA MATRICE GLOBALE
!
    call nmcoma(fonact, &
                modele, carele, &
                ds_material, ds_constitutive, &
                lischa, sddyna, nlDynaDamping, &
                sddisc, numins, iterat, &
                ds_algopara, ds_contact, ds_algorom, &
                ds_print, ds_measure, &
                valinc, solalg, &
                meelem, measse, &
                numedd, numfix, sdnume, &
                solveu, ds_system, &
                maprec, matass, &
                faccvg, ldccvg)
!
! --- ERREUR SANS POSSIBILITE DE CONTINUER
!
    if ((faccvg .eq. 1) .or. (faccvg .eq. 2) .or. (ldccvg .eq. 1)) then
        goto 999
    end if

! - Evaluate second member for correction
    call nmassc(fonact, &
                sddyna, nlDynaDamping, &
                ds_contact, veasse, ds_system, &
                cnpilo, cndonn)
!
! --- ACTUALISATION DES CL CINEMATIQUES
!
    call copisd('CHAMP_GD', 'V', cncine, cncind)
    call nmacin(fonact, matass, depdel, cncind)
!
! --- RESOLUTION
!
    call nmresd(fonact, sddyna, ds_measure, solveu, numedd, &
                r8bid, maprec, matass, cndonn, cnpilo, &
                cncind, solalg, rescvg, ds_algorom)
!
999 continue
!
! --- TRANSFORMATION DES CODES RETOURS EN EVENEMENTS
!
    call nmcret(sderro, 'LDC', ldccvg)
    call nmcret(sderro, 'FAC', faccvg)
    call nmcret(sderro, 'RES', rescvg)
!
! --- EVENEMENT ERREUR ACTIVE ?
!
    call nmltev(sderro, 'ERRI', 'NEWT', lerrit)
!
end subroutine
