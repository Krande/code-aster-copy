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
subroutine nmflma(matrType, mod45, &
                  l_hpp, lModiRigi, &
                  listFuncActi, ds_algopara, &
                  modelZ, caraElem, &
                  ds_material, ds_constitutive, &
                  sddyna, listLoad, &
                  sddisc, numeTime, &
                  ds_posttimestep, nbDofExcl, &
                  hval_incr, hval_algo, &
                  numeDof, ds_system, &
                  ds_measure, hval_meelem, &
                  matrAsse, matrGeom)
!
    use NonLin_Datastructure_type
    use NonLinearElem_module, only: elemSuper
    use NonLinear_module, only: updateLoadBCMatrix, compElemGeom, getOption
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/ascoma.h"
#include "asterfort/asmama.h"
#include "asterfort/asmari.h"
#include "asterfort/asmatr.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/infdbg.h"
#include "asterfort/isfonc.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/matide.h"
#include "asterfort/matr_asse_syme.h"
#include "asterfort/nmcha0.h"
#include "asterfort/nmchai.h"
#include "asterfort/nmchcp.h"
#include "asterfort/nmchex.h"
#include "asterfort/nmrigi.h"
#include "asterfort/NonLinear_type.h"
#include "asterfort/utmess.h"
!
    character(len=16), intent(in) :: matrType
    character(len=4), intent(in) :: mod45
    aster_logical, intent(in) :: l_hpp, lModiRigi
    integer(kind=8), intent(in) :: listFuncActi(*)
    type(NL_DS_AlgoPara), intent(in) :: ds_algopara
    character(len=*), intent(in) :: modelZ
    character(len=24), intent(in) :: caraElem
    type(NL_DS_Material), intent(in) :: ds_material
    type(NL_DS_Constitutive), intent(in) :: ds_constitutive
    character(len=19), intent(in) :: sddyna, listLoad
    character(len=19), intent(in) :: sddisc
    integer(kind=8), intent(in) :: numeTime
    type(NL_DS_PostTimeStep), intent(in) :: ds_posttimestep
    integer(kind=8), intent(in) :: nbDofExcl
    character(len=19), intent(in) :: hval_algo(*), hval_incr(*)
    character(len=24), intent(in) :: numeDof
    type(NL_DS_System), intent(in) :: ds_system
    type(NL_DS_Measure), intent(inout) :: ds_measure
    character(len=19), intent(in) :: hval_meelem(*)
    character(len=19), intent(out) :: matrAsse, matrGeom
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (CALCUL - UTILITAIRE)
!
! CALCUL DE LA MATRICE GLOBALE POUR FLAMBEMENT/MODES VIBRATOIRES
!
! --------------------------------------------------------------------------------------------------
!
! In  matrType         : type of matrix
! IN  MOD45  : TYPE DE CALCUL DE MODES PROPRES
!              'VIBR'     MODES VIBRATOIRES
!              'FLAM'     MODES DE FLAMBEMENT
! IN  l_hpp   : TYPE DE DEFORMATIONS
!                0        PETITES DEFORMATIONS (MATR. GEOM.)
!                1        GRANDES DEFORMATIONS (PAS DE MATR. GEOM.)
! IN  MODELE : MODELE
! IN  NUMEDD : NUME_DDL (VARIABLE AU COURS DU CALCUL)
! IN  NUMFIX : NUME_DDL (FIXE AU COURS DU CALCUL)
! In  ds_material      : datastructure for material parameters
! IN  CARELE : CARACTERISTIQUES DES ELEMENTS DE STRUCTURE
! In  ds_constitutive  : datastructure for constitutive laws management
! IN  LISCHA : LISTE DES CHARGES
! IO  ds_measure       : datastructure for measure and statistics management
! IN  SDDYNA : SD POUR LA DYNAMIQUE
! In  ds_algopara      : datastructure for algorithm parameters
! In  ds_system        : datastructure for non-linear system management
! IN  SDDISC : SD DISC_INST
! IN  PREMIE : SI PREMIER INSTANT DE CALCUL
! IN  NUMINS : NUMERO D'INSTANT
! IN  VALINC : VARIABLE CHAPEAU POUR INCREMENTS VARIABLES
! IN  SOLALG : VARIABLE CHAPEAU POUR INCREMENTS SOLUTIONS
! IN  MEELEM : MATRICES ELEMENTAIRES
! IN  MEASSE : MATRICE ASSEMBLEE
! IN  NDDLE  : NOMBRE DE DDL A EXCLURE
! In  ds_posttimestep  : datastructure for post-treatment at each time step
! IN  MODRIG : MODIFICATION OU NON DE LA RAIDEUR
! OUT MATASS : MATRICE DE RAIDEUR ASSEMBLEE GLOBALE
! OUT MATGEO : DEUXIEME MATRICE ASSEMBLEE POUR LE PROBLEME AUX VP :
!              - MATRICE GEOMETRIQUE GLOBALE (CAS FLAMBEMENT)
!              - MATRICE DE RAIDEUR ASSEMBLEE GLOBALE (CAS FLAMBEMENT)
!              - MATRICE DE MASSE (CAS MODES DYNAMIQUES)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: phaseType = POST_BUCKLING
    integer(kind=8), parameter :: iterNewtPred = 0
    integer(kind=8), parameter :: zvalin = 28
    character(len=16), parameter :: modlag = 'MODI_LAGR_OUI'
    character(len=8), parameter :: tdiag = 'MAX_ABS'
    character(len=19), parameter :: matrRigiSyme = '&&NMFLMA.RIGISYME'
    integer(kind=8) :: ifm, niv
    aster_logical :: l_update_matr
    aster_logical :: lRigiCompute, lSuperElement, lNeumUndead
    character(len=16) :: nonLinearOption
    integer(kind=8) :: reincr
    character(len=8) :: answer
    character(len=19) :: massElem, geomElem
    character(len=19) :: depplu, vitplu, accplu, sigplu, varplu, hval_incrCopy(zvalin)
    integer(kind=8) :: nmax, ldccvg
    character(len=24) :: superElem
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE13_69')
    end if

! - Initializations
    matrAsse = "&&NMFLAM.MATASS"
    matrGeom = "&&NMFLAM.MAGEOM"

! - Copy hat-variable
    call nmchai('VALINC', 'LONMAX', nmax)
    ASSERT(nmax .eq. zvalin)
    call nmcha0('VALINC', 'ALLINI', ' ', hval_incrCopy)
    call nmchcp('VALINC', hval_incr, hval_incrCopy)

! - Active functionnalities
    lNeumUndead = isfonc(listFuncActi, 'NEUM_UNDEAD')
    lSuperElement = isfonc(listFuncActi, 'MACR_ELEM_STAT')

! - Extract hat-variables
    call nmchex(hval_meelem, 'MEELEM', 'MEMASS', massElem)
    call nmchex(hval_meelem, 'MEELEM', 'MEGEOM', geomElem)

! --- ON UTILISE CHAMP EN T+ PAS T- CAR DEPDEL=0 ET
! --- MATRICE = RIGI_MECA_TANG (VOIR FICHE 18779)
! --- SAUF VARI. COMM.
    call nmchex(hval_incr, 'VALINC', 'DEPPLU', depplu)
    call nmchex(hval_incr, 'VALINC', 'SIGPLU', sigplu)
    call nmchex(hval_incr, 'VALINC', 'VARPLU', varplu)
    call nmchex(hval_incr, 'VALINC', 'VITPLU', vitplu)
    call nmchex(hval_incr, 'VALINC', 'ACCPLU', accplu)
    call nmcha0('VALINC', 'DEPMOI', depplu, hval_incrCopy)
    call nmcha0('VALINC', 'VITMOI', vitplu, hval_incrCopy)
    call nmcha0('VALINC', 'VARMOI', varplu, hval_incrCopy)
    call nmcha0('VALINC', 'ACCMOI', accplu, hval_incrCopy)
    call nmcha0('VALINC', 'SIGMOI', sigplu, hval_incrCopy)

! - Get parameter
    reincr = ds_algopara%reac_incr

! - Update global matrix ?
    if ((reincr .eq. 0) .and. (numeTime .ne. 1)) then
        l_update_matr = ASTER_FALSE
    end if
    if (numeTime .eq. 1) then
        l_update_matr = ASTER_TRUE
    end if
    if ((reincr .ne. 0) .and. (numeTime .ne. 1)) then
        l_update_matr = mod(numeTime-1, reincr) .eq. 0
    end if
    lRigiCompute = l_update_matr

! - Select non-linear option for compute matrices
    call getOption(phaseType, listFuncActi, matrType, nonLinearOption)

! - Compute rigidity elementary matrices
    if (lRigiCompute) then
        call nmrigi(modelZ, caraElem, &
                    ds_material, ds_constitutive, &
                    listFuncActi, iterNewtPred, &
                    sddyna, ds_measure, ds_system, &
                    hval_incrCopy, hval_algo, &
                    nonLinearOption, ldccvg)
        ASSERT(ldccvg .eq. 0)
    end if

! - Update elementary matrices for loads and boundary conditions (undead cases)
    call updateLoadBCMatrix(listFuncActi, listLoad, &
                            sddisc, numeTime, &
                            modelZ, caraElem, &
                            ds_material, ds_constitutive, &
                            hval_incrCopy, hval_algo, &
                            hval_meelem)

! - Compute elementary matrices for geometry (buckling)
    if (mod45 .eq. 'FLAM') then
        if (l_hpp) then
            call compElemGeom(modelZ, caraElem, &
                              ds_material, &
                              hval_incrCopy, hval_meelem)
        end if
    end if

! - Compute elementary matrices for super-elements
    if (lSuperElement) then
        call nmchex(hval_meelem, 'MEELEM', 'MESSTR', superElem)
        call elemSuper(modelZ, superElem)
    end if

! - Assemble rigidity matrix
    call asmari(ds_system, hval_meelem, listLoad, matrRigiSyme)
    matrAsse = matrRigiSyme

! - Undead loads
    if (l_update_matr) then
        if (lNeumUndead) then
            call ascoma(hval_meelem, numeDof, listLoad, matrAsse)
        end if
    end if

! - MODIFICATION EVENTUELLE DE LA MATRICE DE RAIDEUR
    if ((nbDofExcl .ne. 0) .and. lModiRigi) then
        call matide(matrAsse, nbDofExcl, ds_posttimestep%stab_para%list_dof_excl, &
                    modlag, tdiag, 10.d0)
    end if

! - Geometry matrix for HPP case
    if (mod45 .eq. 'FLAM') then
        if (l_hpp) then
            call asmatr(1, geomElem, ' ', numeDof, &
                        listLoad, 'ZERO', 'V', 1, matrGeom)
            if ((nbDofExcl .ne. 0) .and. lModiRigi) then
                call matide(matrGeom, nbDofExcl, ds_posttimestep%stab_para%list_dof_excl, &
                            modlag, tdiag, 10.d0)
            end if
        else
            matrGeom = matrAsse
        end if
    else if (mod45 .eq. 'VIBR') then
        call asmama(massElem, ' ', numeDof, listLoad, matrGeom)
    else
        ASSERT(ASTER_FALSE)
    end if

! --- VERIFICATION POUR MODE_VIBR QUE LES DEUX MATRICES SONT SYMETRIQUES
    if (mod45 .eq. 'VIBR') then
        call dismoi('TYPE_MATRICE', matrAsse, 'MATR_ASSE', repk=answer)
        if (answer .eq. 'NON_SYM') then
            call matr_asse_syme(matrAsse)
            call utmess('A', 'MECANONLINE5_56')
        else
            call dismoi('TYPE_MATRICE', matrGeom, 'MATR_ASSE', repk=answer)
            if (answer .eq. 'NON_SYM') then
                call utmess('F', 'MECANONLINE5_56')
            end if
        end if
    end if
!
    call jedema()
!
end subroutine
