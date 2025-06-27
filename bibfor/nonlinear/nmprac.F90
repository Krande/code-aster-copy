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
subroutine nmprac(listFuncActi, listLoad, numeDof, solveu, &
                  sddyna, nlDynaDamping, &
                  ds_measure, ds_contact, &
                  meelem, measse, &
                  maprec, matrAsse, &
                  faccvg)
!
    use NonLin_Datastructure_type
    use NonLinearDyna_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/asmama.h"
#include "asterfort/dismoi.h"
#include "asterfort/infdbg.h"
#include "asterfort/isfonc.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/mtdsc2.h"
#include "asterfort/nmchex.h"
#include "asterfort/nmmatr.h"
#include "asterfort/nmrinc.h"
#include "asterfort/nmtime.h"
#include "asterfort/NonLinear_type.h"
#include "asterfort/preres.h"
#include "asterfort/utmess.h"
!
    integer(kind=8), intent(in) :: listFuncActi(*)
    character(len=19), intent(in) :: listLoad
    character(len=24), intent(in) :: numeDof
    character(len=19), intent(in) :: sddyna
    type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
    type(NL_DS_Measure), intent(inout) :: ds_measure
    character(len=19), intent(in) :: solveu
    character(len=19), intent(in) :: meelem(*), measse(*)
    type(NL_DS_Contact), intent(in) :: ds_contact
    character(len=19), intent(in) :: maprec
    character(len=19), intent(inout) :: matrAsse
    integer(kind=8), intent(out) :: faccvg
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (CALCUL - UTILITAIRE)
!
! CALCUL DE LA MATRICE GLOBALE ACCELERATION INITIALE
!
! --------------------------------------------------------------------------------------------------
!
! IN  NUMEDD : NUME_DDL (VARIABLE AU COURS DU CALCUL)
! IN  LISCHA : LISTE DES CHARGES
! In  ds_contact       : datastructure for contact management
! IO  ds_measure       : datastructure for measure and statistics management
! In  sddyna           : name of datastructure for dynamic parameters
! In  nlDynaDamping    : damping parameters
! IN  SOLVEU : SOLVEUR
! IN  MEELEM : VARIABLE CHAPEAU POUR NOM DES MATR_ELEM
! IN  MEASSE : VARIABLE CHAPEAU POUR NOM DES MATR_ASSE
! OUT MATASS : MATRICE DE RESOLUTION ASSEMBLEE
! OUT MAPREC : MATRICE DE RESOLUTION ASSEMBLEE - PRECONDITIONNEMENT
! OUT FACCVG : CODE RETOUR (INDIQUE SI LA MATRICE EST SINGULIERE)
!                   O -> MATRICE INVERSIBLE
!                   1 -> MATRICE SINGULIERE
!                   2 -> MATRICE PRESQUE SINGULIERE
!                   3 -> ON NE SAIT PAS SI LA MATRICE EST SINGULIERE
!
! ----------------------------------------------------------------------
!
    integer(kind=8), parameter :: phaseType = ACCEL_INIT
    integer(kind=8), parameter :: numeTimeFirst = 1
    aster_logical :: lContContinu
    integer(kind=8) :: iEqua, nbEqua
    integer(kind=8) :: iadia, ibid, lres
    character(len=8) :: answer
    integer(kind=8) :: jvalm, zislv1, zislv3
    integer(kind=8) :: ifm, niv
    character(len=19) :: asseMass, elemMass, elemDiri
    integer(kind=8), pointer :: slvi(:) => null()
    character(len=24), pointer :: slvk(:) => null()
!
! ----------------------------------------------------------------------
!
    call jemarq()
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE13_76')
    end if

! - Initializations
    faccvg = -1

! - Active functionnalites
    lContContinu = isfonc(listFuncActi, 'CONT_CONTINU')

! - Get hat-variables
    call nmchex(meelem, 'MEELEM', 'MEMASS', elemMass)
    call nmchex(meelem, 'MEELEM', 'MEDIRI', elemDiri)
    call nmchex(measse, 'MEASSE', 'MEMASS', asseMass)

! - Assemble mass matrix
    call asmama(elemMass, elemDiri, numeDof, listLoad, asseMass)

! - Compute global matrix of system
    call nmmatr(phaseType, listFuncActi, listLoad, numeDof, &
                sddyna, nlDynaDamping, &
                numeTimeFirst, ds_contact, meelem, measse, &
                matrAsse)
!
! --- SI METHODE CONTINUE ON REMPLACE LES TERMES DIAGONAUX NULS PAR
! --- DES UNS POUR POUVOIR INVERSER LA MATRICE ASSEMBLE MATASS
!
    if (lContContinu) then
        call mtdsc2(matrAsse, 'SXDI', 'L', iadia)
        call dismoi('MATR_DISTR', matrAsse, 'MATR_ASSE', repk=answer)
        if (answer .eq. 'OUI') then
            call jeveuo(matrAsse//'.&INT', 'L', lres)
            nbEqua = zi(lres+5)
        else
            call dismoi('NB_EQUA', numeDof, 'NUME_DDL', repi=nbEqua)
        end if
        call jeveuo(jexnum(matrAsse//'.VALM', 1), 'E', jvalm)
        do iEqua = 1, nbEqua
            if (abs(zr(jvalm-1+zi(iadia-1+iEqua))) .le. r8prem()) then
                zr(jvalm-1+zi(iadia-1+iEqua)) = 1.d0
            end if
        end do
    end if

! --- ON ACTIVE LA DETECTION DE SINGULARITE (NPREC=8)
! --- ON EVITE L'ARRET FATAL LORS DE L'INVERSION DE LA MATRICE
    call jeveuo(solveu//'.SLVK', 'L', vk24=slvk)
    if ((slvk(1) (1:5) .eq. 'MUMPS') .or. (slvk(1) (1:4) .eq. 'LDLT') .or. &
        (slvk(1) (1:10) .eq. 'MULT_FRONT')) then
        call jeveuo(solveu//'.SLVI', 'E', vi=slvi)
        zislv1 = slvi(1)
        zislv3 = slvi(3)
        slvi(1) = 8
        slvi(3) = 2
    end if

! --- FACTORISATION DE LA MATRICE ASSEMBLEE GLOBALE
    call nmtime(ds_measure, 'Init', 'Factor')
    call nmtime(ds_measure, 'Launch', 'Factor')
    call preres(solveu, 'V', faccvg, maprec, matrAsse, ibid, -9999)
    call nmtime(ds_measure, 'Stop', 'Factor')
    call nmrinc(ds_measure, 'Factor')

! - RETABLISSEMENT CODE
    if ((slvk(1) (1:5) .eq. 'MUMPS') .or. (slvk(1) (1:4) .eq. 'LDLT') .or. &
        (slvk(1) (1:10) .eq. 'MULT_FRONT')) then
        slvi(1) = zislv1
        slvi(3) = zislv3
    end if
    call jedema()
!
end subroutine
