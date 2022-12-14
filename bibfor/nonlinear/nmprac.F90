! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
subroutine nmprac(fonact, lischa, numedd, solveu     ,&
                  sddyna, ds_measure, ds_contact,&
                  meelem, measse, maprec, matass    , faccvg)
!
use NonLin_Datastructure_type
!
implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
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
#include "asterfort/preres.h"
#include "asterfort/utmess.h"
#include "asterfort/asmama.h"
!
integer :: fonact(*)
character(len=19) :: sddyna, lischa
type(NL_DS_Measure), intent(inout) :: ds_measure
character(len=24) :: numedd
character(len=19) :: solveu
character(len=19) :: meelem(*), measse(*)
type(NL_DS_Contact), intent(in) :: ds_contact
character(len=19) :: maprec, matass
integer :: faccvg
!
! ----------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (CALCUL - UTILITAIRE)
!
! CALCUL DE LA MATRICE GLOBALE ACCELERATION INITIALE
!
! ----------------------------------------------------------------------
!
! IN  NUMEDD : NUME_DDL (VARIABLE AU COURS DU CALCUL)
! IN  LISCHA : LISTE DES CHARGES
! In  ds_contact       : datastructure for contact management
! IO  ds_measure       : datastructure for measure and statistics management
! IN  SDDYNA : SD POUR LA DYNAMIQUE
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
    aster_logical :: lctcc
    integer :: ieq, ibid, numins
    integer :: iadia, neq, lres, neql
    character(len=8) :: kmatd
    integer :: jvalm, zislv1, zislv3
    integer :: ifm, niv
    character(len=19) :: masse, memass, mediri
    integer, pointer :: slvi(:) => null()
!
! ----------------------------------------------------------------------
!
    call jemarq()
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE13_76')
    endif
!
! --- INITIALISATIONS
!
    faccvg = -1
    numins = 1
    call dismoi('NB_EQUA', numedd, 'NUME_DDL', repi=neq)
!
! --- FONCTIONNALITES ACTIVEES
!
    lctcc = isfonc(fonact,'CONT_CONTINU')
!
! --- DECOMPACTION DES VARIABLES CHAPEAUX
!
    call nmchex(meelem, 'MEELEM', 'MEMASS', memass)
    call nmchex(meelem, 'MEELEM', 'MEDIRI', mediri)
    call nmchex(measse, 'MEASSE', 'MEMASS', masse)
!
! --- ASSEMBLAGE DE LA MATRICE MASSE
!
    call asmama(memass, mediri, numedd, lischa, masse)
!
! --- CALCUL DE LA MATRICE ASSEMBLEE GLOBALE
!
    call nmmatr('ACCEL_INIT', fonact    , lischa, numedd, sddyna,&
                numins      , ds_contact, meelem, measse, matass)
!
! --- SI METHODE CONTINUE ON REMPLACE LES TERMES DIAGONAUX NULS PAR
! --- DES UNS POUR POUVOIR INVERSER LA MATRICE ASSEMBLE MATASS
!
    if (lctcc) then
        call mtdsc2(matass, 'SXDI', 'L', iadia)
        call dismoi('MATR_DISTR', matass, 'MATR_ASSE', repk=kmatd)
        if (kmatd .eq. 'OUI') then
            call jeveuo(matass//'.&INT', 'L', lres)
            neql = zi(lres+5)
        else
            neql = neq
        endif
        call jeveuo(jexnum(matass//'.VALM', 1), 'E', jvalm)
        do ieq = 1, neql
            if (zr(jvalm-1+zi(iadia-1+ieq)) .eq. 0.d0) then
                zr(jvalm-1+zi(iadia-1+ieq)) = 1.d0
            endif
        end do
    endif
!
! --- ON ACTIVE LA DETECTION DE SINGULARITE (NPREC=8)
! --- ON EVITE L'ARRET FATAL LORS DE L'INVERSION DE LA MATRICE
!
    call jeveuo(solveu//'.SLVI', 'E', vi=slvi)
    zislv1 = slvi(1)
    zislv3 = slvi(3)
    slvi(1) = 8
    slvi(3) = 2
!
! --- FACTORISATION DE LA MATRICE ASSEMBLEE GLOBALE
!
    call nmtime(ds_measure, 'Init'  , 'Factor')
    call nmtime(ds_measure, 'Launch', 'Factor')
    call preres(solveu, 'V', faccvg, maprec, matass,&
                ibid, -9999)
    call nmtime(ds_measure, 'Stop', 'Factor')
    call nmrinc(ds_measure, 'Factor')
!
! --- RETABLISSEMENT CODE
!
    slvi(1) = zislv1
    slvi(3) = zislv3
!
! --- LA MATRICE PEUT ETRE QUASI-SINGULIERE PAR EXEMPLE POUR LES DKT
!
    if (faccvg .eq. 1) then
        call utmess('A', 'MECANONLINE_78')
    endif
!
    call jedema()
!
end subroutine
