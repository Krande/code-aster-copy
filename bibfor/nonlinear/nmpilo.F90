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
subroutine nmpilo(sdpilo, deltat, rho, solalg, veasse, &
                  modele, ds_material, ds_constitutive, valinc, &
                  nbatte, numedd, nbeffe, eta, pilcvg, &
                  carele)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmchex.h"
#include "asterfort/nmpial.h"
#include "asterfort/nmpidd.h"
#include "asterfort/nmpila.h"
#include "asterfort/nmpipe.h"
!
    integer(kind=8) :: nbatte, nbeffe
    integer(kind=8) :: pilcvg
    real(kind=8) :: deltat, rho, eta(nbatte)
    character(len=19) :: sdpilo
    type(NL_DS_Constitutive), intent(in) :: ds_constitutive
    character(len=19) :: solalg(*), veasse(*), valinc(*)
    character(len=24) :: modele, carele
    type(NL_DS_Material), intent(in) :: ds_material
    character(len=24) :: numedd
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME - PILOTAGE)
!
! RESOLUTION DE L'EQUATION DE PILOTAGE
!
! --------------------------------------------------------------------------------------------------
!
!       DU = DUN + RHO.DU0 + ETA.RHO.DU1
!       AVEC UN CHARGEMENT F0 + ETA_VRAI.F1
!       ET ETA_VRAI = ETA + ETAN * (1-RHO)
!       ETA_MIN < ETA_VRAI < ETA_MAX
!
!
! IN  SDPILO : SD PILOTAGE
! IN  DELTAT : INCREMENT DE TEMPS
! IN  RHO    : VALEUR DU PAS DE RECHERCHE LINEAIRE
! IN  MODELE : MODELE
! IN  NUMEDD : NUME_DDL
! In  ds_material      : datastructure for material parameters
! In  ds_constitutive  : datastructure for constitutive laws management
! IN  CARELE : CARACTERISTIQUES DES ELEMENTS DE STRUCTURE
! IN  VALINC : VARIABLE CHAPEAU POUR INCREMENTS VARIABLES
! IN  SOLALG : VARIABLE CHAPEAU POUR INCREMENTS SOLUTIONS
! IN  VEASSE : VARIABLE CHAPEAU POUR NOM DES VECT_ASSE
! IN  NBATTE : NOMBRE DE SOLUTIONS ATTENDUES
! OUT NBEFFE : NOMBRE DE SOLUTIONS EFFECTIVES
! OUT ETA    : ETA_PILOTAGE
! OUT PILCVG : CODE DE CONVERGENCE POUR LE PILOTAGE
!                -1 : PAS DE CALCUL DU PILOTAGE
!                 0 : CAS DU FONCTIONNEMENT NORMAL
!                 1 : PAS DE SOLUTION
!                 2 : BORNE ATTEINTE -> FIN DU CALCUL
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: neq, i
    real(kind=8) :: dtau, etrmin, etrmax, coef
    character(len=19) :: ddepl0, ddepl1
    character(len=24) :: typpil
    character(len=19) :: ligrpi, cartyp, careta
    character(len=19) :: depmoi, depdel, deppr1, deppr2
    character(len=19) :: cnfepi
    real(kind=8), pointer :: dep0(:) => null()
    real(kind=8), pointer :: dep1(:) => null()
    real(kind=8), pointer :: du0(:) => null()
    real(kind=8), pointer :: du1(:) => null()
    character(len=24), pointer :: pltk(:) => null()
    real(kind=8), pointer :: plir(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infdbg('PILOTAGE', ifm, niv)
!
! --- AFFICHAGE
!
    if (niv .ge. 2) then
        write (ifm, *) '<PILOTAGE> ... CALCUL DU ETA_PILOTAGE'
    end if
!
! --- INITIALISATIONS
!
    pilcvg = -1
    ddepl0 = '&&CNPART.CHP1'
    ddepl1 = '&&CNPART.CHP2'
!
! --- DECOMPACTION VARIABLES CHAPEAUX
!
    call nmchex(valinc, 'VALINC', 'DEPMOI', depmoi)
    call nmchex(veasse, 'VEASSE', 'CNFEPI', cnfepi)
    call nmchex(solalg, 'SOLALG', 'DEPDEL', depdel)
    call nmchex(solalg, 'SOLALG', 'DEPPR1', deppr1)
    call nmchex(solalg, 'SOLALG', 'DEPPR2', deppr2)
    call dismoi('NB_EQUA', numedd, 'NUME_DDL', repi=neq)
!
! --- LECTURE DONNEES PILOTAGE
!
    call jeveuo(sdpilo(1:19)//'.PLTK', 'L', vk24=pltk)
    call jeveuo(sdpilo(1:19)//'.PLIR', 'L', vr=plir)
    etrmin = plir(5)
    etrmax = plir(4)
    typpil = pltk(1)
    coef = plir(1)
    dtau = deltat/coef
    ligrpi = pltk(2) (1:19)
    cartyp = pltk(3) (1:19)
    careta = pltk(4) (1:19)
!
! --- INCREMENTS DE DEPLACEMENT RHO.DU0 ET RHO.DU1
!
    call jeveuo(deppr1(1:19)//'.VALE', 'L', vr=du0)
    call jeveuo(deppr2(1:19)//'.VALE', 'L', vr=du1)
    call jeveuo(ddepl0(1:19)//'.VALE', 'E', vr=dep0)
    call jeveuo(ddepl1(1:19)//'.VALE', 'E', vr=dep1)
    do i = 1, neq
        dep0(i) = rho*du0(i)
        dep1(i) = du1(i)
    end do
!
    if (niv .ge. 2) then
        write (ifm, *) '<PILOTAGE> ...... SECOND MEMBRE DTAU : ', dtau
    end if
!
! --- PILOTAGE PAR UN DDL IMPOSE
!
    if (typpil .eq. 'DDL_IMPO') then
        call nmpidd(numedd, sdpilo, dtau, depdel, ddepl0, &
                    ddepl1, eta(1), pilcvg, nbeffe)
!
! --- PILOTAGE POUR ANALYSE LIMITE : TRAVAIL UNITAIRE
!
    else if (typpil .eq. 'ANA_LIM') then
        call nmpial(numedd, depdel, depmoi, cnfepi, ddepl0, &
                    ddepl1, eta(1), pilcvg, nbeffe)
!
! --- PILOTAGE PAR LONGUEUR D'ARC
!
    else if (typpil .eq. 'LONG_ARC') then
        call nmpila(numedd, sdpilo, dtau, depdel, &
                    ddepl0, ddepl1, nbeffe, eta, pilcvg)
!
! --- PILOTAGE PAR CRITERE
!
    else if (typpil .eq. 'PRED_ELAS' .or. typpil .eq. 'DEFORMATION') then
        call nmpipe(modele, ligrpi, cartyp, careta, ds_material, &
                    ds_constitutive, valinc, depdel, ddepl0, &
                    ddepl1, dtau, nbeffe, eta, pilcvg, &
                    typpil, carele)
    else
        ASSERT(ASTER_FALSE)
    end if
!
! --- LE CALCUL DE PILOTAGE A FORCEMENT ETE REALISE
!
    ASSERT(pilcvg .ge. 0)
!
! --- RECADRAGE DANS LES BORNES ETA_PILO_R_MIN ET ETA_PILO_R_MAX
!
    if (pilcvg .ne. 1) then
        if (nbeffe .eq. 2) then
            if ((eta(1) .lt. etrmin) .or. (eta(1) .gt. etrmax)) then
                nbeffe = nbeffe-1
                eta(1) = eta(2)
            end if
            if ((eta(2) .lt. etrmin) .or. (eta(2) .gt. etrmax)) then
                nbeffe = nbeffe-1
            end if
        end if
        if (nbeffe .eq. 1) then
            if ((eta(1) .lt. etrmin) .or. (eta(1) .gt. etrmax)) then
                nbeffe = nbeffe-1
            end if
        end if
        if (nbeffe .eq. 0) then
            pilcvg = 1
        end if
    end if
!
! --- LE CALCUL DE PILOTAGE A FORCEMENT ETE REALISE
!
    ASSERT(pilcvg .ge. 0)
!
    call jedema()
end subroutine
