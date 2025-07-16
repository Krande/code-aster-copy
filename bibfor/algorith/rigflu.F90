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
subroutine rigflu(modelZ, matecoZ, &
                  nbLoad, loadNameZ, &
                  solverZ, numeDof, matrAsse)
!
    use listLoad_module
!
    implicit none
!
#include "asterfort/addModelLigrel.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/asmatr.h"
#include "asterfort/mecact.h"
#include "asterfort/merith.h"
#include "asterfort/numero.h"
#include "asterfort/preres.h"
#include "jeveux.h"
!
    integer(kind=8), intent(in) :: nbLoad
    character(len=*), intent(in) :: modelZ, loadNameZ, matecoZ, solverZ
    character(len=14), intent(out) :: numeDof
    character(len=8), intent(out) :: matrAsse
!
! --------------------------------------------------------------------------------------------------
!
! BUT : CETTE ROUTINE CALCULE LA MATRICE ASSEMBLEE DE RIGIDITE
!       FLUIDE S'APPUYANT SUR UN MODELE THERMIQUE
!     IN : MODELE : NOM DU MODELE FLUIDE UTILISE
!        : TIME   : INSTANT DU CALCUL
!        : NBCHAR : NOMBRE DE CHARGE
!        : CHAR   : NOM DE LA CHARGE
!        : MATE   : CHAMP DE MATERIAU
!        : SOLVEZ : METHODE DE RESOLUTION 'MULT_FRONT','LDLT' OU 'GCPC'
!     OUT: MA     : MATRICE ASSEMBLEE DE RIGIDITE FLUIDE
!        : NU     : NUMEROTATION ASSOCIEE
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24), parameter :: listLoad = '&&OP0152.INFCHA', timeMap = '&TIME'
    integer(kind=8) :: ibid, ierr
    character(len=8) :: matrElem, model, loadName
    character(len=24) :: mateco
    character(len=19) :: solver, matrElem19
    integer(kind=8) :: nbLigr
    character(len=24), pointer :: listLigr(:) => null()
    character(len=19), parameter :: maprec = "&&OP0152.MAPREC"
    integer(kind=8), parameter :: nbCmp = 6
    character(len=8), parameter :: cmpName(nbCmp) = (/ &
                                   'INST    ', 'DELTAT  ', 'THETA   ', &
                                   'KHI     ', 'R       ', 'RHO     '/)
    real(kind=8), parameter :: cmpVale(nbCmp) = (/ &
                               0.d0, 1.d0, 1.d0, &
                               0.d0, 0.d0, 0.d0/)
!
! --------------------------------------------------------------------------------------------------
!
    matrAsse = '&MATAS'
    numeDof = '&&RIGFLU.NUM'
    matrElem = '&MATEL'
    solver = solverZ
    model = modelZ
    mateco = matecoZ
    loadName = loadNameZ

! - Create list of loads (based on thermal phenomen)
    call creaListLoadFSIOne(model, nbLoad, loadName, listLoad)

! - Create map for time parameters
    call mecact('V', timeMap, 'MODELE', model, 'INST_R', &
                ncmp=nbCmp, lnomcmp=cmpName, vr=cmpVale)

! - CALCUL DE LA MATRICE ELEMENTAIRE DE RAIDEUR DU FLUIDE
    call merith(model, loadName, mateco, " ", &
                timeMap, matrElem, "V")

! - Add LIGREL from model
    nbLigr = 0
    call addModelLigrel(model, nbLigr, listLigr)

! - Get list of LIGREL from loads
    call getListLoadLigrel(listLoad, nbLigr, listLigr)

! - Create numbering
    call numero(numeDof, 'VV', &
                nbLigr, listLigr)
    AS_DEALLOCATE(vk24=listLigr)

! - Assemblying
    matrElem19 = matrElem
    call asmatr(1, matrElem19, ' ', numeDof, &
                listLoad, 'ZERO', 'V', 1, matrAsse)

! - Factor
    call preres(solver, 'V', ierr, maprec, matrAsse, &
                ibid, -9999)
!
end subroutine
