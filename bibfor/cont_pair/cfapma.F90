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

subroutine cfapma(noma, newgeo, ds_contact, lctfd, &
                  ndimg, izone, posnoe, numnoe, &
                  coorne, posmam, ksipr1, ksipr2, tau1m, &
                  tau2m, iliai)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/r8prem.h"
#include "asterfort/cfaddm.h"
#include "asterfort/cfcoor.h"
#include "asterfort/cfnewj.h"
#include "asterfort/cfnumm.h"
#include "asterfort/cfposn.h"
#include "asterfort/cfreli.h"
#include "asterfort/cftanr.h"
#include "asterfort/mmnorm.h"
#include "asterfort/utmess.h"
#include "asterfort/int_to_char8.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=8) :: noma
    type(NL_DS_Contact), intent(in) :: ds_contact
    character(len=19) :: newgeo
    real(kind=8) :: coorne(3), ksipr1, ksipr2
    real(kind=8) :: tau1m(3), tau2m(3)
    integer(kind=8) :: izone, ndimg
    integer(kind=8) :: posmam
    integer(kind=8) :: posnoe, numnoe
    integer(kind=8) :: iliai
    aster_logical :: lctfd
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (METHODES DISCRETES - APPARIEMENT)
!
! RECOPIE DE LA SD APPARIEMENT - CAS MAIT_ESCL
!
! ----------------------------------------------------------------------
!
!
! IN  NOMA   : NOM DU MAILLAGE
! IN  NEWGEO : NOUVELLE GEOMETRIE (AVEC DEPLACEMENT GEOMETRIQUE)
! In  ds_contact       : datastructure for contact management
! IN  LCTFD  : FROTTEMENT
! IN  NDIMG  : DIMENSION DE L'ESPACE
! IN  IZONE  : NUMERO DE LA ZONE DE CONTACT
! IN  POSNOE : INDICES DANS CONTNO DU NOEUD ESCLAVE
! IN  NUMNOE : NUMERO ABSOLU DU NOEUD ESCLAVE
! IN  COORNE : COORDONNEES DU NOEUD ESCLAVE
! IN  POSMAM : INDICES DANS CONTNO DE LA MAILLE MAITRE
! IN  KSIPR1 : COORDONNEE PARAMETRIQUE SUR MAITRE DU POINT ESCLAVE
!              PROJETE
! IN  KSIPR2 : COORDONNEE PARAMETRIQUE SUR MAITRE DU POINT ESCLAVE
!              PROJETE
! IN  TAU1M  : PREMIERE TANGENTE SUR LA MAILLE MAITRE
! IN  TAU2M  : SECONDE TANGENTE SUR LA MAILLE MAITRE
! IN  KSI2   : SECONDE COORD. DE LA PROJECTION
! IN  ILIAI  : INDICE DE LA LIAISON COURANTE
!
!

    real(kind=8) :: tau1(3), tau2(3)
    real(kind=8) :: norm(3), noor
    real(kind=8) :: jeu
    real(kind=8) :: coornp(3)
    character(len=8) :: nomnoe
    integer(kind=8) :: nbnom, nummam
    integer(kind=8) :: posnsm(9)
    real(kind=8) :: coefno(9)
!
! ----------------------------------------------------------------------

!
! --- NUMERO DE LA MAILLE MAITRE
!
    call cfnumm(ds_contact%sdcont_defi, posmam, nummam)
!
! --- CARACTERISTIQUES DE LA MAILLE MAITRE
!
    call cfposn(ds_contact%sdcont_defi, posmam, posnsm, nbnom)
!
! --- COORDONNEES PROJECTION DU NOEUD ESCLAVE SUR LA MAILLE MAITRE
!
    call cfcoor(noma, ds_contact%sdcont_defi, newgeo, posmam, ksipr1, &
                ksipr2, coornp)
!
! --- RE-DEFINITION BASE TANGENTE SUIVANT OPTIONS
!
    call cftanr(noma, ndimg, ds_contact, izone, &
                posnoe, 'MAIL', posmam, nummam, ksipr1, &
                ksipr2, tau1m, tau2m, tau1, tau2)
!
! --- CALCUL DE LA NORMALE INTERIEURE
!
    call mmnorm(ndimg, tau1, tau2, norm, noor)
    if (noor .le. r8prem()) then
        nomnoe = int_to_char8(numnoe)
        call utmess('F', 'CONTACT3_26', sk=nomnoe)
    end if
!
! --- CALCUL DU JEU
!
    call cfnewj(ndimg, coorne, coornp, norm, jeu)
!
! --- COEFFICIENT DE LA RELATION LINEAIRE SUR NOEUD MAITRE
!
    call cfreli(noma, nummam, nbnom, ksipr1, ksipr2, &
                coefno)
!
! --- AJOUT DE LA LIAISON NODALE
!
    call cfaddm(ds_contact, lctfd, posnoe, iliai, &
                ndimg, nbnom, posnsm, coefno, tau1, &
                tau2, norm, jeu, coornp)
!
end subroutine
