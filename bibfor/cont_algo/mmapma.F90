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

subroutine mmapma(mesh, ds_contact, model_ndim, i_zone, &
                  lexfro, typint, aliase, posmae, node_mast_nume, &
                  nnomae, elem_mast_indx, elem_mast_nume, ksipr1, ksipr2, &
                  tau1m, tau2m, iptm, iptc, norm, &
                  nommam)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/r8prem.h"
#include "asterfort/mmgaus.h"
#include "asterfort/mmnorm.h"
#include "asterfort/mmnumn.h"
#include "asterfort/mmpnoe.h"
#include "asterfort/mmsauv.h"
#include "asterfort/mmtanr.h"
#include "asterfort/utmess.h"
#include "asterfort/int_to_char8.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=8) :: mesh
    character(len=8) :: aliase
    type(NL_DS_Contact), intent(in) :: ds_contact
    real(kind=8) :: ksipr1, ksipr2
    integer(kind=8) :: model_ndim
    integer(kind=8) :: posmae, node_mast_nume
    integer(kind=8) :: elem_mast_indx, elem_mast_nume, nnomae
    integer(kind=8) :: i_zone, iptm, iptc
    integer(kind=8) :: typint
    real(kind=8) :: tau1m(3), tau2m(3), norm(3)
    character(len=8) :: nommam
    aster_logical :: lexfro
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (METHODES CONTINUES - APPARIEMENT)
!
! RECOPIE DE LA SD APPARIEMENT - CAS MAIT_ESCL
!
! ----------------------------------------------------------------------
!
! IN  LSSFRO : IL Y A DES NOEUDS DANS SANS_GROUP_NO_FR
! IN  NOMA   : NOM DU MAILLAGE
! IN  ALIASE : NOM D'ALIAS DE L'ELEMENT ESCLAVE
! In  ds_contact       : datastructure for contact management
! IN  NDIMG  : DIMENSION DE L'ESPACE
! IN  IZONE  : ZONE DE CONTACT ACTIVE
! IN  LEXFRO : LE POINT D'INTEGRATION DOIT-IL ETRE EXCLUS DU FROTTEMENT?
! IN  POSMAM : POSITION DE LA MAILLE MAITRE DANS LES SD CONTACT
! IN  NUMMAM : NUMERO ABSOLU MAILLE MAITRE QUI RECOIT LA PROJECTION
! IN  POSMAE : POSITION DE LA MAILLE ESCLAVE DANS LES SD CONTACT
! IN  NNOMAE : NOMBRE DE NOEUDS DE LA MAILLE ESCLAVE
! IN  TYPINT : TYPE D'INTEGRATION
! IN  IPTM   : NUMERO DU POINT D'INTEGRATION DANS LA MAILLE
! IN  KSIPR1 : PREMIERE COORDONNEE PARAMETRIQUE PT CONTACT PROJETE
!              SUR MAILLE MAITRE
! IN  KSIPR2 : SECONDE COORDONNEE PARAMETRIQUE PT CONTACT PROJETE
!              SUR MAILLE MAITRE
! IN  TAU1M  : PREMIERE TANGENTE SUR LA MAILLE MAITRE AU POINT ESCLAVE
!              PROJETE
! IN  TAU2M  : SECONDE TANGENTE SUR LA MAILLE MAITRE AU POINT ESCLAVE
!              PROJETE
! OUT NORM   : NORMALE FINALE
! OUT NOMMAM : NOM DE LA MAILLE MAITRE
!
!
!
!
    real(kind=8) :: noor
    real(kind=8) :: ksipc1, ksipc2, wpc
    real(kind=8) :: tau1(3), tau2(3)
    integer(kind=8) :: node_slav_indx, node_slav_nume
!
! ----------------------------------------------------------------------
!
!
! --- POSITION DU NOEUD ESCLAVE SI INTEGRATION AUX NOEUDS
!
    call mmpnoe(ds_contact%sdcont_defi, posmae, aliase, typint, iptm, &
                node_slav_indx)
!
! --- NUMERO ABSOLU DU POINT DE CONTACT
!
    call mmnumn(mesh, typint, node_mast_nume, nnomae, iptm, &
                node_slav_nume)
!
! --- RE-DEFINITION BASE TANGENTE SUIVANT OPTIONS
!
    call mmtanr(mesh, model_ndim, ds_contact, i_zone, &
                lexfro, node_slav_indx, ksipr1, ksipr2, elem_mast_indx, &
                elem_mast_nume, tau1m, tau2m, tau1, tau2)
!
! --- CALCUL DE LA NORMALE
!
    call mmnorm(model_ndim, tau1, tau2, norm, noor)
    if (noor .le. r8prem()) then
        nommam = int_to_char8(elem_mast_nume)
        call utmess('F', 'CONTACT3_24', sk=nommam)
    end if
!
! --- POIDS ET COORDONNEES DU POINT DE CONTACT
!
    call mmgaus(aliase, typint, iptm, ksipc1, ksipc2, &
                wpc)
!
! --- SAUVEGARDE APPARIEMENT
!
    call mmsauv(ds_contact, i_zone, iptc, elem_mast_nume, ksipr1, &
                ksipr2, tau1, tau2, node_mast_nume, node_slav_nume, &
                ksipc1, ksipc2, wpc)
!
end subroutine
