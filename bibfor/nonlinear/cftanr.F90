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

subroutine cftanr(noma, ndimg, ds_contact, izone, &
                  posnoe, typenm, posenm, numenm, ksipr1, &
                  ksipr2, tau1m, tau2m, tau1, tau2)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/apvect.h"
#include "asterfort/assert.h"
#include "asterfort/cfchno.h"
#include "asterfort/cfdisl.h"
#include "asterfort/cfinvm.h"
#include "asterfort/cfnben.h"
#include "asterfort/cfnomm.h"
#include "asterfort/cfnors.h"
#include "asterfort/cfnumm.h"
#include "asterfort/mmelty.h"
#include "asterfort/mminfi.h"
#include "asterfort/mminfl.h"
#include "asterfort/mminfr.h"
#include "asterfort/utmess.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=8) :: noma
    integer(kind=8) :: posenm, posnoe, numenm
    integer(kind=8) :: izone
    integer(kind=8) :: ndimg
    real(kind=8) :: ksipr1, ksipr2
    character(len=4) :: typenm
    type(NL_DS_Contact), intent(in) :: ds_contact
    real(kind=8) :: tau1(3), tau2(3)
    real(kind=8) :: tau1m(3), tau2m(3)
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (METHODES MAILLEES - APPARIEMENT)
!
! MOD. LES VECTEURS TANGENTS RESULTANTS SUIVANT OPTIONS
!
! ----------------------------------------------------------------------
!
!  NB: LE REPERE EST ORTHONORME ET TEL QUE LA NORMALE POINTE VERS
!  L'INTERIEUR DE LA MAILLE
!
! IN  NOMA   : NOM DU MAILLAGE
! IN  NDIMG  : DIMENSION DU MODELE
! In  ds_contact       : datastructure for contact management
! IN  IZONE  : ZONE DE CONTACT ACTIVE
! IN  TYPENM : TYPE DE L'ENTITE MAITRE RECEVANT LA PROJECTION
!               'MAIL' UNE MAILLE
!               'NOEU' UN NOEUD
! IN  POSENM : NUMERO ENTITE MAITRE QUI RECOIT LA PROJECTION
! IN  NUMENM : NUMERO ABSOLU ENTITE MAITRE QUI RECOIT LA PROJECTION
! IN  POSNOE : NOEUD ESCLAVE
! IN  KSIPR1 : COORDONNEE PARAMETRIQUE SUR MAITRE DU POINT ESCLAVE
!              PROJETE
! IN  KSIPR2 : COORDONNEE PARAMETRIQUE SUR MAITRE DU POINT ESCLAVE
!              PROJETE
! IN  TAU1M  : PREMIERE TANGENTE SUR LA MAILLE MAITRE AU POINT ESCLAVE
!              PROJETE
! IN  TAU2M  : SECONDE TANGENTE SUR LA MAILLE MAITRE AU POINT ESCLAVE
!              PROJETE
! OUT TAU1   : PREMIER VECTEUR TANGENT LOCAL AU POINT ESCLAVE PROJETE
! OUT TAU2   : SECOND VECTEUR TANGENT LOCAL AU POINT ESCLAVE PROJETE
!
!
!
!
    aster_logical :: lliss, lmfixe, lefixe, lmait, lescl
    aster_logical :: lpoutr, lpoint
    integer(kind=8) :: ima
    integer(kind=8) :: posmam, posmae, nummae, nummam
    integer(kind=8) :: itypem, itypee
    character(len=8) :: aliase, aliasm, nommam, nommae
    real(kind=8) :: r8bid, vector(3)
    real(kind=8) :: tau1e(3), tau2e(3)
    character(len=19) :: sdappa
    integer(kind=8) :: nbma, jdeciv
!
! ----------------------------------------------------------------------
!
!
! --- LECTURE APPARIEMENT
!
    sdappa = ds_contact%sdcont_solv(1:14)//'.APPA'
!
! --- INITIALISATIONS
!
    lmfixe = .false.
    lefixe = .false.
    lmait = .false.
    lescl = .false.
    vector(1) = 0.d0
    vector(2) = 0.d0
    vector(3) = 0.d0
!
! --- LISSAGE OU PAS ?
!
    lliss = cfdisl(ds_contact%sdcont_defi, 'LISSAGE')
!
! --- NORMALES A MODIFIER
!
    if (mminfl(ds_contact%sdcont_defi, 'MAIT', izone)) then
        lmait = .true.
        lescl = .false.
    else if (mminfl(ds_contact%sdcont_defi, 'MAIT_ESCL', izone)) then
        lmait = .true.
        lescl = .true.
    else if (mminfl(ds_contact%sdcont_defi, 'ESCL', izone)) then
        lmait = .false.
        lescl = .true.
    else
        ASSERT(.false.)
    end if
!
! --- INCOMPATIBILITE SCHEMA INTEGRATION GAUSS AVEC OPTION ESCLAVE
!
    if (lescl) then
        if (posnoe .eq. 0) then
            call utmess('F', 'CONTACT_98')
        else if (posnoe .eq. -1) then
            call utmess('F', 'CONTACT_99')
        end if
    end if
!
! --- DIMENSION MAILLE ESCLAVE: ON PREND LA PREMIERE MAILLE ATTACHEE
! --- AU NOEUD ESCLAVE
!
    if (lescl) then
        call cfnben(ds_contact%sdcont_defi, posnoe, 'CONINV', nbma, jdeciv)
        ima = 1
        call cfinvm(ds_contact%sdcont_defi, jdeciv, ima, posmae)
        call cfnumm(ds_contact%sdcont_defi, posmae, nummae)
        call cfnomm(noma, ds_contact%sdcont_defi, 'MAIL', posmae, nommae)
        call mmelty(noma, nummae, aliase)
    end if
!
! --- DIMENSION MAILLE MAITRE
!
    if (lmait) then
!
! --- RECUP. MAILLE SI APPARIEMENT NODAL
!
        if (typenm .eq. 'NOEU') then
            call cfnben(ds_contact%sdcont_defi, posenm, 'CONINV', nbma, jdeciv)
            ima = 1
            call cfinvm(ds_contact%sdcont_defi, jdeciv, ima, posmam)
        else
            posmam = posenm
        end if
        call cfnumm(ds_contact%sdcont_defi, posmam, nummam)
        call cfnomm(noma, ds_contact%sdcont_defi, 'MAIL', posmam, nommam)
        call mmelty(noma, nummam, aliasm)
    end if
!
! --- RECUPERATION TANGENTES ESCLAVES SI NECESSSAIRE
!
    if (lescl) then
        call apvect(sdappa, 'APPARI_NOEUD_TAU1', posnoe, tau1e)
        call apvect(sdappa, 'APPARI_NOEUD_TAU2', posnoe, tau2e)
    end if
!
! --- MODIFICATIONS DES NORMALES MAITRES
!
    if (lmait) then
        itypem = mminfi(ds_contact%sdcont_defi, 'VECT_MAIT', izone)
        if (itypem .ne. 0) then
            vector(1) = mminfr(ds_contact%sdcont_defi, 'VECT_MAIT_DIRX', izone)
            vector(2) = mminfr(ds_contact%sdcont_defi, 'VECT_MAIT_DIRY', izone)
            vector(3) = mminfr(ds_contact%sdcont_defi, 'VECT_MAIT_DIRZ', izone)
        end if
        if ((ndimg .eq. 2) .and. (itypem .eq. 2)) then
            call utmess('F', 'CONTACT3_43', sk=nommam)
        end if
        lpoutr = (ndimg .eq. 3) .and. (aliasm(1:2) .eq. 'SE')
        lpoint = aliasm .eq. 'POI1'
        if (lpoint) then
            call utmess('F', 'CONTACT3_75', sk=nommam)
        end if
        call cfnors(noma, ds_contact, posmam, typenm, &
                    numenm, lpoutr, lpoint, ksipr1, ksipr2, &
                    lliss, itypem, vector, tau1m, tau2m, &
                    lmfixe)
    end if
!
! --- MODIFICATIONS DES NORMALES ESCLAVES
!
    if (lescl) then
        itypee = mminfi(ds_contact%sdcont_defi, 'VECT_ESCL', izone)
        if (itypee .ne. 0) then
            vector(1) = mminfr(ds_contact%sdcont_defi, 'VECT_ESCL_DIRX', izone)
            vector(2) = mminfr(ds_contact%sdcont_defi, 'VECT_ESCL_DIRY', izone)
            vector(3) = mminfr(ds_contact%sdcont_defi, 'VECT_ESCL_DIRZ', izone)
        end if
        if ((ndimg .eq. 2) .and. (itypee .eq. 2)) then
            call utmess('F', 'CONTACT3_43', sk=nommae)
        end if
        lpoutr = (ndimg .eq. 3) .and. (aliase(1:2) .eq. 'SE')
        lpoint = aliase .eq. 'POI1'
        call cfnors(noma, ds_contact, posmae, typenm, &
                    numenm, lpoutr, lpoint, r8bid, r8bid, &
                    .false._1, itypee, vector, tau1e, tau2e, &
                    lefixe)
    end if
!
! --- CHOIX DE LA NORMALE -> CALCUL DES TANGENTES
!
    call cfchno(noma, ds_contact, ndimg, posnoe, typenm, &
                numenm, lmait, lescl, lmfixe, lefixe, &
                tau1m, tau2m, tau1e, tau2e, tau1, &
                tau2)
!
end subroutine
