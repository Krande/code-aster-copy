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

subroutine cfapre(mesh, ds_contact, time_curr)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/apinfi.h"
#include "asterfort/apinfr.h"
#include "asterfort/apvect.h"
#include "asterfort/assert.h"
#include "asterfort/cfapma.h"
#include "asterfort/cfapno.h"
#include "asterfort/cfappi.h"
#include "asterfort/cfcorn.h"
#include "asterfort/cfdisd.h"
#include "asterfort/cfdisi.h"
#include "asterfort/cfdisl.h"
#include "asterfort/cfdist.h"
#include "asterfort/cfecrd.h"
#include "asterfort/cfmmco.h"
#include "asterfort/cfnumn.h"
#include "asterfort/cfparz.h"
#include "asterfort/infdbg.h"
#include "asterfort/mminfi.h"
#include "asterfort/mminfl.h"
#include "asterfort/mminfr.h"
#include "asterfort/int_to_char8.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=8), intent(in) :: mesh
    type(NL_DS_Contact), intent(in) :: ds_contact
    real(kind=8), intent(in) :: time_curr
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Solve
!
! Discrete method - Save pairing in contact datastructures
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  ds_contact       : datastructure for contact management
! In  time_curr        : current time
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    character(len=19) :: sdappa, newgeo
    integer(kind=8) :: izone, i, iliai, ip
    integer(kind=8) :: jdecne
    integer(kind=8) :: inoe
    integer(kind=8) :: posmae, posnoe(1), posmam, posnom(1)
    integer(kind=8) :: numnoe(1)
    integer(kind=8) :: entapp, typapp
    aster_logical :: lctfd
    integer(kind=8) :: nzoco, ndimg, nbpt, nbliai
    integer(kind=8) :: nesmax
    aster_logical :: lveri
    character(len=8) :: nomnoe
    real(kind=8) :: ksipr1, ksipr2, tau1m(3), tau2m(3)
    real(kind=8) :: coorne(3), gap_user
    real(kind=8) :: coefff, coefpn, coefpt, coefte
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('CONTACT', ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<CONTACT> ... Save pairing in contact datastructures'
    end if
!
! - Pairing datastructure
!
    sdappa = ds_contact%sdcont_solv(1:14)//'.APPA'
!
! - New geometry name
!
    newgeo = ds_contact%sdcont_solv(1:14)//'.NEWG'
!
! --- INFOS SUR LA CHARGE DE CONTACT
!
    lctfd = cfdisl(ds_contact%sdcont_defi, 'FROT_DISCRET')
!
! --- NOMBRE TOTAL DE NOEUDS ESCLAVES ET DIMENSION DU PROBLEME
!
    nzoco = cfdisi(ds_contact%sdcont_defi, 'NZOCO')
    ndimg = cfdisi(ds_contact%sdcont_defi, 'NDIM')
!
! --- INITIALISATIONS
!
    ip = 1
    iliai = 0
    posmae = 0
!
! --- BOUCLE SUR LES ZONES
!
    do izone = 1, nzoco
!
! ----- INFORMATION SUR LA ZONE
!
        nbpt = mminfi(ds_contact%sdcont_defi, 'NBPT', izone)
        jdecne = mminfi(ds_contact%sdcont_defi, 'JDECNE', izone)
!
! ----- MODE VERIF: ON SAUTE LES POINTS
!
        lveri = mminfl(ds_contact%sdcont_defi, 'VERIF', izone)
        if (lveri) then
            ip = ip+nbpt
            goto 25
        end if
!
! ----- COEFFICIENTS
!
        call cfmmco(ds_contact, izone, 'E_N', 'L', coefpn)
        call cfmmco(ds_contact, izone, 'E_T', 'L', coefpt)
        coefff = mminfr(ds_contact%sdcont_defi, 'COEF_COULOMB', izone)
        coefte = mminfr(ds_contact%sdcont_defi, 'COEF_MATR_FROT', izone)
!
! ----- BOUCLE SUR LES NOEUDS DE CONTACT
!
        do i = 1, nbpt
!
! ------- NOEUD ESCLAVE COURANT
!
            inoe = i
            posnoe = jdecne+inoe
!
! ------- INDICE ABSOLU DANS LE MAILLAGE DU NOEUD
!
            call cfnumn(ds_contact%sdcont_defi, 1, posnoe(1), numnoe(1))
!
! ------- COORDONNEES DU NOEUD
!
            call cfcorn(newgeo, numnoe(1), coorne)
!
! ------- NOM DU NOEUD
!
            nomnoe = int_to_char8(numnoe(1))
!
! ------- INFOS APPARIEMENT
!
            call apinfi(sdappa, 'APPARI_TYPE', ip, typapp)
            call apinfi(sdappa, 'APPARI_ENTITE', ip, entapp)
            call apinfr(sdappa, 'APPARI_PROJ_KSI1', ip, ksipr1)
            call apinfr(sdappa, 'APPARI_PROJ_KSI2', ip, ksipr2)
            call apvect(sdappa, 'APPARI_TAU1', ip, tau1m)
            call apvect(sdappa, 'APPARI_TAU2', ip, tau2m)
!
! ------- RECOPIE APPARIEMENT
!
            if (typapp .lt. 0) then
                if (niv .ge. 2) then
                    call cfappi(mesh, ds_contact%sdcont_defi, nomnoe, typapp, entapp)
                end if
                goto 35
            else if (typapp .eq. 1) then
! --------- CARAC. MAITRE
                posnom = entapp
! --------- LIAISON DE CONTACT EFFECTIVE
                iliai = iliai+1
! --------- CALCUL LIAISON
                call cfapno(mesh, newgeo, ds_contact, lctfd, &
                            ndimg, izone, posnoe(1), numnoe(1), &
                            coorne, posnom(1), tau1m, tau2m, iliai)
!
            else if (typapp .eq. 2) then
! --------- CARAC. MAITRE
                posmam = entapp
! --------- LIAISON DE CONTACT EFFECTIVE
                iliai = iliai+1
! --------- CALCUL LIAISON
                call cfapma(mesh, newgeo, ds_contact, lctfd, &
                            ndimg, izone, posnoe(1), numnoe(1), &
                            coorne, posmam, ksipr1, ksipr2, tau1m, &
                            tau2m, iliai)
            else
                ASSERT(.false.)
            end if
!
! ------- CALCUL DU JEU FICTIF DE LA ZONE
!
            call cfdist(ds_contact, izone, posmae, coorne, time_curr, &
                        gap_user, node_slav_indx_=posnoe(1))
!
! ------- CARACTERISTIQUES DE LA LIAISON POUR LA ZONE
!
            call cfparz(ds_contact, iliai, coefff, coefpn, coefpt, &
                        coefte, gap_user, izone, ip, numnoe(1), &
                        posnoe(1))
!
35          continue
!
! ------- POINT SUIVANT
!
            ip = ip+1
            ASSERT(iliai .le. ip)
!
        end do
25      continue
    end do
!
! --- NOMBRE DE LIAISONS EFFECTIVES
!
    nbliai = iliai
    call cfecrd(ds_contact%sdcont_solv, 'NBLIAI', nbliai)
    nesmax = cfdisd(ds_contact%sdcont_solv, 'NESMAX')
    ASSERT(nbliai .le. nesmax)
!
end subroutine
