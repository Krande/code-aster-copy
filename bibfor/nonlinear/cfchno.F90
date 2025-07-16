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
subroutine cfchno(noma, ds_contact, ndimg, posnoe, typenm, &
                  numenm, lmait, lescl, lmfixe, lefixe, &
                  tau1m, tau2m, tau1e, tau2e, tau1, &
                  tau2)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/cfnomm.h"
#include "asterfort/cfnorm.h"
#include "asterfort/mmmron.h"
#include "asterfort/mmnorm.h"
#include "asterfort/mmtann.h"
#include "asterfort/utmess.h"
#include "blas/dcopy.h"
#include "asterfort/int_to_char8.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=8) :: noma
    character(len=4) :: typenm
    integer(kind=8) :: ndimg
    integer(kind=8) :: posnoe, numenm
    type(NL_DS_Contact), intent(in) :: ds_contact
    real(kind=8) :: tau1m(3), tau2m(3)
    real(kind=8) :: tau1e(3), tau2e(3)
    real(kind=8) :: tau1(3), tau2(3)
    aster_logical :: lmfixe, lefixe
    aster_logical :: lmait, lescl
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (TOUTES METHODES - APPARIEMENT)
!
! CALCUL DES VECTEURS TANGENTS RESULTANTS SUIVANT OPTION
!
! ----------------------------------------------------------------------
!
! IN  NOMA   : NOM DU MAILLAGE
! IN  NDIMG  : DIMENSION DU MODELE
! In  ds_contact       : datastructure for contact management
! IN  LMFIXE : .TRUE. SI LA NORMALE MAITRE EST FIXE OU VECT_Y
! IN  LEFIXE : .TRUE. SI LA NORMALE ESCLAVE EST FIXE OU VECT_Y
! IN  LMAIT  : .TRUE. SI LA NORMALE = MAITRE  / MAIT_ESCL
! IN  LESCL  : .TRUE. SI LA NORMALE = ESCLAVE / MAIT_ESCL
! IN  POSNOE : NOEUD ESCLAVE (NUMERO DANS SD CONTACT)
! IN  TYPENM : TYPE DE L'ENTITE MAITRE RECEVANT LA PROJECTION
!               'MAIL' UNE MAILLE
!               'NOEU' UN NOEUD
! IN  NUMENM : NUMERO ABSOLU ENTITE MAITRE QUI RECOIT LA PROJECTION
! IN  TAU1M  : PREMIERE TANGENTE SUR LA MAILLE MAITRE AU POINT ESCLAVE
!              PROJETE
! IN  TAU2M  : SECONDE TANGENTE SUR LA MAILLE MAITRE AU POINT ESCLAVE
!              PROJETE
! IN  TAU1E  : PREMIERE TANGENTE AU NOEUD ESCLAVE
! IN  TAU2E  : SECONDE TANGENTE AU NOEUD ESCLAVE
! OUT TAU1   : PREMIERE TANGENTE LOCALE AU POINT ESCLAVE PROJETE
! OUT TAU2   : SECONDE TANGENTE LOCALE AU POINT ESCLAVE PROJETE
!
!
!
!
    integer(kind=8) :: i, niverr
    character(len=24) :: valk(2)
    real(kind=8) :: noor
    real(kind=8) :: enorm(3), mnorm(3), norm(3)
    character(len=8) :: nomnoe, nomenm
    blas_int :: b_incx, b_incy, b_n
!
! ----------------------------------------------------------------------
!
!
! --- NOM DE L'ENTITE (NOEUD OU MAILLE)
!
    if (typenm .eq. 'MAIL') then
        nomenm = int_to_char8(numenm)
    else if (typenm .eq. 'NOEU') then
        nomenm = int_to_char8(numenm)
    else
        ASSERT(.false.)
    end if
!
! --- NOM DU NOEUD ESCLAVE
!
    if (posnoe .le. 0) then
        nomnoe = 'PT CONT.'
    else
        call cfnomm(noma, ds_contact%sdcont_defi, 'NOEU', posnoe, nomnoe)
    end if
    valk(1) = nomnoe
    valk(2) = nomenm
!
! --- NORMALE AU NOEUD ESCLAVE: EXTERIEURE
!
    if (lescl) then
        call cfnorm(ndimg, tau1e, tau2e, enorm, noor)
        if (noor .le. r8prem()) then
            call utmess('F', 'CONTACT3_26', sk=nomnoe)
        end if
    end if
!
! --- NORMALE A LA MAILLE MAITRE: INTERIEURE
!
    if (lmait) then
        call mmnorm(ndimg, tau1m, tau2m, mnorm, noor)
        if (noor .le. r8prem()) then
            if (typenm .eq. 'MAIL') then
                call utmess('F', 'CONTACT3_27', sk=nomenm)
            else if (typenm .eq. 'NOEU') then
                call utmess('F', 'CONTACT3_26', sk=nomenm)
            else
                ASSERT(.false.)
            end if
        end if
    end if
!
! --- CALCUL DE LA NORMALE
!
    if (lmait .and. (.not. lescl)) then
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, mnorm, b_incx, norm, b_incy)
    else if (lmait .and. lescl) then
        do i = 1, 3
            norm(i) = (enorm(i)+mnorm(i))/2.d0
        end do
    else if (lescl) then
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, enorm, b_incx, norm, b_incy)
    else
        ASSERT(.false.)
    end if
!
! --- RECOPIE DES TANGENTES SI NORMALE FIXE
!
    if (lmfixe) then
        if (lmait .and. (.not. lescl)) then
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, tau1m, b_incx, tau1, b_incy)
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, tau2m, b_incx, tau2, b_incy)
        else
            ASSERT(.false.)
        end if
    end if
    if (lefixe) then
        if (lescl .and. (.not. lmait)) then
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, tau1e, b_incx, tau2, b_incy)
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, tau2e, b_incx, tau1, b_incy)
        else
            ASSERT(.false.)
        end if
    end if
!
! --- RE-CALCUL DES TANGENTES SI NORMALE AUTO
!
    if ((.not. lmfixe) .and. (.not. lefixe)) then
        call mmmron(ndimg, norm, tau1, tau2)
    end if
!
! --- NORMALISATION DES TANGENTES
!
    call mmtann(ndimg, tau1, tau2, niverr)
    if (niverr .eq. 1) then
        if (typenm .eq. 'MAIL') then
            call utmess('F', 'CONTACT3_31', nk=2, valk=valk)
        else if (typenm .eq. 'NOEU') then
            call utmess('F', 'CONTACT3_35', nk=2, valk=valk)
        else
            ASSERT(.false.)
        end if
    end if
!
end subroutine
