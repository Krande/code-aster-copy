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
subroutine mmexfr(mesh, ds_contact, i_zone, elem_mast_indx, tau1, &
                  tau2)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterc/r8prem.h"
#include "asterfort/cfdisi.h"
#include "asterfort/cfnomm.h"
#include "asterfort/mminfi.h"
#include "asterfort/mminfr.h"
#include "asterfort/mmnorm.h"
#include "asterfort/normev.h"
#include "asterfort/provec.h"
#include "asterfort/utmess.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/ddot.h"
#include "blas/dscal.h"
!
! person_in_charge: thomas.de-soza at edf.fr
!
    character(len=8), intent(in) :: mesh
    type(NL_DS_Contact), intent(in) :: ds_contact
    integer(kind=8), intent(in) :: i_zone
    integer(kind=8), intent(in) :: elem_mast_indx
    real(kind=8), intent(out) :: tau1(3)
    real(kind=8), intent(out) :: tau2(3)
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (METHODE CONTINUE - APPARIEMENT - UTILITAIRE)
!
! REDEFINIT LA BASE TANGENTE LOCALE POUR SANS_GROUP_NO_FR
!
! ----------------------------------------------------------------------
!
! IN  NOMA   : NOM DU MAILLAGE
! In  ds_contact       : datastructure for contact management
! IN  IZONE  : NUMERO DE LA ZONE DE CONTACT
! IN  POSMAM : POSITION DE LA MAILLE MAITRE DANS LES SD CONTACT
! OUT TAU1   : PREMIER VECTEUR TANGENT
! OUT TAU2   : SECOND VECTEUR TANGENT
!
    integer(kind=8) :: model_ndim, ndirex
    real(kind=8) :: vdirex(3), norm(3), norme
    real(kind=8) :: tau1fr(3), tau2fr(3)
    real(kind=8) :: extau1, extau2
    character(len=8) :: elem_mast_name
    blas_int :: b_incx, b_incy, b_n
!
! ----------------------------------------------------------------------
!
    model_ndim = cfdisi(ds_contact%sdcont_defi, 'NDIM')
!
! --- NOMBRE DE DIRECTIONS A EXCLURE POUR LA ZONE
!
    ndirex = mminfi(ds_contact%sdcont_defi, 'EXCL_DIR', i_zone)
!
! --- REDEFINITION DU REPERE SI NECESSAIRE (UNE DIRECTION EXCLUE EN 3D)
!
    if ((model_ndim .eq. 3) .and. (ndirex .eq. 1)) then
! ----- DIRECTION D'EXCLUSION
        vdirex(1) = mminfr(ds_contact%sdcont_defi, 'EXCL_FROT_DIRX', i_zone)
        vdirex(2) = mminfr(ds_contact%sdcont_defi, 'EXCL_FROT_DIRY', i_zone)
        vdirex(3) = mminfr(ds_contact%sdcont_defi, 'EXCL_FROT_DIRZ', i_zone)
! ----- ON LA PROJETTE SUR LE PLAN TANGENT
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, tau1, b_incx, tau1fr, b_incy)
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, tau2, b_incx, tau2fr, b_incy)
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        extau1 = ddot(b_n, vdirex, b_incx, tau1, b_incy)
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        extau2 = ddot(b_n, vdirex, b_incx, tau2, b_incy)
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        call dscal(b_n, extau1, tau1fr, b_incx)
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, extau2, tau2fr, b_incx, tau1fr, &
                   b_incy)
        call normev(tau1fr, norme)
        if (norme .le. r8prem()) then
            call cfnomm(mesh, ds_contact%sdcont_defi, 'MAIL', elem_mast_indx, elem_mast_name)
            call utmess('F', 'CONTACT3_18', sk=elem_mast_name, si=i_zone, nr=3, &
                        valr=vdirex)
        end if
! ----- ON CALCULE TAU2FR PAR PROD. VECT.
        call mmnorm(model_ndim, tau1, tau2, norm, norme)
        call provec(tau1fr, norm, tau2fr)
! ----- RECOPIE
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, tau1fr, b_incx, tau1, b_incy)
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, tau2fr, b_incx, tau2, b_incy)
    end if
!
end subroutine
