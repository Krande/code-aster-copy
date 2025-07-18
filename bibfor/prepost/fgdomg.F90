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
subroutine fgdomg(method, nommat, nomnap, nomfon, valmin, &
                  valmax, ncyc, dommag)
!       ================================================================
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/fgdoba.h"
#include "asterfort/fgdohs.h"
#include "asterfort/fgdoma.h"
#include "asterfort/fgdowh.h"
#include "asterfort/fgtaep.h"
#include "asterfort/fgtaes.h"
#include "asterfort/rcpare.h"
#include "asterfort/utmess.h"
    character(len=*) :: method
    character(len=8) :: nommat, nomnap, nomfon
    real(kind=8) :: valmin(*), valmax(*), dommag
    integer(kind=8) :: ncyc
!
!       CALCUL DU DOMMAGE PAR DIFFERENTES METHODES
!       ----------------------------------------------------------------
!       IN  METHOD  METHODE DE CALCUL  DE DOMMAGE EMPLOYEE
!                       /WOHLER
!                       /MANSON_COFFIN
!                       /TAHERI_MANSON
!                       /TAHERI_MIXTE
!           NOMMAT  NOM DU CHAM_MATER
!           NOMNAP  NOM DE LA NAPPE POUR LOI DE TAHERI
!           NOMFON  NOM DE LA FONCTION POUR LOI DE TAHERI
!           CYCLE   TABLE DES CYCLES DETECTES (INDICES DEB ET FIN)
!           NCYC    NOMBRE  DE  CYCLE
!       OUT DOMMAG  VALEUR DU DOMMAGE
!       ----------------------------------------------------------------
    integer(kind=8) :: icodwo, icodba, icodhs
    character(len=32) :: pheno
    character(len=16) :: cara, k16b
    aster_logical :: lke, lhaigh
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ivcorr, ivke
    real(kind=8), pointer :: vdommag(:) => null()
!-----------------------------------------------------------------------
    dommag = 0.d0
    pheno = 'FATIGUE'
    AS_ALLOCATE(vr=vdommag, size=ncyc)
    lke = .false.
    lhaigh = .false.
    ivke = 0
    ivcorr = 0
!
    if (method .eq. 'WOHLER') then
        cara = 'WOHLER'
        call rcpare(nommat, pheno, cara, icodwo)
        cara = 'A_BASQUIN'
        call rcpare(nommat, pheno, cara, icodba)
        cara = 'A0'
        call rcpare(nommat, pheno, cara, icodhs)
        if (icodwo .eq. 0) then
            call fgdowh(nommat, ncyc, valmin, valmax, lke, &
                        zr(ivke), lhaigh, zr(ivcorr), vdommag)
        else if (icodba .eq. 0) then
            call fgdoba(nommat, ncyc, valmin, valmax, lke, &
                        zr(ivke), lhaigh, zr(ivcorr), vdommag)
        else if (icodhs .eq. 0) then
            call fgdohs(nommat, ncyc, valmin, valmax, lke, &
                        zr(ivke), lhaigh, zr(ivcorr), vdommag)
        end if
    else if (method .eq. 'MANSON_COFFIN') then
        call fgdoma(nommat, ncyc, valmin, valmax, vdommag)
    else if (method .eq. 'TAHERI_MANSON') then
        call fgtaep(nommat, nomfon, nomnap, ncyc, valmin, &
                    valmax, vdommag)
    else if (method .eq. 'TAHERI_MIXTE') then
        call fgtaes(nommat, nomnap, ncyc, valmin, valmax, &
                    vdommag)
    else
        k16b = method(1:16)
        call utmess('F', 'PREPOST_4', sk=k16b)
    end if
!
    do i = 1, ncyc
        dommag = dommag+vdommag(i)
    end do
!
    AS_DEALLOCATE(vr=vdommag)
end subroutine
