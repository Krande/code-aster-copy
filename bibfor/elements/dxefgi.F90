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

subroutine dxefgi(nomte, pgl, epsini, sigt)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/dxmate.h"
#include "asterfort/utmess.h"
!
    character(len=16) :: nomte
    real(kind=8) :: pgl(3, 3)
    real(kind=8) :: epsini(6)
    real(kind=8) :: sigt(*)
!     ------------------------------------------------------------------
! --- EFFORTS GENERALISES D'ORIGINE THERMIQUE AUX POINTS D'INTEGRATION
! --- POUR LES ELEMENTS COQUES A FACETTES PLANES :
! --- DST, DKT, DSQ, DKQ, Q4G
! --- CALCULES A PARTIR D'UN CHAMP DE DEFORMATIONS INITIALES QUI EST
! --- POUR L'INSTANT CONSTANT PAR ELEMENT ET QUI NE PREND PAS EN
! --- COMPTE LES DEFORMATIONS INITIALES DE CISAILLEMENT TRANSVERSE.
!     ------------------------------------------------------------------
!     IN  NOMTE        : NOM DU TYPE D'ELEMENT
!     IN  PGL(3,3)     : MATRICE DE PASSAGE DU REPERE GLOBAL AU REPERE
!                        LOCAL
!     IN  EPSINI(6)    : DEFORMATIONS INITIALES CONSTANTES SUR L'ELEMENT
!                        DANS L'ORDRE : EPXX, EPYY, EPXY, KXX, KYY, KXY
!     OUT SIGT(1)      : EFFORTS  GENERALISES D'ORIGINE THERMIQUE
!                        AUX POINTS D'INTEGRATION
    integer(kind=8) :: multic
    real(kind=8) :: df(3, 3), dm(3, 3), dmf(3, 3), dc(2, 2), dci(2, 2)
    real(kind=8) :: dmc(3, 2), dfc(3, 2)
    real(kind=8) :: kxx, kyy, kxy, t2iu(4), t2ui(4), t1ve(9)
    aster_logical :: coupmf
!     ------------------------------------------------------------------
!
! --- INITIALISATIONS :
!     -----------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, igau, nno, npg, ncomp
    real(kind=8) :: epxx, epxy, epyy, zero
!-----------------------------------------------------------------------
    zero = 0.0d0
    ncomp = 6
!
    do i = 1, 32
        sigt(i) = zero
    end do
!
    if (nomte .eq. 'MEDKTR3 ' .or. nomte .eq. 'MEDSTR3 ' .or. nomte .eq. 'MEDKTG3 ') then
!
        npg = 3
        nno = 3
!
    else if (nomte .eq. 'MEDKQU4 ' .or. nomte .eq. 'MEDSQU4 ' .or. &
             nomte .eq. 'MEQ4QU4 ' .or. nomte .eq. 'MEDKQG4 ') then
        npg = 4
        nno = 4
!
    else
        call utmess('F', 'ELEMENTS_14', sk=nomte)
    end if
!
! --- CALCUL DES MATRICES DE HOOKE DE FLEXION, MEMBRANE,
! --- MEMBRANE-FLEXION, CISAILLEMENT, CISAILLEMENT INVERSE
!     ----------------------------------------------------
!
    call dxmate('RIGI', df, dm, dmf, dc, &
                dci, dmc, dfc, nno, pgl, &
                multic, coupmf, t2iu, t2ui, t1ve)

!
! --- BOUCLE SUR LES POINTS D'INTEGRATION
!     -----------------------------------
    do igau = 1, npg
!
        epxx = epsini(ncomp*(igau-1)+1)
        epyy = epsini(ncomp*(igau-1)+2)
        epxy = 2.d0*epsini(ncomp*(igau-1)+3)
        kxx = epsini(ncomp*(igau-1)+4)
        kyy = epsini(ncomp*(igau-1)+5)
        kxy = 2.d0*epsini(ncomp*(igau-1)+6)
!
        sigt(1+8*(igau-1)) = dm(1, 1)*epxx+dm(1, 2)*epyy+dm(1, 3)*epxy
        sigt(2+8*(igau-1)) = dm(2, 1)*epxx+dm(2, 2)*epyy+dm(2, 3)*epxy
        sigt(3+8*(igau-1)) = dm(3, 1)*epxx+dm(3, 2)*epyy+dm(3, 3)*epxy
!
        sigt(4+8*(igau-1)) = df(1, 1)*kxx+df(1, 2)*kyy+df(1, 3)*kxy
        sigt(5+8*(igau-1)) = df(2, 1)*kxx+df(2, 2)*kyy+df(2, 3)*kxy
        sigt(6+8*(igau-1)) = df(3, 1)*kxx+df(3, 2)*kyy+df(3, 3)*kxy
!
        sigt(7+8*(igau-1)) = zero
        sigt(8+8*(igau-1)) = zero
    end do
!
end subroutine
