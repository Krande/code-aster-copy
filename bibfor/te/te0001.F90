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

subroutine te0001(option, nomte)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8dgrd.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
#include "asterfort/matrot.h"
#include "asterfort/utmess.h"
    character(len=16) :: option, nomte
!     CALCUL DES TERMES DE FORC_NOD_6DDL, 3DDL, 2DDL
!     -----------------------------------------------------------------
!     EN ENTREE :
!        OPTION : NOM DE L'OPTION A CALCULER
!        NOMTE  : NOM DU TYPE_ELEMENT
!     -----------------------------------------------------------------
!     -----------------------------------------------------------------
    real(kind=8) :: dgrd
    real(kind=8) :: valpar(4), angl(3), mat(3, 3), vect(6)
    character(len=8) :: nompar(4), nomfon
    aster_logical :: langl
!     -----------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ier, j, jdimp, jgeom, jtime, jvec
    integer(kind=8) :: nbpar, nddl, nddl1
!-----------------------------------------------------------------------
    if (nomte .eq. 'FORCE_NOD_6DDL') nddl1 = 6
    if (nomte .eq. 'FORCE_NOD_3DDL') nddl1 = 3
    if (nomte .eq. 'FORCE_NOD_2DDL') nddl1 = 2
    if (nomte .eq. 'FORCE_NOD_COQ2D') nddl1 = 3
    nddl = nddl1
    if (nomte .eq. 'FORCE_NOD_COQ2D') nddl = 2
!
    if (option .eq. 'CHAR_MECA_FORC_R') then
        call jevech('PGEOMER', 'L', jgeom)
        call jevech('PFORNOR', 'L', jdimp)
        call jevech('PVECTUR', 'E', jvec)
        do i = 1, nddl1
            zr(jvec-1+i) = zr(jdimp-1+i)
        end do
        langl = zr(jdimp+nddl1) .lt. 0.d0
        do i = 1, 3
            angl(i) = zr(jdimp+nddl1+i)
        end do
    else if (option .eq. 'CHAR_MECA_FORC_F') then
        nbpar = 4
        nompar(1) = 'X'
        nompar(2) = 'Y'
        nompar(3) = 'Z'
        nompar(4) = 'INST'
        call jevech('PGEOMER', 'L', jgeom)
        call jevech('PINSTR', 'L', jtime)
        call jevech('PVECTUR', 'E', jvec)
        valpar(1) = zr(jgeom-1+1)
        valpar(2) = zr(jgeom-1+2)
        valpar(3) = zr(jgeom-1+3)
        valpar(4) = zr(jtime-1+1)
        call jevech('PFORNOF', 'L', jdimp)
        do i = 1, nddl1
            nomfon = zk8(jdimp-1+i)
            ier = 0
            call fointe('FM', nomfon, nbpar, nompar, valpar, &
                        zr(jvec-1+i), ier)
        end do
        langl = zk8(jdimp+nddl1) .eq. 'UTILISAT'
        if (langl) then
            dgrd = r8dgrd()
            do i = 1, 3
                nomfon = zk8(jdimp+nddl1+i)
                ier = 0
                call fointe('FM', nomfon, nbpar, nompar, valpar, &
                            angl(i), ier)
                angl(i) = angl(i)*dgrd
            end do
        end if
    else
        call utmess('F', 'ELEMENTS2_61', sk=option)
    end if
!
!     --- PROJECTION DANS LE REPERE ABSOLU ---
    if (langl) then
        call matrot(angl, mat)
        do i = 1, min(nddl, 3)
            vect(i) = 0.d0
            do j = 1, min(nddl, 3)
                vect(i) = vect(i)+mat(j, i)*zr(jvec-1+j)
            end do
        end do
        do i = 4, min(nddl, 6)
            vect(i) = 0.d0
            do j = 4, min(nddl, 6)
                vect(i) = vect(i)+mat(j-3, i-3)*zr(jvec-1+j)
            end do
        end do
        do i = 1, nddl
            zr(jvec-1+i) = vect(i)
        end do
    end if
end subroutine
