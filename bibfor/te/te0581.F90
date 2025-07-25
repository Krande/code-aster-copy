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

subroutine te0581(option, nomte)
    implicit none
!  DESCRIPTION : REALISE L'OPTION ADD_SIGM :
!  -----------   ADDITION DE CONTRAINTES AUX ELEMENTS
!
!  IN     : OPTION : CHARACTER*16 , SCALAIRE
!                    OPTION DE CALCUL
!  IN     : NOMTE  : CHARACTER*16 , SCALAIRE
!                    NOM DU TYPE ELEMENT
!
!-------------------   DECLARATION DES VARIABLES   ---------------------
!
!
! ARGUMENTS
! ---------
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
    character(len=16) :: option, nomte
!
! VARIABLES LOCALES
! -----------------
    integer(kind=8) :: jtab1(7), jtab2(7), jtab3(7), iret
    integer(kind=8) :: j1, j2, j3, nbcmp, nbsp, nbsp2, nbsp1, nbpt
    character(len=24) :: valk(2)
    real(kind=8) :: v1, v2
!
!-------------------   DEBUT DU CODE EXECUTABLE    ---------------------
!
    call tecach('ONO', 'PEPCON1', 'L', iret, nval=7, &
                itab=jtab1)
    call tecach('ONO', 'PEPCON2', 'L', iret, nval=7, &
                itab=jtab2)
    call tecach('ONO', 'PEPCON3', 'E', iret, nval=7, &
                itab=jtab3)
!
    if ((jtab1(1) .eq. 0) .or. (jtab2(1) .eq. 0) .or. (jtab3(1) .eq. 0)) then
        valk(1) = option
        valk(2) = nomte
        call utmess('F', 'ELEMENTS4_38', nk=2, valk=valk)
    end if
!
    if ((jtab1(2) .ne. jtab2(2)) .or. (jtab1(3) .ne. jtab2(3))) then
        call utmess('F', 'ELEMENTS4_39', sk=option)
    end if
!
!
    nbpt = jtab3(3)
    ASSERT(nbpt .eq. jtab1(3))
!
    nbsp = jtab3(7)
    nbsp1 = jtab1(7)
    nbsp2 = jtab2(7)
    ASSERT(nbsp .eq. nbsp1)
    ASSERT((nbsp2 .eq. nbsp) .or. (nbsp2 .eq. 1))
!
    nbcmp = jtab3(2)/nbpt
    ASSERT(jtab3(2) .eq. nbcmp*nbpt)
!
    do j1 = 1, nbpt
        do j2 = 1, nbsp
            do j3 = 1, nbcmp
                v1 = zr(jtab1(1)+(j1-1)*nbsp*nbcmp+(j2-1)*nbcmp+j3-1)
                if (nbsp .eq. nbsp2) then
                    v2 = zr(jtab2(1)+(j1-1)*nbsp*nbcmp+(j2-1)*nbcmp+j3-1)
                else
                    ASSERT(nbsp2 .eq. 1)
                    v2 = zr(jtab2(1)+(j1-1)*nbcmp+j3-1)
                end if
!              -- LA PRECONTRAINTE != 0 N'EST AUTORISEE QUE POUR
!              LES ELEMENTS DE BARRE ET LES ELEMENTS CABLE_GAINE
!              (CABLES DE PRECONTRAINTE) :
                if (v2 .ne. 0.d0) then
                    ASSERT(nomte .eq. 'MECA_BARRE' .or. nomte .eq. 'MECGSEG3')
                end if
                zr(jtab3(1)+(j1-1)*nbsp*nbcmp+(j2-1)*nbcmp+j3-1) = v1+v2
            end do
        end do
    end do
!
end subroutine
