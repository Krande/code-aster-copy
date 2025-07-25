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

subroutine te0582(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/mavec.h"
#include "asterfort/pmavec.h"
#include "asterfort/tumass.h"
#include "asterfort/turigi.h"
#include "asterfort/utmess.h"
    character(len=16) :: option, nomte
! ......................................................................
!
!    - FONCTION REALISEE:  CALCUL DES MATRICES ELEMENTAIRES
!                          TUYAU ET LES VECTEURS ELEMENTAIRES DE FORCES
!                          D ACCELERATION
!                          OPTION : RIGI_MECA, MASS_MECA, M_GAMMA
!                          SERT A DIMENSIONNER LES MATRICES
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
    integer(kind=8) :: nbrddm
    parameter(nbrddm=156)
    integer(kind=8) :: npg, ipoids, ivf
    integer(kind=8) :: ndim, nnos, nno, jcoopg, idfdk, jdfd2, jgano
    real(kind=8) :: mass(nbrddm*nbrddm), k(nbrddm*nbrddm)
    integer(kind=8) :: m, nbrddl, nc, iacce, ivect, imass
!
!
!
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, &
                     npg=npg, jpoids=ipoids, jcoopg=jcoopg, jvf=ivf, jdfde=idfdk, &
                     jdfd2=jdfd2, jgano=jgano)
!
    m = 3
    if (nomte .eq. 'MET6SEG3') m = 6
!
!     FORMULE GENERALE
!
    nbrddl = nno*(6+3+6*(m-1))
!
!     VERIFS PRAGMATIQUES
!
    if (nbrddl .gt. nbrddm) then
        call utmess('F', 'ELEMENTS4_40')
    end if
    if (nomte .eq. 'MET3SEG3') then
        if (nbrddl .ne. 63) then
            call utmess('F', 'ELEMENTS4_41')
        end if
    else if (nomte .eq. 'MET6SEG3') then
        if (nbrddl .ne. 117) then
            call utmess('F', 'ELEMENTS4_41')
        end if
    else if (nomte .eq. 'MET3SEG4') then
        if (nbrddl .ne. 84) then
            call utmess('F', 'ELEMENTS4_41')
        end if
    else
        call utmess('F', 'ELEMENTS4_42')
    end if
!
    if (option .eq. 'RIGI_MECA') then
        call turigi(nomte, nbrddl, k)
    else if ((option .eq. 'MASS_MECA') .or. (option .eq. 'M_GAMMA')) then
        call tumass(nomte, nbrddl, mass)
    end if
!
    if (option .eq. 'MASS_MECA') then
        call jevech('PMATUUR', 'E', imass)
!     DIMENSION DE LA MATRICE STOCKEE SOUS FORME VECTEUR
        nc = nbrddl*(nbrddl+1)/2
        call mavec(mass, nbrddl, zr(imass), nc)
    else if (option .eq. 'M_GAMMA') then
        call jevech('PACCELR', 'L', iacce)
        call jevech('PVECTUR', 'E', ivect)
        call pmavec('ZERO', nbrddl, mass, zr(iacce), zr(ivect))
    end if
!
end subroutine
