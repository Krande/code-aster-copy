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

subroutine te0399(option, nomte)

! person_in_charge: mohamed.torkhani at edf.fr

    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elref1.h"
#include "asterfort/arlref.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/arlapl.h"

    character(len=16) :: nomte, option

! ----------------------------------------------------------------------

! CALCUL DES MATRICES DE COUPLAGE ARLEQUIN
! OPTION ARLQ_MATR

! CREATION DES 2 MATRICES DE COUPLAGE POUR LES ELEMENTS 1D ET 3D
!

! ----------------------------------------------------------------------

! IN  OPTION : OPTION DE CALCUL
! IN  NOMTE  : NOM DU TYPE ELEMENT

    integer(kind=8) :: nbnomx
    parameter(nbnomx=27)
    integer(kind=8) :: nbgamx
    parameter(nbgamx=64)

    integer(kind=8) :: nns, nnos
    integer(kind=8) :: npgs, ipoids, ivfs, idfdes, jrefe1, jrefe2
    integer(kind=8) :: jfamil, jinfor, jcoopg, jdfd2, jgano
    integer(kind=8) :: ndim
    character(len=8) :: nomfam, elrfs
    character(len=8) :: elrf1, elrf2
    character(len=16) :: nomte1, nomte2
    integer(kind=8) :: nn1, nn2

! ----------------------------------------------------------------------

! --- FAMILLE D'INTEGRATION

    call jevech('PFAMILK', 'L', jfamil)
    nomfam = zk8(jfamil)

! --- INFORMATIONS SUR MAILLES COUPLEES

    call jevech('PINFORR', 'L', jinfor)
    ndim = nint(zr(jinfor+1-1))
    nn1 = nint(zr(jinfor+2-1))
    nn2 = nint(zr(jinfor+3-1))
    call jevech('PREFE1K', 'L', jrefe1)
    elrf1 = zk8(jrefe1)
    call jevech('PREFE2K', 'L', jrefe2)
    elrf2 = zk8(jrefe2)

! --- SCHEMA INTEGRATION MAILLE SUPPORT

    call elref1(elrfs)

    if (elrf2 == 'SE2') then
        nomte2 = 'MECA_POU_D_T'
    end if
    if (elrf1 == 'H20') then
        nomte1 = 'MECA_HEXA20'
    elseif (elrf1 == 'HE8') then
        nomte1 = 'MECA_HEXA8'
    elseif (elrf1 == 'P15') then
        nomte1 = 'MECA_PENTA15'
    elseif (elrf1 == 'PE6') then
        nomte1 = 'MECA_PENTA6'
    elseif (elrf1 == 'T10') then
        nomte1 = 'MECA_TETRA10'
    elseif (elrf1 == 'TE4') then
        nomte1 = 'MECA_TETRA4'
    end if

    if (nomte1 .ne. nomte) then
        call arlref(elrefe=elrf1, fami=nomfam, nomte=nomte1, ndim=ndim, nno=nns, nnos=nnos, &
                    npg=npgs, jpoids=ipoids, jcoopg=jcoopg, jvf=ivfs, jdfde=idfdes, &
                    jdfd2=jdfd2, jgano=jgano)
        nns = nn1
    else
        call elrefe_info(elrefe=elrfs, fami=nomfam, ndim=ndim, nno=nns, nnos=nnos, &
                         npg=npgs, jpoids=ipoids, jcoopg=jcoopg, jvf=ivfs, jdfde=idfdes, &
                         jdfd2=jdfd2, jgano=jgano)
    end if

! --- VERIFICATIONS

    ASSERT(nn1 .le. nbnomx)
    ASSERT(nn2 .le. nbnomx)
    ASSERT(nns .le. nbnomx)
    ASSERT(npgs .le. nbgamx)

    call arlapl(ndim, nns, nn1, nn2, nomte, npgs, ipoids, ivfs, idfdes)

end subroutine
