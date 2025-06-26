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
subroutine utelvf(elrefa, famil, nomjv, npg, nno)
!
    implicit none
!
#include "MeshTypes_type.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elraca.h"
#include "asterfort/elraga.h"
#include "asterfort/elrfvf.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: npg, nno
    character(len=8) :: elrefa, famil
    character(len=*) :: nomjv
!
! BUT: RECUPERER LES VALEURS DES FONCTIONS DE FORME
! ----------------------------------------------------------------------
!   IN   ELREFA : NOM DE L'ELREFA (K8)
!        FAMIL  : NOM DE LA FAMILLE DE POINTS DE GAUSS :
!                 'FPG1','FPG3',...
!   IN   NOMJV  : NOM JEVEUX POUR STOCKER LES FONCTIONS DE FORME
!   OUT  NPG    : NOMBRE DE POINTS DE GAUSS
!        NNO    : NOMBRE DE NOEUDS DU TYPE_MAILLE
! ----------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbpgmx = 27
    integer(kind=8) :: nbpg(MT_NBFAMX), ndim, nnos, nbfpg
    integer(kind=8) :: ifam, decal, ipg, ino, jvr
    real(kind=8) :: xpg(3*nbpgmx), poipg(nbpgmx), ff(MT_NNOMAX)
    character(len=8) :: nofpg(MT_NBFAMX)
! DEB ------------------------------------------------------------------
!

! - Get list of integration schemes of geometric support
    call elraca(elrefa, &
                nbfpg_=nbfpg, fapg_=nofpg, nbpg_=nbpg, &
                ndim_=ndim, nno_=nno, nnos_=nnos)

    ASSERT((ndim .ge. 0) .and. (ndim .le. 3))
    ASSERT((nno .gt. 0) .and. (nno .le. MT_NNOMAX))
    ASSERT((nbfpg .gt. 0) .and. (nbfpg .le. MT_NBFAMX))
!
    do ifam = 1, nbfpg
        if (nofpg(ifam) .eq. famil) goto 12
    end do
    ASSERT(ASTER_FALSE)
12  continue
!
    npg = nbpg(ifam)
    ASSERT((npg .gt. 0) .and. (npg .le. nbpgmx))
!
    call wkvect(nomjv, 'V V R', npg*nno, jvr)
!
!       -- COORDONNEES ET POIDS DES POINTS DE GAUSS :
!       ------------------------------------------------
    call elraga(elrefa, nofpg(ifam), ndim, npg, xpg, &
                poipg)
!
!     -- VALEURS DES FONCTIONS DE FORME :
!     ------------------------------------------------
    decal = 0
    do ipg = 1, npg
        call elrfvf(elrefa, xpg(ndim*(ipg-1)+1), ff, nno)
        do ino = 1, nno
            decal = decal+1
            zr(jvr-1+decal) = ff(ino)
        end do
    end do
!
end subroutine
