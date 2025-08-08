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

subroutine acevd2(noma, nomo, mcf, lmax, nbocc)
!
!
    implicit none
    integer(kind=8) :: lmax, nbocc
    character(len=8) :: noma, nomo
    character(len=*) :: mcf
!
! --------------------------------------------------------------------------------------------------
!
!        AFFE_CARA_ELEM
!           TEST DES CARACTERISTIQUES POUR LES ELEMENTS DISCRET
!
! --------------------------------------------------------------------------------------------------
!
! IN
!     NOMA     : NOM DU MAILLAGE
!     NOMO     : NOM DU MODELE
!     MCF      :  MOT CLEF
!     LMAX     : NOMBRE MAX DE MAILLE OU GROUPE DE MAILLE
!     NBOCC    : NOMBRE D'OCCURENCES DU MOT CLE DISCRET
! --------------------------------------------------------------------------------------------------
! person_in_charge: jean-luc.flejou at edf.fr
!
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/acevtr.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/getvem.h"
#include "asterfort/getvtx.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/verdis.h"
!
! --------------------------------------------------------------------------------------------------
    integer(kind=8), parameter :: nbcar = 100
    integer(kind=8) :: ier, i3d, i2d, ndim, ioc, ng, nm, ncar, icar
    integer(kind=8) :: ii, nbma, ialima
    character(len=8) :: nomu, car(nbcar)
    character(len=16) :: concep, cmd
    character(len=24) :: tmpdis, grpma
! --------------------------------------------------------------------------------------------------
    character(len=24), pointer :: zjdls(:) => null()
! --------------------------------------------------------------------------------------------------
!
    call getres(nomu, concep, cmd)
    tmpdis = nomu//'.DISCRET'
    grpma = noma//'.GROUPEMA       '
!
!   Vérification des dimensions / modélisations
    ier = 0
    call verdis(nomo, noma, 'F', i3d, i2d, ndim, ier)
    ASSERT((mcf .eq. 'DISCRET_2D') .or. (mcf .eq. 'DISCRET'))
!
    AS_ALLOCATE(vk24=zjdls, size=lmax)
!
!   Boucle sur les occurences
    do ioc = 1, nbocc
        call getvem(noma, 'GROUP_MA', mcf, 'GROUP_MA', ioc, lmax, zjdls, ng)
        call getvem(noma, 'MAILLE', mcf, 'MAILLE', ioc, lmax, zjdls, nm)
        call getvtx(mcf, 'CARA', iocc=ioc, nbval=nbcar, vect=car, nbret=ncar)
!
        if (ncar .gt. ncar) then
            ASSERT(.false.)
        end if
        do icar = 1, ncar
            if (car(icar) (3:4) .eq. 'TR') then
!               GROUP_MA = toutes les mailles de tous les groupes de mailles
                if (ng .gt. 0) then
                    do ii = 1, ng
                        call jelira(jexnom(grpma, zjdls(ii)), 'LONUTI', nbma)
                        call jeveuo(jexnom(grpma, zjdls(ii)), 'L', ialima)
                        call acevtr(nomo, 2, zk24(1), zi(ialima), nbma, ndim)
                    end do
                end if
!               MAILLE = toutes les mailles  de la liste de mailles
                if (nm .gt. 0) then
                    call acevtr(nomo, 1, zjdls, zi(1), nm, ndim)
                end if
            end if
        end do
    end do
!
    AS_DEALLOCATE(vk24=zjdls)
end subroutine
