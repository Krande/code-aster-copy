! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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

subroutine te0008(option, nomte)
!
!
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/bsigmc.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/nbsigm.h"
#include "asterfort/r8inir.h"
#include "asterfort/tecach.h"
#include "asterfort/terefe.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
!
    character(len=16) :: option, nomte
!
! ----------------------------------------------------------------------
!
! ELEMENTS ISOPARAMETRIQUES 2D ET 3D
!
! CALCUL DES OPTIONS FORC_NODA ET REFE_FORC_NODA
!
! ----------------------------------------------------------------------
!
!
! IN  OPTION : OPTION DE CALCUL
! IN  NOMTE  : NOM DU TYPE ELEMENT
!
!
!
!
    real(kind=8) :: zero, sigref
    real(kind=8) :: nharm, bsigm(81), geo(81), sigtmp(162), ftemp(81)
    integer :: nbsig, ndim, nno, nnos, npg1
    integer :: ipoids, ivf, idfde, jgano
    integer :: igeom, ivectu
    integer :: jvDisp, jvSief, jvCompor
    integer :: i, j, nbinco, iretc
!
! ----------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, &
                     npg=npg1, jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
! --- INITIALISATIONS
!
    zero = 0.d0
    nharm = zero
    nbinco = nno*ndim
    ASSERT(nbinco .le. 81)
!
! --- NOMBRE DE CONTRAINTES ASSOCIE A L'ELEMENT
!
    nbsig = nbsigm()
!
! --- PARAMETRE EN SORTIE: VECTEUR DES FORCES INTERNES (BT*SIGMA)
!
    call jevech('PVECTUR', 'E', ivectu)
!
! --- PARAMETRE EN ENTREE: GEOMETRIE
!
    call jevech('PGEOMER', 'L', igeom)
    call dcopy(nbinco, zr(igeom), 1, geo, 1)
!
    if (option .eq. 'FORC_NODA') then
        call jevech('PSIEFR', 'L', jvSief)
        call tecach('ONO', 'PCOMPOR', 'L', iretc, iad=jvCompor)
        if (iretc .eq. 0) then
            if (zk16(jvCompor+2) (1:6) .ne. 'PETIT ') then
                call jevech('PDEPLAR', 'L', jvDisp)
                call daxpy(nbinco, 1.d0, zr(jvDisp), 1, geo, 1)
            end if
        end if
        call bsigmc(nno, ndim, nbsig, npg1, ipoids, &
                    ivf, idfde, geo, nharm, zr(jvSief), &
                    bsigm)
        call dcopy(nbinco, bsigm, 1, zr(ivectu), 1)
!
    else if (option .eq. 'REFE_FORC_NODA') then
!
        call terefe('SIGM_REFE', 'MECA_ISO', sigref)
!
        call r8inir(nbsig*npg1, 0.d0, sigtmp, 1)
        call r8inir(nbinco, 0.d0, ftemp, 1)
!
        do i = 1, nbsig*npg1
!
            sigtmp(i) = sigref
            call bsigmc(nno, ndim, nbsig, npg1, ipoids, &
                        ivf, idfde, geo, nharm, sigtmp, &
                        bsigm)
!
            do j = 1, nbinco
                ftemp(j) = ftemp(j)+abs(bsigm(j))
            end do
!
            sigtmp(i) = 0.d0
!
        end do
!
        call daxpy(nbinco, 1.d0/npg1, ftemp, 1, zr(ivectu), 1)
!
    end if
!
end subroutine
