! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

subroutine te0056(option, nomte)
!
    use FE_topo_module
    use FE_quadrature_module
    use FE_basis_module
!.......................................................................
    implicit none
#include "jeveux.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe_varc.h"
#include "asterfort/jevech.h"
!
    character(len=16) :: option, nomte
!
!     BUT: CALCUL DU SECOND MEMBRE ELEMENTAIRE EN THERMIQUE CORRESPON-
!          DANT A UNE SOURCE VOLUMIQUE FONCTION
!          ELEMENTS ISOPARAMETRIQUES 3D
!
!          OPTION : 'CHAR_THER_SOUR_F'
!
!     ENTREES  ---> OPTION : OPTION DE CALCUL
!          ---> NOMTE  : NOM DU TYPE ELEMENT
!.......................................................................
!
    type(FE_Cell) :: FECell
    type(FE_Quadrature) :: FEQuadCell
    type(FE_basis_cell) :: FEBasis

    character(len=8) :: nompar(4)
    real(kind=8) :: sourc, theta
    real(kind=8) :: valpar(4), funcEF(27)
    integer :: ipoids, ivf, idfde, igeom
    integer :: jgano, nno, kp, npg1, i, ivectt, isour, itemps
!
!
!-----------------------------------------------------------------------
    integer :: ier, ndim, nnos
    real(kind=8) :: soun, sounp1
!-----------------------------------------------------------------------
    call FECell%init()

    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, &
                     npg=npg1, jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)

    call FEQuadCell%initCell(FECell, npg1)
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PSOURCF', 'L', isour)
    call jevech('PTEMPSR', 'L', itemps)
    call jevech('PVECTTR', 'E', ivectt)
!
    theta = zr(itemps+2)
    nompar(1) = 'X'
    nompar(2) = 'Y'
    nompar(3) = 'Z'
    nompar(4) = 'INST'
!
    do i = 1, nno
        zr(ivectt-1+i) = 0.0d0
    end do
!
!    BOUCLE SUR LES POINTS DE GAUSS
!
    do kp = 1, FEQuadCell%nbQuadPoints
!
!    CALCUL DE SOURC
!
        valpar(1:3) = FECell%evalCoor(FEQuadCell%points(1:3, kp))
        valpar(4) = zr(itemps)

!       EC : je voulais mettre fami = RIGI et kpg = kp
!       mais il n'y a que FPG1 dans MATER pour cet élément
        call fointe_varc('FM', 'FPG1', 1, 1, '+', &
                         zk8(isour), 4, nompar, valpar, &
                         sounp1, ier)
!        call fointe('FM', zk8(isour), 4, nompar, valpar,&
!                    sounp1, ier)
        if (theta .ne. 1.0d0) then
            valpar(4) = zr(itemps)-zr(itemps+1)
            call fointe_varc('FM', 'FPG1', 1, 1, '+', &
                             zk8(isour), 4, nompar, valpar, &
                             soun, ier)
!            call fointe('FM', zk8(isour), 4, nompar, valpar,&
!                        soun, ier)
        else
            soun = 0.d0
        end if
        if (theta < -0.5) then
            sourc = sounp1
        else
            sourc = theta*sounp1+(1.0d0-theta)*soun
        end if
!
        call FEBasis%func(FECell, FEQuadCell%points(1:3, kp), funcEF)
!
        do i = 1, nno
            zr(ivectt+i-1) = zr(ivectt+i-1)+FEQuadCell%weights(kp)*sourc*funcEF(i)
        end do
!
    end do
!
end subroutine
