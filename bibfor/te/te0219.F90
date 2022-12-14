! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
subroutine te0219(option, nomte)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/connec.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/rcvalb.h"
#include "asterfort/teattr.h"
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES VECTEURS ELEMENTAIRES
!                          OPTION : 'CHAR_THER_GRAI_R/F  '
!                          EN 2D
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
    integer :: icodre(1), kpg, spt
    character(len=8) :: grxf, gryf, nompar(3), elrefe, alias8, fami, poum
    real(kind=8) :: dfdx(9), dfdy(9), poids, x, y, valres(1)
    real(kind=8) :: coorse(18), vectt(9), grx, gry, valpar(3)
    integer :: ndim, nno, nnos, kp, npg, i, k, ivectt, igrai
    integer :: ipoids, ivf, idfde, igeom, imate, jgano
    integer :: nnop2, c(6, 9), ise, nse, itemps, j, ier, ibid
!
!
    aster_logical :: fonc
!
!
    call elref1(elrefe)
!
    if (lteatt('LUMPE','OUI')) then
        call teattr('S', 'ALIAS8', alias8, ibid)
        if (alias8(6:8) .eq. 'QU9') elrefe='QU4'
        if (alias8(6:8) .eq. 'TR6') elrefe='TR3'
    endif
!
    call elrefe_info(elrefe=elrefe, fami='RIGI', ndim=ndim, nno=nno, nnos=nnos,&
                     npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
    call jevech('PGEOMER', 'L', igeom)
!
    if (option .eq. 'CHAR_THER_GRAI_R') then
        fonc=.false.
        call jevech('PGRAINR', 'L', igrai)
        grx=zr(igrai)
        gry=zr(igrai+1)
    else if (option.eq.'CHAR_THER_GRAI_F') then
        fonc=.true.
        call jevech('PTEMPSR', 'L', itemps)
        call jevech('PGRAINF', 'L', igrai)
        grxf=zk8(igrai)
        gryf=zk8(igrai+1)
        nompar(1)='X'
        nompar(2)='Y'
        nompar(3)='INST'
        valpar(3) = zr(itemps)
    endif
!
    call jevech('PMATERC', 'L', imate)
    call jevech('PVECTTR', 'E', ivectt)
    fami='FPG1'
    kpg=1
    spt=1
    poum='+'
    call rcvalb(fami, kpg, spt, poum, zi(imate),&
                ' ', 'THER', 0, ' ', [0.d0],&
                1, 'LAMBDA', valres, icodre, 1)
!
    call connec(nomte, nse, nnop2, c)
!
    do i = 1, nnop2
        vectt(i)=0.d0
    end do
!
!     BOUCLE SUR LES SOUS-ELEMENTS
    do ise = 1, nse
!
        do i = 1, nno
            do j = 1, 2
                coorse(2*(i-1)+j) = zr(igeom-1+2*(c(ise,i)-1)+j)
            end do
        end do
!
        do kp = 1, npg
            k=(kp-1)*nno
            call dfdm2d(nno, kp, ipoids, idfde, coorse,&
                        poids, dfdx, dfdy)
            x = 0.d0
            y = 0.d0
            do i = 1, nno
                x = x + coorse(2*(i-1)+1) * zr(ivf+k+i-1)
                y = y + coorse(2*(i-1)+2) * zr(ivf+k+i-1)
            end do
!
            if (fonc) then
                valpar(1) = x
                valpar(2) = y
                call fointe('FM', grxf, 3, nompar, valpar,&
                            grx, ier)
                call fointe('FM', gryf, 3, nompar, valpar,&
                            gry, ier)
            endif
!
            if (lteatt('AXIS','OUI')) poids = poids*x
            poids = poids*valres(1)
!
            do i = 1, nno
                vectt(c(ise,i)) = vectt( c(ise,i)) + poids*( dfdx(i)* grx+dfdy(i)*gry)
            end do
        end do
    end do
!
    do i = 1, nnop2
        zr(ivectt-1+i)=vectt(i)
    end do
!
end subroutine
