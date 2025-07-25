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
subroutine te0270(option, nomte)
!    - FONCTION REALISEE:  CALCUL DES VECTEURS ELEMENTAIRES
!                          OPTION : 'CHAR_THER_ECHA_F'
!                          ELEMENTS FOURIER
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
!----------------------------------------------------------------------
! CORPS DU PROGRAMME
    implicit none
!
#include "jeveux.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
#include "asterfort/vff2dn.h"
!
    character(len=16) :: option, nomte
    integer(kind=8) :: nbres
    parameter(nbres=3)
    character(len=8) :: nompar(nbres)
    real(kind=8) :: valpar(nbres), poids, r, z, nx, ny, tpg, theta, coen, coenp1
    real(kind=8) :: texn, texnp1
    integer(kind=8) :: nno, nnos, jgano, ndim, kp, npg, ipoids, ivf, idfde, igeom
    integer(kind=8) :: itemps, ivectt, k, i, itex, icoefh, itemp, icode
!
!
!====
! 1.1 PREALABLES: RECUPERATION ADRESSES FONCTIONS DE FORMES...
!====
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
!====
! 1.2 PREALABLES LIES AUX RECHERCHES DE DONNEES GENERALES
!====
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PINSTR', 'L', itemps)
    call jevech('PTEMPER', 'L', itemp)
    call jevech('PCOEFHF', 'L', icoefh)
    call jevech('PT_EXTF', 'L', itex)
    call jevech('PVECTTR', 'E', ivectt)
!
!====
! 2. CALCULS TERMES DE MASSE
!====
    theta = zr(itemps+2)
    do kp = 1, npg
        k = (kp-1)*nno
        call vff2dn(ndim, nno, kp, ipoids, idfde, &
                    zr(igeom), nx, ny, poids)
        r = 0.d0
        z = 0.d0
        tpg = 0.d0
        do i = 1, nno
            r = r+zr(igeom+2*i-2)*zr(ivf+k+i-1)
            z = z+zr(igeom+2*i-1)*zr(ivf+k+i-1)
            tpg = tpg+zr(itemp+i-1)*zr(ivf+k+i-1)
        end do
        poids = poids*r
        valpar(1) = r
        nompar(1) = 'X'
        valpar(2) = z
        nompar(2) = 'Y'
        nompar(3) = 'INST'
        valpar(3) = zr(itemps)
        call fointe('FM', zk8(icoefh), 3, nompar, valpar, &
                    coenp1, icode)
        valpar(3) = zr(itemps)-zr(itemps+1)
        call fointe('FM', zk8(icoefh), 3, nompar, valpar, &
                    coen, icode)
        valpar(3) = zr(itemps)
        call fointe('FM', zk8(itex), 3, nompar, valpar, &
                    texnp1, icode)
        valpar(3) = zr(itemps)-zr(itemps+1)
        call fointe('FM', zk8(itex), 3, nompar, valpar, &
                    texn, icode)
!
        if (theta < -0.5) then
            do i = 1, nno
                zr(ivectt+i-1) = zr(ivectt+i-1)+poids*zr(ivf+k+i-1)*coenp1*(texnp1-tpg)
            end do
        else
            do i = 1, nno
                zr(ivectt+i-1) = zr(ivectt+i-1)+poids*zr(ivf+k+i-1)*(theta*coenp1*texnp1+(1.0d0-&
                                &theta)*coen*(texn-tpg))
            end do
        end if
    end do
end subroutine
