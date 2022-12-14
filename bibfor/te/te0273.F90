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
subroutine te0273(option, nomte)
!
!     BUT: CALCUL DES VECTEURS ELEMENTAIRES EN THERMIQUE
!          CORRESPONDANT AU FLUX NON-LINEAIRE
!          SUR DES FACES D'ELEMENTS ISOPARAMETRIQUES 3D
!
!          OPTION : 'CHAR_THER_FLUNL'
!
!     ENTREES  ---> OPTION : OPTION DE CALCUL
!              ---> NOMTE  : NOM DU TYPE ELEMENT
!----------------------------------------------------------------------
! CORPS DU PROGRAMME
    implicit none
!
#include "jeveux.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/foderi.h"
#include "asterfort/jevech.h"
!
    character(len=16) :: option, nomte
    real(kind=8) :: nx, ny, nz, sx(9, 9), sy(9, 9), sz(9, 9), jac, theta, tpg
    real(kind=8) :: alpha, rbid
    integer :: ndim, nno, npg1, ipoids, ivf, idfdx, idfdy, igeom, iflux, itempr
    integer :: itemps, ino, jno, ivectt, i, j, kp, kdec, ldec, idec, jdec, nnos
    integer :: jgano
    character(len=8) :: coef
!
!
!====
! 1.1 PREALABLES: RECUPERATION ADRESSES FONCTIONS DE FORMES...
!====
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg1,&
                     jpoids=ipoids, jvf=ivf, jdfde=idfdx, jgano=jgano)
    idfdy = idfdx + 1
!
!
!====
! 1.2 PREALABLES LIES AUX RECHERCHES DE DONNEES GENERALES
!====
!
! RECUPERATION DE T-
    call jevech('PTEMPER', 'L', itempr)
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PTEMPSR', 'L', itemps)
! FLUX NON-LIN
    call jevech('PFLUXNL', 'L', iflux)
    call jevech('PVECTTR', 'E', ivectt)
!
!====
! 1.3 PREALABLES LIES AUX CALCULS
!====
!
    theta = zr(itemps+2)
    coef = zk8(iflux)
    if (coef(1:7) .eq. '&FOZERO') goto 999
!
!    CALCUL DES PRODUITS VECTORIELS OMI   OMJ
!
    do ino = 1, nno
        i = igeom + 3*(ino-1) -1
        do jno = 1, nno
            j = igeom + 3*(jno-1) -1
            sx(ino,jno) = zr(i+2) * zr(j+3) - zr(i+3) * zr(j+2)
            sy(ino,jno) = zr(i+3) * zr(j+1) - zr(i+1) * zr(j+3)
            sz(ino,jno) = zr(i+1) * zr(j+2) - zr(i+2) * zr(j+1)
        end do
    end do
!
!====
! 2. CALCULS TERMES
!====
!    BOUCLE SUR LES POINTS DE GAUSS
    do kp = 1, npg1
!
        kdec = (kp-1)*nno*ndim
        ldec = (kp-1)*nno
        nx = 0.0d0
        ny = 0.0d0
        nz = 0.0d0
!
!   CALCUL DE LA NORMALE AU POINT DE GAUSS KP
!
        do i = 1, nno
            idec = (i-1)*ndim
            do j = 1, nno
                jdec = (j-1)*ndim
                nx = nx+ zr(idfdx+kdec+idec)* zr(idfdy+kdec+jdec)* sx( i,j)
                ny = ny+ zr(idfdx+kdec+idec)* zr(idfdy+kdec+jdec)* sy( i,j)
                nz = nz+ zr(idfdx+kdec+idec)* zr(idfdy+kdec+jdec)* sz( i,j)
            enddo
        enddo
!
!   CALCUL DU JACOBIEN AU POINT DE GAUSS KP
        jac = sqrt(nx*nx + ny*ny + nz*nz)
!
        tpg = 0.d0
        do i = 1, nno
! CALCUL DE T-
            tpg = tpg + zr(itempr+i-1) * zr(ivf+ldec+i-1)
        end do
!
! OBTENTION DES PARAMETRES MATERIAU
        call foderi(coef, tpg, alpha, rbid)
!
        do i = 1, nno
            zr(ivectt+i-1) = zr(ivectt+i-1) + zr(ipoids+kp-1)*jac* (1.d0-theta)*alpha*zr(ivf+ldec&
                             &+i-1)
        end do
!
! FIN BOUCLE SUR LES PT DE GAUSS
    end do
! EXIT SI LE FLUX EST LA FONCTION NULLE
999 continue
! FIN ------------------------------------------------------------------
end subroutine
