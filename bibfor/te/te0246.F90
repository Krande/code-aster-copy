! --------------------------------------------------------------------
! Copyright (C) 2019 Christophe Durand - www.code-aster.org
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
!
subroutine te0246(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterc/r8dgrd.h"
#include "asterfort/connec.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/ntfcma.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcdiff.h"
#include "asterfort/rcfode.h"
#include "asterfort/teattr.h"
!
    character(len=16) :: option, nomte
! ----------------------------------------------------------------------
!    - FONCTION REALISEE:  CALCUL DES MATRICES ELEMENTAIRES
!                          OPTION : 'MASS_THER_TANG'
!                          ELEMENTS 2D LUMPES
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
!
! THERMIQUE NON LINEAIRE
!
    integer :: nbres
    parameter(nbres=3)
    integer :: icodre(nbres)
    character(len=8) :: elrefe, alias8
    character(len=32) :: phenom
    real(kind=8) :: r8bid, rhocp
    real(kind=8) :: dfdx(9), dfdy(9), poids, r, tpgi
    real(kind=8) :: mt(9, 9), coorse(18)
    integer ::  nno, kp, i, j, ij, k, ifon(6)
    integer :: igeom, imate
    integer :: icomp, itempi, imattt, ipoid, npg
    integer :: c(6, 9), ise, nse, nnop2, ivf, idfde
    integer :: ibid
    aster_logical :: aniso
!
!====
! 1.1 PREALABLES: RECUPERATION ADRESSES FONCTIONS DE FORMES...
!====
    call elref1(elrefe)
!
    if (lteatt('LUMPE', 'OUI')) then
        call teattr('S', 'ALIAS8', alias8, ibid)
        if (alias8(6:8) .eq. 'QU9') elrefe = 'QU4'
        if (alias8(6:8) .eq. 'TR6') elrefe = 'TR3'
    end if
    call elrefe_info(elrefe=elrefe, fami='MASS', nno=nno, &
                     npg=npg, jpoids=ipoid, jvf=ivf, jdfde=idfde)
!
!====
! 1.2 PREALABLES LIES AUX RECHERCHES DE DONNEES GENERALES
!====
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
    call jevech('PTEMPEI', 'L', itempi)
    call jevech('PCOMPOR', 'L', icomp)
    call jevech('PMATTTR', 'E', imattt)
!====
! 1.4 PREALABLES LIES A L ANISOTROPIE EN THERMIQUE ET RECUPERATION PARAMETRES MATERIAU
!====
    if (zk16(icomp) (1:5) .eq. 'THER_') then
        call rccoma(zi(imate), 'THER', 1, phenom, icodre(1))
        aniso = .false.
        if (phenom(1:12) .eq. 'THER_NL_ORTH') then
            aniso = .true.
        end if
        call ntfcma(zk16(icomp), zi(imate), aniso, ifon)
    end if
!====
! 1.5 PREALABLES LIES AUX ELEMENTS LUMPES
!====
!  CALCUL ISO-P2 : ELTS P2 DECOMPOSES EN SOUS-ELTS LINEAIRES
!
    call connec(nomte, nse, nnop2, c)
!
    do i = 1, nnop2
        do j = 1, nnop2
            mt(i, j) = 0.d0
        end do
    end do
!
!====
! 2. CALCULS DU TERME DE L'OPTION
!====
! ----- 2EME FAMILLE DE PTS DE GAUSS/BOUCLE SUR LES SOUS-ELEMENTS
!
    do ise = 1, nse
!
        do i = 1, nno
            do j = 1, 2
                coorse(2*(i-1)+j) = zr(igeom-1+2*(c(ise, i)-1)+j)
            end do
        end do
!
        if (zk16(icomp) (1:5) .eq. 'THER_') then
! ------- TERME DE MASSE : 3EME FAMILLE DE PTS DE GAUSS -----------
!
            do i = 1, nno
                do j = 1, 2
                    coorse(2*(i-1)+j) = zr(igeom-1+2*(c(ise, i)-1)+j)
                end do
            end do
!
            do kp = 1, npg
                k = (kp-1)*nno
                call dfdm2d(nno, kp, ipoid, idfde, coorse, &
                            poids, dfdx, dfdy)
                r = 0.d0
                tpgi = 0.d0
                do i = 1, nno
                    r = r+coorse(2*(i-1)+1)*zr(ivf+k+i-1)
                    tpgi = tpgi+zr(itempi-1+c(ise, i))*zr(ivf+k+i-1)
                end do
                if (lteatt('AXIS', 'OUI')) poids = poids*r
                call rcfode(ifon(1), tpgi, r8bid, rhocp)
!
                do i = 1, nno
                    do j = 1, nno
                        mt(c(ise, i), c(ise, j)) = mt(c(ise, i), c(ise, j))+&
                                                & poids*rhocp*&
                                                & zr(ivf+k+i-1)*zr(ivf+k+j-1)
                    end do
                end do
            end do
!
! --- SECHAGE
!
        else if (zk16(icomp) (1:5) .eq. 'SECH_') then
!
! ------- TERME DE MASSE : 3EME FAMILLE DE PTS DE GAUSS -----------
!
            do i = 1, nno
                do j = 1, 2
                    coorse(2*(i-1)+j) = zr(igeom-1+2*(c(ise, i)-1)+j)
                end do
            end do
!
            do kp = 1, npg
                k = (kp-1)*nno
                call dfdm2d(nno, kp, ipoid, idfde, coorse, &
                            poids, dfdx, dfdy)
                r = 0.d0
                do i = 1, nno
                    r = r+coorse(2*(i-1)+1)*zr(ivf+k+i-1)
                end do
                if (lteatt('AXIS', 'OUI')) poids = poids*r
!
                do i = 1, nno
!
                    do j = 1, nno
                        mt(c(ise, i), c(ise, j)) = mt(c(ise, i), c(ise, j))+poids*&
                                                &(zr(ivf+k+i-1)*zr(ivf+k+j-1))
                    end do
                end do
            end do
        end if
!
! FIN DE LA BOUCLE SUR LES SOUS-ELEMENTS
!
    end do
!
! MISE SOUS FORME DE VECTEUR
    ij = imattt-1
    do i = 1, nnop2
        do j = 1, i
            ij = ij+1
            zr(ij) = mt(i, j)
        end do
    end do
! FIN ------------------------------------------------------------------
end subroutine
