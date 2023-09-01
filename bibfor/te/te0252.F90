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
subroutine te0252(option, nomte)
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
#include "asterfort/ppgan2.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcdiff.h"
#include "asterfort/rcfode.h"
#include "asterfort/rcvalb.h"
#include "asterfort/runge6.h"
#include "asterfort/teattr.h"
!
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES VECTEURS RESIDUS
!                          OPTION : 'MASS_THER_RESI'
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
    character(len=32) :: phenom
    real(kind=8) :: beta, deltat, tpg
    real(kind=8) :: dfdx(9), dfdy(9), poids, r, r8bid
    real(kind=8) :: hydrgm(9), hydrgp(9)
    real(kind=8) :: coorse(18), vectt(9), err
    real(kind=8) :: chal(1), tpgm
    character(len=8) :: elrefe, alias8
    integer :: nno, kp, i, j, k, itemps, ifon(6)
    integer :: igeom, imate
    integer :: icomp, itempi, iveres, ipoid2, npg2
    integer :: c(6, 9), ise, nse, nnop2, ivf2, idfde2
    integer :: ibid, jgano2
    integer :: ihydr, ihydrp, itempr
    aster_logical :: aniso
! ----------------------------------------------------------------------
! PARAMETER ASSOCIE AU MATERIAU CODE
!
! --- INDMAT : INDICE SAUVEGARDE POUR LE MATERIAU
!
!C      PARAMETER        ( INDMAT = 8 )
!
! DEB ------------------------------------------------------------------
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
                     npg=npg2, jpoids=ipoid2, jvf=ivf2, jdfde=idfde2, jgano=jgano2)
!
!====
! 1.2 PREALABLES LIES AUX RECHERCHES DE DONNEES GENERALES
!====
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
    call jevech('PTEMPSR', 'L', itemps)
    call jevech('PTEMPEI', 'L', itempi)
    call jevech('PCOMPOR', 'L', icomp)
    call jevech('PRESIDU', 'E', iveres)
!
    deltat = zr(itemps+1)
!
!====
! 1.3 PREALABLES LIES AU SECHAGE
!====
    if (zk16(icomp) (1:5) .eq. 'THER_') then
        call rccoma(zi(imate), 'THER', 1, phenom, icodre(1))
        aniso = .false.
        if (phenom(1:12) .eq. 'THER_NL_ORTH') then
            aniso = .true.
        end if
        call ntfcma(zk16(icomp), zi(imate), aniso, ifon)
!
    end if
!
!====
! 1.5 PREALABLES LIES A L'HYDRATATION
!====
    if (zk16(icomp) (1:9) .eq. 'THER_HYDR') then
        call jevech('PHYDRPM', 'L', ihydr)
        call jevech('PHYDRPP', 'E', ihydrp)
        call jevech('PTEMPER', 'L', itempr)
        call rcvalb('FPG1', 1, 1, '+', zi(imate), &
                    ' ', 'THER_HYDR', 0, ' ', [r8bid], &
                    1, 'CHALHYDR', chal, icodre(1), 1)
        do kp = 1, npg2
            k = nno*(kp-1)
            hydrgm(kp) = 0.d0
            do i = 1, nno
                hydrgm(kp) = hydrgm(kp)+zr(ihydr-1+i)*zr(ivf2+k+i-1)
            end do
        end do
    end if
!====
! 1.6 PREALABLES LIES AUX ELEMENTS LUMPES
!====
!  CALCUL ISO-P2 : ELTS P2 DECOMPOSES EN SOUS-ELTS LINEAIRES
!
    call connec(nomte, nse, nnop2, c)
    do i = 1, nnop2
        vectt(i) = 0.d0
    end do
!
!====
! 2. CALCULS DU TERME DE RIGIDITE DE L'OPTION
!====
! BOUCLE SUR LES SOUS-ELEMENTS
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
!
!====
! 3. CALCULS DU TERME DE MASSE DE L'OPTION
!====
! ------- TERME DE MASSE : 3EME FAMILLE DE PTS DE GAUSS -----------
!
            do i = 1, nno
                do j = 1, 2
                    coorse(2*(i-1)+j) = zr(igeom-1+2*(c(ise, i)-1)+j)
                end do
            end do
!
            do kp = 1, npg2
                k = (kp-1)*nno
                call dfdm2d(nno, kp, ipoid2, idfde2, coorse, &
                            poids, dfdx, dfdy)
                r = 0.d0
                tpg = 0.d0
                do i = 1, nno
                    r = r+coorse(2*(i-1)+1)*zr(ivf2+k+i-1)
                    tpg = tpg+zr(itempi-1+c(ise, i))*zr(ivf2+k+i-1)
                end do
!
! ---  RESOLUTION DE L EQUATION D HYDRATATION
!
                if (zk16(icomp) (1:9) .eq. 'THER_HYDR') then
                    tpgm = 0.d0
                    do i = 1, nno
                        tpgm = tpgm+zr(itempr-1+c(ise, i))*zr(ivf2+k+i-1)
                    end do
                    call runge6(ifon(3), deltat, tpg, tpgm, hydrgm(kp), &
                                hydrgp(kp), err)
                end if
!
                call rcfode(ifon(1), tpg, beta, r8bid)
                if (lteatt('AXIS', 'OUI')) poids = poids*r
                if (zk16(icomp) (1:9) .eq. 'THER_HYDR') then
! --- THERMIQUE NON LINEAIRE AVEC HYDRATATION
                    do i = 1, nno
                        k = (kp-1)*nno
                        vectt(c(ise, i)) = vectt(c(ise, i))+&
                                          & poids*(beta-chal(1)*hydrgp(kp))*&
                                          & zr(ivf2+k+i-1)
                    end do
                else
! --- THERMIQUE NON LINEAIRE SEULE
                    do i = 1, nno
                        vectt(c(ise, i)) = vectt(c(ise, i))+&
                                          & poids*beta*zr(ivf2+k+i-1)
                    end do
                end if
            end do
!
        else if (zk16(icomp) (1:5) .eq. 'SECH_') then
!
! ------- TERME DE MASSE : 3EME FAMILLE DE PTS DE GAUSS -----------
!
            do kp = 1, npg2
                k = (kp-1)*nno
                call dfdm2d(nno, kp, ipoid2, idfde2, coorse, &
                            poids, dfdx, dfdy)
                r = 0.d0
                tpg = 0.d0
                do i = 1, nno
                    r = r+coorse(2*(i-1)+1)*zr(ivf2+k+i-1)
                    tpg = tpg+zr(itempi-1+c(ise, i))*zr(ivf2+k+i-1)
                end do
                if (lteatt('AXIS', 'OUI')) poids = poids*r
!
                do i = 1, nno
                    k = (kp-1)*nno
                    vectt(c(ise, i)) = vectt( &
                                       c(ise, i))+poids*(zr(ivf2+k+i-1)*tpg)
                end do
            end do
!
        end if
!
    end do
!
! MISE SOUS FORME DE VECTEUR
    do i = 1, nnop2
        zr(iveres-1+i) = vectt(i)
    end do
    if (zk16(icomp) (1:9) .eq. 'THER_HYDR') call ppgan2(jgano2, 1, 1, hydrgp, zr(ihydrp))
! FIN ------------------------------------------------------------------
end subroutine
