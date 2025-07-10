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
subroutine te0433(option, nomte)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/cargri.h"
#include "asterfort/dxqpgl.h"
#include "asterfort/dxtpgl.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/getDensity.h"
#include "asterfort/jevech.h"
#include "asterfort/nmgrib.h"
#include "asterfort/lteatt.h"
#include "asterfort/rcvalb.h"
#include "asterfort/verift.h"
!
    character(len=16) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
!    - FONCTION REALISEE:  CALCUL DES OPTIONS DE POST-TRAITEMENT :
!                                  - SIEF_ELGA
!                                  - EPOT_ELEM
!                                  - EPSI_ELGA
!                                  - MASS_INER
!                          POUR LES GRILLES MEMBRANES EXCENTREES OU NON
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: codres(2)
    character(len=4), parameter :: fami = 'RIGI'
    character(len=16) :: nomres(2)
    integer(kind=8) :: nddl, nno, npg, i, j, n, kpg
    integer(kind=8) :: ipoids, ivf, idfde, igeom, imate
    integer(kind=8) :: icontp, imass, idepl, idefo, inr
    real(kind=8) :: dff(2, 8), vff(8), b(6, 8), p(3, 6), jac
    real(kind=8) :: dir11(3), densit, pgl(3, 3), distn, vecn(3)
    real(kind=8) :: epsm, epsg(9), epsthe, sig, sigg(9), rho, valres(2), epot
    real(kind=8) :: x(8), y(8), z(8), volume, cdg(3), ppg, xxi, yyi, zzi
    real(kind=8) :: matine(6), vro
    aster_logical :: lexc
!
! --------------------------------------------------------------------------------------------------
!

    lexc = (lteatt('MODELI', 'GRC'))
!
! - FONCTIONS DE FORMES ET POINTS DE GAUSS
!
    call elrefe_info(fami=fami, nno=nno, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde)
!
! - PARAMETRES EN ENTREE
!
    call jevech('PGEOMER', 'L', igeom)
!
    if ((option .eq. 'SIEF_ELGA') .or. (option .eq. 'EPOT_ELEM')) then
        call jevech('PDEPLAR', 'L', idepl)
        call jevech('PMATERC', 'L', imate)
!
    else if (option .eq. 'EPSI_ELGA') then
        call jevech('PDEPLAR', 'L', idepl)
        epsg = 0.d0
!
    else if (option .eq. 'MASS_INER') then
        call jevech('PMATERC', 'L', imate)
    end if
!
! - PARAMETRES EN SORTIE
!
    if (option .eq. 'SIEF_ELGA') then
        call jevech('PCONTRR', 'E', icontp)
        sigg = 0.d0
!
    else if (option .eq. 'EPOT_ELEM') then
        call jevech('PENERDR', 'E', inr)
        epot = 0.d0
!
    else if (option .eq. 'EPSI_ELGA') then
        call jevech('PDEFOPG', 'E', idefo)
!
    else if (option .eq. 'MASS_INER') then
        call jevech('PMASSINE', 'E', imass)
        call getDensity(zi(imate), rho)
    end if
!
! - LECTURE DES CARACTERISTIQUES DE GRILLE ET
!   CALCUL DE LA DIRECTION D'ARMATURE
!
    call cargri(lexc, densit, distn, dir11)
!
! - SI EXCENTREE : RECUPERATION DE LA NORMALE ET DE L'EXCENTREMENT
!
    if (lexc) then
!
        if (nomte .eq. 'MEGCTR3') then
            call dxtpgl(zr(igeom), pgl)
        else if (nomte .eq. 'MEGCQU4') then
            call dxqpgl(zr(igeom), pgl)
        end if
!
        do i = 1, 3
            vecn(i) = distn*pgl(3, i)
        end do
        nddl = 6
!
    else
        nddl = 3
    end if
!
! - COORDONNEES PHYSIQUES DES NOEUDS
!
    if (option .eq. 'MASS_INER') then
        do i = 1, nno
            x(i) = zr(igeom+3*(i-1))
            y(i) = zr(igeom+3*i-2)
            z(i) = zr(igeom+3*i-1)
        end do
        if (lexc) then
            x(i) = x(i)+vecn(1)
            y(i) = y(i)+vecn(2)
            z(i) = z(i)+vecn(3)
        end if
        cdg = 0.d0
        matine = 0.d0
    end if
!
    volume = 0.d0
!
! - DEBUT DE LA BOUCLE SUR LES POINTS DE GAUSS
!
    do kpg = 1, npg
!
! --- MISE SOUS FORME DE TABLEAU DES VALEURS DES FONCTIONS DE FORME
!     ET DES DERIVEES DE FONCTION DE FORME
!
        do n = 1, nno
            vff(n) = zr(ivf+(kpg-1)*nno+n-1)
            dff(1, n) = zr(idfde+(kpg-1)*nno*2+(n-1)*2)
            dff(2, n) = zr(idfde+(kpg-1)*nno*2+(n-1)*2+1)
        end do
!
! --- CALCUL DE LA MATRICE "B" : DEPL NODAL --> EPS11 ET DU JACOBIEN
!
        call nmgrib(nno, zr(igeom), dff, dir11, lexc, &
                    vecn, b, jac, p)
!
! --- SIEF_ELGA, EPOT_ELEM : ON CALCULE LA CONTRAINTE AU PG
!
        if ((option .eq. 'SIEF_ELGA') .or. (option .eq. 'EPOT_ELEM')) then
!
!         CALCUL DE LA DEFORMATION AU PG
            epsm = 0.d0
            do i = 1, nno
                do j = 1, nddl
                    epsm = epsm+b(j, i)*zr(idepl+(i-1)*nddl+j-1)
                end do
            end do
!
!         CALCUL DE LA CONTRAINTE
            call verift(fami, kpg, 1, '+', zi(imate), &
                        epsth_=epsthe)
            nomres(1) = 'E'
            call rcvalb(fami, kpg, 1, '+', zi(imate), &
                        ' ', 'ELAS', 0, ' ', [0.d0], &
                        1, nomres, valres, codres, 0)
            epsm = epsm-epsthe
            sig = valres(1)*epsm
!
            if (option .eq. 'EPOT_ELEM') then
                epot = epot+(sig*epsm*zr(ipoids+kpg-1)*jac*densit)/2
            else
                sigg(kpg) = sig
            end if
!
! --- EPSI_ELGA : ON CALCULE LA DEFORMATION AU PG
!
        else if (option .eq. 'EPSI_ELGA') then
!
!         CALCUL DE LA DEFORMATION AU PG
            do i = 1, nno
                do j = 1, nddl
                    epsg(kpg) = epsg(kpg)+b(j, i)*zr(idepl+(i-1)*nddl+j-1)
                end do
            end do
!
! --- MASS_INER : ON SOMME LA CONTRIBUTION DU PG A LA MASSE TOTALE
!
        else if (option .eq. 'MASS_INER') then
            volume = volume+zr(ipoids+kpg-1)*densit*jac
            ppg = zr(ipoids+kpg-1)*jac*densit
            do i = 1, nno
                cdg(1) = cdg(1)+ppg*vff(i)*x(i)
                cdg(2) = cdg(2)+ppg*vff(i)*y(i)
                cdg(3) = cdg(3)+ppg*vff(i)*z(i)
                xxi = 0.d0
                yyi = 0.d0
                zzi = 0.d0
                do j = 1, nno
                    xxi = xxi+x(i)*vff(i)*vff(j)*x(j)
                    yyi = yyi+y(i)*vff(i)*vff(j)*y(j)
                    zzi = zzi+z(i)*vff(i)*vff(j)*z(j)
                    matine(2) = matine(2)+x(i)*vff(i)*vff(j)*y(j)*ppg
                    matine(4) = matine(4)+x(i)*vff(i)*vff(j)*z(j)*ppg
                    matine(5) = matine(5)+y(i)*vff(i)*vff(j)*z(j)*ppg
                end do
                matine(1) = matine(1)+ppg*(yyi+zzi)
                matine(3) = matine(3)+ppg*(xxi+zzi)
                matine(6) = matine(6)+ppg*(xxi+yyi)
            end do
        end if
!
! - FIN DE LA BOUCLE SUR LES POINTS DE GAUSS
    end do
!
    if (option .eq. 'SIEF_ELGA') then
        do kpg = 1, npg
            zr(icontp+kpg-1) = sigg(kpg)
        end do
!
    else if (option .eq. 'EPOT_ELEM') then
        zr(inr) = epot
!
    else if (option .eq. 'EPSI_ELGA') then
        do kpg = 1, npg
            zr(idefo+kpg-1) = epsg(kpg)
        end do
!
    else if (option .eq. 'MASS_INER') then
        vro = rho/volume
        zr(imass) = rho*volume
        zr(imass+1) = cdg(1)/volume
        zr(imass+2) = cdg(2)/volume
        zr(imass+3) = cdg(3)/volume
        zr(imass+4) = matine(1)*rho-vro*(cdg(2)*cdg(2)+cdg(3)*cdg(3))
        zr(imass+5) = matine(3)*rho-vro*(cdg(1)*cdg(1)+cdg(3)*cdg(3))
        zr(imass+6) = matine(6)*rho-vro*(cdg(1)*cdg(1)+cdg(2)*cdg(2))
        zr(imass+7) = matine(2)*rho-vro*(cdg(1)*cdg(2))
        zr(imass+8) = matine(4)*rho-vro*(cdg(1)*cdg(3))
        zr(imass+9) = matine(5)*rho-vro*(cdg(2)*cdg(3))
    end if
!
end subroutine
