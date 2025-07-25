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
subroutine te0436(option, nomte)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/r8dgrd.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/mbcine.h"
#include "asterfort/mbrigi.h"
#include "asterfort/getDensity.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/verift.h"
!
    character(len=16) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
!    - FONCTION REALISEE:  CALCUL DES OPTIONS DE POST-TRAITEMENT :
!                                  - SIEF_ELGA
!                                  - EFGE_ELGA
!                                  - EPOT_ELEM
!                                  - EPSI_ELGA
!                                  - MASS_INER
!                          POUR LES MEMBRANES
!    - ARGUMENTS :
!        DONNEES :      OPTION       -->  OPTION DE CALCUL
!                       NOMTE        -->  NOM DU TYPE ELEMENT
!
! --------------------------------------------------------------------------------------------------
!
    character(len=4), parameter :: fami = 'RIGI'
    integer(kind=8) :: nddl, nno, npg, ncomp
    integer(kind=8) :: i, j, n, c, cc, kpg, j1, j2, k
    integer(kind=8) :: ipoids, ivf, idfde, iret
    integer(kind=8) :: igeom, icacoq, imate, idepl, icontp, inr, idefo, imass, icompo
    real(kind=8) :: dff(2, 9), vff(9), b(3, 3, 9), jac
    real(kind=8) :: alpha, beta, epot
    real(kind=8) :: epsm(3), epsg(3, 9), epsthe, sig(3), sigg(3, 9), rig(3, 3)
    real(kind=8) :: rho, rhog
    real(kind=8) :: x(9), y(9), z(9), surfac, cdg(3), ppg, xxi, yyi, zzi
    real(kind=8) :: matine(6)
    real(kind=8) :: vro
!
! --------------------------------------------------------------------------------------------------
!
    ncomp = 3
    nddl = 3
!
! - FONCTIONS DE FORMES ET POINTS DE GAUSS
!
    call elrefe_info(fami=fami, nno=nno, &
                     npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfde)

    if (option .eq. 'EFGE_ELGA') then
! ---   c'est facile : il n'y a qu'a recopier :
        call jevech('PSIEFR', 'L', j1)
        call jevech('PEFGER', 'E', j2)
        do k = 1, 3*npg
            zr(j2-1+k) = zr(j1-1+k)
        end do
        goto 999
    end if
!
! - PARAMETRES EN ENTREE
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PCACOQU', 'L', icacoq)
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
    end if

! - ON INTERDIT CERTAINES OPTIONS POUR LES GRANDES DEFORMATIONS
    call tecach('NNO', 'PCOMPOR', 'L', iret, iad=icompo)
    if (((option .eq. 'EPSI_ELGA') .or. (option .eq. 'EPOT_ELEM')) .and. &
        (iret .eq. 0) .and. (zk16(icompo+2) (1:9) .eq. 'GROT_GDEP')) then
        call utmess('F', 'MEMBRANE_8', sk=option)
    end if

!
! - LE VECTEUR NORME QUI DETERMINE LE REPERE LOCAL DE LA MEMBRANE
!   (COMPORTEMENT ANISOTROPE)
!
    alpha = zr(icacoq+1)*r8dgrd()
    beta = zr(icacoq+2)*r8dgrd()
!
! - COORDONNEES PHYSIQUES DES NOEUDS
!
    if (option .eq. 'MASS_INER') then
        do i = 1, nno
            x(i) = zr(igeom+3*(i-1))
            y(i) = zr(igeom+3*i-2)
            z(i) = zr(igeom+3*i-1)
        end do
        cdg = 0.d0
        matine = 0.d0
        surfac = 0.d0
        rhog = 0.d0
    end if
!
!
! - DEBUT DE LA BOUCLE SUR LES POINTS DE GAUSS
!
    do kpg = 1, npg
!
! --- MISE SOUS FORME DE TABLEAU DES VALEURS ET DES DERIVEES
!     DES FONCTIONS DE FORME
!
        do n = 1, nno
            vff(n) = zr(ivf+(kpg-1)*nno+n-1)
            dff(1, n) = zr(idfde+(kpg-1)*nno*2+(n-1)*2)
            dff(2, n) = zr(idfde+(kpg-1)*nno*2+(n-1)*2+1)
        end do
!
! --- CALCUL DE LA MATRICE "B" :
!              DEPL NODAL --> DEFORMATIONS MEMBRANAIRES ET JACOBIEN
!
        call mbcine(nno, zr(igeom), dff, alpha, beta, &
                    b, jac)
!
! --- SIEF_ELGA, EPOT_ELEM : ON CALCULE LA CONTRAINTE AU PG
!
        if ((option .eq. 'SIEF_ELGA') .or. (option .eq. 'EPOT_ELEM')) then
!
! ------    CALCUL DE LA DEFORMATION MEMBRANAIRE DANS LE REPERE LOCAL
            epsm = 0.d0
            do n = 1, nno
                do i = 1, nddl
                    do c = 1, ncomp
                        epsm(c) = epsm(c)+b(c, i, n)*zr(idepl+(n-1)*nddl+i-1)
                    end do
                end do
            end do
!
! ------    RETRAIT DE LA DEFORMATION THERMIQUE
            call verift(fami, kpg, 1, '+', zi(imate), &
                        epsth_=epsthe)
            epsm(1) = epsm(1)-epsthe
            epsm(2) = epsm(2)-epsthe
!
! ------    CALCUL DE LA CONTRAINTE AU PG
            call mbrigi(fami, kpg, imate, rig)

            sig = 0.d0
            do c = 1, ncomp
                do cc = 1, ncomp
                    sig(c) = sig(c)+epsm(cc)*rig(cc, c)
                end do
            end do
!
            if (option .eq. 'EPOT_ELEM') then
                do c = 1, ncomp
                    epot = epot+(sig(c)*epsm(c)*zr(ipoids+kpg-1)*jac)/2
                end do
            else
                do c = 1, ncomp
                    sigg(c, kpg) = sig(c)
                end do
            end if
!
! --- EPSI_ELGA : ON CALCULE LA DEFORMATION AU PG
!
        else if (option .eq. 'EPSI_ELGA') then
!
! ------    CALCUL DE LA DEFORMATION MEMBRANAIRE DANS LE REPERE LOCAL
            do n = 1, nno
                do i = 1, nddl
                    do c = 1, ncomp
                        epsg(c, kpg) = epsg(c, kpg)+b(c, i, n)*zr(idepl+( &
                                                                  n-1)*nddl+i-1)
                    end do
                end do
            end do
!
! --- MASS_INER : ON SOMME LA CONTRIBUTION DU PG A LA MASSE TOTALE
!
        else if (option .eq. 'MASS_INER') then
            call getDensity(zi(imate), rho, 'ELAS_MEMBRANE')
            surfac = surfac+zr(ipoids+kpg-1)*jac
            rhog = rhog+rho*zr(ipoids+kpg-1)*jac
            ppg = zr(ipoids+kpg-1)*jac
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
            do c = 1, ncomp
                zr(icontp+(kpg-1)*ncomp+c-1) = sigg(c, kpg)
            end do
        end do
!
    else if (option .eq. 'EPOT_ELEM') then
        zr(inr) = epot
!
    else if (option .eq. 'EPSI_ELGA') then
        do kpg = 1, npg
            do c = 1, ncomp
                zr(idefo+(kpg-1)*ncomp+c-1) = epsg(c, kpg)
            end do
        end do
!
    else if (option .eq. 'MASS_INER') then
        vro = rhog/surfac
        zr(imass) = rhog*surfac
        zr(imass+1) = cdg(1)/surfac
        zr(imass+2) = cdg(2)/surfac
        zr(imass+3) = cdg(3)/surfac
        zr(imass+4) = matine(1)*rhog-vro*(cdg(2)*cdg(2)+cdg(3)*cdg(3))
        zr(imass+5) = matine(3)*rhog-vro*(cdg(1)*cdg(1)+cdg(3)*cdg(3))
        zr(imass+6) = matine(6)*rhog-vro*(cdg(1)*cdg(1)+cdg(2)*cdg(2))
        zr(imass+7) = matine(2)*rhog-vro*(cdg(1)*cdg(2))
        zr(imass+8) = matine(4)*rhog-vro*(cdg(1)*cdg(3))
        zr(imass+9) = matine(5)*rhog-vro*(cdg(2)*cdg(3))
    end if

999 continue

end subroutine
