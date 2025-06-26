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
subroutine te0326(option, nomte)
    implicit none
!
!...................................................................
!
! BUT: CALCUL DES VECTEURS ELEMENTAIRES POUR CALCULER LE 2 IEME MEMBRE
!                          DE LA
! FORMULATION VARIATIONNELLE PERMETTANT D'OBTENIR LE POTENTIEL  PHI2
!          DANS LA DECOMPOSITION DU POTENTIEL FLUCTUANT
!          ELEMENTS ISOPARAMETRIQUES 2D
!
!          OPTION : CHAR_THER_PHID_R
!     ENTREES  ---> OPTION : OPTION DE CALCUL
!          ---> NOMTE  : NOM DU TYPE ELEMENT
!..................................................................
!
#include "jeveux.h"
#include "asterfort/divgra.h"
#include "asterfort/e1e2nn.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/rcvalb.h"
#include "asterfort/subacv.h"
#include "asterfort/sumetr.h"
!
    integer(kind=8) :: icodre(1)
    character(len=8) :: fami, poum
    character(len=16) :: nomte, option
    real(kind=8) :: jac(9), nx(9), ny(9), nz(9)
    real(kind=8) :: sx(9, 9), sy(9, 9), sz(9, 9)
    real(kind=8) :: norm(3, 9), rho(1)
    real(kind=8) :: acloc(3, 9), acc(3, 9), flufn(9)
    real(kind=8) :: vibar(2, 9), e1(3, 9), e2(3, 9)
    real(kind=8) :: divsig(9), xin(9), cova(3, 3), metr(2, 2), a(2, 2)
    real(kind=8) :: jc, cnva(3, 3), e1n(3, 9), e2n(3, 9)
    real(kind=8) :: san(9), can(9), dxn(2)
    real(kind=8) :: dxinde(9), dxindk(9)
    real(kind=8) :: dfde(9, 9), dfdk(9, 9), nxn(9), nyn(9), nzn(9)
    real(kind=8) :: normn(3, 9), j1n(9), j2n(9), vibarn(2, 9), gphgxn(9)
    integer(kind=8) :: ipoids, ivf, idfdx, idfdy, igeom
    integer(kind=8) :: ndim, nno, ipg, npg1, ivectt, imate
    integer(kind=8) :: idec, jdec, kdec, ldec, spt
    integer(kind=8) :: i, j, k, iacce, idim, ii, ino, itemp, jno, nnos, jgano
!-----------------------------------------------------------------------
!
!     CALCUL DES DERIVEES PREMIERES DES FONCTIONS DE FORME
!     POUR LES ELEMENTS QUAD4 ET QUAD8
!
    call elrefe_info(fami='NOEU', ndim=ndim, nno=nno, nnos=nnos, npg=npg1, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfdx, jgano=jgano)
    idfdy = idfdx+1
    do ii = 1, nno
        kdec = (ii-1)*nno*ndim
        do j = 1, nno
            idec = (j-1)*ndim
            dfde(j, ii) = zr(idfdx+kdec+idec)
            dfdk(j, ii) = zr(idfdy+kdec+idec)
        end do
    end do
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg1, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfdx, jgano=jgano)
    idfdy = idfdx+1
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
    call jevech('PVECTTR', 'E', ivectt)
    call jevech('PACCELR', 'L', iacce)
    call jevech('PTEMPER', 'L', itemp)
!
    fami = 'RIGI'
    spt = 1
    poum = '+'
!
    do i = 1, nno
        acloc(1, i) = 0.0d0
        acloc(2, i) = 0.0d0
        acloc(3, i) = 0.0d0
    end do
    k = 0
    do i = 1, nno
        do idim = 1, 3
            k = k+1
            acloc(idim, i) = zr(iacce+k-1)
        end do
    end do
    do ipg = 1, npg1
        acc(1, ipg) = 0.0d0
        acc(2, ipg) = 0.0d0
        acc(3, ipg) = 0.0d0
    end do
    do ipg = 1, npg1
        ldec = (ipg-1)*nno
        do i = 1, nno
            acc(1, ipg) = acc(1, ipg)+acloc(1, i)*zr(ivf+ldec+i-1)
            acc(2, ipg) = acc(2, ipg)+acloc(2, i)*zr(ivf+ldec+i-1)
            acc(3, ipg) = acc(3, ipg)+acloc(3, i)*zr(ivf+ldec+i-1)
        end do
    end do
!
    do i = 1, nno
        zr(ivectt+i-1) = 0.d0
    end do
!
! --- CALCUL DES VECTEURS E1, E2 TANGENTS A L'ELEMENT NON NORMALISES
!     ET DES VECTEURS UNITAIRES NORMES
!
    do ipg = 1, npg1
        kdec = (ipg-1)*nno*ndim
!
        e1(1, ipg) = 0.0d0
        e1(2, ipg) = 0.0d0
        e1(3, ipg) = 0.0d0
!
        e2(1, ipg) = 0.0d0
        e2(2, ipg) = 0.0d0
        e2(3, ipg) = 0.0d0
!
        do j = 1, nno
            idec = (j-1)*ndim
!
            e1(1, ipg) = e1(1, ipg)+zr(igeom+3*(j-1)-1+1)*zr(idfdx+ &
                                                             kdec+idec)
            e1(2, ipg) = e1(2, ipg)+zr(igeom+3*(j-1)-1+2)*zr(idfdx+ &
                                                             kdec+idec)
            e1(3, ipg) = e1(3, ipg)+zr(igeom+3*(j-1)-1+3)*zr(idfdx+ &
                                                             kdec+idec)
!
            e2(1, ipg) = e2(1, ipg)+zr(igeom+3*(j-1)-1+1)*zr(idfdy+ &
                                                             kdec+idec)
            e2(2, ipg) = e2(2, ipg)+zr(igeom+3*(j-1)-1+2)*zr(idfdy+ &
                                                             kdec+idec)
            e2(3, ipg) = e2(3, ipg)+zr(igeom+3*(j-1)-1+3)*zr(idfdy+ &
                                                             kdec+idec)
        end do
    end do
!
! --- CALCUL DES PRODUITS VECTORIELS OMI X OMJ
!
    do ino = 1, nno
        i = igeom+3*(ino-1)-1
        do jno = 1, nno
            j = igeom+3*(jno-1)-1
            sx(ino, jno) = zr(i+2)*zr(j+3)-zr(i+3)*zr(j+2)
            sy(ino, jno) = zr(i+3)*zr(j+1)-zr(i+1)*zr(j+3)
            sz(ino, jno) = zr(i+1)*zr(j+2)-zr(i+2)*zr(j+1)
        end do
    end do
!
! --- BOUCLE SUR LES POINTS DE GAUSS
!
    do ipg = 1, npg1
!
        kdec = (ipg-1)*nno*ndim
        ldec = (ipg-1)*nno
!
        nx(ipg) = 0.0d0
        ny(ipg) = 0.0d0
        nz(ipg) = 0.0d0
!
        do i = 1, nno
            idec = (i-1)*ndim
            do j = 1, nno
                jdec = (j-1)*ndim
!
                nx(ipg) = nx(ipg)+zr(idfdx+kdec+idec)*zr(idfdy+kdec+jdec)*sx(i, j)
                ny(ipg) = ny(ipg)+zr(idfdx+kdec+idec)*zr(idfdy+kdec+jdec)*sy(i, j)
                nz(ipg) = nz(ipg)+zr(idfdx+kdec+idec)*zr(idfdy+kdec+jdec)*sz(i, j)
            end do
        end do
!
! ------ CALCUL DU JACOBIEN AU POINT DE GAUSS IPG
!
        jac(ipg) = sqrt(nx(ipg)*nx(ipg)+ny(ipg)*ny(ipg)+nz(ipg)*nz(ipg))
!
! ------ CALCUL DE LA NORMALE UNITAIRE
!
        norm(1, ipg) = nx(ipg)/jac(ipg)
        norm(2, ipg) = ny(ipg)/jac(ipg)
        norm(3, ipg) = nz(ipg)/jac(ipg)
!
    end do
!
! --- CALCUL DU PRODUIT (XI.N) AUX NOEUDS
!
    call e1e2nn(nno, dfde, dfdk, e1n, e2n, &
                nxn, nyn, nzn, normn, j1n, &
                j2n, san, can)
!
    do i = 1, nno
        xin(i) = acloc(1, i)*normn(1, i)+acloc(2, i)*normn(2, i)+acloc( &
                 3, i)*normn(3, i)
    end do
!
    do ipg = 1, npg1
!
! CALCUL DES DERIVEES DES (XIN) (GRADSIGMA(XI.N)) AUX POINTS DE GAUSS
!
        dxinde(ipg) = 0.0d0
        dxindk(ipg) = 0.0d0
!
        kdec = (ipg-1)*nno*ndim
        do i = 1, nno
            idec = (i-1)*ndim
            dxinde(ipg) = dxinde(ipg)+zr(idfdx+kdec+idec)*xin(i)
            dxindk(ipg) = dxindk(ipg)+zr(idfdy+kdec+idec)*xin(i)
        end do
!
        dxn(1) = dxinde(ipg)
        dxn(2) = dxindk(ipg)
!
        vibar(1, ipg) = 0.0d0
        vibar(2, ipg) = 0.0d0
!
! ------ CALCUL DE GRAD(PHIBARRE):VITESSE FLUIDE PERMANENTE
!        VIBAR  AUX POINTS DE GAUSS
!
        do i = 1, nno
            idec = (i-1)*ndim
            vibar(1, ipg) = vibar(1, ipg)+zr(itemp+i-1)*zr(idfdx+kdec+ &
                                                           idec)
            vibar(2, ipg) = vibar(2, ipg)+zr(itemp+i-1)*zr(idfdy+kdec+ &
                                                           idec)
        end do
!
! ------ CONSTITUTION DE LA BASE DES DEUX VECTEURS COVARIANTS
!        AUX POINTS DE GAUSS
!
        do i = 1, 3
            cova(i, 1) = e1(i, ipg)
            cova(i, 2) = e2(i, ipg)
        end do
!
! ------ ON CALCULE LE TENSEUR METRIQUE
!
        call sumetr(cova, metr, jc)
!
! ------ CALCUL DE LA BASE CONTRAVARIANTE
!
        call subacv(cova, metr, jc, cnva, a)
!
! ------ CALCUL DU PRODUIT SCALAIRE GRAD(POTENTIEL PERMANENT)*
!        GRAD(DEPLACEMENT NORMAL) RELATIVEMENT A LA BASE CONTRAVARIANTE
!
        gphgxn(ipg) = 0.d0
        do i = 1, 2
            do j = 1, 2
                gphgxn(ipg) = gphgxn(ipg)+a(i, j)*vibar(i, ipg)*dxn(j)
            end do
        end do
!
    end do
!
!
! --- CALCUL DE LA DIVSIGMA(GRAD(PHIBAR)) AUX POINTS DE GAUSS
!
    call divgra(e1, e2, dfde, dfdk, vibarn, &
                divsig)
!
! --- CALCUL DU FLUX FLUIDE NORMAL AUX POINTS DE GAUSS
!
    do ipg = 1, npg1
        flufn(ipg) = 0.0d0
        flufn(ipg) = acc(1, ipg)*norm(1, ipg)+acc(2, ipg)*norm(2, ipg)+acc(3, ipg)*norm(3, ipg)
    end do
!
! --- CALCUL DU VECTEUR ASSEMBLE SECOND MEMBRE
!     POUR LE CALCUL DU SECOND POTENTIEL INSTATIONNAIRE PHI2
!
    do ipg = 1, npg1
        ldec = (ipg-1)*nno

        call rcvalb(fami, ipg, spt, poum, zi(imate), &
                    ' ', 'THER', 0, ' ', [0.d0], &
                    1, 'RHO_CP', rho, icodre, 1)

        do i = 1, nno
            zr(ivectt+i-1) = zr(ivectt+i-1)+rho(1)*jac(ipg)*zr(ipoids+ipg-1)*zr(ivf+ldec+i-1)&
                            & *(flufn(ipg)*divsig(ipg)+gphgxn(ipg))
!
        end do
    end do
!
end subroutine
