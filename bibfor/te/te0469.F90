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
subroutine te0469(option, nomte)
!.......................................................................
    implicit none
!
!     BUT: CONSTRUCTION DU VECTEUR DES FORCES CORRESPONDANT A UN
!          CHARGEMENT FORCE_ARETE POUR LES ELEMENTS ISOPARAMETRIQUES 3D.
!          (I.E. CALCUL DU VECTEUR DES FORCES NODALES EQUIVALENTES
!                A UNE FORCE LINEIQUE LE LONG D'UNE ARETE POUR
!                LES ELEMENTS 3D).
!
!          OPTIONS : 'CHAR_MECA_FR1D3D'
!                    'CHAR_MECA_FF1D3D'
!
!     ENTREES  ---> OPTION : OPTION DE CALCUL
!              ---> NOMTE  : NOM DU TYPE ELEMENT
!.......................................................................
!
#include "jeveux.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
#include "asterfort/vff3d.h"
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, idecno, idecpg, idfdk, idflin, ier, igau
    integer(kind=8) :: igeom, ino, ipoids, itemps, ivectu, ivf, jgano
    integer(kind=8) :: nbnomx, ndim, nno, nnos, npg
    real(kind=8) :: zero
!-----------------------------------------------------------------------
    parameter(nbnomx=27)
    character(len=8) :: nompar(4)
    character(len=16) :: nomte, option
    real(kind=8) :: fx(nbnomx), fy(nbnomx), fz(nbnomx), valpar(4)
    real(kind=8) :: xyzgau(3), forlin(3), jacob
    real(kind=8) :: fxlin(5), fylin(5), fzlin(5)
!
!
!
! --- CARACTERISTIQUES DU TYPE D'ELEMENT :
! --- INTEGRATION ET INTERPOLATION
!      ----------------------------
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfdk, jgano=jgano)
!
! --- INITIALISATIONS :
!     -----------------
    zero = 0.0d0
!
    do i = 1, nbnomx
        fx(i) = zero
        fy(i) = zero
        fz(i) = zero
    end do
!
    do i = 1, npg
        fxlin(i) = zero
        fylin(i) = zero
        fzlin(i) = zero
    end do
!
! --- RECUPERATION DES COORDONNEES DES CONNECTIVITES :
!     ----------------------------------------------
    call jevech('PGEOMER', 'L', igeom)
!
! --- OPTION 'CHAR_MECA_FR1D3D'
! --- CAS OU LES DONNEES DES FORCES LINEIQUES SONT DES REELS :
!     ------------------------------------------------------
    if (option .eq. 'CHAR_MECA_FR1D3D') then
!
! ---  RECUPERATION DES VALEURS DE LA FORCE LINEIQUE A APPLIQUER
! ---  SUR L'ELEMENT D'ARETE :
!      ---------------------
        call jevech('PFR1D3D', 'L', idflin)
!
! ---  BOUCLE SUR LES POINTS D'INTEGRATION
!      -----------------------------------
        do igau = 1, npg
!
            idecpg = nno*(igau-1)
!
!
! ---    CALCUL DE LA FORCE LINEIQUE AUX POINTS D'INTEGRATION :
!        -----------------------------------------------------
            do ino = 1, nno
!
                fxlin(igau) = fxlin(igau)+zr(ivf+idecpg+ino-1)*zr(idflin+1-1)
                fylin(igau) = fylin(igau)+zr(ivf+idecpg+ino-1)*zr(idflin+2-1)
                fzlin(igau) = fzlin(igau)+zr(ivf+idecpg+ino-1)*zr(idflin+3-1)
            end do
        end do
!
! ---    BOUCLE SUR LES POINTS D'INTEGRATION
!        -----------------------------------
        do igau = 1, npg
!
            idecpg = nno*(igau-1)
!
! ---    CALCUL DU PRODUIT JACOBIEN*POIDS AU POINT D'INTEGRATION
! ---    COURANT :
!        -------
            call vff3d(nno, zr(ipoids+igau-1), zr(idfdk+idecpg), zr(igeom), jacob)
!
! ---    CALCUL DE LA CONTRIBUTION AU VECTEUR DES FORCES NODALES
! ---    DU CHARGEMENT LINEIQUE AU POINT D'INTEGRATION COURANT :
!        -----------------------------------------------------
            do ino = 1, nno
!
                fx(ino) = fx(ino)+zr(ivf+idecpg+ino-1)*fxlin(igau)*jacob
                fy(ino) = fy(ino)+zr(ivf+idecpg+ino-1)*fylin(igau)*jacob
                fz(ino) = fz(ino)+zr(ivf+idecpg+ino-1)*fzlin(igau)*jacob
            end do
!
        end do
!
! --- OPTION 'CHAR_MECA_FF1D3D'
! --- CAS OU LES DONNEES DES FORCES LINEIQUES SONT DES FONCTIONS
! --- DES COORDONNEES ET DU TEMPS :
!     ---------------------------
    else if (option .eq. 'CHAR_MECA_FF1D3D') then
!
! ---  RECUPERATION DES NOMS DES FONCTIONS REPRESENTANT LA FORCE
! ---  LINEIQUE A APPLIQUER SUR L'ELEMENT D'ARETE :
!      ------------------------------------------
        call jevech('PFF1D3D', 'L', idflin)
!
! ---  RECUPERATION DE L'INSTANT D'INTERPOLATION :
!      -----------------------------------------
        call jevech('PINSTR', 'L', itemps)
!
! ---  AFFECTATION DES VARIABLES POUR L'INTERPOLATION :
!      ----------------------------------------------
        nompar(1) = 'X'
        nompar(2) = 'Y'
        nompar(3) = 'Z'
        nompar(4) = 'INST'
!
        valpar(4) = zr(itemps)
!
! ---  BOUCLE SUR LES POINTS D'INTEGRATION
!      -----------------------------------
        do igau = 1, npg
!
            idecpg = nno*(igau-1)
!
            do i = 1, 3
                xyzgau(i) = zero
                forlin(i) = zero
            end do
!
! ---    CALCUL DES COORDONNEES DU POINT D'INTEGRATION COURANT :
!        -----------------------------------------------------
            do ino = 1, nno
!
                idecno = 3*(ino-1)-1
!
                xyzgau(1) = xyzgau(1)+zr(ivf+idecpg+ino-1)*zr(igeom+1+idecno)
                xyzgau(2) = xyzgau(2)+zr(ivf+idecpg+ino-1)*zr(igeom+2+idecno)
                xyzgau(3) = xyzgau(3)+zr(ivf+idecpg+ino-1)*zr(igeom+3+idecno)
            end do
!
! ---    INTERPOLATION DES FORCES LINEIQUES EN FONCTION DES
! ---    COORDONNEES ET DU TEMPS :
!        -----------------------
            valpar(1) = xyzgau(1)
            valpar(2) = xyzgau(2)
            valpar(3) = xyzgau(3)
!
            call fointe('FM', zk8(idflin+1-1), 4, nompar, valpar, &
                        forlin(1), ier)
            call fointe('FM', zk8(idflin+2-1), 4, nompar, valpar, &
                        forlin(2), ier)
            call fointe('FM', zk8(idflin+3-1), 4, nompar, valpar, &
                        forlin(3), ier)
!
! ---    CALCUL DU PRODUIT JACOBIEN*POIDS AU POINT D'INTEGRATION
! ---    COURANT :
!        -------
            call vff3d(nno, zr(ipoids+igau-1), zr(idfdk+idecpg), zr(igeom), jacob)
!
! ---    CALCUL DE LA CONTRIBUTION AU VECTEUR DES FORCES NODALES
! ---    DU CHARGEMENT LINEIQUE AU POINT D'INTEGRATION COURANT :
!        -----------------------------------------------------
            do ino = 1, nno
!
                fx(ino) = fx(ino)+zr(ivf+idecpg+ino-1)*forlin(1)*jacob
                fy(ino) = fy(ino)+zr(ivf+idecpg+ino-1)*forlin(2)*jacob
                fz(ino) = fz(ino)+zr(ivf+idecpg+ino-1)*forlin(3)*jacob
            end do
!
        end do
!
    end if
!
! ---- RECUPERATION ET AFFECTATION DU VECTEUR FORCES NODALES EN SORTIE :
!      ---------------------------------------------------------------
    call jevech('PVECTUR', 'E', ivectu)
!
    do ino = 1, nno
        zr(ivectu+3*(ino-1)+1-1) = fx(ino)
        zr(ivectu+3*(ino-1)+2-1) = fy(ino)
        zr(ivectu+3*(ino-1)+3-1) = fz(ino)
    end do
!
end subroutine
