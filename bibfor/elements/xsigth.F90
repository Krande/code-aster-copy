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

subroutine xsigth(ndim, lonch, inst, nbsig, sigth)
!
!
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dmatmc.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/epstmc.h"
#include "asterfort/iselli.h"
#include "asterfort/jevech.h"
#include "asterfort/rccoma.h"
#include "asterfort/tecach.h"
!
    integer(kind=8) :: ndim, nbsig, lonch(10)
    real(kind=8) :: sigth(*), inst
!
! ----------------------------------------------------------------------
! FONCTION REALISEE:
!
!      CALCUL DES CONTRAINTES THERMIQUES POUR LES ELEMENTS X-FEM
!
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  PINTT   : COORDONNÉES DES POINTS D'INTERSECTION
! IN  LONCH   : LONGUEURS DES CHAMPS UTILISÉES
! IN  SIGMA   : CONTRAINTES DE CAUCHY AUX POINTS DE GAUSS DES SOUS-ÉLTS
! IN  NBSIG   : NOMBRE DE CONTRAINTES ASSOCIE A L'ELEMENT
! IN  PMILT   : COORDONNEES DES POINTS MILIEUX
! IN  INST    : INSTANT
! IN  NBSIG   : DIMENSION DU TENSEUR DES CONTRAINTES
!
! OUT SIGTH   : CONTRAINTES THERMIQUES
!
! ----------------------------------------------------------------------
!
    real(kind=8) :: r8bi3(3), epsth(6), d(36)
    integer(kind=8) :: nse, idecpg, idebs, iret, ipg, i, ise, npg, j
    integer(kind=8) :: imate, irese, nno, ibid, kpg
    character(len=8) :: elrefp, elrese(6), fami(6)
    character(len=16) :: phenom, option
!
    parameter(option='CHAR_MECA_TEMP_R')
    data elrese/'SE2', 'TR3', 'TE4', 'SE3', 'TR6', 'T10'/
    data fami/'BID', 'XINT', 'XINT', 'BID', 'XINT', 'XINT'/
!
! ----------------------------------------------------------------------
!
    call elref1(elrefp)
!
!     INITIALISATION DES VECTEURS BIDONS
    r8bi3(:) = 0.d0
!
!     ON AUTORISE UNIQUEMENT L'ISOTROPIE
    call jevech('PMATERC', 'L', imate)
    call rccoma(zi(imate), 'ELAS', 1, phenom, iret)
    ASSERT(iret .eq. 0 .and. phenom .eq. 'ELAS')
    call tecach('ONO', 'PCAMASS', 'L', iret, iad=ibid)
    ASSERT(iret .ne. 0)
!
!     SOUS-ELEMENT DE REFERENCE : RECUP DE NNO ET NPG
    if (.not. iselli(elrefp)) then
        irese = 3
    else
        irese = 0
    end if
    call elrefe_info(elrefe=elrese(ndim+irese), fami=fami(ndim+irese), nno=nno, &
                     npg=npg)
!
!     RÉCUPÉRATION DE LA SUBDIVISION DE L'ÉLÉMENT EN NSE SOUS ELEMENT
    nse = lonch(1)
!
!       BOUCLE SUR LES NSE SOUS-ELEMENTS
    do ise = 1, nse
!
!       DEBUT DE LA ZONE MÉMOIRE DE SIGMA CORRESPONDANTE
        idecpg = npg*(ise-1)
        idebs = nbsig*idecpg
        if (ndim .eq. 3) then
            ASSERT(nbsig .eq. 6)
        else if (ndim .eq. 2) then
            ASSERT(nbsig .eq. 4)
        end if
!
!       BOUCLE SUR LES POINTS DE GAUSS DU SOUS-ELEMENT
        do kpg = 1, npg
!
!         NUMERO DU PG DANS LE FAMILLE 'XFEM'
            ipg = idecpg+kpg
!
!         CALCUL DES DEFORMATIONS THERMIQUES EPSTH
            epsth(:) = 0.d0
            call epstmc('XFEM', ndim, inst, '+', ipg, &
                        1, r8bi3, zi(imate), option, &
                        epsth)
!
!         CALCUL DE LA MATRICE DE HOOKE (MATERIAU ISOTROPE)
            d(:) = 0.d0
            call dmatmc('XFEM', zi(imate), inst, '+', &
                        ipg, 1, r8bi3, nbsig, &
                        d)
!
!         CONTRAINTES THERMIQUES AU PG COURANT
            do i = 1, nbsig
                do j = 1, nbsig
                    sigth(idebs+nbsig*(kpg-1)+i) = sigth( &
                                                   idebs+nbsig*(kpg-1)+i)+d(j+(i-1)*nbsig)*epsth&
                                                   &(j &
                                                   )
                end do
            end do
        end do
    end do
!
end subroutine
