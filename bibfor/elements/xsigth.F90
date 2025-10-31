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
subroutine xsigth(ndim, lonch, time, nbsig, sigth)
!
    use BehaviourStrain_type
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/dmatmc.h"
#include "asterfort/ElasticityMaterial_type.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/epstmc.h"
#include "asterfort/iselli.h"
#include "asterfort/jevech.h"
#include "asterfort/get_elas_id.h"
#include "asterfort/tecach.h"
#include "jeveux.h"
!
    integer(kind=8) :: ndim, nbsig, lonch(10)
    real(kind=8) :: sigth(*), time
!
! --------------------------------------------------------------------------------------------------
!
!      CALCUL DES CONTRAINTES THERMIQUES POUR LES ELEMENTS X-FEM
!
! --------------------------------------------------------------------------------------------------
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
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter :: anglNaut(3) = 0.d0
    integer(kind=8), parameter :: ksp = 1
    real(kind=8) :: epsiTher(6), d(36)
    integer(kind=8) :: nse, idecpg, idebs, iret, kpg, i, ise, npg, j
    integer(kind=8) :: jvMater, irese, nno, ibid, kpgXFEM
    character(len=8) :: elrefp
    integer(kind=8) :: elasID
    type(All_Varc_Strain) :: allVarcStrain
    character(len=8), parameter :: elrese(6) = (/'SE2', 'TR3', 'TE4', 'SE3', 'TR6', 'T10'/)
    character(len=8), parameter :: fami(6) = (/'BID ', 'XINT', 'XINT', 'BID ', 'XINT', 'XINT'/)
!
! --------------------------------------------------------------------------------------------------
!
    call elref1(elrefp)

!   ON AUTORISE UNIQUEMENT L'ISOTROPIE
    call jevech('PMATERC', 'L', jvMater)
    call get_elas_id(zi(jvMater), elasID)
    ASSERT(elasID .eq. ELAS_ISOT)
    call tecach('ONO', 'PCAMASS', 'L', iret, iad=ibid)
    ASSERT(iret .ne. 0)

!   SOUS-ELEMENT DE REFERENCE : RECUP DE NNO ET NPG
    if (.not. iselli(elrefp)) then
        irese = 3
    else
        irese = 0
    end if
    call elrefe_info(elrefe=elrese(ndim+irese), fami=fami(ndim+irese), nno=nno, &
                     npg=npg)

!   RÉCUPÉRATION DE LA SUBDIVISION DE L'ÉLÉMENT EN NSE SOUS ELEMENT
    nse = lonch(1)

!   BOUCLE SUR LES NSE SOUS-ELEMENTS
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

! ----- Loop on XFEM Gauss points
        do kpgXFEM = 1, npg
            kpg = idecpg+kpgXFEM
!
!         CALCUL DES DEFORMATIONS THERMIQUES EPSTH
            epsiTher = 0.d0
            call epstmc('XFEM', '+', kpg, ksp, ndim, &
                        time, anglNaut, zi(jvMater), &
                        VARC_STRAIN_TEMP, allVarcStrain, &
                        epsiTher)
!
!         CALCUL DE LA MATRICE DE HOOKE (MATERIAU ISOTROPE)
            d = 0.d0
            call dmatmc('XFEM', zi(jvMater), time, '+', &
                        kpg, ksp, anglNaut, nbsig, &
                        d)
!
!         CONTRAINTES THERMIQUES AU PG COURANT
            do i = 1, nbsig
                do j = 1, nbsig
                    sigth(idebs+nbsig*(kpgXFEM-1)+i) = sigth(idebs+nbsig*(kpgXFEM-1)+i)+ &
                                                       d(j+(i-1)*nbsig)*epsiTher(j)
                end do
            end do
        end do
    end do
!
end subroutine
