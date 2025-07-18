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
subroutine te0007(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/bsigmc.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/nbsigm.h"
#include "asterfort/tecach.h"
!
    character(len=16) :: option, nomte
! ----------------------------------------------------------------------
! FONCTION REALISEE:  CALCUL DE L'OPTION FORC_NODA
!                       EN 2D POUR ELEMENTS NON LOCAUX A GRAD. DE DEF.
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
    real(kind=8) :: nharm, bsigm(18), geo(18)
! DEB ------------------------------------------------------------------
! ---- CARACTERISTIQUES DU TYPE D'ELEMENT :
! ---- GEOMETRIE ET INTEGRATION
!      ------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, icomp, jvSief, jvDisp, idfde, igeom, ipoids
    integer(kind=8) :: iretc, ivectu, ivf, jgano, kp, ku
    integer(kind=8) :: n, nbsig, ndim, ndimsi, nno, nnos, npg
!
    real(kind=8) :: zero
!-----------------------------------------------------------------------
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
! --- INITIALISATIONS :
!     -----------------
    zero = 0.0d0
    nharm = zero
!
! - SPECIFICATION DE LA DIMENSION
!
    if (lteatt('AXIS', 'OUI')) then
        ndim = 2
    else if (lteatt('C_PLAN', 'OUI')) then
        ndim = 2
    else if (lteatt('D_PLAN', 'OUI')) then
        ndim = 2
    end if
!
    ndimsi = ndim*2
!
! ---- NOMBRE DE CONTRAINTES ASSOCIE A L'ELEMENT
!      -----------------------------------------
    nbsig = nbsigm()
!
! ---- PARAMETRES EN ENTREE
!      --------------------
! ----     COORDONNEES DES CONNECTIVITES
    call jevech('PGEOMER', 'L', igeom)
! ----     CONTRAINTES AUX POINTS D'INTEGRATION
    call jevech('PSIEFR', 'L', jvSief)
!
!         CHAMPS POUR LA REACTUALISATION DE LA GEOMETRIE
    do i = 1, ndim*nno
        geo(i) = zr(igeom-1+i)
    end do
    call jevech('PDEPLAR', 'L', jvDisp)
    call tecach('ONO', 'PCOMPOR', 'L', iretc, iad=icomp)
    if ((iretc .eq. 0)) then
        if (zk16(icomp+2) (1:6) .ne. 'PETIT ') then
            do i = 1, ndim*nno
                geo(i) = geo(i)+zr(jvDisp-1+i)
            end do
        end if
    end if
! ---- PARAMETRES EN SORTIE
!      --------------------
! ----     VECTEUR DES FORCES INTERNES (BT*SIGMA)
    call jevech('PVECTUR', 'E', ivectu)
!
!
! ---- CALCUL DU VECTEUR DES FORCES INTERNES (BT*SIGMA) :
!      --------------------------------------------------
    call bsigmc(nno, ndim, nbsig, npg, ipoids, &
                ivf, idfde, zr(igeom), nharm, zr(jvSief), &
                bsigm)
!
! ---- AFFECTATION DU VECTEUR EN SORTIE :
!      ----------------------------------
    do n = 1, nnos
        do i = 1, ndim
            ku = (ndimsi+ndim)*(n-1)+i
            kp = ndim*(n-1)+i
            zr(ivectu+ku-1) = bsigm(kp)
        end do
        do i = 1, ndimsi
            ku = (ndimsi+ndim)*(n-1)+i+ndim
            zr(ivectu+ku-1) = 0.d0
        end do
    end do
    do n = nnos+1, nno
        do i = 1, ndim
            ku = (ndimsi+ndim)*nnos+ndim*(n-nnos-1)+i
            kp = ndim*(n-1)+i
            zr(ivectu+ku-1) = bsigm(kp)
        end do
    end do
!
! FIN ------------------------------------------------------------------
end subroutine
