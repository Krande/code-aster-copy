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
subroutine xrmes2(ndim, nbnase, cpt, in, ivois, &
                  jsigse, nno, nbcmp, jcnset, dsg11, &
                  dsg22, dsg12)
! aslint: disable=W1306
    implicit none
#include "jeveux.h"
    integer(kind=8) :: ndim, nbnase, cpt, in, ivois, jsigse, nno, nbcmp, jcnset
    real(kind=8) :: dsg11(nbnase), dsg22(nbnase), dsg12(nbnase)
!  XFEM - ERREUR EN MECANIQUE - TERME DE SAUT - DIMENSION 2
!  *       *        **                   *                *
! =====================================================================
!
!     BUT:
!         CALCUL DU SAUT DE CONTRAINTE (AUX NOEUDS) ENTRE LE
!         SOUS-ELEMENT XFEM COURANT ET LE SOUS-ELEMENT VOISIN PAR
!         RAPPORT AU BORD (ARETE) COURANT.
!
!
!     ARGUMENTS:
!     ----------
!
!      ENTREE :
!-------------
! IN   NDIM   : DIMENSION
! IN   NBNASE : NOMBRE DE NOEUDS PAR ARRÊTE DU SOUS-ELEMENT
! IN   CPT    : NUMERO DU SOUS-ELEMENT
! IN   IN     : NUMERO DU NOEUD (I.E. DE L'ARETE)
! IN   IVOIS  : ADRESSE DANS ZI DES VOISINS DU SOUS-ELEMENT
! IN   JSIGSE : ADRESSE DES CONTRAINTES AUX NOEUDS DU SOUS-ELEMENT
! IN   NBCMP  : NOMBRE DE COMPOSANTE DE SIGMA
! IN   JCNSET : ADRESSE DANS ZI DE LA CONNECTIVITE
!
!      SORTIE :
!-------------
! OUT  DSG11  : SAUT DE CONTRAINTE AUX NOEUDS COMPOSANTE 11
! OUT  DSG22  : SAUT DE CONTRAINTE AUX NOEUDS COMPOSANTE 22
! OUT  DSG12  : SAUT DE CONTRAINTE AUX NOEUDS COMPOSANTE 12
!
! ......................................................................
!
!
!
!
    real(kind=8) :: sig11(nbnase), sig22(nbnase), sig12(nbnase)
    real(kind=8) :: sigv11(nbnase), sigv22(nbnase), sigv12(nbnase)
    integer(kind=8) :: j, k, insui, inv, invsui, inx, inxv, insux
!
! ----------------------------------------------------------------------
!
    if (in .eq. nno) then
        insui = 1
    else
        insui = in+1
    end if
!
! --- ON ETABLIT LA CORRESPONDANCE ENTRE LA NUMEROTATION DES NOEUDS
! --- LOCALE À LA MAILLE COURANTE ET LA NUMEROTATION DES NOEUDS LOCALE
! --- À LA MAILLE VOISINE. POUR CELA ON UTILISE LA NUMEROTATION XFEM
!
    inx = zi(jcnset-1+(ndim+1)*(cpt-1)+in)
    insux = zi(jcnset-1+(ndim+1)*(cpt-1)+insui)
!
    do j = 1, nno
        inxv = zi(jcnset-1+(ndim+1)*(ivois-1)+j)
        if (inxv .eq. inx) then
            inv = j
        end if
        if (inxv .eq. insux) then
            invsui = j
        end if
    end do
!
! --- LE CAS QUADRATIQUE N'EST PAS PRÉVU CAR LES ELEMENTS SOUS DECOUPE
! --- SONT TOUJOURS LINEAIRE !
!
! --- RECUPERATION DES CONTRAINTES AUX NOEUDS DE L'ARRÊTE COURANTE
! --- POUR LE SOUS-ELEMENT COURANT
!
    sig11(1) = zr(jsigse-1+nbcmp*nno*(cpt-1)+nbcmp*(in-1)+1)
    sig11(2) = zr(jsigse-1+nbcmp*nno*(cpt-1)+nbcmp*(insui-1)+1)
!
    sig22(1) = zr(jsigse-1+nbcmp*nno*(cpt-1)+nbcmp*(in-1)+2)
    sig22(2) = zr(jsigse-1+nbcmp*nno*(cpt-1)+nbcmp*(insui-1)+2)
!
    sig12(1) = zr(jsigse-1+nbcmp*nno*(cpt-1)+nbcmp*(in-1)+4)
    sig12(2) = zr(jsigse-1+nbcmp*nno*(cpt-1)+nbcmp*(insui-1)+4)
!
! --- RECUPERATION DES CONTRAINTES AUX NOEUDS DE L'ARRÊTE COURANTE
! --- POUR LE SOUS-ELEMENT VOISIN PAR RAPPORT À CETTE ARRÊTE
!
    sigv11(1) = zr(jsigse-1+nbcmp*nno*(ivois-1)+nbcmp*(inv-1)+1)
    sigv11(2) = zr(jsigse-1+nbcmp*nno*(ivois-1)+nbcmp*(invsui-1)+1)
!
    sigv22(1) = zr(jsigse-1+nbcmp*nno*(ivois-1)+nbcmp*(inv-1)+2)
    sigv22(2) = zr(jsigse-1+nbcmp*nno*(ivois-1)+nbcmp*(invsui-1)+2)
!
    sigv12(1) = zr(jsigse-1+nbcmp*nno*(ivois-1)+nbcmp*(inv-1)+4)
    sigv12(2) = zr(jsigse-1+nbcmp*nno*(ivois-1)+nbcmp*(invsui-1)+4)
!
! --- CALCUL DU SAUT DE CONTRAINTES
!
    do k = 1, nbnase
        dsg11(k) = sig11(k)-sigv11(k)
        dsg22(k) = sig22(k)-sigv22(k)
        dsg12(k) = sig12(k)-sigv12(k)
    end do
!
end subroutine
