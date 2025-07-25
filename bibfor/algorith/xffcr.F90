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
subroutine xffcr(nfon, jfono, jbaso, jtailo, jindpt, &
                 typfon, jfon, jnofaf, jbas, jtail)
!
!
    implicit none
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/padist.h"
    integer(kind=8) :: nfon, jfono, jbaso, jtailo, jindpt, jfon, jbas, jtail
    integer(kind=8) :: jnofaf
    character(len=19) :: typfon
!
! ----------------------------------------------------------------------
!
! ROUTINE XFEM
!
!              ORDONNANCEMENT DES VECTEURS BASEFOND, FONDFISS ET
!              FOND.TAILLE_R
!
!              REMPLISSAGE ET ORDONNANCEMENT DU VECTEUR NOFACPTFON
!
! ----------------------------------------------------------------------
!
!
! IN  NFON  :  NOMBRE DE POINTS AU FOND DE FISSURE
!     JFONO :  ADRESSE DES POINTS DU FOND DE FISSURE DÉSORDONNÉS
!     JBASO :  ADRESSE DES DIRECTIONS DE PROPAGATION DÉSORDONNÉES
!     JTAILO:  ADRESSE DES TAILLES MAXIMALES DE MAILLES DÉSORDONNÉES
!     JINDPT:  ADRESSE DES INDICES DES POINTS ORDONNES
!     TYPFON:  TYPE DU FOND DE FISSURE (OUVERT OU FERME)
!
! OUT JFON  :  ADRESSE DES POINTS DU FOND DE FISSURE ORDONNÉS
!     JNOFAF:  ADRESSE DES NUMERO DES NOEUDS DES FACES DES ELEMENTS
!              PARENTS QUI CONTIENNENT LES POINTS DU FOND DE FISSURE
!              ORDONNES (VECTEUR NOFACPTFON)
!     JBAS  :  ADRESSE DES DIRECTIONS DE PROPAGATION ORDONNÉES
!     JTAIL :  ADRESSE DES TAILLES MAXIMALES DE MAILLES ORDONNÉES
!
!
    integer(kind=8) :: indipt, ipt, k
    real(kind=8) :: m(3), p(3)
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
    do ipt = 1, nfon
!
        indipt = zi(jindpt-1+ipt)
!
        do k = 1, 3
!
            zr(jfon-1+4*(ipt-1)+k) = zr(jfono-1+11*(indipt-1)+k)
            zi(jnofaf-1+4*(ipt-1)+k) = int(zr(jfono-1+11*(indipt-1)+4+k))
            zr(jbas-1+6*(ipt-1)+k) = zr(jbaso-1+6*(indipt-1)+k)
            zr(jbas-1+6*(ipt-1)+k+3) = zr(jbaso-1+6*(indipt-1)+3+k)
!
        end do
!
        zr(jfon-1+4*(ipt-1)+4) = zr(jfono-1+11*(indipt-1)+4)
        zi(jnofaf-1+4*(ipt-1)+4) = int(zr(jfono-1+11*(indipt-1)+8))
        zr(jtail-1+ipt) = zr(jtailo-1+indipt)
!
    end do
!
!     CAS D'UN FOND FERME: PREMIER POINT DU FOND = DERNIER POINT
    if (typfon .eq. 'FERME') then
!
        nfon = nfon+1
!
        do k = 1, 3
!
            zr(jfon-1+4*(nfon-1)+k) = zr(jfon-1+4*(1-1)+k)
            zi(jnofaf-1+4*(nfon-1)+k) = zi(jnofaf-1+4*(1-1)+k)
            zr(jbas-1+6*(nfon-1)+k) = zr(jbas-1+6*(1-1)+k)
            zr(jbas-1+6*(nfon-1)+k+3) = zr(jbas-1+6*(1-1)+3+k)
!
            p(k) = zr(jfon-1+4*(nfon-1)+k)
            m(k) = zr(jfon-1+4*(nfon-2)+k)
!
        end do
!
        zr(jfon-1+4*(nfon-1)+4) = zr(jfon-1+4*(nfon-2)+4)+padist(3, m, &
                                                                 p)
        zi(jnofaf-1+4*(nfon-1)+4) = zi(jnofaf-1+4*(1-1)+4)
        zr(jtail-1+nfon) = zr(jtail-1+1)
!
    end if
!
    call jedema()
end subroutine
