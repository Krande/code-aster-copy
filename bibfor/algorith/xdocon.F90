! --------------------------------------------------------------------
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
subroutine xdocon(algocr, algofr, cface, contac, coefcp, &
                  coeffp, coefcr, coeffr, elc, fpg, &
                  ifiss, ivff, jcface, jdonco, jlonch, &
                  mu, nspfis, ncompd, ndim, nface, &
                  ninter, nnof, nomte, npgf, nptf, &
                  rela)
! aslint: disable=W1504
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/teattr.h"
#include "asterfort/xminte.h"
! RECUPERATION DIVERSES INFOS LIEES AU CONTACT XFEM DANS UN TE
!
! OUT ALGOCR : ALGO DE CONTACT (1:LAG, 2:PENA, 3:COHESIF)
! OUT ALGOFR : ALGO DE FROTTEMENT (1: LAGRANGIEN, 2: PENALISATION)
! OUT CFACE  : TABLEAU CONNECTIVITE FACETTES DE CONTACT
! IN  CONTAC : TYPE DE CONTACT (P1P1, P1P2)
! OUT COEFCP : COEFF PENALISATION CONTACT
! OUT COEFFP : COEFF PENALISATION FROTTEMENT
! OUT COEFCR : COEFF AUGM CONTACT
! OUT COEFFR : COEFF AUGM FROTTEMENT
! OUT COHES  : VECTEUR VARIABLES INTERNES LOI COHESIVE
! OUT IAINT  : ADRESSE TOPOFAC.AI FISSURE COURANTE
! IN  IFISS  : NUMERO FISSURE
! OUT INDCO  : TABLEAU DES STATUTS FISSURE COURANTE
! OUT IVFF   : ADRESSE TABLEAU DES FONCTIONS DE FORME
! IN  JAINT  : ADRESSE TOPOFAC.AI INITIAL
! IN  JBASEC : ADRESSE TOPOFAC.BA INITIAL
! IN  JCFACE : ADRESSE TOPOFAC.CF INITIAL
! IN  JDONCO : ADRESSE CHAM_ELEM DONNEES DU CONTACT
! IN  JINDCO : ADRESSE STOCKAGE CHAMP STATUT DE CONTACT
! IN  JLONCH : ADRESSE STOCKAGE CHAMP TOPOFAC.LO
! IN  JPTINT : ADRESSE STOCKAGE CHAMP TOPOFAC.PI
! IN  JSEUIL : ADRESSE STOCKAGE CHAMP SEUIL FROTTEMENT
! OUT MU     : COEFFICIENT DE COULOMB
! IN  NCOMP[X] : NOMBRE DE COMPOSANTES CHAMP LOCAL DS LE CATALOGUE
! IN NDIM    : DIMENSION MODELE
! OUT NFACE  : NOMBRE DE FACETTES POUR L ELEMENT ET LA FISSURE
! OUT NINTER : NOMBRE DE POINTS D INTERSECTION
! OUT NNOF   : NOMBRE DE NOEUDS TOTAL
!               D UN EL DE CONTACT (UNE FACETTE)
! IN  NOMTE  : TYPE ELEMENT
! OUT NPGF   : NOMBRE DE POINTS DE GAUSS PAR FACETTE
! OUT SEUIL  : TABLEAU DES SEUILS DE FROTTEMENT
    integer :: algocr, algofr, cface(30, 6), contac
    integer :: i, ibid, idfdef, ifiss
    integer :: ipoidf, ivff, j, jcface
    integer :: jdonco, jlonch
    integer :: nspfis, ncompd, ndim
    integer :: nface, nint, ninteg, ninter, nnof
    integer :: npgf, nptf
    real(kind=8) :: coefcp, coeffp, coefcr, coeffr, mu, rela
    character(len=8) :: elc, fpg
    character(len=16) :: enr, nomte
!
!
! --- RECUPERATION TOPOLOGIE FACETTES
! --- DECALAGE D INDICE POUR TOPOFAC (MULTIFISSURATION)
!
!   INITIALISATION
    rela = 0.d0
!
    ninter = zi(jlonch+3*(ifiss-1)-1+1)
    if (ninter .eq. 0) goto 99
!
! --- SCHEMA D'INTEGRATION NUMERIQUE ET ELEMENT
!     DISCUSSION VOIR BOOK IV 18/10/2004 ET BOOK VI 06/07/2005
    ninteg = nint(zr(jdonco-1+(ifiss-1)*ncompd+4))
    call xminte(ndim, ninteg, fpg)
!
    if (ndim .eq. 3) then
        if (contac .le. 2) then
            elc = 'TR3'
        else
            elc = 'TR3'
        end if
    else if (ndim .eq. 2) then
        if (contac .le. 2) then
            elc = 'SE2'
        else
            elc = 'SE3'
        end if
    end if
    call elrefe_info(elrefe=elc, fami=fpg, nno=nnof, npg=npgf, jpoids=ipoidf, &
                     jvf=ivff, jdfde=idfdef)
!
! --- RECUPERATIONS DONNEES ET COEFFS DU CONTACT
!
    algocr = nint(zr(jdonco-1+(ifiss-1)*ncompd+6))
    algofr = nint(zr(jdonco-1+(ifiss-1)*ncompd+7))
    coefcr = zr(jdonco-1+(ifiss-1)*ncompd+1)
    mu = zr(jdonco-1+(ifiss-1)*ncompd+2)
    coeffr = zr(jdonco-1+(ifiss-1)*ncompd+3)
    coefcp = zr(jdonco-1+(ifiss-1)*ncompd+8)
    coeffp = zr(jdonco-1+(ifiss-1)*ncompd+9)
!     VERIFICATIONS
    rela = 0.d0
    if (algocr .eq. 1) then
        ASSERT(coefcr .ne. 0.d0 .and. coefcp .eq. 0.d0)
    else if (algocr .eq. 2) then
        ASSERT(coefcr .eq. 0.d0 .and. coefcp .ne. 0.d0)
    else if (algocr .eq. 3) then
        call teattr('S', 'XFEM', enr, ibid, typel=nomte)
        ASSERT(enr(1:3) .eq. 'XHC')
        rela = zr(jdonco-1+(ifiss-1)*ncompd+10)
        ASSERT(rela .gt. 0.d0)
        coefcr = 1.d0
    else
        ASSERT(.false.)
    end if
!
    if (algofr .eq. 1) then
        ASSERT(coeffr .ne. 0.d0 .and. coeffp .eq. 0.d0)
    else if (algofr .eq. 2) then
        ASSERT(coeffr .eq. 0.d0 .and. coeffp .ne. 0.d0)
    else
        ASSERT(algofr .eq. 0)
    end if
!
! --- DECALAGE INDICE POUR LA MULTIFISSURATION
!
    nface = zi(jlonch+3*(ifiss-1)-1+2)
    nptf = zi(jlonch+3*(ifiss-1)-1+3)
    do i = 1, nface
        do j = 1, nptf
            cface(i, j) = zi(jcface-1+nptf*(i-1)+j)
        end do
    end do
!
    nspfis = npgf*nface
!
99  continue
!
end subroutine
