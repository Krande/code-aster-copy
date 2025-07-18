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

subroutine te0144(option, nomte)
!
!
! --------------------------------------------------------------------------------------------------
!     CALCUL
!       - DU VECTEUR ELEMENTAIRE EFFORT GENERALISE,
!       - DU VECTEUR ELEMENTAIRE CONTRAINTE
!     POUR LES ELEMENTS DE POUTRE D'EULER ET DE TIMOSHENKO.
!
! --------------------------------------------------------------------------------------------------
!
! IN  OPTION : K16 : NOM DE L'OPTION A CALCULER
!        'SIEF_ELGA'
! IN  NOMTE  : K16 : NOM DU TYPE ELEMENT
!        'MECA_POU_D_E' : POUTRE DROITE D'EULER       (SECTION VARIABLE)
!        'MECA_POU_D_T' : POUTRE DROITE DE TIMOSHENKO (SECTION VARIABLE)
!
! --------------------------------------------------------------------------------------------------
!
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lonele.h"
#include "asterfort/matrot.h"
#include "asterfort/moytem.h"
#include "asterfort/pmavec.h"
#include "asterfort/porigi.h"
#include "asterfort/poutre_modloc.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utpvgl.h"
#include "asterfort/vecma.h"
#include "asterfort/verifm.h"
!
    character(len=*) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
    integer(kind=8) :: nbres, npg, nno, nc, nnoc, ncc, jeffo, lmater, iret
    integer(kind=8) :: iret1, itype, lorien, jdepl, i, j
!
    real(kind=8) :: ul(12), ug(12), pgl(3, 3), klc(12, 12), klv(78)
    real(kind=8) :: fl(12), epsith
    real(kind=8) :: temp, tvar
    real(kind=8) :: e, xnu, xl
!
! --------------------------------------------------------------------------------------------------
    parameter(nbres=2)
    real(kind=8) :: valres(nbres)
    integer(kind=8) :: codres(nbres)
    character(len=16) :: nomres(nbres)
    data nomres/'E', 'NU'/
! --------------------------------------------------------------------------------------------------
!
    nno = 2
    nc = 6
    nnoc = 1
    ncc = 6
!
    if (option .eq. 'SIEF_ELGA') then
        call jevech('PCONTRR', 'E', jeffo)
    else
!       option non programmee
        ASSERT(.false.)
    end if
!
!   point de gauss de l'element
    call elrefe_info(fami='RIGI', npg=npg)
    ASSERT((npg .eq. 2) .or. (npg .eq. 3))
!   recuperation des caracteristiques materiaux
    call jevech('PMATERC', 'L', lmater)
!
!   recuperation de la temperature
    call verifm('RIGI', npg, 1, '+', zi(lmater), epsith, iret)
    call moytem('RIGI', npg, 1, '+', temp, iret1)
!
    call rcvalb('RIGI', 1, 1, '+', zi(lmater), ' ', 'ELAS', 1, 'TEMP', [temp], &
                2, nomres, valres, codres, 1)
    e = valres(1)
    xnu = valres(2)
!
!   calcul de la matrice de rigidite locale
    call porigi(nomte, e, xnu, -1.d0, klv)
!
!   matrice rigidite ligne > matrice rigidite carre
    call vecma(klv, 78, klc, 12)
!
!   recuperation des caracteristiques generales des sections
    call poutre_modloc('CAGNPO', ['TVAR'], 1, valeur=tvar)
    itype = nint(tvar)
!
!   recuperation des coordonnees des noeuds
    xl = lonele()
!
!   matrice de rotation pgl
    call jevech('PCAORIE', 'L', lorien)
    call matrot(zr(lorien), pgl)
!
    call jevech('PDEPLAR', 'L', jdepl)
    do i = 1, 12
        ug(i) = zr(jdepl+i-1)
    end do
!
!   vecteur deplacement local  UL = PGL * UG
    call utpvgl(nno, nc, pgl, ug, ul)
!
!   vecteur effort local  FL = KLC * UL
    call pmavec('ZERO', 12, klc, ul, fl)
!
!   tenir compte des efforts dus a la dilatation
    if (epsith .ne. 0.0d0) then
        do i = 1, 12
            ug(i) = 0.0d0
        end do
        ug(1) = -epsith*xl
        ug(7) = -ug(1)
!       calcul des forces induites
        do i = 1, 6
            do j = 1, 6
                fl(i) = fl(i)-klc(i, j)*ug(j)
                fl(i+6) = fl(i+6)-klc(i+6, j+6)*ug(j+6)
            end do
        end do
    end if
!   archivage
    if (npg .eq. 2) then
        ASSERT(.false.)
!        do i = 1, 6
!            zr(jeffo+i-1) = -fl(i)
!            zr(jeffo+i+6-1) = fl(i+6)
!        enddo
    else
!       dans le cas 3 points de gauss il ne faut pas mettre 0 sur un des points de gauss
!        cf doc aster pour la numerotation ou elrega
!              noeud        n1       n2       n3
!              position   (-0.7777 , 0.0000 , 0.7777)
!           c'est lineaire : (n2) = ( (n1) + (n3) )/2
        do i = 1, 6
            zr(jeffo+i-1) = -fl(i)
            zr(jeffo+i+6-1) = (-fl(i)+fl(i+6))*0.5d0
            zr(jeffo+i+12-1) = fl(i+6)
        end do
    end if
end subroutine
