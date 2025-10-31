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
subroutine corddl(admodl, lcmodl, idprn1, idprn2, ili, &
                  mode, nec, ncmp, n, k, &
                  nddloc, pos)
! aslint: disable=
    implicit none
!
!
!     IN
!     --
!     IDPRN1 IDPRN2 : ADRESSES DE PRNO ( OBJET ET POINTEUR DE LONGUEUR)
!     ILI : NUMERO DE LIGREL
!     MODE : MODE
!     NEC  : NBEC(GD)
!
!     N    : NUMERO GLOBAL DU NOEUD
!     K    : NUMERO LOCAL DU NOEUD ( DANS ELEMENT )
!     OUT
!     ---
!
!     NDDLOC : NBRE DE DDL SUPPORTES PAR CE NOEUD SUR L ELEMENT
!     POS    : TABLEAU DE CORRESPONDANCE AVEC LES DDL SUR LE NOEUD
!          EN TANT QUE NOEUD GLOBAL
! ----------------------------------------------------------------------
#include "jeveux.h"
#include "asterfort/entcod.h"
!-----------------------------------------------------------------------
    integer(kind=8) :: i, idprn1, idprn2, iec, iecdg, iecdl, ili
    integer(kind=8) :: in, k, mode, n, nbecmx, ncmp
    integer(kind=8) :: nddloc, nec
!-----------------------------------------------------------------------
    parameter(nbecmx=11)
    integer(kind=8) :: ifin(nbecmx)
    integer(kind=8) :: admodl, lcmodl, pos(1)
    integer(kind=8) :: ecodg, ecodl
!
!     ECODG : ENTIER CODE DE N
!     ECODL : ENTIER CODE DE K
!
!
!
!
!
!     FONCTION D ACCES A PRNO
!
#define prno(ili,nunoel,l) zi(idprn1-1+zi(idprn2+ili-1)+ (nunoel-1)* (nec+2)+l-1)
! - DEB ----------------------------------------------------------------
!
! --- IFIN DONNE POUR CHAQUE ENTIER CODE LE NOMBRE MAX DE DDLS
! --- QUE L'ON PEUT TROUVER SUR CET ENTIER :
!     ------------------------------------
    do iec = 1, nec-1
        ifin(iec) = 30
    end do
    ifin(nec) = ncmp-30*(nec-1)
!
    in = 0
    nddloc = 0
!
!
    do iec = 1, nec
        ecodg = prno(ili, n, 2+iec)
        if (ecodg .eq. 0) goto 20
!
        ecodl = entcod(admodl, lcmodl, nec, mode, k, iec)
!
        do i = 1, ifin(iec)
            ecodg = ecodg/2
            if (ecodg .eq. 0) goto 20
            ecodl = ecodl/2
            iecdg = iand(1, ecodg)
            if (iecdg .gt. 0) then
                in = in+1
                iecdl = iand(1, ecodl)
                if (iecdl .gt. 0) then
                    nddloc = nddloc+1
                    pos(nddloc) = in
                end if
            end if
        end do
20      continue
    end do
!
end subroutine
