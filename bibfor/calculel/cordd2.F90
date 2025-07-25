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
subroutine cordd2(jprn1, jprn2, ili, ecodl, nec, &
                  ncmp, n, nddloc, pos)
! aslint: disable=
    implicit none
!
!     BUT:
!     ----
!     ROUTINE PARALLELE A CORDDL POUR LES SOUS-STRUCTURES.
!
!     IN
!     --
!     JPRN1,JPRN2 : ADRESSES DE PRNO ( OBJET ET POINTEUR DE LONGUEUR)
!     NEC  : NBEC(GD) (SUPPOSE = 1!)
!     ILI  : NUMERO DU LIGREL (ICI TOUJOURS 1).
!     N    : NUMERO GLOBAL DU NOEUD
!     ECODL(*) : ENTIER CODE DU NUMERO LOCAL DU NOEUD
!
!     OUT
!     ---
!     NDDLOC : NBRE DE DDL SUPPORTES PAR CE NOEUD SUR L ELEMENT
!     POS    : TABLEAU DE CORRESPONDANCE AVEC LES DDL SUR LE NOEUD
!              EN TANT QUE NOEUD GLOBAL
! ----------------------------------------------------------------------
#include "jeveux.h"
#include "asterfort/assert.h"
    integer(kind=8) :: nbecmx, ncmp
    parameter(nbecmx=10)
    integer(kind=8) :: ifin(nbecmx)
    integer(kind=8) :: pos(1)
    integer(kind=8) :: ecodg, ecodl(*)
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
    integer(kind=8) :: ili, jprn1, jprn2, nec, iec, in
    integer(kind=8) :: nddloc, n, i, iecdg, iecdl
!
!     FONCTION D ACCES A PRNO
#define prno(ili,nunoel,l) zi(jprn1-1+zi(jprn2+ili-1)+ (nunoel-1)* (nec+2)+l-1)
!
! - DEB ----------------------------------------------------------------
!
! --- IFIN DONNE POUR CHAQUE ENTIER CODE LE NOMBRE MAX DE DDLS
! --- QUE L'ON PEUT TROUVER SUR CET ENTIER :
!     ------------------------------------
    ASSERT(nec .gt. 0 .and. nec .le. nbecmx)
    ASSERT(ncmp .gt. 0 .and. ncmp .le. 30*nbecmx)
    do iec = 1, nec-1
        ifin(iec) = 30
    end do
    ifin(nec) = ncmp-30*(nec-1)
!
    in = 0
    nddloc = 0
!
! --- AJOUT DE LA BOUCLE 20 SUR LE NOMBRE D'ENTIERS CODES
! --- PAR RAPPORT A LA VERSION NORMALE DE CORDD2 , LES INSTRUCTIONS
! --- NE CHANGENT PAS, EXCEPTEE LA DEFINITION DE ECODL ET ECODG
! --- OU INTERVIENT L'INDICE D'ENTIER CODE :
!     ------------------------------------
    do iec = 1, nec
!
!      ECODG = PRNO(ILI,N,3)
        ecodg = prno(ili, n, 2+iec)
!
        do i = 1, ifin(iec)
            ecodg = ecodg/2
            ecodl(iec) = ecodl(iec)/2
            iecdg = iand(1, ecodg)
            if (iecdg .gt. 0) then
                in = in+1
                iecdl = iand(1, ecodl(iec))
                if (iecdl .gt. 0) then
                    nddloc = nddloc+1
                    pos(nddloc) = in
                end if
            end if
        end do
!
    end do
!
end subroutine
