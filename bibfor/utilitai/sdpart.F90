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
subroutine sdpart(nbsd, nbsdp0, sdloc)
!    - FONCTION REALISEE : REPARTITION DE SOUS-DOMAINES PAR PROCESSEUR
!
! ARGUMENTS D'APPELS
! IN  NBSD   : (I) : NOMBRE DE SOUS-DOMAINES A REPARTIR
! IN  NBSDP0 : (I) : NOMBRE DE SOUS-DOMAINES A DONNER AU PROCESSEUR 0
! OUT SDLOC  : (I) : TABLEAU DE TAILLE NBSD
!                    SDLOC (I) = 1 SI LE SOUS-DOMAINE I EST TRAITE
!                                  LOCALEMENT
!----------------------------------------------------------------------
! person_in_charge: thomas.de-soza at edf.fr
! CORPS DU PROGRAMME
    implicit none
! aslint: disable=W1306
#include "asterf_types.h"
#include "asterfort/asmpi_info.h"
! DECLARATION PARAMETRES D'APPELS
    integer(kind=8) :: nbsd, sdloc(nbsd), nbsdp0
!
! DECLARATION VARIABLES LOCALES
    integer(kind=8) :: nbproc, rang, i
    integer(kind=8) :: nbsdpp, sdrest, npdeb, nsddeb, nsdfin
    integer(kind=8) :: iproc, iproc1, decal
    mpi_int :: mrank, msize
!----------------------------------------------------------------------
!
! --- INITIALISATIONS
    call asmpi_info(rank=mrank, size=msize)
    rang = to_aster_int(mrank)
    nbproc = to_aster_int(msize)
!
! --- EN SEQUENTIEL ON GAGNE DU TEMPS
    if (nbproc .eq. 1) then
        do i = 1, nbsd
            sdloc(i) = 1
        end do
        goto 999
    end if
!
    do i = 1, nbsd
        sdloc(i) = 0
    end do
!
! --- PAS DE TRAITEMENT PARTICULIER DU PROC. 0
    if (nbsdp0 .eq. 0) then
        nbsdpp = nbsd/nbproc
        sdrest = nbsd-(nbproc*nbsdpp)
        npdeb = 0
    else
!
! --- DELESTAGE DU PROC. 0
        if (rang .eq. 0) then
            do i = 1, nbsdp0
                sdloc(i) = 1
            end do
        end if
!
! ----- RESTE REPARTI ENTRE LES PROC. RESTANTS
        nbsdpp = (nbsd-nbsdp0)/(nbproc-1)
        sdrest = (nbsd-nbsdp0)-((nbproc-1)*nbsdpp)
        npdeb = 1
    end if
!
    do iproc = npdeb, nbproc-1
        if (iproc .eq. rang) then
! --------- INDICE RELATIF DU PROCESSEUR A EQUILIBRER
            iproc1 = iproc-npdeb
! --------- BORNES DES SDS A LUI ATTRIBUER
            nsddeb = 1+nbsdp0+iproc1*nbsdpp
            nsdfin = nbsdp0+(iproc1+1)*nbsdpp
! --------- REPARTITION DES SD RESTANTS (AU PLUS NBPROC-1 OU -2)
! --------- PARMI LES PROC. DE NUMERO (1 OU 2) A NBPROC-1
            if (iproc1 .gt. sdrest) then
! ----------- TOUS LES SD RESTANTS ONT ETE TRAITES
                decal = sdrest
            else
                if (iproc1 .eq. 0) then
! ------------- LE PROC. 0 OU 1 NE RECOIVENT PAS DE SD. SUPPLEMENTAIRES
                    decal = 0
                else
! ------------- LE PROC. IPROC1 RECOIT 1 SD SUPPLEMENTAIRE
                    decal = iproc1-1
                    nsdfin = nsdfin+1
                end if
            end if
! --------- ATTRIBUTION DES SD AUX PROC.
            do i = nsddeb, nsdfin
                sdloc(decal+i) = 1
            end do
        end if
    end do
!
999 continue
!
end subroutine
