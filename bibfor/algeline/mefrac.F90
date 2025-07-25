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
subroutine mefrac(mailla, nbgrmx, nomrac, nbgrma, nomcyl)
    implicit none
!
#include "jeveux.h"
#include "asterfort/iunifi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
!
    integer(kind=8) :: nbgrmx, nbgrma
    character(len=8) :: mailla
    character(len=24) :: nomcyl(*), nomrac
!     RECHERCHE DES GROUPES DE MAILLES PRESENTANT UNE RACINE COMMUNE
!     SOUS LA FORME *RACINE*, OU *RACINE, OU RACINE*.
!     OPERATEUR APPELANT : OP0144 , FLUST3
! ----------------------------------------------------------------------
!     OPTION DE CALCUL   : CALC_FLUI_STRU , CALCUL DES PARAMETRES DE
!     COUPLAGE FLUIDE-STRUCTURE POUR UNE CONFIGURATION DE TYPE "FAISCEAU
!     DE TUBES SOUS ECOULEMENT AXIAL"
! ----------------------------------------------------------------------
! IN  : MAILLA : NOM DU MAILLAGE
! IN  : NBGRMX : NOMBRE DE GROUPES DE MAILLES DU MAILLAGE
! IN  : NOMRAC : NOM DE LA RACINE COMMUNE
! OUT : NBGRMA : NOMBRE DE GROUPES DE MAILLES AVEC LA RACINE COMMUNE
!                NOMRAC
! OUT : NOMCYL : NOMS DES GROUPES DE MAILLES AVEC LA RACINE COMMUNE
!                NOMRAC
! ----------------------------------------------------------------------
!     ------------------------------------------------------------------
    character(len=24) :: nomgri
!     ------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, icar, ifm, ipre, isuf, j
    integer(kind=8) :: ndeb, nfin, nt
!-----------------------------------------------------------------------
    call jemarq()
!
    nbgrma = 0
    ipre = 0
    isuf = 0
    icar = 0
    if (nomrac(1:1) .eq. '*') then
        ipre = 1
    end if
    do i = 2, 8
        if (nomrac(i:i) .eq. '*') then
            isuf = 1
            icar = i-1-ipre
            goto 20
        else if (nomrac(i:i) .eq. ' ') then
            icar = i-1-ipre
            goto 20
        end if
    end do
20  continue
!
    if (isuf .eq. 0 .and. ipre .eq. 0) then
        call utmess('F', 'ALGELINE_86', sk=nomrac)
    end if
!
    if (ipre .eq. 0) then
        do i = 1, nbgrmx
            call jenuno(jexnum(mailla//'.GROUPEMA', i), nomgri)
            if (nomrac(1:icar) .eq. nomgri(1:icar)) then
                nbgrma = nbgrma+1
                nomcyl(nbgrma) = nomgri
            end if
        end do
    else if (ipre .eq. 1 .and. isuf .eq. 0) then
        do i = 1, nbgrmx
            call jenuno(jexnum(mailla//'.GROUPEMA', i), nomgri)
            do j = 2, 8-icar
                if (nomrac(2:2) .eq. nomgri(j:j)) then
                    if (nomrac(2:icar+1) .eq. nomgri(j:(j+icar-1)) .and. &
                        nomgri((j+icar):(j+icar)) .eq. ' ') then
                        nbgrma = nbgrma+1
                        nomcyl(nbgrma) = nomgri
                        goto 60
                    end if
                end if
            end do
            j = 8-icar+1
            if (nomrac(2:2) .eq. nomgri(j:j)) then
                if (nomrac(2:icar+1) .eq. nomgri(j:(j+icar-1))) then
                    nbgrma = nbgrma+1
                    nomcyl(nbgrma) = nomgri
                end if
            end if
60          continue
        end do
    else if (ipre .eq. 1 .and. isuf .ne. 0) then
        do i = 1, nbgrmx
            call jenuno(jexnum(mailla//'.GROUPEMA', i), nomgri)
            do j = 1, 8-icar
                if (nomrac(2:2) .eq. nomgri(j:j)) then
                    if (nomrac(2:icar+1) .eq. nomgri(j:(j+icar-1))) then
                        nbgrma = nbgrma+1
                        nomcyl(nbgrma) = nomgri
                        goto 90
                    end if
                end if
            end do
90          continue
        end do
    end if
    if (nbgrma .eq. 0) then
        call utmess('F', 'ALGELINE_87')
    end if
!
    ifm = iunifi('MESSAGE')
    write (ifm, *) '==============================================='&
     &   , '================================='
    write (ifm, *) '           GROUPES DE MAILLES SELECTIONNES '&
     &   , 'POUR LA RACINE COMMUNE'
    write (ifm, *) '==============================================='&
     &   , '================================='
    nt = int(nbgrma/8)
    do i = 1, nt
        ndeb = nt*(i-1)+1
        nfin = nt*(i-1)+8
        write (ifm, 6001) (nomcyl(j), j=ndeb, nfin)
    end do
    if ((nt*8) .lt. nbgrma) then
        ndeb = nt*8+1
        nfin = nbgrma
        write (ifm, 6001) (nomcyl(j), j=ndeb, nfin)
    end if
    write (ifm, *) '==============================================='&
     &   , '================================='
!
6001 format(1x, 6(2x, a24))
!
!
    call jedema()
end subroutine
