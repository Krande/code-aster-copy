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
subroutine acyel1(nmcolz, nomobz, nobl, nobc, okpart, &
                  lilig, nblig, licol, nbcol, cmat, &
                  ndim, ideb, jdeb, x)
    implicit none
!
!***********************************************************************
!    P. RICHARD     DATE 10/04/91
!-----------------------------------------------------------------------
!  BUT:  ASSEMBLER SI ELLE EXISTE LA SOUS-MATRICE  CORRESPONDANT
!  A UN NOM OBJET DE COLLECTION DANS UNE MATRICE COMPLEXE AVEC
!   UN ASSEMBLAGE EN UN TEMPS (ADAPTE AU CYCLIQUE)
!
!  SI OKPART VRAI ON ASSEMBLE QUE LA LISTE DE LIGNE ET DE COLONNES
!                 DONNEES
!
!  SI OKPART FAUX ON ASSEMBLE TOUT
!
!-----------------------------------------------------------------------
!
! NMCOLZ   /I/: NOM K24 DE LA COLLECTION
! NOMOBZ   /I/: NOM K8 DE L'OBJET DE COLLECTION
! NOBL     /I/: NOMBRE DE LIGNE DE LA MATRICE ELEMENTAIRE
! NOBC     /I/: NOMBRE DE COLONNES DE LA MATRICE ELEMENTAIRE
! OKPART   /I/: INDICATEUR SI ASSEMBLAGE PARTIEL
! LILIG    /I/: LISTE DES INDICE DE LIGNE A ASSEMBLER (SI OKPART=TRUE)
! NBLIG    /I/: NOMBRE DE LIGNES DE LA LISTE
! LICOL    /I/: LISTE INDICES DE COLONNES A ASSEMBLER (SI OKPART=TRUE)
! NBLIG    /I/: NOMBRE DE COLONNES DE LA LISTE
! CMAT     /M/: MATRICE RECEPTRICE COMPLEXE
! NDIM     /I/: DIMENSION DE LA MATRICE RECEPTRICE CARREE
! IDEB     /I/: INDICE DE PREMIERE LIGNE RECEPTRICE
! JDEB     /I/: INDICE DE PREMIERE COLONNE RECEPTRICE
! X        /I/: COEFFICIENT ASSEMBLAGE
!
!
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/ampcpr.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iad, ibid, ideb, iret, j, jdeb
    integer(kind=8) :: llob, nbcol, nblig, ndim, nobc, nobl
    real(kind=8) :: x
    character(len=8) :: nomob
    character(len=24) :: nomcol
    complex(kind=8) :: cmat(ndim, ndim)
    character(len=*) :: nmcolz, nomobz
    integer(kind=8) :: lilig(nblig), licol(nbcol)
    aster_logical :: okpart
!-----------------------------------------------------------------------
    call jemarq()
    nomob = nomobz
    nomcol = nmcolz
!
    call jenonu(jexnom(nomcol(1:15)//'.REPE.MAT', nomob), iret)
    if (iret .eq. 0) goto 999
!
    call jenonu(jexnom(nomcol(1:15)//'.REPE.MAT', nomob), ibid)
    call jeveuo(jexnum(nomcol, ibid), 'L', llob)
!
    if (okpart) then
!
! SI ASSEMBLAGE PARTIEL ON TRAITE LIGNE PAR LIGNE
!        ET COLONNE PAR COLONNE
!
        do j = 1, nbcol
            do i = 1, nblig
                iad = llob+(licol(j)-1)*nobl+lilig(i)-1
                call ampcpr(cmat, ndim, ndim, zr(iad), 1, &
                            1, ideb+i-1, jdeb+j-1, x, 1, &
                            1)
            end do
        end do
!
    else
!
!  SI ASSEMBLAGE COMPLET ON TRAITE TOUT D'UN COUP
!
        call ampcpr(cmat, ndim, ndim, zr(llob), nobl, &
                    nobc, ideb, jdeb, x, 1, &
                    1)
!
    end if
!
!
999 continue
    call jedema()
end subroutine
