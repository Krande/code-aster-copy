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

subroutine mag152(n9, n10, nomres, nugene, modmec, &
                  modgen, nbloc, indice)
    implicit none
! AUTEUR : G.ROUSSEAU
! CREATION DE LA MATRICE ASSEMBLEE GENERALISEE :
!      - OBJET    .UALF
!      - STOCKAGE .SMOS
! ET REMPLISSAGE DE SES OBJETS AUTRES QUE LE .UALF
!---------------------------------------------------------------------
#include "jeveux.h"
#include "asterfort/getvid.h"
#include "asterfort/jecrec.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: indice
    integer(kind=8) :: jrefa, i, iaconl, iadesc, iblo
    integer(kind=8) :: somme
    integer(kind=8) :: jsmde, n1bloc, n2bloc
    integer(kind=8) :: nbid, nbloc, ntbloc, nueq, nhmax
    integer(kind=8) :: n9, n10, hc
    character(len=8) :: nomres, modmec, nummod
    character(len=8) :: modgen
    character(len=14) :: num14, nugene
    character(len=19) :: nomsto
    integer(kind=4), pointer :: smhc(:) => null()
    integer(kind=8), pointer :: smdi(:) => null()
    character(len=24), pointer :: refn(:) => null()
! -----------------------------------------------------------------
!
!        CAS NUME_DDL_GENE  PRESENT
!
    call jemarq()
!
    call wkvect(nomres//'           .REFA', 'G V K24', 20, jrefa)
    zk24(jrefa-1+11) = 'MPI_COMPLET'
    nomsto = nugene//'.SMOS'
!
    if ((n9 .gt. 0)) then
        call jeveuo(nomsto//'.SMDE', 'L', jsmde)
        nueq = zi(jsmde-1+1)
        ntbloc = zi(jsmde-1+2)
        nbloc = zi(jsmde-1+3)
!       nhmax = zi(jscde-1+4)
        call jeveuo(nomsto//'.SMDI', 'L', vi=smdi)
        nhmax = 0
        do i = 1, nueq
            hc = smdi(i)
            if (i .gt. 1) hc = hc-smdi(i-1)
            nhmax = max(nhmax, hc)
        end do
!
!
! TEST SUR LE MODE DE STOCKAGE : SI ON N EST PAS EN STOCKAGE
! LIGNE DE CIEL PLEIN ON PLANTE
!
        if (nueq .ne. nhmax) then
            call utmess('A', 'ALGORITH5_16')
        end if
!
        if ((nueq*(nueq+1)/2) .gt. (nbloc*ntbloc)) then
            call utmess('F', 'ALGORITH5_17')
        end if
!
! CALCUL DU NOMBRE DE TERME PAR BLOC ET TOTAL
!
        call jeveuo(nomsto//'.SMHC', 'L', vi4=smhc)
!
        somme = 0
!
        do iblo = 1, nbloc
!
!----------------------------------------------------------------
!
!         BOUCLE SUR LES COLONNES DE LA MATRICE ASSEMBLEE
!
            n1bloc = 1
            n2bloc = nueq
!
!
            do i = n1bloc, n2bloc
                hc = smdi(i)
                if (i .gt. 1) hc = hc-smdi(i-1)
                somme = somme+hc
            end do
        end do
!
        if ((nueq*(nueq+1)/2) .ne. somme) then
            call utmess('F', 'ALGORITH5_18')
        end if
!
!
!
        call jecrec(nomres//'           .UALF', 'G V R', 'NU', 'DISPERSE', 'CONSTANT', &
                    nbloc)
        call jeecra(nomres//'           .UALF', 'LONMAX', ntbloc)
!
!
        call wkvect(nomres//'           .CONL', 'G V R', nueq, iaconl)
!
!       CAS DU CHAM_NO
!
    else
!
        call jeveuo(nomsto//'.SMDE', 'L', jsmde)
        nueq = zi(jsmde-1+1)
        nbloc = 1
        ntbloc = nueq*(nueq+1)/2
!
        call jecrec(nomres//'           .UALF', 'G V R', 'NU', 'DISPERSE', 'CONSTANT', &
                    nbloc)
        call jeecra(nomres//'           .UALF', 'LONMAX', ntbloc)
        call wkvect(nomres//'           .CONL', 'G V R', nueq, iaconl)
!
    end if
!
! ----------- CREATION ET REMPLISSAGE DU .DESC ---------------
    call wkvect(nomres//'           .DESC', 'G V I', 3, iadesc)
    zi(iadesc) = 2
    zi(iadesc+1) = ntbloc
    zi(iadesc+2) = 2
!
! ----------- REMPLISSAGE DU .REFA
!---------------------ET DU .CONL ---------------------------
!
!
    if (n10 .gt. 0) then
        zk24(jrefa-1+1) = ' '
!
    else if (indice .eq. 1) then
        call getvid(' ', 'NUME_DDL_GENE', scal=nummod, nbret=nbid)
        num14 = nummod
        call jeveuo(num14//'.NUME.REFN', 'L', vk24=refn)
        zk24(jrefa-1+1) = refn(1)
!
    else
        zk24(jrefa-1+1) = modmec
    end if
!
    zk24(jrefa-1+2) = nugene
    zk24(jrefa-1+9) = 'MS'
    zk24(jrefa-1+10) = 'GENE'
!
    do i = 1, nueq
        zr(iaconl+i-1) = 1.0d0
    end do
!
    call jedema()
end subroutine
