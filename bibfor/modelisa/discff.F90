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
subroutine discff(nbfonc, nomfon, nbp1, nbp2, disc, &
                  vale)
    implicit none
!     DISCRETISATION DES FONCTIONS DE FORME
!     APPELANT : SPECFF
!-----------------------------------------------------------------------
! IN  : NBFONC : NOMBRE DE TABLE_FONCTIONS
! IN  : NOMFON : LISTE DES NOMS DES CONCEPTS DE TYPE TABLE_FONCTION
! IN  : NBP1   : NOMBRE DE POINTS DE DISCRETISATION DES FONCTIONS DE
!                FORME SUR LA PREMIERE DIRECTION D'ESPACE DU REPERE
! IN  : NBP2   : NOMBRE DE POINTS DE DISCRETISATION DES FONCTIONS DE
!                FORME SUR LA DEUXIEME DIRECTION D'ESPACE DU REPERE
! OUT : DISC   : VALEURS DU PARAMETRE SUR 0,2L   - DIM : NBP1+NBP2
! OUT : VALE   : TABLEAU DES VALEURS DES FONCTIONS SUR 0,2L
!                - DIM : (NBP1+NBP2,NBFONC)
!
! REMARQUE :
! ----------
! LES DEUX FONCTIONS DE FORME SONT CONCATENEES ET LES VALEURS DES
! DISCRETISATIONS PASSENT DE (0,L) POUR LA PREMIERE FONCTION ET (0,L)
! POUR LA DEUXIEME A UNE SEULE FONCTION DEFINIE SUR (0,2L): (0,L),(L,2L)
!
#include "jeveux.h"
#include "asterfort/fointr.h"
#include "asterfort/jedema.h"
#include "asterfort/jelibe.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: nbfonc, nbp1, nbp2
    character(len=8) :: nomfon(nbfonc)
    character(len=19) :: tbfonc, fonc1, fonc2
    real(kind=8) :: disc(nbp1+nbp2), vale(nbp1+nbp2, nbfonc)
!
    character(len=24) :: fvale1, fprol1, fvale2, fprol2, lstfon, tblp
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: ier1, ier2, ifo, ip, itblp, ilfon
    integer(kind=8) :: ifabs1, iford1, ifpro1, ifabs2, iford2, ifpro2
    real(kind=8) :: long
!-----------------------------------------------------------------------
    call jemarq()
!
!
!-----1.VALEURS DU PARAMETRE : ON RETIENT LA DISCRETISATION DE LA
!-----  PREMIERE FONCTION DE FORME
!
    tbfonc = nomfon(1)
    tblp = tbfonc//'.TBLP'
!
    call jeveuo(tblp, 'L', itblp)
    lstfon = zk24(itblp+2)
    call jeveuo(lstfon, 'L', ilfon)
    fonc1 = zk8(ilfon)//'           '
    fonc2 = zk8(ilfon+1)//'           '
!
    call jelibe(lstfon)
    call jelibe(tblp)
!
    fvale1 = fonc1//'.VALE'
    fvale2 = fonc2//'.VALE'
!
    call jeveuo(fvale1, 'L', ifabs1)
    call jeveuo(fvale2, 'L', ifabs2)
    do ip = 1, nbp1
        disc(ip) = zr(ifabs1+ip-1)
    end do
    long = disc(nbp1)
    do ip = 1, nbp2
        disc(nbp1+ip) = zr(ifabs2+ip-1)
    end do
!
!-----2.VALEURS DES FONCTIONS
!
!-----2.1.ON RECOPIE LES VALEURS DE LA PREMIERE FONCTION
!
    iford1 = ifabs1+nbp1
    iford2 = ifabs2+nbp2
    do ip = 1, nbp1
        vale(ip, 1) = zr(iford1+ip-1)
    end do
    do ip = 1, nbp2
        vale(nbp1+ip, 1) = zr(iford2+ip-1)
    end do
!
!-----2.2.ON INTERPOLE LES AUTRES FONCTIONS LE CAS ECHEANT
!
    if (nbfonc .gt. 1) then
!
        do ifo = 2, nbfonc
!
            tbfonc = nomfon(ifo)
            tblp = tbfonc//'.TBLP'
!
            call jeveuo(tblp, 'L', itblp)
            lstfon = zk24(itblp+2)
            call jeveuo(lstfon, 'L', ilfon)
            fonc1 = zk8(ilfon)//'           '
            fonc2 = zk8(ilfon+1)//'           '
!
            call jelibe(lstfon)
            call jelibe(tblp)
            call jelibe(fvale1)
            call jelibe(fvale2)
!
            fvale1 = fonc1//'.VALE'
            fprol1 = fonc1//'.PROL'
!
            fvale2 = fonc2//'.VALE'
            fprol2 = fonc2//'.PROL'
!
            call jeveuo(fvale1, 'L', ifabs1)
            call jeveuo(fprol1, 'L', ifpro1)
!
            call jeveuo(fvale2, 'L', ifabs2)
            call jeveuo(fprol2, 'L', ifpro2)
!
            iford1 = ifabs1+nbp1
            iford2 = ifabs2+nbp2
!
            call fointr(fonc1, zk24(ifpro1), nbp1, zr(ifabs1), zr(iford1), &
                        nbp1, disc(1), vale(1, ifo), ier1)
            call fointr(fonc2, zk24(ifpro2), nbp2, zr(ifabs2), zr(iford2), &
                        nbp2, disc(nbp1+1), vale(nbp1+1, ifo), ier2)
!
!
            if (ier1 .ne. 0 .or. ier2 .ne. 0) then
                call utmess('F', 'MODELISA4_40')
            end if
!
            call jelibe(fprol1)
            call jelibe(fprol2)
!
        end do
!
!     --- CONCATENATION DES DISCRETISATIONS DES FONCTIONS DE FORME
!       - SUR (0,2L)
        long = disc(nbp1)
        do ip = 1, nbp2
            disc(nbp1+ip) = zr(ifabs2+ip-1)+long
        end do
!
    end if
!
    call jelibe(fvale1)
    call jelibe(fvale2)
!
    call jedema()
end subroutine
