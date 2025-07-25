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
subroutine irgmec(numold, ima, connex, nbord2, tabd, &
                  tabl, tabv, partie, jtype, nbno, &
                  listno, icmp, ifi, iwri, iadmax, &
                  ordr, chamsy, nomcon, lresu)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/cesexi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/rsexch.h"
#include "asterfort/utmess.h"
!
    integer(kind=8) :: numold(*), tabd(*), tabl(*), tabv(*), nbno
    integer(kind=8) :: listno(*), icmp, ifi, ima, nbord2, iadmax, jtype, ordr(nbord2)
    aster_logical :: iwri, lresu
    character(len=*) :: partie, chamsy, nomcon
    character(len=24) :: connex
!
!     BUT: ECRITURE D'UNE CMP D'UN CHAMP PAR ELEMENT
!     POUR UN TYPE D'ELEMENT AU FORMAT GMSH
!
!     ENTREE:
!     NUMOLD : I   : TABLEAU DE CORRESPONDANCE NOUV MAILLE ANC. MAILLE
!     IMA    : I   : NUMERO NOUVELLE MAILLE
!     CONNEX : I   : CONNECTIVITE ANCIEN MAILLAGE
!     NBORD2 : I   : NOMBRE DE NUM D'ORDRE
!     TABD   : I   : DECRIPTEURS DU CHMAP SIMPLE A IMMRIMER
!     TABL   : I   : DECRIPTEURS DU CHMAP SIMPLE A IMMRIMER
!     TABV   : I   : DECRIPTEURS DU CHMAP SIMPLE A IMMRIMER
!     PARTIE : K4  : IMPRESSION DE LA PARTIE COMPLEXE OU REELLE DU CHAMP
!     JTYPE  : I   : ADRESSE DU TYPE DU CHAMP ( REEL OU COMPLEXE )
!     NBNO   : I   : NOMBRE NOEUD DE LA NOUVELLE MAILLE
!     LISTNO : I   : LISTE DES NOEUDS DE LA NOUVELLE MAILLE
!     ICMP   : I   : NUMERO COMPOSANTE CHAMP
!     IFI    : I   : NUMERO D'UNITE LOGIQUE DU FICHIER GMSH
!     IWRI   : L   : INDIQUE SI ON DOIT ECRIRE
!     ORDR   : I   : LISTE DES NUMEROS D'ORDRE A IMPRIMER
!     CHAMSY : K16 : NOM DU CHAM_ELEM A ECRIRE
!     NOMCON : K19 : NOM DU CONCEPT A IMPRIMER
!     LRESU  : L   : LOGIQUE INDIQUANT SI NOMCON EST UNE SD RESULTAT
!     SORTIE
!     IADMAX  : I   : MAX DES IAD SI >0 LE CHAMP EXISTE POUR LA MAILLE
!
!     ------------------------------------------------------------------
    integer(kind=8) :: imaold, jcnold, ior, jcesd, jcesl, jcesv, nbpt, nbsp, j, ino
    integer(kind=8) :: itrou, ipt, inold, isp, iad, iret
    real(kind=8) :: vale
    character(len=16) :: chams2
    character(len=19) :: ch19, nomco2
!
!     ------------------------------------------------------------------
!
    call jemarq()
!
    chams2 = chamsy
    nomco2 = nomcon
!
    imaold = numold(ima)
    call jeveuo(jexnum(connex, imaold), 'L', jcnold)
!
! --- ON NE TRAITE QUE LES CHAMPS A 1 SOUS-POINT,
!     ET UNE SEULE VALEUR SCALAIRE (COMPOSANTE K DE LA BOUCLE 51)
!
    isp = 1
    iadmax = 0
    ch19 = ' '
    do ior = 1, nbord2
        if (lresu) then
            call rsexch(' ', nomco2, chams2, ordr(ior), ch19, &
                        iret)
            if (iret .ne. 0) goto 11
        end if
        jcesd = tabd(ior)
        jcesl = tabl(ior)
        jcesv = tabv(ior)
        nbpt = zi(jcesd-1+5+4*(imaold-1)+1)
        nbsp = zi(jcesd-1+5+4*(imaold-1)+2)
        if (nbsp .ne. 1) then
            call utmess('F', 'PREPOST2_57')
        end if
        itrou = 0
        if (zk8(jtype-1+ior) .eq. 'R') then
            do j = 1, nbno
                ino = listno(j)
                itrou = 0
                do ipt = 1, nbpt
                    inold = zi(jcnold-1+ipt)
                    if (ino .eq. inold) then
                        itrou = 1
                        call cesexi('C', jcesd, jcesl, imaold, ipt, &
                                    isp, icmp, iad)
                        if (iad .gt. 0) then
                            vale = zr(jcesv-1+iad)
                            if (abs(vale) .le. 1.d-99) vale = 0.d0
                            iadmax = iad
                        else
                            vale = 0.d0
                        end if
                        goto 15
                    end if
                end do
15              continue
                if (iwri) write (ifi, 1000) vale
            end do
            if (itrou .eq. 0) then
                call utmess('F', 'PREPOST2_58')
            end if
        else if (zk8(jtype-1+ior) .eq. 'C') then
            do j = 1, nbno
                ino = listno(j)
                itrou = 0
                do ipt = 1, nbpt
                    inold = zi(jcnold-1+ipt)
                    if (ino .eq. inold) then
                        itrou = 1
                        call cesexi('C', jcesd, jcesl, imaold, ipt, &
                                    isp, icmp, iad)
                        if (iad .gt. 0) then
                            if (partie .eq. 'REEL') then
                                vale = dble(zr(jcesv-1+iad))
                            else if (partie .eq. 'IMAG') then
                                vale = dimag(zc(jcesv-1+iad))
                            end if
                            if (abs(vale) .le. 1.d-99) vale = 0.d0
                            iadmax = iad
                        else
                            vale = 0.d0
                        end if
                        goto 25
                    end if
                end do
                if (itrou .eq. 0) then
                    call utmess('F', 'PREPOST2_58')
                end if
25              continue
                if (iwri) write (ifi, 1000) vale
            end do
            if (itrou .eq. 0) then
                call utmess('F', 'PREPOST2_58')
            end if
        end if
11      continue
    end do
!
    call jedema()
!
1000 format(1p, e15.7e3)
!
end subroutine
