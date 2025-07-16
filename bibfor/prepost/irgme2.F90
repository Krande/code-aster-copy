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
subroutine irgme2(numold, ima, connex, nbord2, tabd, &
                  tabl, tabv, partie, jtype, nbno, &
                  listno, nbcmp, ifi, iadmax)
    implicit none
#include "jeveux.h"
#include "asterfort/cesexi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/r8inir.h"
#include "asterfort/utmess.h"
!
    integer(kind=8) :: numold(*), tabd(*), tabl(*), tabv(*), nbno
    integer(kind=8) :: listno(*), nbcmp, ifi, ima, nbord2, iadmax, jtype
    character(len=24) :: connex
    character(len=*) :: partie
!
!     BUT: ECRITURE DES CMP D'UN CHAMP TENSORIEL PAR ELEMENT
!     POUR UN TYPE D'ELEMENT AU FORMAT GMSH
!
!     ENTREE:
!     NUMOLD : I   : TABLEAU DE CORRESPONDANCE NOUV MAILLE ANC. MAILLE
!     IMA    : I   : NUMERO NOUVELLE MAILLE
!     CONNEX : I   : CONNECTIVITE ANCIEN MAILLAGE
!     NBORD2 : I   : NOMBRE DE NUM D'ORDRE
!     TABD   : I   : DESCRIPTEURS DU CHAMP SIMPLE A IMPRIMER
!     TABL   : I   : DESCRIPTEURS DU CHAMP SIMPLE A IMPRIMER
!     TABV   : I   : DESCRIPTEURS DU CHAMP SIMPLE A IMPRIMER
!     PARTIE : K4  : IMPRESSION DE LA PARTIE COMPLEXE OU REELLE DU CHAMP
!     JTYPE  : I   : ADRESSE DU TYPE DU CHAMP ( REEL OU COMPLEXE )
!     NBNO   : I   : NOMBRE NOEUD DE LA NOUVELLE MAILLE
!     LISTNO : I   : LISTE DES NOEUDS DE LA NOUVELLE MAILLE
!     NBCMP  : I   : NOMBRE DE COMPOSANTES DU CHAMP
!     IFI    : I   : NUMERO D'UNITE LOGIQUE DU FICHIER GMSH
!     SORTIE
!     IADMAX  : I   : MAX DES IAD SI >0 LE CHAMP EXISTE POUR LA MAILLE
!
!     ------------------------------------------------------------------
    integer(kind=8) :: imaold, jcnold, ior, jcesd, jcesl, jcesv, nbpt, nbsp, j, ino
    integer(kind=8) :: itrou, ipt, inold, isp, iad, k
    real(kind=8) :: vale
    real(kind=8) :: val2(6)
!
!     ------------------------------------------------------------------
!
    call jemarq()
!
    imaold = numold(ima)
    call jeveuo(jexnum(connex, imaold), 'L', jcnold)
!
! --- ON NE TRAITE QUE LES CHAMPS A 1 SOUS-POINT,
!     ET UNE SEULE VALEUR SCALAIRE (COMPOSANTE K DE LA BOUCLE 51)
!
    isp = 1
    iadmax = 0
    k = 0
    call r8inir(6, 0.d0, val2, 1)
    do ior = 1, nbord2
        jcesd = tabd(ior)
        jcesl = tabl(ior)
        jcesv = tabv(ior)
        nbpt = zi(jcesd-1+5+4*(imaold-1)+1)
        nbsp = zi(jcesd-1+5+4*(imaold-1)+2)
        if (nbsp .ne. 1) then
            call utmess('F', 'PREPOST2_57')
        end if
!
        itrou = 0
        if (zk8(jtype-1+ior) .eq. 'R') then
            do j = 1, nbno
                ino = listno(j)
                itrou = 0
                do ipt = 1, nbpt
                    inold = zi(jcnold-1+ipt)
                    if (ino .eq. inold) then
                        itrou = 1
                        do k = 1, nbcmp
                            call cesexi('C', jcesd, jcesl, imaold, ipt, &
                                        isp, k, iad)
                            if (iad .gt. 0) then
                                vale = zr(jcesv-1+iad)
                                if (abs(vale) .le. 1.d-99) vale = 0.d0
                                val2(k) = vale
                                iadmax = iad
                            else
                                vale = 0.d0
                                val2(k) = vale
                            end if
                        end do
                    end if
                    write (ifi, 1010) val2(1), val2(4), val2(5), val2(4), &
                        val2(2), val2(6), val2(5), val2(6), val2(3)
                    goto 15
                end do
                if (itrou .eq. 0) then
                    call utmess('F', 'PREPOST2_58')
                end if
15              continue
            end do
        else if (zk8(jtype-1+ior) .eq. 'C') then
            do j = 1, nbno
                ino = listno(j)
                itrou = 0
                do ipt = 1, nbpt
                    inold = zi(jcnold-1+ipt)
                    if (ino .eq. inold) then
                        itrou = 1
                        do k = 1, nbcmp
                            call cesexi('C', jcesd, jcesl, imaold, ipt, &
                                        isp, k, iad)
                            if (iad .gt. 0) then
                                if (partie .eq. 'REEL') then
                                    vale = dble(zc(jcesv-1+iad))
                                else if (partie .eq. 'IMAG') then
                                    vale = dimag(zc(jcesv-1+iad))
                                end if
                                if (abs(vale) .le. 1.d-99) vale = 0.d0
                                val2(k) = vale
                                iadmax = iad
                            else
                                vale = 0.d0
                                val2(k) = vale
                            end if
                        end do
                    end if
                    write (ifi, 1010) val2(1), val2(4), val2(5), val2(4), &
                        val2(2), val2(6), val2(5), val2(6), val2(3)
                    goto 25
                end do
                if (itrou .eq. 0) then
                    call utmess('F', 'PREPOST2_58')
                end if
25              continue
            end do
        end if
    end do
!
    call jedema()
!
1010 format(1p, 9(e15.7e3))
!
end subroutine
