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
subroutine mltfc1(nbloc, ncbloc, decal, supnd, fils, &
                  frere, seq, lgsn, lfront, adress, &
                  local, adpile, nbass, pile, lgpile, &
                  adper, t1, t2, factol, factou, &
                  typsym, ad, eps, ier, nbb, &
                  cl, cu, diag)
!
! person_in_charge: olivier.boiteau at edf.fr
!     VERSION MODIFIEE POUR L' APPEL A DGEMV (PRODUITS MATRICE-VECTEUR)
!     LE STOCKAGE DES COLONNES DE LA FACTORISEE EST MODIFIE, ET AINSI
!      ADPER LES COLONNES FORMENT UN BLOC RECTANGULAIRE
!
! aslint: disable=W1504
    implicit none
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jelibe.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/mlnflm.h"
#include "asterfort/mlnfmj.h"
#include "asterfort/mltaff.h"
#include "asterfort/mltafp.h"
#include "asterfort/mltf21.h"
#include "asterfort/mltflm.h"
#include "asterfort/mltfmj.h"
!
    integer(kind=8) :: pmin, nbb
    parameter(pmin=10)
    integer(kind=8) :: nbloc, ncbloc(*), decal(*)
    integer(kind=8) :: lgsn(*), lfront(*), lgpile, typsym
    integer(kind=4) :: local(*)
    integer(kind=8) :: nbass(*), adpile(*), fils(*), supnd(*)
    integer(kind=8) :: adress(*), frere(*), seq(*), ad(*), ier
    real(kind=8) :: pile(*), eps, cl(nbb, nbb, *), cu(nbb, nbb, *), diag(*)
    character(len=24) :: factol, factou
!
    real(kind=8) :: t1(*), t2(*)
    integer(kind=8) :: adper(*), ifacl, ifacu, lmatf
    integer(kind=8) :: itemp, i, j, isnd, sni, sn, n, m, p, nl, nc, mem, adfacl, ib, nb
    integer(kind=8) :: lm1, lm2, adfacu, long, iad, adfac0, adfac
!
    call jemarq()
    itemp = 1
    mem = 0
    isnd = 0
    do i = 1, lgpile
        pile(i) = 0.d0
    end do
!
    do ib = 1, nbloc
        call jeveuo(jexnum(factol, ib), 'E', ifacl)
        if (typsym .eq. 0) then
            call jeveuo(jexnum(factou, ib), 'E', ifacu)
        end if
        adfac0 = ifacl-1
        do nc = 1, ncbloc(ib)
            isnd = isnd+1
            sni = seq(isnd)
            long = adress(sni+1)-adress(sni)
            iad = supnd(sni)
            p = lgsn(sni)
!
!
            m = lfront(sni)
            n = m+p
            lmatf = (m*(m+1))/2
!
            lm1 = lmatf
            if (typsym .eq. 0) lmatf = 2*lmatf
!         CHANGTPOUR L' APPEL A DGEMV
            do i = 1, p
                adper(i) = (i-1)*n+i
            end do
            do i = p, n-1
                adper(i+1) = 1+(n+(n-i+1))*i/2
            end do
            sn = fils(sni)
            do j = 1, lmatf
                pile(itemp+j-1) = 0.d0
            end do
            adfacl = ifacl-1+decal(sni)
            if (typsym .eq. 0) adfacu = ifacu-1+decal(sni)
40          continue
!     DO WHILE (SN.NE.0)
            if (sn .ne. 0) then
                nl = lgsn(sn)
                nb = nbass(sn)
                lm2 = (lfront(sn)*(lfront(sn)+1))/2
                call mltafp(lfront(sn), nb, adper, zr(adfacl), pile(adpile(sn)), &
                            local(adress(sn)+nl))
                call mltaff(lfront(sn), nb, adper, pile(itemp), pile(adpile(sn)), &
                            local(adress(sn)+nl), p)
                if (typsym .eq. 0) then
                    call mltafp(lfront(sn), nb, adper, zr(adfacu), pile(adpile(sn)+lm2), &
                                local(adress(sn)+nl))
                    call mltaff(lfront(sn), nb, adper, pile(itemp+lm1), pile(adpile(sn)+lm2), &
                                local(adress(sn)+nl), p)
                end if
                sn = frere(sn)
                goto 40
!     FIN DO WHILE
            end if
            if (p .le. pmin .and. typsym .ne. 0) then
                call mltf21(p, zr(adfacl), pile(itemp), n, t1, &
                            t2, eps, ier)
                if (ier .ne. 0) goto 999
            else
                if (typsym .eq. 0) then
                    call mlnflm(nbb, n, p, zr(adfacl), zr(adfacu), &
                                adper, t1, t2, ad, eps, &
                                ier, cl, cu)
                    if (ier .ne. 0) goto 999
                    call mlnfmj(nbb, n, p, zr(adfacl), zr(adfacu), &
                                pile(itemp), pile(itemp+lm1), adper, t1, t2, &
                                cl, cu)
                else
                    call mltflm(nbb, n, p, zr(adfacl), adper, &
                                t1, ad, eps, ier, cl)
                    if (ier .ne. 0) goto 999
                    call mltfmj(nbb, n, p, zr(adfacl), pile(itemp), &
                                adper, t1, cl)
                end if
            end if
            if (fils(sni) .ne. 0) then
                mem = max(mem, (itemp+lmatf-1))
                do j = 1, lmatf
                    pile(adpile(fils(sni))+j-1) = pile(itemp+j-1)
                end do
                adpile(sni) = adpile(fils(sni))
                itemp = adpile(sni)+lmatf
            else
                adpile(sni) = itemp
                itemp = itemp+lmatf
            end if
            mem = max(mem, itemp)
!
            do i = 1, p
                adfac = adfac0+long*(i-1)+i
                diag(iad+i-1) = zr(adfac)
            end do
            adfac0 = adfac0+long*lgsn(sni)
        end do
!
        call jelibe(jexnum(factol, ib))
        if (typsym .eq. 0) call jelibe(jexnum(factou, ib))
    end do
999 continue
    if (ier .ne. 0) ier = ier+supnd(sni)-1
    call jedema()
end subroutine
