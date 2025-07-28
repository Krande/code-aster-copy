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

subroutine asasve(vechar, numedd, typres, vachar)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/assvec.h"
#include "asterfort/corich.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/gcnco2.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/vemare.h"
#include "asterfort/reajre.h"
#include "asterfort/utmess.h"
#include "asterfort/vtcopy.h"
#include "asterfort/vtcreb.h"
#include "asterfort/wkvect.h"
!
    character(len=*) :: numedd, typres, vechar
    character(len=24) :: vachar
!
! BUT : ASSEMBLER UN VECT_ELEM RESPECTANT CERTAINES CONVENTIONS
!  =============================================================
!
! IN/JXVAR  VECHAR  (K19) : VECT_ELEM A ASSEMBLER.
! IN/JXIN   NUMEDD  (K14): NUME_DDL DU SYSTEME ASSEMBLE
! IN        TYPRES  (K1) : 'R' OU 'C' (REELS OU COMPLEXES)
! OUT/JXOU VACHAR  (K24): OBJET JEVEUX CONTENANT LA LISTE DES CHAM_NO
!                    RESULTAT DE L'ASSEMBLAGE DES DIFFERENTS RESUELEM
!                    DU VECT_ELEM (VECHAR).
!
!  REMARQUES IMPORTANTES :
!  ======================
!
!   - LE NOM DU VACHAR EST OBTENU EN PRENANT LE NOM DU VECT_ELEM
!     ET EN METTANT UN "A" EN 4EME POSITION
!
!   - POUR CHAQUE RESUELEM DU VECT_ELEM (VECHAR), ON CREE UN CHAM_NO
!     PAR 1 APPEL A ASSVEC.
!
!   - SI LE VECT_ELEM EST TRUANDE (PAR EXEMPLE S'IL VIENT DE VECHME)
!     CERTAINS DES RESUELEM N'EN SONT PAS : CE SONT DEJA DES CHAM_NO
!     DANS CE CAS, ON NE L'ASSEMBLE PAS, MAIS ON LE RECOPIE.
!
!   - SI LE VECT_ELEM EST BIDON ON REND UN VACHAR BIDON CONTENANT
!     1 SEUL CHAM_NO NUL.
!     1 VECT_ELEM EST BIDON SI IL NE CONTIENT AUCUN CHAMP (LONUTI=0)
!
!   - ATTENTION : LE VECT_ELEM EST DETRUIT A LA FIN DE LA ROUTINE
!
!
!
!
!
    integer(kind=8) :: nbvec, ityp, jass, i, iret, icha
    integer(kind=8) :: n1, jvacha
    aster_logical :: bidon
    character(len=4) :: tych
    character(len=8) :: modele, newnom, vacha8
    character(len=19) :: chamno, resuElem, vecele
    character(len=24), pointer :: relr(:) => null()
!
! DEB ------------------------------------------------------------------
    call jemarq()
!
    vecele = vechar
    vacha8 = vecele(1:3)//'A'//vecele(5:8)
    chamno = vacha8//'.???????'
    newnom = '.0000000'
!
!
!     1. SI LE VECT_ELEM N'EXISTE PAS : ERREUR FATALE
!     --------------------------------------------------------
    call jeexin(vecele//'.RELR', iret)
    if (iret .eq. 0) then
        call utmess('F', 'ALGORITH_13', sk=vecele)
    end if
    call jelira(vecele//'.RELR', 'LONUTI', nbvec)
    call jeveuo(vecele//'.RELR', 'E', vk24=relr)
!
!
!     2. DESTRUCTION ET RE-ALLOCATION DE VACHAR :
!     --------------------------------------------------------
    call jeexin(vacha8, iret)
    if (iret .gt. 0) then
        call jeveuo(vacha8, 'L', jvacha)
        call jelira(vacha8, 'LONMAX', n1)
        do i = 1, n1
            call detrsd('CHAMP_GD', zk24(jvacha-1+i) (1:19))
        end do
        call jedetr(vacha8)
    end if
    call wkvect(vacha8, 'V V K24', max(nbvec, 1), jass)
!
!
!     2. SI IL N'Y A RIEN A FAIRE (VECT_ELEM BIDON):
!     --------------------------------------------------------
    bidon = .false.
    if (nbvec .eq. 0) bidon = .true.
!
    if (bidon) then
        call gcnco2(newnom)
        chamno(10:16) = newnom(2:8)
        call corich('E', chamno, ichin_=-2)
        call vtcreb(chamno, 'V', typres, &
                    nume_ddlz=numedd)
        zk24(jass-1+1) = chamno
        goto 30
    end if
!
!
!     3. SI IL FAUT FAIRE QUELQUE CHOSE :
!     --------------------------------------------------------
    call dismoi('NOM_MODELE', numedd, 'NUME_DDL', repk=modele)
    call vemare('V', '&&ASASVE', modele)
    call reajre('&&ASASVE', ' ', 'V')
!
    ityp = 1
    if (typres .eq. 'C') ityp = 2
!
    do i = 1, nbvec
        resuElem = relr(i) (1:19)
        call corich('L', resuElem, ichout_=icha)
        ASSERT((icha .ne. 0) .and. (icha .ge. -2))
!
        call gcnco2(newnom)
        chamno(10:16) = newnom(2:8)
        call corich('E', chamno, ichin_=icha)
        zk24(jass+i-1) = chamno
!
!       -- SI LE RESUELEM EST UN RESUELEM !
        call dismoi('TYPE_CHAMP', resuElem, 'CHAMP', repk=tych)
        if (tych .eq. 'RESL') then
            call jedetr('&&ASASVE           .RELR')
            call reajre('&&ASASVE', resuElem, 'V')
            call assvec('V', chamno, 1, '&&ASASVE           .RELR', [1.d0], &
                        numedd, vectScalType_=ityp)
!
!
!       -- SI LE RESUELEM N'EST PAS UN RESUELEM !(CHAM_NO)
        else if (tych .eq. 'NOEU') then
            call vtcreb(chamno, 'V', typres, nume_ddlz=numedd)
            call vtcopy(resuElem, chamno, iret)
            if (iret .ne. 0) then
                call utmess("A", "FIELD0_16")
            end if
!
        else
            ASSERT(.false.)
        end if
!
    end do
    call jedetr('&&ASASVE           .RELR')
    call jedetr('&&ASASVE           .RERR')
!
!
30  continue
!
!
!
!     DESTRUCTION DU VECT_ELEM :
!     -----------------------------------
    do i = 1, nbvec
        call corich('S', relr(i))
        call detrsd('CHAMP_GD', relr(i))
    end do
    call jedetr(vecele//'.RELR')
    call jedetr(vecele//'.RERR')
!
    vachar = vacha8
!
!
    call jedema()
end subroutine
