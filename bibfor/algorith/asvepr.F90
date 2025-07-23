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

subroutine asvepr(lischa, vecelz, typres, numedd)
!
!
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/assvec.h"
#include "asterfort/chor2c.h"
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
#include "asterfort/lisltc.h"
#include "asterfort/vemare.h"
#include "asterfort/reajre.h"
#include "asterfort/utmess.h"
#include "asterfort/vtcopy.h"
#include "asterfort/vtcreb.h"
#include "asterfort/wkvect.h"
    character(len=19) :: lischa
    character(len=*) :: vecelz, numedd
    character(len=1) :: typres
!
! ----------------------------------------------------------------------
!
! PREPARATION DE L'ASSEMBLAGE D'UN VECT_ELEM
!
! ----------------------------------------------------------------------
!
!
! IN  LISCHA : SD LISTE DES CHARGES
! IN  VECELE : NOM DU VECT_ELEM
! IN  TYPRES : TYPE DU RESULTAT 'R' OU 'C' (REELS OU COMPLEXES)
! IN  NUMEDD : NUME_DDL DU SYSTEME ASSEMBLE
!
! LE VACHAR CONTIENT LA LISTE DES CHAM_NO RESULTAT DE L'ASSEMBLAGE DES
! DIFFERENTS RESU_ELEM DU VECT_ELEM
!
! * POUR CHAQUE RESU_ELEM DU VECT_ELEM, ON CREE UN CHAM_NO PAR UN APPEL
!   A ASSVEC
! * SI LE VECT_ELEM EST TRUANDE (VOIR VECHME), CERTAINS DES RESU_ELEM
!   N'EN SONT PAS : CE SONT DEJA DES CHAM_NO (CHARGEMENT VECT_ASSE DANS
!   AFFE_CHAR_MECA). DANS CE CAS, ON NE L'ASSEMBLE PAS, MAIS ON LE
!   RECOPIE.
! * SI LE VECT_ELEM EST BIDON ON REND UN VACHAR BIDON CONTENANT
!     1 SEUL CHAM_NO NUL.
!     1 VECT_ELEM EST BIDON SI IL NE CONTIENT AUCUN CHAMP (LONUTI=0)
!
! ATTENTION : LE VECT_ELEM EST DETRUIT A LA FIN DE LA ROUTINE
!
!
!
!
    character(len=19) :: vecele, chamno
    character(len=24) :: vachar
    integer(kind=8) :: jvacha
    character(len=24) :: resuElem
    character(len=8) :: newnom, modele, typech, typsca
    integer(kind=8) :: ivach, nbvach
    integer(kind=8) :: nbvec
    integer(kind=8) :: iret, ivec, ichar, ityprs
    character(len=4) :: tyresl
    character(len=1) :: typchn
    character(len=24), pointer :: relr(:) => null()
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATIONS
!
    newnom = '.0000000'
    vecele = vecelz
    ASSERT(typres .eq. 'R' .or. typres .eq. 'C')
!
! --- LE VECT_ELEM EXISTE-IL ?
!
    call jeexin(vecele//'.RELR', iret)
    if (iret .eq. 0) then
        nbvec = 0
    else
        call jelira(vecele//'.RELR', 'LONUTI', nbvec)
    end if
!
! --- NOM DU CHAMNO
!
    chamno = vecele(1:8)//'.???????'
!
! --- NOM DU VACHAR
!
    vachar = vecele(1:19)//'.CHNO'
!
! --- DESTRUCTION DU VACHAR
!
    call jeexin(vachar, iret)
    if (iret .gt. 0) then
        call jeveuo(vachar, 'L', jvacha)
        call jelira(vachar, 'LONMAX', nbvach)
        do ivach = 1, nbvach
            call detrsd('CHAMP_GD', zk24(jvacha-1+ivach) (1:19))
        end do
        call jedetr(vachar)
    end if
!
! --- CREATION DU VACHAR
!
    call wkvect(vachar, 'V V K24', max(nbvec, 1), jvacha)
!
! --- SI IL N'Y A RIEN A FAIRE, ON CREE UN CHAM_NO BIDON
!
    if (nbvec .eq. 0) then
        call gcnco2(newnom)
        chamno(10:16) = newnom(2:8)
        call corich('E', chamno, ichin_=-2)
        call vtcreb(chamno, 'V', typres, &
                    nume_ddlz=numedd)
        zk24(jvacha-1+1) = chamno
        goto 99
    end if
!
! --- CREER L'OBJET .RERR DU VECT_ELEM
!
    call dismoi('NOM_MODELE', numedd, 'NUME_DDL', repk=modele)
    call vemare('V', '&&ASVEPR', modele)
!
! --- INITIALISER L'OBJET .RELR DU VECT_ELEM
!
    call reajre('&&ASVEPR', ' ', 'V')
    call jeveuo(vecele//'.RELR', 'E', vk24=relr)
!
! --- ASSEMBLAGE DES VECT_ELEM
!
    do ivec = 1, nbvec
!
! ----- NOM DU RESU_ELEM
!
        resuElem = relr(ivec)
!
! ----- PREPARATION DU NOM DU CHAM_NO
!
        call gcnco2(newnom)
        chamno(10:16) = newnom(2:8)
        zk24(jvacha-1+ivec) = chamno
!
! ----- ENREGISTREMENT DU NUMERO DE LA CHARGE DANS LE CHAM_NO
!
        call corich('L', resuElem, ichout_=ichar)
        call corich('E', chamno, ichin_=ichar)
!
! ----- TYPE DU RESU_ELEM
!
        call dismoi('TYPE_CHAMP', resuElem, 'CHAMP', repk=tyresl)
        ASSERT(tyresl .eq. 'RESL' .or. tyresl .eq. 'NOEU')
!
! ----- SI LE RESU_ELEM EST UN VRAI RESU_ELEM (ISSU DE CALCUL)
!
        if (tyresl .eq. 'RESL') then
            call jedetr('&&ASVEPR           .RELR')
            call reajre('&&ASVEPR', resuElem, 'V')
            call dismoi('TYPE_SCA', resuElem, 'RESUELEM', repk=typsca)
            if (typsca .eq. 'R') then
                ityprs = 1
            else if (typsca .eq. 'C') then
                ityprs = 2
            else
                ASSERT(.false.)
            end if
            call assvec('V', chamno, 1, '&&ASVEPR           .RELR', [1.d0], &
                        numedd, vectScalType_=ityprs)
        end if
!
! ----- SI LE RESU_ELEM EST UN FAUX RESU_ELEM (CHAM_NO)
!
        if (tyresl .eq. 'NOEU') then
            call lisltc(lischa, ichar, typech)
            typchn = 'R'
            if (typech .eq. 'COMP') typchn = 'C'
            call vtcreb(chamno, 'V', typchn, nume_ddlz=numedd)

            call vtcopy(resuElem, chamno, iret)
            if (iret .ne. 0) then
                call utmess("F", "FIELD0_16")
            end if
        end if
!
    end do
    call jedetr('&&ASVEPR           .RELR')
    call jedetr('&&ASVEPR           .RERR')
!
! --- TRANSFORMATION EN COMPLEXE SI NECESSAIRE
!
    if (typres .eq. 'C') then
        call chor2c(lischa, vecele)
    end if
!
! --- DESTRUCTION DU VECT_ELEM
!
    call jeexin(vecele//'.RELR', iret)
    do ivec = 1, nbvec
        resuElem = relr(ivec)
        call corich('S', resuElem)
        call detrsd('CHAMP_GD', resuElem)
    end do
    call jedetr(vecele//'.RELR')
    call jedetr(vecele//'.RERR')
!
99  continue
!
    call jedema()
end subroutine
