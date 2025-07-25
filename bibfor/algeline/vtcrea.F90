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

subroutine vtcrea(champ, crefe, base, typc, neq)
    implicit none
#include "jeveux.h"
#include "asterfort/dismoi.h"
#include "asterfort/copisd.h"
#include "asterfort/exisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/codent.h"
#include "asterfort/idensd.h"
#include "asterfort/gnomsd.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jemarq.h"
#include "asterfort/sdchgd.h"
#include "asterfort/wkvect.h"
    character(len=*) :: champ, base, typc
    character(len=24) :: crefe(*)
!     CREATION D'UNE STRUCTURE CHAM_NO A PARTIR D'UN MODELE : CREFE
!     LE CHAM_NO MODELE NE DOIT PAS ETRE A REPRESENTATION CONSTANTE.
!     ------------------------------------------------------------------
!     IN  CHAMP  : K19 : NOM DU CHAM_NO A CREER
!     IN  CREFE  : K24 : CONTENU DE L'OBJET .REFE D'UN CHAM_NO MODELE
!                (1) :  K8  : MODELE
!                (2) :  K19 : NUME_EQUA
! IN  BASE   : CH1 : NOM DE LA BASE SUR LAQUELLE LE CHAM_NO DOIT ETRE
!                    CREER
!     IN  TYPC   :     : TYPE DES VALEURS DU CHAM_NO A CREER
!              'R'  ==> COEFFICIENTS REELS
!              'C'  ==> COEFFICIENTS COMPLEXES
!              'K8' ==> COEFFICIENTS CARACTERE*8
!     REMARQUE:  AUCUN CONTROLE SUR LE "TYPC" QUE L'ON PASSE TEL QUEL
!                A JEVEUX
!     ------------------------------------------------------------------
!     PRECAUTIONS D'EMPLOI :
!       1) LE CHAM_NO "CHAMP" NE DOIT PAS EXISTER
!       2) LES COEFFICIENTS DU CHAM_NO "CHAMP" NE SONT PAS AFFECTES
!                 (I.E.  LE .VALE EST VIERGE)
!     ------------------------------------------------------------------
!
!
!     ------------------------------------------------------------------
    integer(kind=8) :: lchamp, jrefn, prev, iexi
    character(len=1) :: classe
    character(len=1) :: type, type2
    character(len=8) :: nomgd
    character(len=19) :: nume_equa, nume_equa_tmp
    character(len=24) :: vale, refe, noojb
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: ibid, neq
!-----------------------------------------------------------------------
    data vale/'                   .VALE'/
    data refe/'                   .REFE'/
!     DEB --------------------------------------------------------------
    call jemarq()
    classe = base(1:1)
    if (typc(1:1) .eq. 'K') then
        type = 'F'
    else
        type = typc(1:1)
    end if
!
    nume_equa = crefe(2) (1:19)
    call dismoi('TYPE_SCA', nume_equa, 'NUME_EQUA', repk=type2)
    if (type .ne. type2) then
        call dismoi('NOM_GD', nume_equa, 'NUME_EQUA', repk=nomgd)
        nomgd(6:6) = type(1:1)
        noojb = '12345678.NUME000000.PRNO'
        call gnomsd(champ, noojb, 14, 19)
        noojb(1:8) = champ(1:8)
        nume_equa_tmp = noojb(1:19)
        call copisd("NUME_EQUA", classe, nume_equa, nume_equa_tmp)
        nume_equa = nume_equa_tmp
        call jeveuo(nume_equa//".REFN", 'E', jrefn)
        zk24(jrefn+1) = nomgd
        read (nume_equa(14:19), '(i6)') prev
        if (prev > 0) then
            prev = max(0, prev-1)
            nume_equa_tmp = nume_equa
            call codent(prev, "D0", nume_equa_tmp(14:19))
            call exisd('NUME_EQUA', nume_equa_tmp, iexi)
            if (iexi .gt. 0) then
                if (idensd('NUME_EQUA', nume_equa, nume_equa_tmp)) then
                    call detrsd("NUME_EQUA", nume_equa)
                    nume_equa = nume_equa_tmp
                end if
            end if
        end if
    end if
!
!     --- RECOPIE DE L'OBJET .REFE MODELE :
    refe(1:19) = champ
    call wkvect(refe, classe//' V K24', 4, lchamp)
    call jeecra(refe, 'DOCU', ibid, 'CHNO')
    zk24(lchamp-1+2) = nume_equa
!
!     -- CREATION DE L'OBJET .VALE :
    vale(1:19) = champ
    call wkvect(vale, classe//' V '//type, neq, lchamp)
!
!     -- CHANGER LE TYPE SCALAIRE DE LA GRANDEUR ---
    call sdchgd(champ, type)
    call jedema()
end subroutine
