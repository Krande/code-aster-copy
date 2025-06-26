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

subroutine vtcrec(champ, chmod, base, typc, neq)
    implicit none
#include "jeveux.h"
#include "asterfort/dismoi.h"
#include "asterfort/gnomsd.h"
#include "asterfort/copisd.h"
#include "asterfort/jedema.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jedupo.h"
#include "asterfort/jeecra.h"
#include "asterfort/jemarq.h"
#include "asterfort/sdchgd.h"
#include "asterfort/wkvect.h"
    character(len=*) :: champ, base, typc, chmod
! person_in_charge: jacques.pellet at edf.fr
!     ------------------------------------------------------------------
!     CREATION D'UNE STRUCTURE CHAM_NO A PARTIR D'UN MODELE : CHMOD
!     ------------------------------------------------------------------
!     IN  CHAMP  : K19 : NOM DU CHAM_NO A CREER
!     IN  CHMOD  : K29 : NOM DU CHAMP MODELE
!     IN  BASE   : CH1 : NOM DE LA BASE SUR LAQUELLE LE CHAM_NO DOIT
!                        ETRE CREE
!     IN  TYPC   :     : TYPE DES VALEURS DU CHAM_NO A CREER
!              'R'  ==> COEFFICIENTS REELS
!              'C'  ==> COEFFICIENTS COMPLEXES
!              'K8' ==> COEFFICIENTS CARACTERE*8
!     REMARQUE:  AUCUN CONTROLE SUR LE "TYPC" QUE L'ON PASSE TEL QUEL
!                A JEVEUX
!     ------------------------------------------------------------------
!     DETAILS :
!       1) cette routine ne fonctione pas avec les cham_ni a
!          representation constante
!       2) LES COEFFICIENTS DU CHAM_NO "CHAMP" NE SONT PAS AFFECTES
!                 (I.E.  LE .VALE EST VIERGE)
!     ------------------------------------------------------------------
!
!
!     ------------------------------------------------------------------
    integer(kind=8) :: lchamp, jrefn
    character(len=1) :: classe
    character(len=1) :: type, type2
    character(len=8) :: nomgd
    character(len=19) :: nume_equa, nume_equa_new
    character(len=24) :: vale, refe, noojb
!     ------------------------------------------------------------------
    integer(kind=8) :: ibid, neq
    character(len=19) :: chmod2
!     ------------------------------------------------------------------
    data vale/'                   .VALE'/
    data refe/'                   .REFE'/
!     DEB --------------------------------------------------------------
    call jemarq()
!
    refe(1:19) = champ
    vale(1:19) = champ
    chmod2 = chmod
!
    classe = base(1:1)
    if (typc(1:1) .eq. 'K') then
        type = 'F'
    else
        type = typc(1:1)
    end if
!
!     -- RECOPIE DE L'OBJET .REFE MODELE :
    call jedupo(chmod2//'.REFE', classe, refe, .false._1)
    call jeecra(refe, 'DOCU', ibid, 'CHNO')
!
    call dismoi("TYPE_SCA", chmod2, "CHAM_NO", repk=type2)
    if (type .ne. type2) then
        noojb = '12345678.NUME000000.PRNO'
        call gnomsd(champ, noojb, 14, 19)
        noojb(1:8) = champ(1:8)
        nume_equa_new = noojb(1:19)
        call dismoi("NUME_EQUA", chmod2, "CHAM_NO", repk=nume_equa)
        call copisd("NUME_EQUA", classe, nume_equa, nume_equa_new)
        call jeveuo(nume_equa_new//".REFN", "E", jrefn)
        call dismoi("NOM_GD", chmod2, "CHAM_NO", repk=nomgd)
        zk24(jrefn-1+2) = nomgd(1:5)//type
        call jeveuo(champ(1:19)//".REFE", "E", jrefn)
        zk24(jrefn-1+2) = nume_equa_new

    end if
!
!     -- CREATION DE L'OBJET .VALE :
    call wkvect(vale, classe//' V '//type, neq, lchamp)
!
!     -- CHANGER LE TYPE SCALAIRE DE LA GRANDEUR ---
    call sdchgd(champ, type)
!
    call jedema()
end subroutine
