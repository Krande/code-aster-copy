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

subroutine exisd(typesd, nomsd, iret)
    implicit none
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: iret
    character(len=*) :: typesd, nomsd
! ----------------------------------------------------------------------
!  BUT : DETERMINER SI UNE SD EXISTE
!  IN   TYPESD : TYPE DE LA STRUCTURE DE DONNEE A TESTER
!         / 'CARTE'        /'CHAM_NO'      /'CHAM_ELEM'   /'RESUELEM'
!         / 'CHAM_ELEM_S'  /'CHAM_NO_S'
!         / 'CHAMP' (CHAPEAU AUX CHAM_NO/CHAM_ELEM/CARTE/RESUELEM)
!         / 'CHAMP_GD' (CHAPEAU DESUET AUX CHAM_NO/CHAM_ELEM/...)
!         / 'TABLE'
!         / 'RESULTAT'
!         / 'FONCTION'
!         / 'MODELE'
!         / 'PARTITION'
!         / 'MAILLAGE'
!         / 'NUME_DDL'
!         / 'NUME_EQUA'
!         / 'MATR_ASSE'
!       NOMSD   : NOM DE LA STRUCTURE DE DONNEES A TESTER
!
!  OUT:  IRET   : 0 -> LA SD N'EXISTE PAS
!                 1 -> LA SD EXISTE
! ----------------------------------------------------------------------
!
    integer(kind=8) :: i1, i2, i3, i4, i5
    character(len=8) :: ch8
    character(len=16) :: typ2sd
    character(len=19) :: ch
! -DEB------------------------------------------------------------------
!
    call jemarq()
    typ2sd = typesd
!
    iret = 0
!
    if (typ2sd .eq. 'MAILLAGE') then
!     ------------------------------
        ch8 = nomsd
        call jeexin(ch8//'.DIME', i1)
        call jeexin(ch8//'.COORDO', i2)
        if (i1*i2 .ne. 0) iret = 1
!
!
    else if (typ2sd .eq. 'MODELE') then
!     ------------------------------
        ch8 = nomsd
        call jeexin(ch8//'.MODELE    .LIEL', i1)
        if (i1 .ne. 0) iret = 1
!
!
    else if (typ2sd .eq. 'PARTITION') then
!     ------------------------------
        ch8 = nomsd
        call jeexin(ch8//'.PRTI', i1)
        call jeexin(ch8//'.PRTK', i2)
        if (i1*i2 .ne. 0) iret = 1
!
    else if (typ2sd .eq. 'CARTE') then
!     ------------------------------
        ch = nomsd
        call jeexin(ch//'.NOMA', i1)
        call jeexin(ch//'.DESC', i2)
        call jeexin(ch//'.VALE', i3)
        if (i1*i2*i3 .ne. 0) iret = 1
!
!
    else if (typ2sd .eq. 'CHAM_NO') then
!     ------------------------------
        ch = nomsd
        call jeexin(ch//'.REFE', i1)
        call jeexin(ch//'.VALE', i2)
        if (i1*i2 .ne. 0) iret = 1

    else if (typ2sd .eq. 'CHAM_GEOM') then
!     ------------------------------
        ch = nomsd
        call jeexin(ch//'.VALE', i1)
        if (i1 .ne. 0) iret = 1
!
!
    else if (typ2sd .eq. 'CHAM_ELEM') then
!     ------------------------------
        ch = nomsd
        call jeexin(ch//'.CELD', i1)
        call jeexin(ch//'.CELV', i2)
        if (i1*i2 .ne. 0) iret = 1
!
!
    else if (typ2sd .eq. 'RESUELEM') then
!     ------------------------------
        ch = nomsd
        call jeexin(ch//'.DESC', i1)
        call jeexin(ch//'.RESL', i2)
        call jeexin(ch//'.NOLI', i3)
        if (i1*i2*i3 .ne. 0) iret = 1
!
!
    else if ((typ2sd .eq. 'CHAMP') .or. (typ2sd .eq. 'CHAMP_GD')) then
!     -------------------------------------------------------------
        ch = nomsd
!
!       -- CHAM_ELEM ?
        call jeexin(ch//'.CELD', i1)
        call jeexin(ch//'.CELV', i2)
        if (i1*i2 .ne. 0) iret = 1
!
!       -- CHAM_NO OU CARTE ?
        call jeexin(ch//'.VALE', i1)
        if (i1 .ne. 0) iret = 1
!
!       -- RESUELEM ?
        call jeexin(ch//'.DESC', i1)
        call jeexin(ch//'.RESL', i2)
        call jeexin(ch//'.NOLI', i3)
        if (i1*i2*i3 .ne. 0) iret = 1
!
!
    else if (typ2sd .eq. 'CHAM_NO_S') then
!     ------------------------------------
        ch = nomsd
        call jeexin(ch//'.CNSD', i1)
        call jeexin(ch//'.CNSV', i2)
        call jeexin(ch//'.CNSL', i3)
        if (i1*i2*i3 .ne. 0) iret = 1
!
!
    else if (typ2sd .eq. 'CHAM_ELEM_S') then
!     --------------------------------------
        ch = nomsd
        call jeexin(ch//'.CESD', i1)
        call jeexin(ch//'.CESV', i2)
        call jeexin(ch//'.CESL', i3)
        if (i1*i2*i3 .ne. 0) iret = 1
!
!
    else if (typ2sd .eq. 'TABLE') then
!     --------------------------------
        ch = nomsd
        call jeexin(ch//'.TBBA', i1)
        call jeexin(ch//'.TBNP', i2)
        call jeexin(ch//'.TBLP', i3)
        if (i1*i2*i3 .ne. 0) iret = 1
!
!
    else if (typ2sd .eq. 'RESULTAT') then
!     -----------------------------------
        ch = nomsd
        call jeexin(ch//'.DESC', i1)
        call jeexin(ch//'.NOVA', i2)
        call jeexin(ch//'.TAVA', i3)
        call jeexin(ch//'.ORDR', i4)
        call jeexin(ch//'.TACH', i5)
        if (i1*i2*i3*i4*i5 .ne. 0) iret = 1
!
    else if (typ2sd .eq. 'LIGREL') then
!     -----------------------------------
        ch = nomsd
        call jeexin(ch//'.LGRF', i1)
        call jeexin(ch//'.NBNO', i2)
        if (i1*i2 .ne. 0) iret = 1
!
    else if (typ2sd .eq. 'FONCTION') then
!     -----------------------------------
        ch = nomsd
        call jeexin(ch//'.PROL', i1)
        if (i1 .ne. 0) iret = 1
!
    else if (typ2sd .eq. 'MATR_ASSE') then
!     -----------------------------------
        ch = nomsd
        call jeexin(ch//'.REFA', i2)
        call jeexin(ch//'.VALM', i3)
        if (i2*i3 .ne. 0) iret = 1
!
    else if (typ2sd .eq. 'NUME_DDL') then
!     -----------------------------------
        ch = nomsd

        call jeexin(ch(1:14)//'.NUME.DEEQ', i1)
        call jeexin(ch(1:14)//'.NUME.DELG', i2)
        call jeexin(ch(1:14)//'.NUME.LILI', i3)
        call jeexin(ch(1:14)//'.NUME.NUEQ', i4)
        if (i1*i2*i3*i4 .ne. 0) iret = 1
!
    else if (typ2sd .eq. 'NUME_EQUA') then
        !     -----------------------------------
        ch = nomsd
        call jeexin(ch(1:19)//'.PRNO', i1)
        call jeexin(ch(1:19)//'.DEEQ', i2)
        call jeexin(ch(1:19)//'.LILI', i3)
        call jeexin(ch(1:19)//'.NUEQ', i4)
        call jeexin(ch(1:19)//'.REFN', i5)
        if (i1*i2*i3*i4*i5 .ne. 0) iret = 1
!
    else if (typ2sd .eq. 'DOMJOINTS') then
        !     -----------------------------------
        ch = nomsd
        call jeexin(ch(1:19)//'.DOMJ', i1)
        call jeexin(ch(1:19)//'.SEND', i2)
        call jeexin(ch(1:19)//'.RECV', i3)
        call jeexin(ch(1:19)//'.GCOM', i4)
        call jeexin(ch(1:19)//'.PGID', i5)
        if (i1*i2*i3*i4*i5 .ne. 0) iret = 1
!
    else
        call utmess('F', 'UTILITAI_47', sk=typ2sd)
    end if
!
    call jedema()
end subroutine
