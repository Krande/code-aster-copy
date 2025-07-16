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

subroutine copich(base, ch1z, ch2z)
!
! person_in_charge: jacques.pellet at edf.fr
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/copisd.h"
#include "asterfort/dismoi.h"
#include "asterfort/gnomsd.h"
#include "asterfort/jedema.h"
#include "asterfort/jedup1.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
    character(len=1) :: base
    character(len=*) :: ch1z, ch2z
! ----------------------------------------------------------------------
!
!   BUT:
!   DUPLIQUER UN CHAMP_GD SOUS UN AUTRE NOM.
!    L'EXISTENCE DE CH1Z EST OBLIGATOIRE
!   (SI CH2Z EXISTE DEJA, ON L'ECRASE)
!
!     IN       BASE        'G' , 'V' , ... : BASE DE CREATION DE CH2
!     IN       CH1Z    K19  NOM DU CHAMP_GD A DUPLIQUER
!     IN/JXOUT CH2Z    K19  NOM DU CHAMP_GD A CREER
!
!-----------------------------------------------------------------------
!
    character(len=4) :: docu
    character(len=8) :: nomu
    character(len=16) :: concep, cmd
    character(len=19) :: prno, prno2, ch1, ch2, ligr, ligr2
    character(len=24) :: noojb
    integer(kind=8) :: iret
    character(len=24), pointer :: refe(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
!
    ch1 = ch1z
    ch2 = ch2z
!
    call dismoi("DOCU", ch1, "CHAMP", repk=docu, ier=iret, arret="C")
    if (iret > 0) go to 999
!
!     -- CAS DES CHAM_NO :
!     ----------------------
    if (docu .eq. 'CHNO') then
        call jedup1(ch1//'.REFE', base, ch2//'.REFE')
        call jedup1(ch1//'.VALE', base, ch2//'.VALE')
!
!
!       SI LE NOUVEAU CHAM_NO DOIT ETRE CREE SUR 'G', ON IMPOSE
!       QUE LE NOM DU NUME_EQUA DE CE CHAMP  COMMENCE PAR LE NOM
!       UTILISATEUR DU RESULTAT DE LA COMMANDE EN COURS :
!       --------------------------------------------------------------
        if (base .eq. 'G') then
            call getres(nomu, concep, cmd)
            ! nomu peut être vide si appel depuis le c++
            if (nomu == ' ') then
                nomu = ch2(1:8)
            end if
            call dismoi('NUME_EQUA', ch2, 'CHAM_NO', repk=prno)
!         -- REMARQUE : UN CHAM_NO PEUT NE PAS AVOIR DE NUME_EQUA (' '):
            if (prno .ne. ' ') then
                if (prno(1:8) .ne. nomu) then
                    noojb = '12345678.NUMEQ00000.PRNO'
                    call gnomsd(nomu, noojb, 15, 19)
                    prno2 = noojb(1:19)
                    call jeveuo(ch2//'.REFE', 'E', vk24=refe)
                    call copisd('NUME_EQUA', base, prno, prno2)
                    refe(2) = prno2
                end if
            end if
        end if
!
!     -- CAS DES CHAM_GEOM :
!     ----------------------
    else if (docu .eq. 'CHGO') then
        call jedup1(ch1//'.DESC', base, ch2//'.DESC')
        call jedup1(ch1//'.VALE', base, ch2//'.VALE')
!
!
!     -- CAS DES CARTES :
!     ----------------------
    else if (docu .eq. 'CART') then
        call jedup1(ch1//'.DESC', base, ch2//'.DESC')
        call jedup1(ch1//'.LIMA', base, ch2//'.LIMA')
        call jedup1(ch1//'.NOLI', base, ch2//'.NOLI')
        call jedup1(ch1//'.NOMA', base, ch2//'.NOMA')
        call jedup1(ch1//'.VALE', base, ch2//'.VALE')
!
!
!     -- CAS DES CHAM_ELEM :
!     ----------------------
    else if (docu .eq. 'CHML') then
        call jedup1(ch1//'.CELD', base, ch2//'.CELD')
        call jedup1(ch1//'.CELK', base, ch2//'.CELK')
        call jedup1(ch1//'.CELV', base, ch2//'.CELV')
!       dupliquer le LIGREL
        if (base .eq. 'G') then
            call getres(nomu, concep, cmd)
            call dismoi('NOM_LIGREL', ch2, 'CHAM_ELEM', repk=ligr)
            if (ligr(1:8) .ne. nomu .and. ligr(9:15) .ne. '.MODELE') then
                noojb = '12345678.LIGR000000.LIEL'
                call gnomsd(' ', noojb, 14, 19)
                ligr2 = noojb(1:19)
                call copisd('LIGREL', base, ligr, ligr2)
                call jeveuo(ch2//'.CELK', 'E', vk24=refe)
                refe(1) = ligr2
            end if
        end if
!
!     -- CAS DES RESUELEM :
!     ----------------------
    else if (docu .eq. 'RESL') then
        call jedup1(ch1//'.DESC', base, ch2//'.DESC')
        call jedup1(ch1//'.NOLI', base, ch2//'.NOLI')
        call jedup1(ch1//'.RESL', base, ch2//'.RESL')
!
!
    else
        print *, ch1, ", ", docu
        call utmess('F', 'CALCULEL_17')
    end if
!
999 continue
    call jedema()
end subroutine
