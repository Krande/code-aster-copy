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

subroutine modirepcham(resuou, resuin)
!
!
! --------------------------------------------------------------------------------------------------
!
!     COMMANDE : MODI_REPERE / CHAM_GD
!
!   in
!       resuin  : Nom du champ en entrée
!   out
!       resuou  : Nom du champ en sortie
! --------------------------------------------------------------------------------------------------
!
    implicit none
    character(len=19) :: resuou, resuin
!
#include "jeveux.h"
#include "asterfort/calcul.h"
#include "asterfort/cesvar.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/megeom.h"
#include "asterfort/utmess.h"
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv, nret, iret
    character(len=8) :: maillage, modele, carelem, caramail, caramodel
    character(len=16) :: repere
    character(len=19) :: chpass
    character(len=24) :: ligrel, option
!   Pour calcul
    character(len=8) :: lpain(4), lpaou(4)
    character(len=24) :: lchin(4), lchou(4), chgeom
!   Pour les messages
!     integer ::  vali(2)
    character(len=80) :: valk(3)
!
    aster_logical :: lreuse
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    call infmaj()
    call infniv(ifm, niv)

    lreuse = .false.
    if (resuin .eq. resuou) then
        lreuse = .true.
    end if
!
!   Définition du repère utilisé
    call getvtx(' ', 'REPERE', scal=repere, nbret=nret)
    if (nret .eq. 0 .and. .not. lreuse) then
        call utmess('F', 'MODELISA3_2')
    else if (repere .ne. 'GLOBAL_UTIL' .and. .not. lreuse) then
        call utmess('F', 'MODELISA3_3')
    end if
!   Lecture du concept CARA_ELEM
    call getvid(' ', 'CARA_ELEM', scal=carelem, nbret=nret)
    if (nret .eq. 0 .and. .not. lreuse) then
        call utmess('F', 'MODELISA3_7')
    end if
!
!   Informations sur le champ en entrée.
!    call dismoi('TYPE_CHAMP', resuin, 'CHAMP', repk=tychamp)
!    call dismoi('NOM_PARAM',  resuin, 'CHAMP', repk=nompar)
    call dismoi('NOM_OPTION', resuin, 'CHAMP', repk=option)
    call dismoi('NOM_LIGREL', resuin, 'CHAMP', repk=ligrel)
    call dismoi('NOM_MAILLA', resuin, 'CHAMP', repk=maillage)
    call dismoi('NOM_MODELE', resuin, 'CHAMP', repk=modele)
!   Le champ doit être crée avec INI_SP_RIGI
    if (option .ne. 'INI_SP_RIGI') then
        call utmess('F', 'MODELISA3_1')
    end if
!
! --------------------------------------------------------------------------------------------------
!   Vérification que CARCOQUE existe
    call exisd('CARTE', carelem//'.CARCOQUE', iret)
    if (iret .eq. 0) then
        valk(1) = carelem
        valk(2) = 'EP, ALPHA, BETA'
        call utmess('F', 'MODELISA3_4', nk=2, valk=valk)
    end if
!   Vérification que CANBSP existe
    call exisd('CHAM_ELEM', carelem//'.CANBSP', iret)
    if (iret .eq. 0) then
        valk(1) = carelem
        valk(2) = 'COQ_NCOU'
        call utmess('F', 'MODELISA3_4', nk=2, valk=valk)
    end if
!   Nom du maillage sous-jacent à la carte. Le même que celui du champ.
    call dismoi('NOM_MAILLA', carelem, 'CARA_ELEM', repk=caramail)
    if (maillage .ne. caramail) then
        valk(1) = maillage
        valk(2) = caramail
        call utmess('F', 'MODELISA3_5', nk=2, valk=valk)
    end if
!   Nom du modèle sous-jacent à CARA_ELEM. Le même que celui du champ.
    call dismoi('NOM_MODELE', carelem, 'CARA_ELEM', repk=caramodel)
    if (modele .ne. caramodel) then
        valk(1) = modele
        valk(2) = caramodel
        call utmess('F', 'MODELISA3_6', nk=2, valk=valk)
    end if
!
! --------------------------------------------------------------------------------------------------
!   Matrice de passage du repère global vers le repère utilisateur
    chpass = '&&REPCHA.MATPASS'
!
    call megeom(modele, chgeom)
    lchin(1) = chgeom
    lpain(1) = 'PGEOMER'
    lchin(2) = carelem//'.CARCOQUE'
    lpain(2) = 'PCACOQU'
!
    lchou(1) = chpass
    lpaou(1) = 'PMATPASS'
!
    call calcul('C', 'REPERE_LOCAL', ligrel, 2, lchin, &
                lpain, 1, lchou, lpaou, 'V', 'NON')
! --------------------------------------------------------------------------------------------------
!   Changement de repère
    lchin(1) = chpass
    lpain(1) = 'PMATPASS'
    lchin(2) = resuin(1:8)
    lpain(2) = 'PSIEFR'
!
    lchou(1) = resuou
    lpaou(1) = 'PCONTPR'
!
    call cesvar(carelem, ' ', ligrel, lchou(1))
    call calcul('C', 'MODI_REPERE', ligrel, 2, lchin, &
                lpain, 1, lchou, lpaou, 'G', 'NON')
!
! --------------------------------------------------------------------------------------------------
!   Ménage
    call jedetr(chpass)
!
    call jedema()
end subroutine
