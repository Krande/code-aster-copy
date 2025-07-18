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

subroutine ccvepo(modele, resuin, typesd, lisord, nbordr, &
                  option, &
                  nbchre, ioccur, suropt, ligrel, exipou)
    implicit none
!     --- ARGUMENTS ---
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getexm.h"
#include "asterfort/assert.h"
#include "asterfort/cochre.h"
#include "asterfort/copich.h"
#include "asterfort/dismoi.h"
#include "asterfort/exlim1.h"
#include "asterfort/exlima.h"
#include "asterfort/gnomsd.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/medom1.h"
#include "asterfort/utmamo.h"
#include "asterfort/lisnch.h"
#include "asterfort/medome_once.h"
#include "asterfort/utmess.h"
    aster_logical :: exipou
    integer(kind=8) :: nbchre, ioccur
    character(len=8) :: modele, resuin
    character(len=16) :: typesd, option
    character(len=24) :: suropt, ligrel
    integer(kind=8) :: nbordr
    character(len=19) :: lisord
!  CALC_CHAMP - VERIFICATION POUR LES POUTRES
!  -    -       --                    --
! ----------------------------------------------------------------------
!
!  ROUTINE PERMETTANT DE SAVOIR SI DES POUTRES SONT DANS LE LIGREL
!   REDUIT ET DE VERIFIER LES CHARGES REPARTIES
!
! IN  :
!   MODELE  K8   NOM DU MODELE
!   RESUIN  K8   NOM DE LA STRUCTURE DE DONNEES RESULTAT IN
!   TYPESD  K16  TYPE DE LA STRUCTURE DE DONNEES RESULTAT
!
! OUT :
!   NBCHRE  I    NOMBRE DE CHARGES REPARTIES (POUTRES)
!   IOCCUR  I    NUMERO D'OCCURENCE OU SE TROUVE LE CHARGE REPARTIE
!   SUROPT  K24  SUROPTION
!   LIGREL  K24  NOM DU LIGREL A CREER
!   EXIPOU  L    LOGIQUE INDIQUANT LE PRESENCE DE POUTRES
! ----------------------------------------------------------------------
!
    integer(kind=8) :: ierd, ltymo, nbmaal
    integer(kind=8) :: n1, n2
!
    character(len=8) :: k8b, model, cara_elem
    character(len=24) :: mater, mateco
    character(len=16) :: typemo
    character(len=19) :: refe, masse, chdynr, chdepl
    character(len=24) :: noojb
    integer(kind=8) :: nbchar
    character(len=19) :: lischa
    integer(kind=8), pointer :: liste_mailles(:) => null()
    integer(kind=8), pointer :: v_list_store(:) => null()
    character(len=8), pointer :: lcha(:) => null()
!
    call jemarq()
!
    typemo = ' '
    suropt = ' '
    nbchar = 0
    ioccur = 0
    nbchre = 0
    lischa = '&&CCVEPO.LISCHA'
    if (typesd .eq. 'MODE_MECA') then
        call rsadpa(resuin, 'L', 1, 'TYPE_MODE', 1, &
                    0, sjv=ltymo)
        typemo = zk16(ltymo)
    end if
!
!     REDUCTION DU LIGREL SUR LA BASE DE GROUP_MA, GROUP_NO, ETC.
    n1 = getexm(' ', 'GROUP_MA')
    n2 = getexm(' ', 'MAILLE')
    if (n1+n2 .ne. 0) then
        call exlima(' ', 0, 'V', modele, ligrel)
    else
        call utmamo(modele, nbmaal, '&&CCVEPO.LISTE_MAILLES')
        call jeveuo('&&CCVEPO.LISTE_MAILLES', 'L', vi=liste_mailles)
        noojb = '12345678.LIGR000000.LIEL'
        call gnomsd(' ', noojb, 14, 19)
        ligrel = noojb(1:19)
        ASSERT(ligrel .ne. ' ')
        call exlim1(liste_mailles, nbmaal, modele, 'V', ligrel)
        call jedetr('&&CCVEPO.LISTE_MAILLES')
    end if
!
    call dismoi('EXI_POUX', ligrel, 'LIGREL', repk=k8b)
    exipou = .false.
!     SPECIAL POUTRE A LA POUX...
    if (k8b(1:3) .eq. 'OUI') then
        exipou = .true.
!       ON VERIFIE SI DERIERE LE TYPESD MODE_MECA ON TROUVE UN MODE_DYN
        if ((typesd .eq. 'MODE_MECA' .and. typemo(1:8) .eq. 'MODE_DYN') .or. typesd .eq. &
            'DYNA_TRANS' .or. typesd .eq. 'MODE_ACOU' .or. typesd .eq. 'DYNA_HARMO') then
            refe = resuin
            suropt = 'MASS_MECA'
            call dismoi('REF_MASS_PREM', refe, 'RESU_DYNA', repk=masse, arret='C', &
                        ier=ierd)
            if (masse .ne. ' ') then
                call dismoi('SUR_OPTION', masse, 'MATR_ASSE', repk=suropt, arret='C', &
                            ier=ierd)
            end if
            call rsexch('F', resuin, 'DEPL', 1, chdepl, &
                        ierd)
            chdynr = '&&MECALM.M.GAMMA'
            call copich('V', chdepl(1:19), chdynr)
        end if
        call jeveuo(lisord, 'L', vi=v_list_store)
        if (option .eq. 'REAC_NODA') then
            call medome_once(resuin, v_list_store, nbordr, &
                             list_load_=lischa)
            call lisnch(lischa, nbchar)
        else
            call medom1(model, mater, mateco, cara_elem, lischa, nbchar, &
                        resuin, v_list_store(1))
        end if
!       VERIFIE L'UNICITE DE LA CHARGE REPARTIE
        if (nbchar .ne. 0) then
            call jeveuo(lischa//'.LCHA', 'L', vk8=lcha)
            ioccur = 0
            call cochre(lcha, nbchar, nbchre, ioccur)
            if (nbchre .gt. 1) then
                call utmess('F', 'CALCULEL2_92')
            end if
        end if
    end if
!
    call jedema()
!
end subroutine
