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

subroutine op0195()
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/cheksd.h"
#include "asterc/getres.h"
#include "asterc/getfac.h"
#include "asterfort/alchml.h"
#include "asterfort/assert.h"
#include "asterfort/caraff.h"
#include "asterfort/celver.h"
#include "asterfort/chcore.h"
#include "asterfort/chpass.h"
#include "asterfort/chpchd.h"
#include "asterfort/chpeva.h"
#include "asterfort/chprec.h"
#include "asterfort/chreco.h"
#include "asterfort/cnoaff.h"
#include "asterfort/cnocns.h"
#include "asterfort/cnonor.h"
#include "asterfort/cnscno.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/duplisp.h"
#include "asterfort/exisd.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/idensd.h"
#include "asterfort/imprsd.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nopar2.h"
#include "asterfort/titre.h"
#include "asterfort/u195tb.h"
#include "asterfort/utmess.h"
#include "asterfort/varaff.h"
#include "asterfort/x195cb.h"
#include "asterfort/xasdpl.h"
!
!
! --------------------------------------------------------------------------------------------------
!
!   COMMANDE CREA_CHAMP
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: n1, ifm, niv, iret, i11, i12, test, ibid, nocc
    character(len=3) :: prol0
    character(len=4) :: tychr, tych
    character(len=8) :: kbid, model, mesh, chou, nomgd, nomgd2, carel
    character(len=8) :: tsca, nogd, nomgd1, nompar, ma2, ta, ma3
    character(len=16) :: tychr1, opera, optio2, typco, option
    character(len=19) :: modelLigrel, chatmp, celmod, numeq1, cns1, ch1, numeq2, chin, chou2
    character(len=8) :: nu1
    character(len=24), pointer :: v_celmod_celk(:) => null()
    aster_logical :: dbg
    character(len=24) :: valk(4)
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infmaj()
    call infniv(ifm, niv)
!
    dbg = .true.
    dbg = .false.
!
!
! 1- CALCUL DE:
!      OPERA: OPERATION A EFFECTUER
!      MO: MODELE (OU ' ')
!      MA: MAILLAGE (OU ' ')
!      CHOU  : CHAMP RESULTAT
!      TYCHR : TYPE DU CHAMP RESULTAT (CART/NOEU/ELNO/ELGA/ELEM)
!      NOMGD : GRANDEUR ASSOCIEE A CHOU
!      PROL0 :/'OUI' POUR PROLONGER PAR ZERO LE CHAMP RESULTAT
!             (POUR LES CHAM_ELEM ET LES CHAM_NO POUR LESQUELS ON
!              VEUT IMPOSER LA NUMEROTATION DES DDLS).
!      OPTION: OPTION PERMETTANT D'ALLOUER UN CHAM_ELEM "MODELE"
!
!     ------------------------------------------------------------------
!
    call getvtx(' ', 'OPERATION', scal=opera)
!
    call getvid(' ', 'MODELE', scal=model, nbret=n1)
    if (n1 .eq. 0) model = ' '
    call getvid(' ', 'MAILLAGE', scal=mesh, nbret=n1)
    if (n1 .eq. 0) mesh = ' '
    if (model .ne. ' ') then
        call dismoi('NOM_MAILLA', model, 'MODELE', repk=ma2)
        if ((mesh .ne. ' ') .and. (mesh .ne. ma2)) then
            call utmess('F', 'UTILITAI3_21')
        end if
        mesh = ma2
    end if
!
    call getres(chou, typco, kbid)
    call exisd('CHAMP', chou, test)
    if (test .eq. 1) then
        if (.not. ((opera .eq. 'ASSE') .or. (opera .eq. 'COMB'))) then
            call utmess('F', 'UTILITAI3_43')
        end if
    end if
    call getvtx(' ', 'TYPE_CHAM', scal=tychr1)
    tychr = tychr1(1:4)
    nomgd = tychr1(6:13)
    call dismoi('TYPE_SCA', nomgd, 'GRANDEUR', repk=tsca)
!
    call getvtx(' ', 'PROL_ZERO', scal=prol0, nbret=iret)
    if (iret .eq. 0) then
        prol0 = ' '
    end if
    if ((prol0 .eq. 'NON') .and. (tsca .eq. 'R')) prol0 = 'NAN'
!
    call getvtx(' ', 'OPTION', scal=option, nbret=n1)
    if (n1 .eq. 0) option = ' '
!
! 2.  CALCUL DE LIGREL,CELMOD  + QUELQUES VERIFICATIONS   :
!     -------------------------------------------------------------
    modelLigrel = ' '
    celmod = ' '
!
    if (tychr(1:2) .eq. 'EL') then
        if ((opera .eq. 'AFFE') .or. (opera .eq. 'ASSE') .or. (opera .eq. 'ASSE_DEPL') .or. &
            (opera .eq. 'DISC')) then
            if (model .eq. ' ') then
                call utmess('F', 'UTILITAI3_22')
            end if
            call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
!
!         -- CALCUL D'UN CHAM_ELEM "MODELE" : CELMOD
!         ---------------------------------------------------
            if (option .eq. ' ') then
                optio2 = 'TOU_INI_'//tychr
!           -- SI OPERATION 'ASSE', IL Y A PEUT-ETRE UNE MEILLEURE
!              OPTION A CHOISIR PAR DEFAUT QUE TOU_INI_ELXX
                if (opera .eq. 'ASSE') then
                    call getvid('ASSE', 'CHAM_GD', iocc=1, scal=chin)
                    call dismoi('NOM_MAILLA', chin, 'CHAMP', repk=ma2)
                    if (ma2 .ne. mesh) then
                        valk(1) = chin
                        valk(2) = model
                        valk(3) = ma2
                        valk(4) = mesh
                        call utmess('F', 'CALCULEL4_59', nk=4, valk=valk)
                    end if
!
                    call jeexin(chin//'.CELK', iret)
                    call dismoi('NOM_GD', chin, 'CHAMP', repk=nomgd2)
                    if (iret .ne. 0 .and. nomgd .eq. nomgd2) then
                        call dismoi('NOM_OPTION', chin, 'CHAM_ELEM', repk=optio2)
                    end if
                end if
!
            else
                optio2 = option
            end if
            call nopar2(optio2, nomgd, 'INOUT', nompar, 0, iret)
            if (iret .ne. 0) then
                call utmess('F', 'CALCULEL3_86', sk=tychr1)
            end if
            celmod = '&&OP0195.CELMOD'
            call alchml(modelLigrel, optio2, nompar, 'V', celmod, &
                        iret, ' ')
            if (iret .ne. 0) then
                valk(1) = modelLigrel
                valk(2) = nompar
                valk(3) = optio2
                call utmess('E', 'UTILITAI3_23', nk=3, valk=valk)
                call utmess('F', 'CALCULEL3_87', sk=tychr1)
            end if
!
!         VERIFICATION DU TYPE DE CELMOD : ELGA/ELNO/ELEM :
            call jeveuo(celmod//'.CELK', 'L', vk24=v_celmod_celk)
            if (v_celmod_celk(3) .ne. tychr) then
                valk(1) = optio2
                valk(2) = tychr
                call utmess('F', 'UTILITAI3_24', nk=2, valk=valk)
            end if
        end if
    end if
!
!
!
! 3.  TRAITEMENT DU MOT CLE OPERATION :
!     -------------------------------------------------------------
!
!
    if (opera .eq. 'NORMALE') then
!     -----------------------------------------
        if (tychr .eq. 'NOEU') then
            if (nomgd .ne. 'GEOM_R') then
                call utmess('F', 'UTILITAI3_25', sk=opera)
            end if
            call cnonor(model, nomgd, 'G', chou)
!
        else
            valk(1) = opera
            valk(2) = tychr
            call utmess('F', 'UTILITAI3_26', nk=2, valk=valk)
        end if
!
!
!
    else if (opera .eq. 'AFFE') then
!     -----------------------------------------
        if (tychr .eq. 'NOEU') then
            call cnoaff(mesh, nomgd, 'G', chou)
!
        else if (tychr .eq. 'CART') then
            call caraff(mesh, nomgd, 'G', chou)
!
        else if (tychr(1:2) .eq. 'EL') then
            chatmp = '&&OP0195.CHATMP'
            if (nomgd .eq. 'VARI_R') then
                call varaff(mesh, nomgd, 'V', chatmp)
                call chpchd(chatmp, tychr, celmod, prol0, 'G', chou, model)
                call detrsd('CHAM_ELEM_S', chatmp)
            else
                call caraff(mesh, nomgd, 'V', chatmp)
                call chpchd(chatmp, tychr, celmod, prol0, 'G', chou, model)
                call detrsd('CARTE', chatmp)
            end if
        end if
!
!
    else if (opera .eq. 'ASSE') then
!     -----------------------------------------
        call chpass(tychr, mesh, celmod, nomgd, prol0, &
                    chou)
!
!
    else if (opera .eq. 'COMB') then
!     -----------------------------------------
        call x195cb(tychr, nomgd, chou)
!
!
    else if (opera .eq. 'EVAL') then
!     -----------------------------------------
        call chpeva(chou)
!
!
    else if (opera .eq. 'ASSE_DEPL') then
!     -----------------------------------------
        call xasdpl(model, celmod, prol0, chou)
!
!
    else if (opera(1:3) .eq. 'R2C') then
!     -----------------------------------------
        call chcore(chou)
!
!
    else if (opera(1:3) .eq. 'C2R') then
!     -----------------------------------------
        call chreco(chou)
!
!
    else if (opera .eq. 'DISC') then
!     -----------------------------------------
        call getvid(' ', 'CHAM_GD', scal=chin)
        call dismoi('NOM_GD', chin, 'CHAMP', repk=nomgd2)
        if (nomgd .ne. nomgd2) then
!          -- EXCEPTION : NOMGD='VARI_R' ET NOMGD2='VAR2_R'
            if (nomgd .eq. 'VARI_R' .and. nomgd2 .eq. 'VAR2_R') then
            elseif (nomgd .eq. 'VAR2_R' .and. nomgd2 .eq. 'VARI_R') then
            else
                valk(1) = chin
                valk(2) = tychr1
                call utmess('F', 'UTILITAI3_27', nk=2, valk=valk)
            end if
        end if
        call chpchd(chin, tychr, celmod, prol0, 'G', chou, model)
!
!
    else if (opera .eq. 'EXTR') then
!     -----------------------------------------
        call getvid(' ', 'TABLE', scal=ta, nbret=n1)
        if (n1 .eq. 0) then
            call chprec(chou)
!
        else
            call u195tb(chou)
        end if
    end if
!
!
! 4.1  SI ON A CREE UN CHAM_NO, ON PEUT IMPOSER SA NUMEROTATION :
! --------------------------------------------------------------
    if (tychr .eq. 'NOEU') then
        call getvid(' ', 'CHAM_NO', scal=ch1, nbret=i11)
        call getvid(' ', 'NUME_DDL', scal=nu1, nbret=i12)
        if (i12 .eq. 1) then
            call dismoi('NOM_MAILLA', nu1, 'NUME_DDL', repk=ma3)
            if (mesh .ne. ' ') then
                if (ma3 .ne. mesh) then
                    valk(1) = nu1
                    valk(2) = ma3
                    valk(3) = mesh
                    call utmess('F', 'CALCULEL4_69', nk=3, valk=valk)
                end if
            end if
        end if
        if ((i11+i12) .gt. 0) then
            call dismoi('NOM_GD', chou, 'CHAMP', repk=nogd)
            numeq1 = ' '
            if (i11 .gt. 0) then
                call dismoi('NUME_EQUA', ch1, 'CHAM_NO', repk=numeq1)
                call dismoi('NOM_GD', ch1, 'CHAM_NO', repk=nomgd1)
            end if
            if (i12 .gt. 0) then
                call dismoi('NUME_EQUA', nu1, 'NUME_DDL', repk=numeq1)
                call dismoi('NOM_GD', nu1, 'NUME_DDL', repk=nomgd1)
            end if
!
            if (nomgd1 .ne. nogd) then
!           -- ON ACCEPTE LES COUPLES XXXX_R / XXXX_C :
                if (nomgd1(5:7) .eq. '_R ' .or. nomgd1(5:7) .eq. '_C ') then
                    if (nogd(5:7) .eq. '_R ' .or. nogd(5:7) .eq. '_C ') then
                        if (nogd(1:4) .eq. nomgd1(1:4)) goto 1
                    end if
                end if
                valk(1) = nomgd1
                valk(2) = nogd
                call utmess('F', 'UTILITAI6_5', nk=2, valk=valk)
            end if
1           continue
!
            call dismoi('NUME_EQUA', chou, 'CHAM_NO', repk=numeq2)
            if (.not. idensd('NUME_EQUA', numeq1, numeq2)) then
                call getfac('COMB', nocc)
                if (nocc .ne. 0) then
                    call utmess('A', 'CREACHAMP1_14')
                end if
                call getvtx(' ', 'PROL_ZERO', scal=prol0, nbret=iret)
                if (iret .eq. 0) then
                    prol0 = 'NON'
                end if
!               chacun son NUME_EQUA pour éviter les problèmes en cas de suppression
                cns1 = '&&OP0195.CNS1'
                call cnocns(chou, 'V', cns1)
                call detrsd('NUME_EQUA', numeq2)
!               si NUME_DDL pas de problèmes, donc on partage le nume_EQUA
                if (i12 .eq. 1) then
                    call cnscno(cns1, numeq1, prol0, 'G', chou, &
                                'F', ibid)
                else
                    call copisd('NUME_EQUA', 'G', numeq1, numeq2)
                    call cnscno(cns1, numeq2, prol0, 'G', chou, &
                                'F', ibid)
                end if
                call detrsd('CHAM_NO_S', cns1)
            end if
        end if
    end if
!
!
! 4.2  SI ON A CREE UN CHAM_ELEM, ON PEUT VOULOIR AFFECTER TOUS LES SOUS-POINTS :
! -------------------------------------------------------------------------------
    if (tychr(1:2) .eq. 'EL') then
        call getfac('AFFE_SP', nocc)
        if (nocc .eq. 1) then
            call getvid('AFFE_SP', 'CARA_ELEM', iocc=1, scal=carel)
            chou2 = '&&OP0195.CHOU2'
            call duplisp(chou, chou2, carel, 'V')
            call detrsd('CHAMP', chou)
            call copisd('CHAMP', 'G', chou2, chou)
            call detrsd('CHAMP', chou2)
        end if
    end if
!
!
! 5.  AJOUT DU TITRE :
! -----------------------------------------------------
    call titre()
!
!
! 6.  SI INFO:2    ON IMPRIME LE CHAMP RESULTAT :
! ----------------------------------------------
    if (niv .eq. 2) call imprsd('CHAMP', chou, ifm, 'CHAMP RESULTAT DE LA COMMANDE CREA_CHAMP :')
!
!
! 7.  VERIFICATION PROL_ZERO='NON' POUR LES CHAM_ELEM :
! ------------------------------------------------------
    if ((tychr(1:2) .eq. 'EL') .and. (prol0 .eq. 'NAN')) then
        call celver(chou, 'PAS_NAN', 'STOP', iret)
    end if
    if (dbg .and. tychr(1:2) .eq. 'EL') then
        call cheksd(chou, 'SD_CHAM_ELEM', iret)
        ASSERT(iret .eq. 0)
    end if
!
!
! 8.  VERIFICATION DE LA COHERENCE DU MAILLAGE SOUS-JACENT :
! ---------------------------------------------------------
    if (mesh .ne. ' ') then
        call dismoi('NOM_MAILLA', chou, 'CHAMP', repk=ma2)
        valk(1) = ma2
        valk(2) = mesh
        if (mesh .ne. ma2) then
            call utmess('F', 'CALCULEL4_78', nk=2, valk=valk)
        end if
    end if
!
!
! 9.  VERIFICATION DE LA COHERENCE DU CHAMP CREE AVEC TYPE_CHAM :
! ---------------------------------------------------------------
    call dismoi('TYPE_CHAMP', chou, 'CHAMP', repk=tych)
    if (tych .ne. tychr) then
        valk(1) = tychr
        valk(2) = tych
        call utmess('F', 'CALCULEL4_70', nk=2, valk=valk)
    end if
!
    call dismoi('NOM_GD', chou, 'CHAMP', repk=nogd)
    if (nogd .ne. nomgd) then
        valk(1) = nomgd
        valk(2) = nogd
        call utmess('F', 'CALCULEL4_71', nk=2, valk=valk)
    end if
!
!
    call jedema()
!
end subroutine
