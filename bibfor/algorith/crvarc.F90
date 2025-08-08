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

subroutine crvarc()
    implicit none
!
! ------------------------------------------------------------------------------
!
!                           COMMANDE    :   CREA_RESU
!
! ------------------------------------------------------------------------------
!
!       Création d'une SD de type "EVOL_THER" contenant la température sur les
!       couches des coques multicouche à partir :
!           d'un CHAM_GD de fonctions [ INST, [EPAIS|X,Y,Z] ]
!           d'un EVOL_THER contenant TEMP/TEMP_INF/TEMP_SUP
!
! ------------------------------------------------------------------------------
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/cesvar.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exlima.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mecact.h"
#include "asterfort/mecara.h"
#include "asterfort/megeom.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rscrsd.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnoch.h"
#include "asterfort/rsorac.h"
#include "asterfort/utmess.h"
!
! ------------------------------------------------------------------------------
!
    integer(kind=8)         :: nbfac, jinst, nbin, iret, iexi, jcara
    integer(kind=8)         :: n1, n2, inst, nbinst, ordr, nbordr
    real(kind=8)    :: vinst, prec
!
    character(len=5)    :: LesInst, LeCas
    character(len=8)    :: resu, carele, paout, lpain(10), tempef
    character(len=8)    :: model2, modele
    character(len=16)   :: typeres, oper, crit
    character(len=19)   :: ligrel, chout, chinst, listr, temper
    character(len=24)   :: chcara(18), lchin(10), chtemp, chgeom
!
    integer(kind=8)             :: vali, lesordres(1)
    real(kind=8)        :: valr
    character(len=24)   :: valk
!
    integer(kind=8)             :: ibib
    character(len=8)    :: k8bid
    complex(kind=8)     :: c16bid
!
    integer(kind=8), pointer :: lordr(:) => null()
    real(kind=8), pointer :: linst(:) => null()
    character(len=24), pointer :: celk(:) => null()
!
! ------------------------------------------------------------------------------
    call jemarq()
!
    call getfac('PREP_VARC', nbfac)
    if (nbfac .eq. 0) goto 20
    ASSERT(nbfac .eq. 1)
!
    call getres(resu, typeres, oper)
!   C'est CHAM_GD ou EVOL_THER
    LeCas = ' '
    call getvid('PREP_VARC', 'CHAM_GD', iocc=1, nbval=0, nbret=n1)
    call getvid('PREP_VARC', 'EVOL_THER', iocc=1, nbval=0, nbret=n2)
    if (n1 .lt. 0) then
        call getvid('PREP_VARC', 'CHAM_GD', iocc=1, scal=tempef, nbret=n1)
        LeCas = 'CHAMP'
    else if (n2 .lt. 0) then
        call getvid('PREP_VARC', 'EVOL_THER', iocc=1, scal=temper, nbret=n2)
        LeCas = 'THERM'
    else
        ASSERT(.FALSE.)
    end if
!
    call getvid('PREP_VARC', 'MODELE', iocc=1, scal=modele, nbret=n1)
    call getvid('PREP_VARC', 'CARA_ELEM', iocc=1, scal=carele, nbret=n1)
!
!   On vérifie que le cara_elem s'appuie bien sur le modèle
    call jeexin(carele//'.CANBSP    .CELK', iexi)
    if (iexi .eq. 0) then
        call utmess('F', 'CALCULEL4_14', sk=carele)
    end if
    call jeveuo(carele//'.CANBSP    .CELK', 'L', vk24=celk)
    model2 = celk(1) (1:8)
    if (model2 .ne. modele) then
        call utmess('F', 'CALCULEL4_15', sk=carele)
    end if
!
!   Les instants
    LesInst = ' '
    call getvr8('PREP_VARC', 'INST', iocc=1, nbval=0, nbret=n1)
    call getvr8('PREP_VARC', 'LIST_INST', iocc=1, nbval=0, nbret=n2)
    if (n1 .lt. 0) then
        nbinst = -n1
        AS_ALLOCATE(vr=linst, size=nbinst)
        call getvr8('PREP_VARC', 'INST', iocc=1, nbval=nbinst, vect=linst, nbret=n1)
        ASSERT(nbinst .eq. n1)
        LesInst = 'LINST'
    else if (n2 .lt. 0) then
        call getvid('PREP_VARC', 'LIST_INST', iocc=1, scal=listr, nbret=n2)
        ASSERT(n2 .eq. 1)
        call jeveuo(listr//'.VALE', 'L', jinst)
        call jelira(listr//'.VALE', 'LONMAX', nbinst)
        ! Dimensionnement de la liste des instants
        AS_ALLOCATE(vr=linst, size=nbinst)
        do inst = 1, nbinst
            linst(inst) = zr(jinst+inst-1)
        end do
        LesInst = 'LINST'
    end if
!
!   Si EVOL_THER
    prec = 0.0; crit = ' '
    if (LeCas .eq. 'THERM') then
        ! Si INST ou LIST_INST : pour recherche de INST dans la SD
        if (LesInst .eq. 'LINST') then
            call getvr8('PREP_VARC', 'PRECISION', iocc=1, scal=prec, nbret=n1)
            call getvtx('PREP_VARC', 'CRITERE', iocc=1, scal=crit, nbret=n2)
        else
            LesInst = 'TINST'
        end if
    end if
!
!   Si EVOL_THER, on vérifie qu'il y a des NUME_ORDRE
    if (LeCas .eq. 'THERM') then
        call jelira(temper//'.ORDR', 'LONUTI', nbordr)
        call jeveuo(temper//'.ORDR', 'L', vi=lordr)
        ASSERT(nbordr .gt. 0)
        if (LesInst .eq. 'TINST') then
            nbinst = nbordr
        end if
    end if
!
!   Champ en sortie
    call jeexin(resu//'           .DESC', iret)
    if (iret .ne. 0) then
        call utmess('F', 'CALCULEL7_6', sk=resu)
    else
        call rscrsd('G', resu, 'EVOL_THER', nbinst)
    end if
!
!   Paramètres pour CALCUL
    call mecara(carele, chcara)
    paout = 'PTEMPCR'
    lpain(1) = 'PNBSP_I'
    lchin(1) = chcara(1) (1:8)//'.CANBSP'
    lpain(2) = 'PCACOQU'
    lchin(2) = chcara(7)
!
    if (LeCas .eq. 'CHAMP') then
        call dismoi('NOM_LIGREL', modele, 'MODELE', repk=ligrel)
        chinst = '&&CRVRC1.CHINST'
        !
        lpain(3) = 'PTEMPEF'
        lchin(3) = tempef
        lpain(4) = 'PINST_R'
        lchin(4) = chinst
        ! On va chercher la géométrie
        call megeom(modele, chgeom)
        lpain(5) = 'PGEOMER'
        lchin(5) = chgeom
        nbin = 5
    else if (LeCas .eq. 'THERM') then
        ! Examine : GROUP_MA, TOU, MAILLE, ....
        call exlima('PREP_VARC', 1, 'G', modele, ligrel)
        !
        lpain(3) = 'PTEMPER'
        nbin = 3
    end if
!
    if (LeCas .eq. 'CHAMP') then
        ! Boucle sur les instants
        do inst = 1, nbinst
            vinst = linst(inst)
            call mecact('V', chinst, 'MODELE', ligrel, 'INST_R', &
                        ncmp=1, nomcmp='INST', sr=vinst)
            call rsexch(' ', resu, 'TEMP', inst, chout, iret)
            !
            call cesvar(carele, ' ', ligrel, chout)
            call calcul('S', 'PREP_VRC', ligrel, nbin, lchin, &
                        lpain, 1, chout, paout, 'G', 'OUI')
            call detrsd('CHAM_ELEM_S', chout)
            call rsnoch(resu, 'TEMP', inst)
            ! On enregistre dans "resu" : l'instant
            call rsadpa(resu, 'E', 1, 'INST', inst, 0, sjv=jinst)
            zr(jinst) = vinst
            ! On enregistre dans "resu" : le cara_elem
            call rsadpa(resu, 'E', 1, 'CARAELEM', inst, 0, sjv=jcara)
            zk8(jcara) = carele
            !
            call detrsd('CHAMP', chinst)
        end do
    else if (LeCas .eq. 'THERM') then
        ! Boucle sur les instants
        do inst = 1, nbinst
            if (LesInst .eq. 'TINST') then
                ! Si TOUT_INST, c'est pareil que TOUT_ORDRE
                ordr = lordr(inst)
                ! On va chercher l'instant correspondant
                call rsadpa(temper, 'L', 1, 'INST', ordr, 0, sjv=jinst)
                vinst = zr(jinst)
            else
                ! Si INST ou LIST_INST, on cherche le NUME_ORDRE dans la SD
                vinst = linst(inst)
                call rsorac(temper, 'INST', ibib, vinst, k8bid, c16bid, &
                            prec, crit, lesordres, 1, n1)
                if (n1 .lt. 0) then
                    valk = temper
                    valr = vinst
                    vali = -n1
                    call utmess('F', 'ALGORITH12_83', sk=valk, si=vali, sr=valr)
                else if (n1 .eq. 0) then
                    valk = temper
                    valr = vinst
                    call utmess('F', 'ALGORITH12_84', sk=valk, sr=valr)
                end if
                ordr = lesordres(1)
            end if
            call rsexch('F', temper, 'TEMP', ordr, chtemp, iret)
            lchin(3) = chtemp
            call rsexch(' ', resu, 'TEMP', ordr, chout, iret)
            !
            call cesvar(carele, ' ', ligrel, chout)
            call calcul('S', 'PREP_VRC', ligrel, nbin, lchin, &
                        lpain, 1, chout, paout, 'G', 'OUI')
            call detrsd('CHAM_ELEM_S', chout)
            call rsnoch(resu, 'TEMP', ordr)
            ! On enregistre dans "resu" : l'instant
            call rsadpa(resu, 'E', 1, 'INST', ordr, 0, sjv=jinst)
            zr(jinst) = vinst
            ! On enregistre dans "resu" : le cara_elem
            call rsadpa(resu, 'E', 1, 'CARAELEM', ordr, 0, sjv=jcara)
            zk8(jcara) = carele
        end do
    end if
!
    AS_DEALLOCATE(vr=linst)
!
20  continue
    call jedema()
end subroutine
