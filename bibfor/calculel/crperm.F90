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

subroutine crperm()
    implicit none
!
!     COMMANDE:  CREA_RESU
!     TRAITEMENT DU MOT CLE FACTEUR "PERM_CHAM"
!
! ----------------------------------------------------------------------
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterfort/celces.h"
#include "asterfort/cescel.h"
#include "asterfort/cetran.h"
#include "asterfort/cnocns.h"
#include "asterfort/cnscno.h"
#include "asterfort/cntran.h"
#include "asterfort/crpcvg.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvem.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsorac.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
!
!
    integer(kind=8) :: n1, nbcham, iord1(1), iord2, nbperm, nbtrou, ip, ibid, ic
    integer(kind=8) :: iret, jlim1, jlim2, nbma, jlino, nbno2, nncp
    real(kind=8) :: inst1, tran(3), prec
    real(kind=8) :: valr
    complex(kind=8) :: cbid
    character(len=8) :: k8b, crit, resu1, resu2, resu3, ma1, ma2
    character(len=16) :: typres, nomcmd, cham(4), option
    character(len=24) :: valk(2)
    character(len=19) :: numeq
    character(len=24) :: ch1, ch2, chs1, chs2, linoeu, gma1, gma2, lima1, lima2
    character(len=24) :: ligrel, chsi1(4), chsi2(4)
    integer(kind=8), pointer :: ordr(:) => null()
! DEB ------------------------------------------------------------------
    call jemarq()
!
    call getres(resu3, typres, nomcmd)
!
! --- RECUPERATION DES DONNEES UTILISATEUR :
!     ------------------------------------
!
    call getvid(' ', 'RESU_INIT', scal=resu1, nbret=n1)
    call getvr8(' ', 'INST_INIT', scal=inst1, nbret=n1)
    if (n1 .eq. 0) then
        call jelira(resu1//'           .ORDR', 'LONUTI', ibid)
        call jeveuo(resu1//'           .ORDR', 'L', vi=ordr)
        iord1(1) = ordr(ibid)
    else
        call getvr8(' ', 'PRECISION', scal=prec, nbret=n1)
        call getvtx(' ', 'CRITERE', scal=crit, nbret=n1)
        call rsorac(resu1, 'INST', ibid, inst1, k8b, &
                    cbid, prec, crit, iord1, 1, &
                    nbtrou)
        if (nbtrou .eq. 0) then
            valr = inst1
            valk(1) = resu1
            call utmess('F', 'CALCULEL5_70', sk=valk(1), sr=valr)
        else if (nbtrou .ne. 1) then
            valr = inst1
            call utmess('F', 'CALCULEL5_71', sr=valr)
        end if
    end if
    call getvid(' ', 'MAILLAGE_INIT', scal=ma1, nbret=n1)
    call getvid(' ', 'RESU_FINAL', scal=resu2, nbret=n1)
    call getvid(' ', 'MAILLAGE_FINAL', scal=ma2, nbret=n1)
    call getvtx(' ', 'NOM_CHAM', nbval=0, nbret=n1)
    if (n1 .eq. 0) then
        nbcham = 4
        cham(1) = 'DEPL'
        cham(2) = 'SIEF_ELGA'
        cham(3) = 'VARI_ELGA'
        cham(4) = 'STRX_ELGA'
    else
        nbcham = -n1
        call getvtx(' ', 'NOM_CHAM', nbval=nbcham, vect=cham, nbret=n1)
!
    end if
!
    call dismoi('NB_NO_MAILLA', ma2, 'MAILLAGE', repi=nbno2)
    iord2 = 1
!
! --- VERIFICATIONS SUPPLEMENTAIRES :
!     -----------------------------
!
    if (resu2 .ne. resu3) then
        valk(1) = resu3
        valk(2) = resu2
        call utmess('F', 'CALCULEL5_72', nk=2, valk=valk)
    end if
!
    call jelira(resu2//'           .ORDR', 'LONUTI', ibid)
    if (ibid .ne. 1) then
        valk(1) = resu2
        valk(2) = k8b
        call utmess('F', 'CALCULEL5_73', nk=2, valk=valk)
    end if
!
    do ic = 1, nbcham
        call rsexch('F', resu1, cham(ic), iord1(1), ch1, &
                    iret)
        call rsexch('F', resu2, cham(ic), iord2, ch2, &
                    iret)
!
        if (cham(ic) .eq. 'DEPL') then
            chs1 = '&&CRPERM.DEPL_1'
            call cnocns(ch1, 'V', chs1)
            chsi1(ic) = chs1
            chs2 = '&&CRPERM.DEPL_2'
            call cnocns(ch2, 'V', chs2)
            chsi2(ic) = chs2
        else if (cham(ic) .eq. 'SIEF_ELGA') then
            chs1 = '&&CRPERM.SIEF_1'
            call celces(ch1, 'V', chs1)
            chsi1(ic) = chs1
            chs2 = '&&CRPERM.SIEF_2'
            call celces(ch2, 'V', chs2)
            chsi2(ic) = chs2
        else if (cham(ic) .eq. 'VARI_ELGA') then
            chs1 = '&&CRPERM.VARI_1'
            call celces(ch1, 'V', chs1)
            chsi1(ic) = chs1
            chs2 = '&&CRPERM.VARI_2'
            call celces(ch2, 'V', chs2)
            chsi2(ic) = chs2
        else if (cham(ic) .eq. 'STRX_ELGA') then
            chs1 = '&&CRPERM.STRX_1'
            call celces(ch1, 'V', chs1)
            chsi1(ic) = chs1
            chs2 = '&&CRPERM.STRX_2'
            call celces(ch2, 'V', chs2)
            chsi2(ic) = chs2
        end if
!
    end do
!
!
    linoeu = '&&CRPERM.LISTE_NOEU'
    lima1 = '&&CRPERM.LISTE_MA_1'
    lima2 = '&&CRPERM.LISTE_MA_2'
!
    call getfac('PERM_CHAM', nbperm)
!
! --- BOUCLE SUR LES TRANSLATIONS A EFFECTUER :
!     ---------------------------------------
!
    do ip = 1, nbperm
!
        call getvem(ma1, 'GROUP_MA', 'PERM_CHAM', 'GROUP_MA_INIT', ip, &
                    1, gma1, n1)
        call getvem(ma2, 'GROUP_MA', 'PERM_CHAM', 'GROUP_MA_FINAL', ip, &
                    1, gma2, n1)
!
        call getvr8('PERM_CHAM', 'TRAN', iocc=ip, nbval=3, vect=tran, &
                    nbret=n1)
        call getvr8('PERM_CHAM', 'PRECISION', iocc=ip, scal=prec, nbret=n1)
!
! ------ VERIFICATION DES GROUPES DE MAILLES FOURNIES :
!        --------------------------------------------
!
        call wkvect(linoeu, 'V V I', nbno2, jlino)
!
        call crpcvg(ma1, ma2, gma1, gma2, tran, &
                    prec, lima1, lima2, zi(jlino))
!
        call jelira(lima1, 'LONMAX', nbma)
        call jeveuo(lima1, 'L', jlim1)
        call jeveuo(lima2, 'L', jlim2)
!
        do ic = 1, nbcham
!
            chs1 = chsi1(ic)
            chs2 = chsi2(ic)
!
! --------- ON TRANSFERE LES VALEURS DE 1 VERS 2 :
!           ------------------------------------
!
            if (cham(ic) .eq. 'DEPL') then
                call cntran(zi(jlino), nbno2, chs1, chs2)
            else
                call cetran(zi(jlim1), zi(jlim2), nbma, chs1, chs2)
!
            end if
!
        end do
!
        call jedetr(lima1)
        call jedetr(lima2)
        call jedetr(linoeu)
!
    end do
!
    do ic = 1, nbcham
        call rsexch('F', resu2, cham(ic), iord2, ch2, &
                    iret)
        chs1 = chsi1(ic)
        chs2 = chsi2(ic)
        if (cham(ic) .eq. 'DEPL') then
            call dismoi('NUME_EQUA', ch2, 'CHAM_NO', repk=numeq)
            call cnscno(chs2, numeq, 'NON', 'G', ch2, &
                        'F', ibid)
            call detrsd('CHAM_NO_S', chs1)
            call detrsd('CHAM_NO_S', chs2)
        else
            call dismoi('NOM_LIGREL', ch2, 'CHAM_ELEM', repk=ligrel)
            call dismoi('NOM_OPTION', ch2, 'CHAM_ELEM', repk=option)
            call cescel(chs2, ligrel, option, ' ', 'OUI', &
                        nncp, 'G', ch2, 'F', ibid)
            call detrsd('CHAM_ELEM_S', chs1)
            call detrsd('CHAM_ELEM_S', chs2)
        end if
    end do
!
    call jedema()
end subroutine
