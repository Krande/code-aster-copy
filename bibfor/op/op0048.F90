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

subroutine op0048()
    implicit none
#include "asterf_types.h"
#include "asterc/getres.h"
#include "asterc/indik8.h"
#include "asterfort/assert.h"
#include "asterfort/carces.h"
#include "asterfort/cesexi.h"
#include "asterfort/chpchd.h"
#include "asterfort/cnocns.h"
#include "asterfort/cnscno.h"
#include "asterfort/cnscre.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/vrcins.h"
! person_in_charge: sam.cuvilliez at edf.fr
!-----------------------------------------------------------------------
!
!     Operateur cache POST_K_VARC, appele depuis la macro POST_K1_K2_K3
!
!     Produit le cham_no de variables de commande pour une sd_resultat
!     donnee a l'instant considere. Le champ est type de la maniere
!     suivante, en fonction du nom de la variable de commande :
!         TEMP  -> TEMP_R
!         NEUT1 -> NEUT_R
!
!     Les seules varc autorisees sont : TEMP, NEUT1
!
!-----------------------------------------------------------------------
!
    integer(kind=8) :: nvrcmx
    parameter(nvrcmx=2)
    character(len=8) :: tab_nomvrc(nvrcmx), tab_nomgrd(nvrcmx)
    character(len=8) :: tab_nomcmp(nvrcmx)
    data tab_nomvrc/'TEMP', 'NEUT1'/
    data tab_nomgrd/'TEMP_R', 'NEUT_R'/
    data tab_nomcmp/'TEMP', 'X1'/
!
    character(len=6) :: nompro
    parameter(nompro='OP0048')
!
    integer(kind=8) :: ier, iret, ibid, nbno, ino, k, l, nbvarc, jcesd1, jcesl1
    integer(kind=8) :: nbma, ima, iad, itrou, nbcmp, itest, nucmp, cpt_dbl, ind
    real(kind=8) :: inst
    character(len=2) :: codret
    character(len=8) :: chnout, resu, model, chmat, carael, tych, mesh
    character(len=8) :: nomvarc, varc, k8_k, k8_l
    character(len=16) :: nomsym, k16bi1, k16bi2
    character(len=19) :: celvrc, cnovrc, cnsvrc, cnsvr2, ces1
    character(len=19) :: cart2
    integer(kind=8), pointer :: cnsvrc_d(:) => null()
    real(kind=8), pointer :: cnsvrc_v(:) => null()
    real(kind=8), pointer :: cnsvr2_v(:) => null()
    character(len=8), pointer :: cvrcvarc(:) => null()
    character(len=8), pointer :: cnsvrc_c(:) => null()
    character(len=16), pointer :: cesv(:) => null()
    aster_logical, pointer :: cnsvrc_l(:) => null()
    aster_logical, pointer :: cnsvr2_l(:) => null()

!-----------------------------------------------------------------------
!
    call jemarq()
!
!-----------------------------------------------------------------------
!    prealables
!-----------------------------------------------------------------------
!
!   recup des entrees / sortie de l'operateur
    call getvid(' ', 'RESULTAT', scal=resu, nbret=ier)
    call getvr8(' ', 'INST', scal=inst, nbret=ier)
    call getvtx(' ', 'NOM_VARC', scal=nomvarc, nbret=ier)
    call getres(chnout, k16bi1, k16bi2)
!
!   recup du modele
    call dismoi('NOM_MODELE', resu, 'RESULTAT', repk=model)
    ASSERT((model .ne. '#AUCUN') .or. (model .ne. '#PLUSIEURS'))
!
!   recup du maillage
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
    call dismoi('NB_NO_MAILLA', mesh, 'MAILLAGE', repi=nbno)
    call dismoi('NB_MA_MAILLA', mesh, 'MAILLAGE', repi=nbma)
!
!   recup du cham_mater
    call dismoi('CHAM_MATER', resu, 'RESULTAT', repk=chmat)
    ASSERT((chmat .ne. '#AUCUN') .or. (chmat .ne. '#PLUSIEURS'))
!
!   presence de variables de commande obligatoire
    call jeexin(chmat//'.CVRCVARC', iret)
    ASSERT(iret .gt. 0)
!
!   verif sur la presence de CARA_ELEM
    call dismoi('CARA_ELEM', resu, 'RESULTAT', repk=carael)
    ASSERT(carael .ne. '#PLUSIEURS')
    if (carael .eq. '#AUCUN') then
        carael = ' '
    end if
!
!   recup de chmat//'.CVRCVARC'
    call jelira(chmat//'.CVRCVARC', 'LONMAX', nbvarc)
    call jeveuo(chmat//'.CVRCVARC', 'L', vk8=cvrcvarc)
    ASSERT(nbvarc .le. nvrcmx)
!
!   pas de doublons dans '.CVRCVARC' (est-ce d'ailleurs possible?)
    cpt_dbl = 0
    do k = 1, nbvarc
        k8_k = cvrcvarc(k)
        do l = k+1, nbvarc
            k8_l = cvrcvarc(l)
            if (k8_k .eq. k8_l) then
                cpt_dbl = cpt_dbl+1
            end if
        end do
    end do
    ASSERT(cpt_dbl .eq. 0)
!
!   recup de nucmp et itrou
!     nucmp  : numero de la cmp qui correspond a nomvarc dans le champ
!             out de vrcins, l'ordre est celui de '.CVRCVARC'
!     itrou : position de nomvarc dans tab_nomvrc
    nucmp = 0
    itrou = 0
    do k = 1, nbvarc
        itest = indik8(tab_nomvrc, cvrcvarc(k), 1, nvrcmx)
        ASSERT(itest .gt. 0)
        if (cvrcvarc(k) .eq. nomvarc) then
            nucmp = k
            itrou = itest
        end if
    end do
    ASSERT(nucmp .gt. 0)
    ASSERT(itrou .gt. 0)
!
!   verif pour interdire le cas particulier de la VARC TEMP aux points
!   de Gauss xfem (cas du chainage thermo mecanique complet avec xfem)
    do k = 1, nbvarc
        varc = cvrcvarc(k)
        cart2 = chmat//'.'//varc//'.2'
        ces1 = '&&'//nompro//'.CES1'
        call carces(cart2, 'ELEM', ' ', 'V', ces1, 'A', iret)
        ASSERT(iret .eq. 0)
        call jeveuo(ces1//'.CESD', 'L', jcesd1)
        call jeveuo(ces1//'.CESL', 'L', jcesl1)
        call jeveuo(ces1//'.CESV', 'L', vk16=cesv)
        do ima = 1, nbma
            call cesexi('C', jcesd1, jcesl1, ima, 1, 1, 1, iad)
            if (iad .le. 0) then
                cycle
            end if
            iad = iad-1
            nomsym = cesv(iad+4)
            if (nomsym .eq. 'TEMP_ELGA') then
                call utmess('F', 'XFEM_1', sk=chmat)
            end if
        end do
        call detrsd('CHAM_ELEM_S', ces1)
    end do
!
!-----------------------------------------------------------------------
!   recup du cham_elem (ELNO) de varc a l'instant inst
!-----------------------------------------------------------------------
!
    celvrc = '&&'//nompro//'.CELVRC'
    call vrcins(model, chmat, carael, inst, celvrc, codret, nompaz='PVARCNO')
    ASSERT(codret .eq. 'OK')
    call exisd('CHAM_ELEM', celvrc, iret)
    ASSERT(iret .eq. 1)
    call dismoi('TYPE_CHAMP', celvrc, 'CHAMP', repk=tych)
    ASSERT(tych .eq. 'ELNO')
!
!-----------------------------------------------------------------------
!   transformation en cham_no. Le nom de la grandeur depend de nomvarc :
!   TEMP  -> TEMP_R
!   NEUT1 -> NEUT_R
!-----------------------------------------------------------------------
!
!   transformation cham_elem (ELNO) -> cham_no : celvrc -> cnovrc
    cnovrc = '&&'//nompro//'.CNOVRC'
    call chpchd(celvrc, 'NOEU', ' ', 'NON', 'V', cnovrc, model)
    call detrsd('CHAM_ELEM', celvrc)
!
!   transformation cham_no -> cham_no_s : cnovrc -> cnsvrc
    cnsvrc = '&&'//nompro//'.CNSVRC'
    call cnocns(cnovrc, 'V', cnsvrc)
    call detrsd('CHAM_NO', cnovrc)
!
!   recup du nombre de cmp dans cnsvrc
    call jeveuo(cnsvrc//'.CNSD', 'L', vi=cnsvrc_d)
    call jeveuo(cnsvrc//'.CNSC', 'L', vk8=cnsvrc_c)
    nbcmp = cnsvrc_d(2)
    ASSERT(nbcmp .eq. nbvarc)
    ASSERT(nbcmp .ge. nucmp)
!
!   creation d'un cham_no_s cnsvr2 dont le type de la grandeur depend
!   de nomvarc
    cnsvr2 = '&&'//nompro//'.CNSVR2'
    call cnscre(mesh, tab_nomgrd(itrou), 1, [tab_nomcmp(itrou)], 'V', cnsvr2)
!
!   copie de la nucmp-ieme cmp de cnsvrc dans cnsvr2
    call jeveuo(cnsvrc//'.CNSL', 'L', vl=cnsvrc_l)
    call jeveuo(cnsvrc//'.CNSV', 'L', vr=cnsvrc_v)
    call jeveuo(cnsvr2//'.CNSL', 'E', vl=cnsvr2_l)
    call jeveuo(cnsvr2//'.CNSV', 'E', vr=cnsvr2_v)
    do ino = 1, nbno
        ind = nbcmp*(ino-1)+nucmp
        if (cnsvrc_l(ind)) then
            cnsvr2_l(ino) = .true.
            cnsvr2_v(ino) = cnsvrc_v(ind)
        end if
    end do
    call detrsd('CHAM_NO_S', cnsvrc)
!
!   transformation cham_no_s ->  cham_no : cnsvr2 -> chnout
    call cnscno(cnsvr2, ' ', 'NON', 'G', chnout, 'F', ibid)
    call detrsd('CHAM_NO_S', cnsvr2)
!
!-----------------------------------------------------------------------
!   fin
!-----------------------------------------------------------------------
!
    call jedema()
!
end subroutine
