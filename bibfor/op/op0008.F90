! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
!
subroutine op0008()
!
implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/getres.h"
#include "asterfort/assert.h"
#include "asterfort/cvePost.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/jecreo.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/me2mac.h"
#include "asterfort/me2mme_2.h"
#include "asterfort/me2mth.h"
#include "asterfort/mecact.h"
#include "asterfort/rcmfmc.h"
#include "asterfort/ss2mme.h"
#include "asterfort/utmess.h"
!
! --------------------------------------------------------------------------------------------------
!
!                       COMMANDE:  CALC_VECT_ELEM
!
! --------------------------------------------------------------------------------------------------
!
    integer :: ibid, ich, icha, ncha, nh
    integer :: n1, n3, n4, n5, n7, n9
    real(kind=8) :: time, tps(6), vcmpth(4)
    character(len=8) :: vectElem, modele, cara, k8bid
    character(len=8) :: nomcmp(6), mo1, ncmpth(4)
    character(len=16) :: type, oper, suropt
    character(len=24) :: time2, mateco, materi
    aster_logical :: l_ther
    data nomcmp/'INST    ','DELTAT  ','THETA   ','KHI     ',&
     &     'R       ','RHO     '/
    data ncmpth/'TEMP','TEMP_MIL','TEMP_INF','TEMP_SUP'/
    data vcmpth/4*0.d0/
    data tps/0,2*1.0d0,3*0/
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infmaj()
!
    call getres(vectElem, type, oper)
!
    call getvtx(' ', 'OPTION', scal=suropt, nbret=n3)
    l_ther = ASTER_FALSE
    if (suropt .eq. 'CHAR_THER') then
        l_ther = ASTER_TRUE
    endif
!
!     - ON VERIFIE LE NOM DU MODELE:
!     -------------------------------
    modele = ' '
    call getvid(' ', 'MODELE', scal=modele, nbret=n1)
    call getvid(' ', 'CHARGE', nbval=0, nbret=ncha)
!
!
    if (ncha .lt. 0) then
        ncha = -ncha
        call jecreo(vectElem(1:8)//'.CHARGES', 'V V K8')
        n3=max(1,ncha)
        call jeecra(vectElem(1:8)//'.CHARGES', 'LONMAX', n3)
        call jeveuo(vectElem(1:8)//'.CHARGES', 'E', icha)
        call getvid(' ', 'CHARGE', nbval=ncha, vect=zk8(icha), nbret=ibid)
!
        call dismoi('NOM_MODELE', zk8(icha), 'CHARGE', repk=mo1)
        if ((n1.eq.1) .and. (modele.ne.mo1)) then
            call utmess('F', 'CALCULEL3_88')
        endif
!
        modele = mo1
        do ich = 1, ncha
            call dismoi('NOM_MODELE', zk8(icha-1+ich), 'CHARGE', repk=k8bid)
            if (k8bid .ne. modele) then
                call utmess('F', 'CALCULEL3_89')
            endif
        end do
    endif
!
!
    cara = ' '
    materi = ' '
    call getvid(' ', 'CARA_ELEM', scal=cara, nbret=n5)
    call getvid(' ', 'CHAM_MATER', scal=materi, nbret=n4)
    if (n4 .ne. 0) then
        call rcmfmc(materi, mateco, l_ther_ = l_ther)
    else
        mateco = ' '
    endif
!
    call getvr8(' ', 'INST', scal=time, nbret=n7)
    call getvis(' ', 'MODE_FOURIER', scal=nh, nbret=n9)
    if (n9 .eq. 0) nh = 0
!
!     -- VERIFICATION DES CHARGES:
    if ((suropt.eq.'CHAR_MECA') .or. (suropt.eq.'CHAR_MECA_LAGR')) then
        do ich = 1, ncha
            call dismoi('TYPE_CHARGE', zk8(icha-1+ich), 'CHARGE', repk=k8bid)
            if (k8bid(1:5) .ne. 'MECA_') then
                call utmess('F', 'CALCULEL3_91')
            endif
        end do
    endif
!
    if ((suropt.eq.'CHAR_THER')) then
        do ich = 1, ncha
            call dismoi('TYPE_CHARGE', zk8(icha-1+ich), 'CHARGE', repk=k8bid)
            if (k8bid(1:5) .ne. 'THER_') then
                call utmess('F', 'CALCULEL3_92')
            endif
        end do
    endif
!
    if ((suropt.eq.'CHAR_ACOU')) then
        do ich = 1, ncha
            call dismoi('TYPE_CHARGE', zk8(icha-1+ich), 'CHARGE', repk=k8bid)
            if (k8bid(1:5) .ne. 'ACOU_') then
                call utmess('F', 'CALCULEL3_93')
            endif
        end do
    endif
!
!
!
    if (suropt .eq. 'CHAR_MECA') then
!     ----------------------------------
!        -- TRAITEMENT DES ELEMENTS FINIS CLASSIQUES (.RELR)
!           (ET CREATION DE L'OBJET .RERR).
        call me2mme_2(modele, ncha, zk8(icha), materi, mateco, cara,&
                    time, vectElem, nh, 'G')
!
!        -- TRAITEMENT DES SOUS-STRUCTURES EVENTUELLES. (.RELC):
        call ss2mme(modele, vectElem, 'G')
!
!
    else if (suropt.eq.'CHAR_THER') then
!     ----------------------------------
        tps(1) = time
        time2 = '&TIME'
        call mecact('V', time2, 'MODELE', modele//'.MODELE', 'INST_R  ',&
                    ncmp=6, lnomcmp=nomcmp, vr=tps)
        call mecact('V', '&&OP0008.PTEMPER', 'MODELE', modele//'.MODELE', 'TEMP_R',&
                    ncmp=4, lnomcmp=ncmpth, vr=vcmpth)
        call me2mth(modele, ncha, zk8(icha), cara,&
                    time2, '&&OP0008.PTEMPER', vectElem)
    else if (suropt.eq.'CHAR_ACOU') then
        call me2mac(modele, ncha, zk8(icha), materi, mateco, vectElem)
!
    endif

! - Post-treatment
    call cvePost(vectElem)

!
    call jedema()
end subroutine
