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

subroutine pjcorr(nomo2, chbid, cns1z, ces2z, ligrel, &
                  corres, option, nompar, iret)
! person_in_charge: jacques.pellet at edf.fr
!
!
    implicit none
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/alchml.h"
#include "asterfort/assert.h"
#include "asterfort/celces.h"
#include "asterfort/cescre.h"
#include "asterfort/cesexi.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
    character(len=*) :: cns1z, ces2z
    character(len=8) :: nomo2, nompar
    character(len=16) :: corres, option
    character(len=19) :: chbid, ligrel
    integer(kind=8) :: iret
!
! COMMANDE PROJ_CHAMP
!   RECOPIE DES VALEURS PROJETEES AUX NOEUDS DU MAILLAGE 2 PRIME
!     SUR LES POINTS DE GAUSS DU MAILLAGE 2
!   UTILISATION DE LA SD CORRES (TABLEAU AUXILIAIRE .PJEF_EL)
!
!
!
!     ------------------------------------------------------------------
!     VARIABLES LOCALES:
!     ------------------
    character(len=8) :: nomgd, ma, licmp(2)
    character(len=19) :: cns1, ces2, cel2
    character(len=19) :: cham1s, dcel
    character(len=24) :: valk(5)
    integer(kind=8) ::  ipo, nbma
    integer(kind=8) ::  jce2l, jce2d, jce2k, icmp2
!
    integer(kind=8) :: jcesk, jcesd, jcesv, jcesl
!
    integer(kind=8) ::  jcns1l
    integer(kind=8) :: nbno1, nbmax, ncmp1, ncmp2
    integer(kind=8) :: iad2, nval
    integer(kind=8) :: iad, nbpt, nbsp, icmp1
!
    integer(kind=8) :: ima, ipt, isp, jcesc
    integer(kind=8), pointer :: cnsd(:) => null()
    character(len=8), pointer :: cnsk(:) => null()
    character(len=8), pointer :: lgrf(:) => null()
    real(kind=8), pointer :: cnsv(:) => null()
    integer(kind=8), pointer :: pjef_el(:) => null()
    character(len=8), pointer :: ce2c(:) => null()
    character(len=8), pointer :: cnsc(:) => null()
    real(kind=8), pointer :: ce2v(:) => null()
!
!     ------------------------------------------------------------------
!
    call jemarq()
    iret = 0
    cns1 = cns1z
    ces2 = ces2z
!
    cel2 = '&&PJCORR.CEL2'
!
!     1- RECUPERATION D'INFORMATIONS DANS CNS1 :
!     ------------------------------------------
    call jeveuo(cns1//'.CNSK', 'L', vk8=cnsk)
    call jeveuo(cns1//'.CNSD', 'L', vi=cnsd)
    call jeveuo(cns1//'.CNSC', 'L', vk8=cnsc)
    call jeveuo(cns1//'.CNSV', 'L', vr=cnsv)
    call jeveuo(cns1//'.CNSL', 'L', jcns1l)
    call jelira(cns1//'.CNSC', 'LONMAX', ncmp1)
!
    nomgd = cnsk(2)
    nbno1 = cnsd(1)
    nbmax = cnsd(2)
!
!
!------------------------------------------------------------------
!     2- ALLOCATION DE CES2 (ELGA):
!
!
!     -----------------------------
    call detrsd('CHAM_ELEM_S', ces2)
!
    cham1s = 'CHAM1S'
!
    call celces(chbid, 'V', cham1s)
!
!
    call jeveuo(cham1s//'.CESK', 'L', jcesk)
    call jeveuo(cham1s//'.CESD', 'L', jcesd)
    call jeveuo(cham1s//'.CESC', 'L', jcesc)
    call jeveuo(cham1s//'.CESV', 'L', jcesv)
    call jeveuo(cham1s//'.CESL', 'L', jcesl)
!
    call jeveuo(ligrel//'.LGRF', 'L', vk8=lgrf)
    ma = lgrf(1)
!
!
!
!
    if (nomgd .eq. 'VAR2_R') then
!
        dcel = '&&PJCORR'
!
        licmp(1) = 'NPG_DYN'
        licmp(2) = 'NCMP_DYN'
        call cescre('V', dcel, 'ELEM', ma, 'DCEL_I', &
                    2, licmp, [-1], [-1], [-2])
!
        call jeveuo(dcel//'.CESD', 'E', jcesd)
        call jeveuo(dcel//'.CESL', 'E', jcesl)
        call jeveuo(dcel//'.CESV', 'E', jcesv)
!
!
        call dismoi('NB_MA_MAILLA', nomo2, 'MODELE', repi=nbma)
!
        do ima = 1, nbma
            nbpt = zi(jcesd-1+5+4*(ima-1)+1)
            ASSERT(nbpt .eq. 1)
            nbsp = zi(jcesd-1+5+4*(ima-1)+2)
            ASSERT(nbsp .eq. 1)
!
            do ipt = 1, nbpt
                do isp = 1, nbsp
                    call cesexi('C', jcesd, jcesl, ima, ipt, &
                                isp, 1, iad)
                    zi(jcesv-1-iad) = 0
                    zl(jcesl-1-iad) = .true.
                    call cesexi('C', jcesd, jcesl, ima, ipt, &
                                isp, 2, iad)
                    zi(jcesv-1-iad) = nbmax
                    zl(jcesl-1-iad) = .true.
                end do
            end do
        end do
!
        call alchml(ligrel, option, nompar, 'V', cel2, &
                    iret, dcel)
!
        if (iret .eq. 1) then
            valk(1) = nompar
            valk(2) = option
            valk(3) = ligrel
            valk(4) = cel2
            call utmess('F', 'CALCULEL_50', nk=4, valk=valk)
        end if
!
    else
        call alchml(ligrel, option, nompar, 'V', cel2, &
                    iret, ' ')
!
        if (iret .eq. 1) then
            valk(1) = nompar
            valk(2) = option
            valk(3) = ligrel
            valk(4) = cel2
            call utmess('F', 'CALCULEL_50', nk=4, valk=valk)
        end if
!
    end if
!
    call celces(cel2, 'V', ces2)
    call detrsd('CHAM_ELEM', cel2)
!
    call jeveuo(ces2//'.CESD', 'L', jce2d)
    call jeveuo(ces2//'.CESC', 'L', vk8=ce2c)
    call jeveuo(ces2//'.CESV', 'E', vr=ce2v)
    call jeveuo(ces2//'.CESL', 'E', jce2l)
    call jeveuo(ces2//'.CESK', 'L', jce2k)
    call jelira(ces2//'.CESC', 'LONMAX', ncmp2)
!   -- on met les booleens a .false. :
    call jelira(ces2//'.CESL', 'LONMAX', nval)
    zl(jce2l-1+1:jce2l-1+nval) = .false.
!
!
!
!------------------------------------------------------------------
!     3- REMPLISSAGE DES VALEURS DE CES2 :
!     -------------------------------
    call jeveuo(corres//'.PJEF_EL', 'L', vi=pjef_el)
!
    do icmp1 = 1, ncmp1
        icmp2 = indik8(ce2c, cnsc(icmp1), 1, ncmp2)
        if (icmp2 .eq. 0) goto 92
        ASSERT(ce2c(icmp2) .eq. cnsc(icmp1))
!       -- nbno1 est le nombre de pseudo-noeuds du maillage 2
        do ipo = 1, nbno1
            ima = pjef_el(2*ipo-1)
            ipt = pjef_el(2*ipo)
            call cesexi('C', jce2d, jce2l, ima, ipt, &
                        1, icmp2, iad2)
            ASSERT(iad2 .le. 0)
            iad2 = -iad2
            if (iad2 .eq. 0) goto 98
!
            ce2v(iad2) = cnsv(1+(ipo-1)*ncmp1+icmp1-1)
            zl(jce2l-1+iad2) = .true.
98          continue
        end do
92      continue
    end do
!
!
    call jedema()
end subroutine
