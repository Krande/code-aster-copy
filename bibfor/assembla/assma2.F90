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

subroutine assma2(ldistme, lmasym, tt, nu14, ncmp, matel, &
                  c1, jvalm, jtmp2, lgtmp2)
! person_in_charge: jacques.pellet at edf.fr
!
    implicit none

!-----------------------------------------------------------------------
! but : assembler les macro-elements dans une matr_asse
!-----------------------------------------------------------------------
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/indik8.h"
#include "asterfort/ascopr.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/asretm.h"
#include "asterfort/cordd2.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/nbec.h"
#include "asterfort/ssvalm.h"
!-----------------------------------------------------------------------
    aster_logical, intent(in) :: ldistme, lmasym
    character(len=2), intent(in) :: tt
    character(len=14), intent(in) :: nu14
    integer(kind=8), intent(in) :: ncmp
    character(len=19), intent(in) :: matel
    real(kind=8), intent(in) :: c1
    integer(kind=8) :: jvalm(2)
    integer(kind=8) :: jtmp2
    integer(kind=8) :: lgtmp2
!-----------------------------------------------------------------------
    character(len=16) :: optio
    character(len=19) :: ligrel
    aster_logical :: lmesym
    character(len=8) :: mo, ma, nogdco, nogdsi, nomacr
    character(len=14) :: num2
    integer(kind=8) :: nbecmx
    parameter(nbecmx=11)
    integer(kind=8) :: icodla(nbecmx), icodge(nbecmx)
    integer(kind=8) :: i1, i2, iad1, iad11, iad2, iad21
    integer(kind=8) :: jsupma, jnulo1, jprno, jposd1
    integer(kind=8) :: iec, ima, inold, nbterm, jprn1, jprn2
    integer(kind=8) :: jresl, jsmdi, jsmhc, k1
    integer(kind=8) :: k2, n1, nugd, iancmp, lgncmp, icmp
    integer(kind=8) :: nbsma, nbssa, nbvel, nddl1, nddl2
    integer(kind=8) :: nec, nm, nmxcmp, nnoe, i, jec, rang
    integer(kind=8) :: lshift
    integer(kind=8), pointer :: conx(:) => null()
    character(len=8), pointer :: vnomacr(:) => null()
    integer(kind=8), pointer :: sssa(:) => null()
    integer(kind=8), pointer :: nueq(:) => null()
    mpi_int :: mrank, msize

!-----------------------------------------------------------------------
!     FONCTIONS FORMULES :
!-----------------------------------------------------------------------
!
#define zzprno(ili,nunoel,l) zi(jprn1-1+zi(jprn2+ili-1)+ \
    (nunoel-1)*(nec+2)+l-1)
#define numlo1(kno,k) zi(jnulo1-1+2*(kno-1)+k)
#define posdd1(kno,kddl) zi(jposd1-1+nmxcmp*(kno-1)+kddl)
!-----------------------------------------------------------------------
!
!
    call jemarq()
!
    call dismoi('NB_SS_ACTI', matel, 'MATR_ELEM', repi=nbssa)
    if (nbssa .eq. 0) goto 100

!   -- si ldistme, il ne faut pas assembler plusieurs fois les macro-elements.
!      Seul le proc 0 le fait
    call asmpi_info(rank=mrank, size=msize)
    rang = to_aster_int(mrank)
    if (ldistme .and. rang .gt. 0) goto 100

    lmesym = .true.
    do i = 1, nbecmx
        icodla(i) = 0
        icodge(i) = 0
    end do
!
    call dismoi('NOM_MODELE', nu14, 'NUME_DDL', repk=mo)
    call dismoi('NOM_MAILLA', mo, 'MODELE', repk=ma)
    call dismoi('NB_NO_MAILLA', mo, 'MODELE', repi=nm)
    call dismoi('NB_SM_MAILLA', mo, 'MODELE', repi=nbsma)
    call jeveuo(ma//'.NOMACR', 'L', vk8=vnomacr)
    call dismoi('NOM_LIGREL', mo, 'MODELE', repk=ligrel)
    call jeveuo(ligrel//'.SSSA', 'L', vi=sssa)
!
    call jeveuo(nu14//'.SMOS.SMDI', 'L', jsmdi)
    call jeveuo(nu14//'.SMOS.SMHC', 'L', jsmhc)
    call jeveuo(nu14//'.NUME.NUEQ', 'L', vi=nueq)
    call jeveuo(nu14//'.NUME.PRNO', 'L', jprn1)
    call jeveuo(jexatr(nu14//'.NUME.PRNO', 'LONCUM'), 'L', jprn2)
!
    call dismoi('SUR_OPTION', matel, 'MATR_ELEM', repk=optio)
    call dismoi('NOM_GD', nu14, 'NUME_DDL', repk=nogdco)
    call dismoi('NOM_GD_SI', nogdco, 'GRANDEUR', repk=nogdsi)
    call dismoi('NB_CMP_MAX', nogdsi, 'GRANDEUR', repi=nmxcmp)
    call dismoi('NUM_GD_SI', nogdsi, 'GRANDEUR', repi=nugd)
    nec = nbec(nugd)
    call jeveuo(jexnom('&CATA.GD.NOMCMP', nogdsi), 'L', iancmp)
    call jelira(jexnom('&CATA.GD.NOMCMP', nogdsi), 'LONMAX', lgncmp)
    icmp = indik8(zk8(iancmp), 'LAGR', 1, lgncmp)
    if (icmp .gt. 0) then
        jec = (icmp-1)/30+1
        icodla(jec) = lshift(1, icmp-(jec-1)*30)
    end if
!
    call jeveuo('&&ASSMAM.NUMLO1', 'E', jnulo1)
    call jeveuo('&&ASSMAM.POSDD1', 'E', jposd1)
!
!
!
    call ssvalm('DEBUT', optio, mo, ma, 0, &
                jresl, nbvel)
!
!   -- boucle sur les macro-elements :
!   ----------------------------------
    do ima = 1, nbsma
        if (sssa(ima) .eq. 0) goto 90
!
        call jeveuo(jexnum(ma//'.SUPMAIL', ima), 'L', jsupma)
        call jelira(jexnum(ma//'.SUPMAIL', ima), 'LONMAX', nnoe)
!
        nbterm = 0
!
        call ssvalm(' ', optio, mo, ma, ima, &
                    jresl, nbvel)
!
        nomacr = vnomacr(ima)
        call dismoi('NOM_NUME_DDL', nomacr, 'MACR_ELEM_STAT', repk=num2)
        call jeveuo(nomacr//'.CONX', 'L', vi=conx)
        call jeveuo(jexnum(num2//'.NUME.PRNO', 1), 'L', jprno)
!
        do k1 = 1, nnoe
            n1 = zi(jsupma-1+k1)
            if (n1 .gt. nm) then
                do iec = 1, nbecmx
                    icodge(iec) = icodla(iec)
                end do
!
            else
                inold = conx(3*(k1-1)+2)
                do iec = 1, nec
                    icodge(iec) = zi(jprno-1+(nec+2)*(inold-1)+2+iec)
                end do
            end if
!
            iad1 = zzprno(1, n1, 1)
            call cordd2(jprn1, jprn2, 1, icodge, nec, &
                        ncmp, n1, nddl1, zi(jposd1-1+nmxcmp*(k1-1)+1))
            zi(jnulo1-1+2*(k1-1)+1) = iad1
            zi(jnulo1-1+2*(k1-1)+2) = nddl1
            do i1 = 1, nddl1
                do k2 = 1, k1-1
                    iad2 = numlo1(k2, 1)
                    nddl2 = numlo1(k2, 2)
                    do i2 = 1, nddl2
                        iad11 = nueq(iad1+posdd1(k1, i1)-1)
                        iad21 = nueq(iad2+posdd1(k2, i2)-1)
                        call asretm(lmasym, jtmp2, lgtmp2, nbterm, jsmhc, &
                                    jsmdi, iad11, iad21)
                    end do
                end do
                k2 = k1
                iad2 = numlo1(k2, 1)
                nddl2 = numlo1(k2, 2)
                do i2 = 1, i1
                    iad11 = nueq(iad1+posdd1(k1, i1)-1)
                    iad21 = nueq(iad2+posdd1(k2, i2)-1)
                    call asretm(lmasym, jtmp2, lgtmp2, nbterm, jsmhc, &
                                jsmdi, iad11, iad21)
                end do
            end do
        end do
!
!
!       -- pour finir, on recopie effectivement les termes:
        call ascopr(lmasym, lmesym, 'R'//tt(2:2), jtmp2, nbterm, &
                    jresl, c1, jvalm)
90      continue
    end do
    call ssvalm('FIN', optio, mo, ma, ima, &
                jresl, nbvel)
!
!
100 continue
    call jedema()
end subroutine
