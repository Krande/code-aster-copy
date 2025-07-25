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
!
subroutine engtrs(ific, nomsd, typtes, preci, formr)
    implicit none
#include "jeveux.h"
#include "asterc/ismaem.h"
#include "asterc/r8maem.h"
#include "asterc/r8prem.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/lxlgut.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsorac.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: ific
    character(len=8) :: typtes
    character(len=10) :: preci, formr
    character(len=19) :: nomsd
!     COMMANDE:  ENGENDRE_TEST
!                TRAITEMENT DES SD RESULTAT
!
! IN  : IFIC   : NUMERO D'UNITE IMPRESSION
! IN  : NOMSD : NOM D'UNE SD RESULTAT
! IN  : TYPTES : TYPE DU TEST = SOMM_ABS, SOMM
! IN  : PRECI  : PRECISION POUR LE TEST_RESU
! ----------------------------------------------------------------------
!
    integer(kind=8) :: ibid, nbordt(1), vali, jordr, nbnosy, isy, iatach, lg, i, j, iord
    integer(kind=8) :: iret, jvale, long, lg1, lg2
    real(kind=8) :: r8b, valr
    complex(kind=8) :: c16b
    character(len=3) :: type
    character(len=8) :: k8b
    character(len=16) :: nomsym
    character(len=19) :: chextr
    character(len=90) :: form1, form2, form3
!     ------------------------------------------------------------------
!
    call jemarq()
!
    lg1 = lxlgut(formr)
    lg2 = lxlgut(typtes)
    form1 = '(&
&            '' TYPE_TEST= '''''//typtes(1:lg2)//''''', VALE_CALC= '', '//formr(1:lg1)//&
&           ', '' ), '' )'
    form2 = '( '' TYPE_TEST= '''''//typtes(1:lg2)//''''', VALE_CALC_I = '', I9, '' ), '' )'
!
    write (ific, 100)
!
! --- NUMEROS D'ORDRE
!
    call rsorac(nomsd, 'LONUTI', 0, r8b, k8b, &
                c16b, r8b, k8b, nbordt, 1, &
                ibid)
    call wkvect('&&ENGTRS.NUME_ORDRE', 'V V I', nbordt(1), jordr)
    call rsorac(nomsd, 'TOUT_ORDRE', 0, r8b, k8b, &
                c16b, r8b, k8b, zi(jordr), nbordt(1), &
                ibid)
!
! --- NOMS SYMBOLIQUES
!
    call jelira(nomsd//'.DESC', 'NOMMAX', nbnosy)
    do isy = 1, nbnosy
        call jenuno(jexnum(nomsd//'.DESC', isy), nomsym)
        call jenonu(jexnom(nomsd//'.DESC', nomsym), ibid)
        call jeveuo(jexnum(nomsd//'.TACH', ibid), 'L', iatach)
        lg = lxlgut(nomsym)
!
        form3 = '(&
&                '' _F(RESULTAT= '',A8,'', NOM_CHAM= '''''//nomsym(1:lg) &
                //''''', NUME_ORDRE= '',I6,'','' )'
!
        do j = 1, nbordt(1)
            iord = zi(jordr+j-1)
            if (zk24(iatach-1+j) (1:1) .ne. ' ') then
                call rsexch(' ', nomsd, nomsym, iord, chextr, &
                            ibid)
!
                call jeexin(chextr//'.VALE', iret)
                if (iret .ne. 0) then
                    call jeveuo(chextr//'.VALE', 'L', jvale)
                    call jelira(chextr//'.VALE', 'LONMAX', long)
                    if (long .eq. 0) goto 110
                    call jelira(chextr//'.VALE', 'TYPE', cval=type)
                    goto 120
                end if
                call jeexin(chextr//'.CELV', iret)
                if (iret .ne. 0) then
                    call jeveuo(chextr//'.CELV', 'L', jvale)
                    call jelira(chextr//'.CELV', 'LONMAX', long)
                    if (long .eq. 0) goto 110
                    call jelira(chextr//'.CELV', 'TYPE', cval=type)
                    goto 120
                end if
                goto 110
120             continue
!
                write (ific, form3) nomsd(1:8), iord
                write (ific, 102) preci
!
                if (type .eq. 'I') then
                    if (typtes .eq. 'SOMM_ABS') then
                        vali = 0
                        do i = 1, long
                            vali = vali+abs(zi(jvale+i-1))
                        end do
                    else if (typtes .eq. 'SOMM') then
                        vali = 0
                        do i = 1, long
                            vali = vali+zi(jvale+i-1)
                        end do
                    else if (typtes .eq. 'MAX') then
                        vali = -ismaem()
                        do i = 1, long
                            vali = max(vali, zi(jvale+i-1))
                        end do
                    else if (typtes .eq. 'MIN') then
                        vali = ismaem()
                        do i = 1, long
                            vali = min(vali, zi(jvale+i-1))
                        end do
                    end if
                    if (vali .eq. 0) write (ific, 101)
                    write (ific, form2) vali
!
                else if (type .eq. 'R') then
                    if (typtes .eq. 'SOMM_ABS') then
                        valr = 0.d0
                        do i = 1, long
                            valr = valr+abs(zr(jvale+i-1))
                        end do
                    else if (typtes .eq. 'SOMM') then
                        valr = 0.d0
                        do i = 1, long
                            valr = valr+zr(jvale+i-1)
                        end do
                    else if (typtes .eq. 'MAX') then
                        valr = -r8maem()
                        do i = 1, long
                            valr = max(valr, zr(jvale+i-1))
                        end do
                    else if (typtes .eq. 'MIN') then
                        valr = r8maem()
                        do i = 1, long
                            valr = min(valr, zr(jvale+i-1))
                        end do
                    end if
                    if (abs(valr) .le. r8prem()) write (ific, 101)
                    write (ific, form1) valr
                end if
!
            end if
110         continue
        end do
    end do
!
    write (ific, 103)
!
    call jedetr('&&ENGTRS.NUME_ORDRE')
!
    call jedema()
!
100 format('TEST_RESU(RESU=( ')
101 format('              CRITERE= ''ABSOLU'', ')
102 format('              TOLE_MACHINE= ', a10, ',')
103 format('          ),)')
!
end subroutine
