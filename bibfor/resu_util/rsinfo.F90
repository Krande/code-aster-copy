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

subroutine rsinfo(nomcon, ifi)
!
!     RESULTAT - INFORMATION
!
!     ECRITURE DE LA STRUCTURE D'UN CONCEPT RESULTAT
!
! IN  NOMCON : NOM DU CONCEPT A IMPRIMER
! IN  IFI    : UNITE LOGIQUE D'IMPRESSION
!     ------------------------------------------------------------------
! aslint: disable=W1303
    implicit none
!
! 0.1. ==> ARGUMENTS
!
#include "jeveux.h"
#include "asterc/r8vide.h"
#include "asterfort/codent.h"
#include "asterfort/dismoi.h"
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
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnopa.h"
#include "asterfort/rsorac.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=*) :: nomcon
! 0.2. ==> COMMUNS
!
!
! 0.3. ==> VARIABLES LOCALES
!
    character(len=6) :: nompro
    parameter(nompro='RSINFO')
!
    integer(kind=8) :: ifi, ibid, nbnosy, lres, inomsy, i, j, k, isy, ii
    integer(kind=8) :: iatach, lnosy, lnopa, lnupa, latac, lg, lb
    integer(kind=8) :: iret, ltirt, nbac, nbpa, iac, ipar, iad, jpa
    integer(kind=8) :: nbordt, ipcd, ipcf, lpoin, longt, tord(1)
    real(kind=8) :: r8b, rundf
    complex(kind=8) :: c16b
    character(len=8) :: k8b, nomb1, nomgd, ctype
    character(len=16) :: nomsym, blanc, nopara, nopar2
    character(len=19) :: nomd2, noch19
    character(len=80) :: form1, form2, form3
    character(len=2000) :: chain1, chain2, chain3, chain4
!
!     ------------------------------------------------------------------
!====
! 1. PREALABLES
!====
!
    call jemarq()
!
! 1.1. ==> EXISTENCE DE LA STRUCTURE
!
    nomd2 = nomcon
    call jelira(nomd2//'.DESC', 'NOMMAX', nbnosy)
    if (nbnosy .eq. 0) then
        call utmess('A', 'UTILITAI4_34', sk=nomd2)
        goto 999
    end if
!
! 1.2. ==> CONSTANTES
!
!              1234567890123456
    blanc = '                '
    rundf = r8vide()
!
!====
! 2. ON PARCOURT LA STRUCTURE DE RESULTAT
!====
!
! 2.1. ==> UNE STRUCTURE DE RESULTAT A-T-ELLE ETE CONSTRUITE ?
!          SI NON, ON PASSE A LA SUIVANTE.
!
!                                 9012345678901234
    call jeexin(nomd2(1:8)//'           .DESC', j)
    if (j .eq. 0) then
        goto 999
    end if
!
! 2.2. ==> NUMEROS D'ORDRE
!
    call rsorac(nomd2, 'LONUTI', 0, r8b, k8b, &
                c16b, r8b, k8b, tord, 1, &
                ibid)
    nbordt = tord(1)
    if (nbordt .eq. 0) then
        write (ifi, 103) nomd2(1:8)
        goto 999
    end if
!
    if (nbordt .eq. 1) then
        write (ifi, 101) nomd2(1:8)
    else
        write (ifi, 102) nomd2(1:8), nbordt
    end if
!
101 format(/, 1x, 'STRUCTURE DU CONCEPT ', a, ' CALCULE POUR 1', &
            ' NUMERO D''ORDRE')
102 format(/, 1x, 'STRUCTURE DU CONCEPT ', a, ' CALCULE POUR ', i10, &
            ' NUMEROS D''ORDRE')
103 format(/, 1x, 'STRUCTURE DU CONCEPT ', a, ' VIDE, NE CONTENANT AUCUN NUMERO D''ORDRE')
!
    call wkvect('&&'//nompro//'.NUME_ORDRE', 'V V I', nbordt, lres)
    call rsorac(nomd2, 'TOUT_ORDRE', 0, r8b, k8b, &
                c16b, r8b, k8b, zi(lres), nbordt, &
                ibid)
!
! 2.3. ==> NOMS SYMBOLIQUES
!
    call wkvect('&&'//nompro//'.NOM_SYMBOL', 'V V K16', nbnosy, lnosy)
    call wkvect('&&'//nompro//'.NUM_SYMBOL', 'V V K16', nbnosy, latac)
    inomsy = 0
    do isy = 1, nbnosy
        call jenuno(jexnum(nomd2//'.DESC', isy), nomsym)
        call jenonu(jexnom(nomd2//'.DESC', nomsym), ibid)
        call jeveuo(jexnum(nomd2//'.TACH', ibid), 'L', iatach)
        do i = 1, nbordt
            if (zk24(iatach-1+i) (1:1) .ne. ' ') then
                inomsy = inomsy+1
                lg = lxlgut(nomsym)
                lb = (16-lg)/2
                zk16(lnosy+inomsy-1) = blanc(1:lb)//nomsym
                zk16(latac+inomsy-1) = nomsym
                goto 241
            end if
        end do
241     continue
    end do
!
    if (inomsy .eq. 0) then
        write (ifi, '(/,1X,A)') 'LISTE DES NOMS SYMBOLIQUES: AUCUN'
        goto 243
    end if
!
    call codent(inomsy, 'D', nomb1)
    form1 = '(1X,''!'',1X,A10,1X,'//nomb1//'(''!'',A16),''!'')'
    longt = 17*inomsy
    if (longt .gt. 2000) then
        call utmess('A', 'UTILITAI4_36')
        goto 999
    end if
    call codent(longt, 'G', nomb1)
    form2 = '(1X,''!'',1X,I10,1X,''!'',A'//nomb1//')'
    form3 = '(1X,''!'',1X,A10,1X,''!'',A'//nomb1//')'
!
    write (ifi, '(/,1X,A)') 'LISTE DES NOMS SYMBOLIQUES:'
!
    call wkvect('&&'//nompro//'.TIRET', 'V V K16', max(inomsy, 1), ltirt)
    call wkvect('&&'//nompro//'.POINT', 'V V K16', max(inomsy, 1), lpoin)
    do i = 1, inomsy
!                          1234567890123456
        zk16(ltirt+i-1) = '----------------'
        zk16(lpoin+i-1) = '      ...       '
    end do
    write (ifi, form1) '----------', (zk16(ltirt+j-1), j=1, inomsy)
    write (ifi, form1) 'NUME_ORDRE', (zk16(lnosy+j-1), j=1, inomsy)
    write (ifi, form1) '----------', (zk16(ltirt+j-1), j=1, inomsy)
!
    chain1 = ' '
    chain3 = ' '
!
    do i = 1, nbordt
        chain2 = ' '
        chain4 = ' '
        ipcd = 1
!
! RECHERCHE DES NOMS SYMBOLIQUES POUR LE NUMERO D'ORDRE COURANT, I
!
        do isy = 1, inomsy
            ipcf = ipcd+16-1
            nomsym = zk16(latac+isy-1)
            call rsexch(' ', nomd2, nomsym, zi(lres+i-1), noch19, &
                        iret)
            if (iret .eq. 0) then
                call dismoi('NOM_GD', noch19, 'CHAMP', repk=nomgd)
                lg = lxlgut(nomgd)
                lb = (16-lg)/2
                chain2(ipcd:ipcf) = blanc(1:lb)//nomgd(1:lg)//blanc
                chain4(ipcd:ipcf) = zk16(lpoin+isy-1)
            else
                chain2(ipcd:ipcf) = blanc
                chain4(ipcd:ipcf) = blanc
            end if
            ipcd = ipcf+1
            chain2(ipcd:ipcd) = '!'
            chain4(ipcd:ipcd) = '!'
            ipcd = ipcd+1
        end do
!
! ECRITURE : ON ECRIT TOUJOURS LA PREMIERE ET LA DERNIERE LIGNE. AU
!            MILIEU, ON N'ECRIT QUE SI LE TEXTE A CHANGE.
! APRES CHAQUE ECRITURE, ON MEMORISE II, NUMERO D'ORDRE QUI A ETE ECRIT
! . 1ERE LIGNE : LA CHAINE COMPLETE
!
        if (i .eq. 1) then
            write (ifi, form2) zi(lres+i-1), chain2(1:longt)
            ii = 1
!
! . SI LE NOUVEAU TEXTE, CHAIN2, EST DIFFERENT DE CELUI DE LA LIGNE
!   PRECEDENTE, CHAIN1
! . OU SI C'EST LA DERNIERE LIGNE
!
        else if (chain1 .ne. chain2 .or. i .eq. nbordt) then
!          . ON VIENT JUSTE D'ECRIRE CHAIN1
            if (ii .eq. (i-1)) then
                write (ifi, form2) zi(lres+i-1), chain2(1:longt)
!          . ON A ECRIT CHAIN1 DEUX NUMEROS AVANT : ON ECRIT LE
!          NUMERO PRECEDENT, I-1, ET LE COURANT, I.
            else if (ii .eq. (i-2)) then
                write (ifi, form2) zi(lres+i-2), chain1(1:longt)
                write (ifi, form2) zi(lres+i-1), chain2(1:longt)
!          . ON A ECRIT CHAIN1 PLUS DE DEUX NUMEROS AVANT : ON ECRIT
!          UNE LIGNE DE POINTILLES, LE NUMERO PRECEDENT, I-1, ET LE
!          NUMERO COURANT, I.
            else
                if (chain1 .ne. chain2) then
                    if (ii .eq. (i-3)) then
                        write (ifi, form2) zi(lres+i-3), chain1(1: &
                                                                longt)
                    else
                        write (ifi, form3) '       ...', chain3(1:longt)
                    end if
                    write (ifi, form2) zi(lres+i-2), chain1(1:longt)
                else
                    write (ifi, form3) '       ...', chain3(1:longt)
                end if
                write (ifi, form2) zi(lres+i-1), chain2(1:longt)
            end if
            ii = i
!
        end if
        chain1(1:longt) = chain2(1:longt)
        chain3(1:longt) = chain4(1:longt)
!
    end do
!
    write (ifi, form1) '----------', (zk16(ltirt+j-1), j=1, inomsy)
!
243 continue
!
    call jedetr('&&'//nompro//'.POINTEUR')
    call jedetr('&&'//nompro//'.COMPT')
    call jedetr('&&'//nompro//'.TIRET')
    call jedetr('&&'//nompro//'.POINT')
    call jedetr('&&'//nompro//'.NUM_SYMBOL')
    call jedetr('&&'//nompro//'.NOM_SYMBOL')
!
!     ------------------------------------------------------------------
!
! 2.4. ==>      --- LES NOMS DES VARIABLES D'ACCES ---
!
    call rsnopa(nomd2, 2, '&&'//nompro//'.NOMS_PARA', nbac, nbpa)
    call jeveuo('&&'//nompro//'.NOMS_PARA', 'L', jpa)
    if (nbac .ne. 0) then
        write (ifi, '(/,1X,A)') 'LISTE DES NOMS DE VARIABLES D''ACCES:'
        do iac = 1, nbac
            call rsadpa(nomd2, 'L', 1, zk16(jpa-1+iac), zi(lres), &
                        1, sjv=iad, styp=ctype, istop=0)
            if (ctype(1:1) .eq. 'I') then
                write (ifi, '(38X,A,A)') zk16(jpa-1+iac), ' DE TYPE  I'
            else if (ctype(1:1) .eq. 'R') then
                if (zr(iad) .ne. rundf) then
                    write (ifi, '(38X,A,A)') zk16(jpa-1+iac), ' DE TYPE  R'
                end if
            else if (ctype(1:3) .eq. 'K80') then
                write (ifi, '(38X,A,A)') zk16(jpa-1+iac), ' DE TYPE  K80'
            else if (ctype(1:3) .eq. 'K32') then
                write (ifi, '(38X,A,A)') zk16(jpa-1+iac), ' DE TYPE  K32'
            else if (ctype(1:3) .eq. 'K24') then
                write (ifi, '(38X,A,A)') zk16(jpa-1+iac), ' DE TYPE  K24'
            else if (ctype(1:3) .eq. 'K16') then
                write (ifi, '(38X,A,A)') zk16(jpa-1+iac), ' DE TYPE  K16'
            else if (ctype(1:2) .eq. 'K8') then
                write (ifi, '(38X,A,A)') zk16(jpa-1+iac), ' DE TYPE  K8'
            end if
        end do
    end if
!
!     ------------------------------------------------------------------
!
! 2.5. ==>          --- LES NOMS DES PARAMETRES ---
!
    if (nbpa .ne. 0) then
        call wkvect('&&'//nompro//'.NOM_PARA', 'V V K16', nbpa, lnopa)
        call wkvect('&&'//nompro//'.NUM_PARA', 'V V K16', nbpa, lnupa)
        ipar = 0
        do i = 1, nbordt
            do j = 1, nbpa
                nopara = zk16(jpa-1+j+nbac)
                call rsadpa(nomd2, 'L', 1, nopara, zi(lres+i-1), &
                            1, sjv=iad, styp=ctype, istop=0)
                if (ctype(1:1) .eq. 'I') then
                else if (ctype(1:1) .eq. 'R') then
                    if (zr(iad) .eq. rundf) goto 261
                else if (ctype(1:1) .eq. 'K') then
                else
                    goto 261
                end if
                lg = lxlgut(nopara)
                lb = (16-lg)/2
                nopar2 = blanc(1:lb)//nopara
                do k = 1, ipar
                    if (zk16(lnopa+k-1) .eq. nopar2) goto 261
                end do
                ipar = ipar+1
                zk16(lnopa+ipar-1) = nopar2
                zk16(lnupa+ipar-1) = nopara
261             continue
            end do
        end do
!
        call codent(ipar, 'D', nomb1)
        form1 = '(1X,''!'',1X,A10,1X,'//nomb1//'(''!'',A16),''!'')'
        ipcd = 17*ipar
        call codent(ipcd, 'G', nomb1)
        form2 = '(1X,''!'',1X,I10,1X,''!'',A'//nomb1//')'
!
        call wkvect('&&'//nompro//'.TIRET', 'V V K16', max(ipar, 1), ltirt)
        call wkvect('&&'//nompro//'.POINT', 'V V K16', max(ipar, 1), lpoin)
        do i = 1, ipar
            zk16(ltirt+i-1) = '----------------'
            zk16(lpoin+i-1) = '      ...       '
        end do
!
        write (ifi, '(/,1X,A)') 'LISTE DES NOMS DE PARAMETRES:'
        write (ifi, form1) '----------', (zk16(ltirt+j-1), j=1, ipar)
        write (ifi, form1) 'NUME_ORDRE', (zk16(lnopa+j-1), j=1, ipar)
        write (ifi, form1) '----------', (zk16(ltirt+j-1), j=1, ipar)
!
        chain1 = ' '
        do i = 1, nbordt
! RECHERCHE DES NOMS DES PARAMETRES POUR LE NUMERO D'ORDRE COURANT, I
!
            chain2 = ' '
            ipcd = 1
            do j = 1, ipar
                ipcf = ipcd+15
                nopara = zk16(lnupa+j-1)
                call rsadpa(nomd2, 'L', 1, nopara, zi(lres+i-1), &
                            1, sjv=iad, styp=ctype, istop=0)
                if (ctype(1:1) .eq. 'I') then
                    chain2(ipcd:ipcf) = '       I        '
                else if (ctype(1:1) .eq. 'R') then
                    if (zr(iad) .ne. rundf) then
                        chain2(ipcd:ipcf) = '       R        '
                    end if
                else if (ctype(1:3) .eq. 'K80') then
                    chain2(ipcd:ipcf) = '      K80        '
                else if (ctype(1:3) .eq. 'K32') then
                    chain2(ipcd:ipcf) = '      K32        '
                else if (ctype(1:3) .eq. 'K24') then
                    chain2(ipcd:ipcf) = '      K24        '
                else if (ctype(1:3) .eq. 'K16') then
                    chain2(ipcd:ipcf) = '      K16        '
                else if (ctype(1:2) .eq. 'K8') then
                    chain2(ipcd:ipcf) = '       K8        '
                else
                    chain2(ipcd:ipcf) = blanc
                end if
                ipcd = ipcf+1
                chain2(ipcd:ipcd) = '!'
                ipcd = ipcd+1
            end do
!
! ECRITURE : ON ECRIT TOUJOURS LA PREMIERE ET LA DERNIERE LIGNE. AU
!            MILIEU, ON N'ECRIT QUE SI LE TEXTE A CHANGE.
! APRES CHAQUE ECRITURE, ON MEMORISE II, NUMERO D'ORDRE QUI A ETE ECRIT
! . 1ERE LIGNE : LA CHAINE COMPLETE
!
            if (i .eq. 1) then
                write (ifi, form2) zi(lres+i-1), chain2
                ii = 1
!
! . SI LE NOUVEAU TEXTE, CHAIN2, EST DIFFERENT DE CELUI DE LA LIGNE
!   PRECEDENTE, CHAIN1
! . OU SI C'EST LA DERNIERE LIGNE
!
            else if (chain1 .ne. chain2 .or. i .eq. nbordt) then
!          . ON VIENT JUSTE D'ECRIRE CHAIN1
                if (ii .eq. (i-1)) then
                    write (ifi, form2) zi(lres+i-1), chain2
!          . ON A ECRIT CHAIN1 DEUX NUMEROS AVANT : ON ECRIT LE
!          NUMERO PRECEDENT, I-1, ET LE COURANT, I.
                else if (ii .eq. (i-2)) then
                    write (ifi, form2) zi(lres+i-2), chain1
                    write (ifi, form2) zi(lres+i-1), chain2
!          . ON A ECRIT CHAIN1 PLUS DE DEUX NUMEROS AVANT : ON ECRIT
!          UNE LIGNE DE POINTILLES, LE NUMERO PRECEDENT, I-1, ET LE
!          NUMERO COURANT, I.
                else
                    if (chain1 .ne. chain2) then
                        if (ii .eq. (i-3)) then
                            write (ifi, form2) zi(lres+i-3), chain1
                        else
                            write (ifi, form1) '       ...', (zk16(lpoin+k-1), k=1, ipar)
                        end if
                        write (ifi, form2) zi(lres+i-2), chain1
                    else
                        write (ifi, form1) '       ...', (zk16(lpoin+k-1), k=1, ipar)
                    end if
                    write (ifi, form2) zi(lres+i-1), chain2
                end if
                ii = i
            end if
            chain1 = chain2
        end do
        write (ifi, form1) '----------', (zk16(ltirt+k-1), k=1, ipar)
        call jedetr('&&'//nompro//'.TIRET')
        call jedetr('&&'//nompro//'.POINT')
        call jedetr('&&'//nompro//'.NOM_PARA')
        call jedetr('&&'//nompro//'.NUM_PARA')
        call jedetr('&&'//nompro//'.PARA_EXIS')
    end if
!
! 2.6. ==> MENAGE PARTIEL
!
    call jedetr('&&'//nompro//'.NOMS_PARA')
    call jedetr('&&'//nompro//'.NUME_ORDRE')
!
! 3. LA FIN
!
999 continue
    call jedema()
end subroutine
