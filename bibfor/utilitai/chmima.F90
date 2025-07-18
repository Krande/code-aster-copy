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

subroutine chmima(nomsd, nomsy, typcha, typmax, nocham, typresu, &
                  mcfz)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cnocns.h"
#include "asterfort/cnscno.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/posddl.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsutnu.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/idensd.h"
#include "asterfort/int_to_char8.h"
!
    integer(kind=8) :: nbordr
    character(len=*) :: nomsd, nomsy, typmax, nocham, typcha
    character(len=8), optional :: typresu
    character(len=*), optional :: mcfz
!      AFFECTATION DU CHAMP-GD DE NOM NOCHAM  AVEC LES
!      VALEURS MINMAX EN TOUT POINT DES CHAMPS-GD DE TYPE
!      NOMSY DU RESULTAT DE NOM NOMSD
! ----------------------------------------------------------------------
! IN  : NOMSD  : NOM DE LA STRUCTURE "RESULTAT"
! IN  : NOMSY  : NOM SYMBOLIQUE DU CHAMP A CHERCHER.
! IN  : TYPCHA : TYPE DU CHAMP A CHERCHER
! IN  : TYPMAX : TYPE D'OPERATION A EFFECTUER
! VAR : NOCHAM : NOM DU CHAMP CONTENANT LES MINMAX DES
!                CHAMPS DE TYPE NOMSY DU RESULTAT NOMSD.
!
! ----------------------------------------------------------------------
    character(len=4) :: ctyp
    character(len=8) :: typma, crit, noma, nomn, valeur
    character(len=19) :: prno, prn2
    character(len=16) :: noms2
    character(len=19) :: nocha2, chextr, knum, chs1, chs2
    character(len=19) :: mcf
    character(len=5) :: sufv, sufsl
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iad, in, inoe, inumer, iocc
    integer(kind=8) :: iret, ivale, j, jddlx, jddly, jddlz, jdlrx
    integer(kind=8) :: jdlry, jdlrz, jordr, jvpnt, n2, nbnoe, nc
    integer(kind=8) :: neq, np, nvale, neq2, icsl
    real(kind=8) :: epsi, rs1, x, y, z
    logical :: verif
!-----------------------------------------------------------------------
    call jemarq()
    knum = '&&CHMIMA.NUME_ORDRE'
    noms2 = nomsy
    nocha2 = nocham
    typma = typmax
    chs1 = '&&CHMIMA.CHS1'
    chs2 = '&&CHMIMA.CHS2'
!
!     --- LECTURE DU MOT-CLE TYPE_RESU ---
!
    if (present(typresu)) then
        valeur = typresu
    else
        call getvtx(' ', 'TYPE_RESU', scal=valeur, nbret=n2)
    end if
!
!     --- RECUPERATION DES NUMEROS D'ORDRE ---
!
    if (present(mcfz)) then
        mcf = mcfz
        iocc = 1
    else
        mcf = ' '
        iocc = 0
    end if

    call getvr8(mcf, 'PRECISION', iocc=iocc, scal=epsi, nbret=np)
    call getvtx(mcf, 'CRITERE', iocc=iocc, scal=crit, nbret=nc)
!
    call rsutnu(nomsd, mcf, 1, knum, nbordr, &
                epsi, crit, iret)
    if (nbordr .eq. 0) then
        call utmess('F', 'UTILITAI_23')
    end if
    call jeveuo(knum, 'L', jordr)
!
!     --- TRAITEMENT DU PREMIER NUMERO D'ORDRE ---
!
    call rsexch('F', nomsd, noms2, zi(jordr), chextr, &
                iret)

    if (typcha(1:4) .eq. 'NOEU' .and. typma .ne. 'NORM_TRA') then
        call cnocns(chextr, 'V', chs1)
        sufv = '.CNSV'
    else
        call copisd('CHAMP_GD', 'G', chextr(1:19), nocha2(1:19))
        chs1 = nocha2
        call jeexin(nocha2(1:19)//'.VALE', iret)
        if (iret .gt. 0) then
            sufv = '.VALE'
        else
            sufv = '.CELV'
        end if
    end if

    call jelira(chs1(1:19)//sufv, 'LONMAX', neq)
    call jeveuo(chs1(1:19)//sufv, 'E', nvale)
!
    call wkvect('&&CHMIMA.INST', 'V V I', neq, inumer)
    do i = 1, neq
        zi(inumer+i-1) = zi(jordr)
    end do
!
!     --- BOUCLE SUR LES NUMEROS D'ORDRE ---
!
    if (typma .eq. 'MAXI') then
!
        do i = 2, nbordr
!
!         - RECUPERATION DU CHAMP DE TYPE NOMSY
!           CORRESPONDANT AU NUMERO D'ORDRE COURANT
!
            call rsexch('F', nomsd, noms2, zi(jordr+i-1), chextr, &
                        iret)
            if (typcha(1:4) .eq. 'NOEU') then
                call cnocns(chextr, 'V', chs2)
                sufsl = '.CNSL'
            else
                chs2 = chextr
                sufsl = '.CESL'
            end if
!           verification pour les cham_elem
            call jelira(chs2(1:19)//sufv, 'LONMAX', neq2)
            if (neq2 .ne. neq) call utmess('F', 'UTILITAI_25')
!
!         - RECUPERATION DU VALE DU CHAMP EXTRAIT
!
            call jeveuo(chs2//sufv, 'L', ivale)
            call jeexin(chs2//sufsl, iret)
            if (iret .eq. 0) then
                verif = .false.
            else
                verif = .true.
                call jeveuo(chs2//sufsl, 'L', icsl)
            end if
!
            if (verif) then
                do j = 1, neq
                    if (zl(icsl+j-1)) then
                        if (zr(ivale+j-1) .gt. zr(nvale+j-1)) then
                            zr(nvale+j-1) = zr(ivale+j-1)
                            zi(inumer+j-1) = zi(jordr+i-1)
                        end if
                    end if
                end do
            else
                do j = 1, neq
                    if (zr(ivale+j-1) .gt. zr(nvale+j-1)) then
                        zr(nvale+j-1) = zr(ivale+j-1)
                        zi(inumer+j-1) = zi(jordr+i-1)
                    end if
                end do
            end if

            if (typcha(1:4) .eq. 'NOEU') call detrsd('CHAM_NO_S', chs2)
!
        end do
!
    else if (typma .eq. 'MAXI_ABS') then
!
        do i = 2, nbordr
!
!         - RECUPERATION DU CHAMP DE TYPE NOMSY
!           CORRESPONDANT AU NUMERO D'ORDRE COURANT
!
            call rsexch('F', nomsd, noms2, zi(jordr+i-1), chextr, &
                        iret)
            if (typcha(1:4) .eq. 'NOEU') then
                call cnocns(chextr, 'V', chs2)
                sufsl = '.CNSL'
            else
                chs2 = chextr
                sufsl = '.CESL'
            end if
!           verification pour les cham_elem
            call jelira(chs2(1:19)//sufv, 'LONMAX', neq2)
            if (neq2 .ne. neq) call utmess('F', 'UTILITAI_25')
!
!         - RECUPERATION DU VALE DU CHAMP EXTRAIT
!
            call jeveuo(chs2//sufv, 'L', ivale)
            call jeexin(chs2//sufsl, iret)
            if (iret .eq. 0) then
                verif = .false.
            else
                verif = .true.
                call jeveuo(chs2//sufsl, 'L', icsl)
            end if
!
!         - L'OBJET .CNSL EST UTILISE POUR EVITER DE TESTER DES VALEURS NaN
!
            if (valeur .eq. 'VALE_ABS') then
                if (verif) then
                    do j = 1, neq
                        if (i .eq. 2) zr(nvale+j-1) = abs(zr(nvale+j-1))
                        if (zl(icsl+j-1)) then
                            if (abs(zr(ivale+j-1)) .gt. zr(nvale+j-1)) then
                                zr(nvale+j-1) = abs(zr(ivale+j-1))
                                zi(inumer+j-1) = zi(jordr+i-1)
                            end if
                        end if
                    end do
                else
                    do j = 1, neq
                        if (i .eq. 2) zr(nvale+j-1) = abs(zr(nvale+j-1))
                        if (abs(zr(ivale+j-1)) .gt. zr(nvale+j-1)) then
                            zr(nvale+j-1) = abs(zr(ivale+j-1))
                            zi(inumer+j-1) = zi(jordr+i-1)
                        end if
                    end do
                end if
            else
                if (verif) then
                    do j = 1, neq
                        if (zl(icsl+j-1)) then
                            if (abs(zr(ivale+j-1)) .gt. abs(zr(nvale+j-1))) then
                                zr(nvale+j-1) = zr(ivale+j-1)
                                zi(inumer+j-1) = zi(jordr+i-1)
                            end if
                        end if
                    end do
                else
                    do j = 1, neq
                        if (abs(zr(ivale+j-1)) .gt. abs(zr(nvale+j-1))) then
                            zr(nvale+j-1) = zr(ivale+j-1)
                            zi(inumer+j-1) = zi(jordr+i-1)
                        end if
                    end do
                end if
            end if

            if (typcha(1:4) .eq. 'NOEU') call detrsd('CHAM_NO_S', chs2)
!
        end do
!
    else if (typma .eq. 'MINI') then
!
        do i = 2, nbordr
!
!         - RECUPERATION DU CHAMP DE TYPE NOMSY
!           CORRESPONDANT AU NUMERO D'ORDRE COURANT
!if (typcha(1:4).eq.'NOEU' .and. typma.eq.'NORM_TRA')then
            call rsexch('F', nomsd, noms2, zi(jordr+i-1), chextr, &
                        iret)
            if (typcha(1:4) .eq. 'NOEU') then
                call cnocns(chextr, 'V', chs2)
                sufsl = '.CNSL'
            else
                chs2 = chextr
                sufsl = '.CESL'
            end if
!           verification pour les cham_elem
            call jelira(chs2(1:19)//sufv, 'LONMAX', neq2)
            if (neq2 .ne. neq) call utmess('F', 'UTILITAI_25')
!
!         - RECUPERATION DU VALE DU CHAMP EXTRAIT
!
            call jeveuo(chs2//sufv, 'L', ivale)
            call jeexin(chs2//sufsl, iret)
            if (iret .eq. 0) then
                verif = .false.
            else
                verif = .true.
                call jeveuo(chs2//sufsl, 'L', icsl)
            end if
!
            if (verif) then
                do j = 1, neq
                    if (zl(icsl+j-1)) then
                        if (zr(ivale+j-1) .lt. zr(nvale+j-1)) then
                            zr(nvale+j-1) = zr(ivale+j-1)
                            zi(inumer+j-1) = zi(jordr+i-1)
                        end if
                    end if
                end do
            else
                do j = 1, neq
                    if (zr(ivale+j-1) .lt. zr(nvale+j-1)) then
                        zr(nvale+j-1) = zr(ivale+j-1)
                        zi(inumer+j-1) = zi(jordr+i-1)
                    end if
                end do
            end if
            if (typcha(1:4) .eq. 'NOEU') call detrsd('CHAM_NO_S', chs2)
!
        end do
!
    else if (typma .eq. 'MINI_ABS') then
!
        do i = 2, nbordr
!
!         - RECUPERATION DU CHAMP DE TYPE NOMSY
!           CORRESPONDANT AU NUMERO D'ORDRE COURANT
!
            call rsexch('F', nomsd, noms2, zi(jordr+i-1), chextr, &
                        iret)
            if (typcha(1:4) .eq. 'NOEU') then
                call cnocns(chextr, 'V', chs2)
                sufsl = '.CNSL'
            else
                chs2 = chextr
                sufsl = '.CESL'
            end if
!           verification pour les cham_elem
            call jelira(chs2(1:19)//sufv, 'LONMAX', neq2)
            if (neq2 .ne. neq) call utmess('F', 'UTILITAI_25')
!
!         - RECUPERATION DU VALE DU CHAMP EXTRAIT
!
            call jeveuo(chs2//sufv, 'L', ivale)
            call jeexin(chs2//sufsl, iret)
            if (iret .eq. 0) then
                verif = .false.
            else
                verif = .true.
                call jeveuo(chs2//sufsl, 'L', icsl)
            end if
!
            if (valeur .eq. 'VALE_ABS') then
                if (verif) then
                    do j = 1, neq
                        if (i .eq. 2) zr(nvale+j-1) = abs(zr(nvale+j-1))
                        if (zl(icsl+j-1)) then
                            if (abs(zr(ivale+j-1)) .lt. zr(nvale+j-1)) then
                                zr(nvale+j-1) = abs(zr(ivale+j-1))
                                zi(inumer+j-1) = zi(jordr+i-1)
                            end if
                        end if
                    end do
                else
                    do j = 1, neq
                        if (i .eq. 2) zr(nvale+j-1) = abs(zr(nvale+j-1))
                        if (abs(zr(ivale+j-1)) .lt. zr(nvale+j-1)) then
                            zr(nvale+j-1) = abs(zr(ivale+j-1))
                            zi(inumer+j-1) = zi(jordr+i-1)
                        end if
                    end do
                end if
            else
                if (verif) then
                    do j = 1, neq
                        if (zl(icsl+j-1)) then
                            if (abs(zr(ivale+j-1)) .lt. abs(zr(nvale+j-1))) then
                                zr(nvale+j-1) = zr(ivale+j-1)
                                zi(inumer+j-1) = zi(jordr+i-1)
                            end if
                        end if
                    end do
                else
                    do j = 1, neq
                        if (abs(zr(ivale+j-1)) .lt. abs(zr(nvale+j-1))) then
                            zr(nvale+j-1) = zr(ivale+j-1)
                            zi(inumer+j-1) = zi(jordr+i-1)
                        end if
                    end do
                end if
            end if

            if (typcha(1:4) .eq. 'NOEU') call detrsd('CHAM_NO_S', chs2)
!
        end do
!
    else if (typma .eq. 'NORM_TRA') then
        call rsexch('F', nomsd, noms2, zi(jordr), chextr, &
                    iret)
        call jeveuo(chextr//'.VALE', 'L', ivale)
        call dismoi('NUME_EQUA', chextr, 'CHAM_NO', repk=prno)
        call dismoi('NOM_MAILLA', chextr, 'CHAM_NO', repk=noma)
        call dismoi('NB_NO_MAILLA', noma, 'MAILLAGE', repi=nbnoe)
!
        do j = 0, neq-1
            zr(nvale+j) = zr(ivale+j)
        end do
        if (nbordr .eq. 1) goto 58
!
        call wkvect('&&CHMIMA.DDL.DX', 'V V I', nbnoe, jddlx)
        call wkvect('&&CHMIMA.DDL.DY', 'V V I', nbnoe, jddly)
        call wkvect('&&CHMIMA.DDL.DZ', 'V V I', nbnoe, jddlz)
        call wkvect('&&CHMIMA.DDL.DRX', 'V V I', nbnoe, jdlrx)
        call wkvect('&&CHMIMA.DDL.DRY', 'V V I', nbnoe, jdlry)
        call wkvect('&&CHMIMA.DDL.DRZ', 'V V I', nbnoe, jdlrz)
        call wkvect('&&CHMIMA.VALE_P.NT', 'V V R', nbnoe, jvpnt)
!
        do in = 0, nbnoe-1
            nomn = int_to_char8(in+1)
            call posddl('CHAM_NO', chextr, nomn, 'DX', inoe, &
                        zi(jddlx+in))
            call posddl('CHAM_NO', chextr, nomn, 'DY', inoe, &
                        zi(jddly+in))
            call posddl('CHAM_NO', chextr, nomn, 'DZ', inoe, &
                        zi(jddlz+in))
            call posddl('CHAM_NO', chextr, nomn, 'DRX', inoe, &
                        zi(jdlrx+in))
            call posddl('CHAM_NO', chextr, nomn, 'DRY', inoe, &
                        zi(jdlry+in))
            call posddl('CHAM_NO', chextr, nomn, 'DRZ', inoe, &
                        zi(jdlrz+in))
            x = zr(ivale+zi(jddlx+in)-1)
            y = zr(ivale+zi(jddly+in)-1)
            if (zi(jddlz+in) .ne. 0) then
                z = zr(ivale+zi(jddlz+in)-1)
            else
                z = 0.d0
            end if
            zr(jvpnt+in) = sqrt(x**2+y**2+z**2)
            zr(nvale+zi(jddlx+in)-1) = x
            zr(nvale+zi(jddly+in)-1) = y
            if (zi(jddlz+in) .ne. 0) zr(nvale+zi(jddlz+in)-1) = z
            zi(inumer+zi(jddlx+in)-1) = zi(jordr)
            zi(inumer+zi(jddly+in)-1) = zi(jordr)
            if (zi(jddlz+in) .ne. 0) zi(inumer+zi(jddlz+in)-1) = zi(jordr)
            if (zi(jdlrx+in) .ne. 0) zr(nvale+zi(jdlrx+in)-1) = zr(ivale+zi(jdlrx+in)-1)
            if (zi(jdlry+in) .ne. 0) zr(nvale+zi(jdlry+in)-1) = zr(ivale+zi(jdlry+in)-1)
            if (zi(jdlrz+in) .ne. 0) zr(nvale+zi(jdlrz+in)-1) = zr(ivale+zi(jdlrz+in)-1)
        end do
!
        do i = 2, nbordr
            call rsexch('F', nomsd, noms2, zi(jordr+i-1), chextr, &
                        iret)
            call dismoi('NUME_EQUA', chextr, 'CHAM_NO', repk=prn2)
            if (prn2 .ne. prno) then
                if (.not. idensd('NUME_EQUA', prn2, prno)) then
                    call utmess('F', 'UTILITAI_26')
                end if
            end if

            call jeveuo(chextr//'.VALE', 'L', ivale)
!
            do in = 0, nbnoe-1
                x = zr(ivale+zi(jddlx+in)-1)
                y = zr(ivale+zi(jddly+in)-1)
                if (zi(jddlz+in) .ne. 0) then
                    z = zr(ivale+zi(jddlz+in)-1)
                else
                    z = 0.d0
                end if
                rs1 = sqrt(x**2+y**2+z**2)
                if (rs1 .gt. zr(jvpnt+in)) then
                    zr(jvpnt+in) = rs1
                    zi(inumer+zi(jddlx+in)-1) = zi(jordr+i-1)
                    zi(inumer+zi(jddly+in)-1) = zi(jordr+i-1)
                    if (zi(jddlz+in) .ne. 0) zi(inumer+zi(jddlz+in)-1) = zi(jordr+i-1)
                    zr(nvale+zi(jddlx+in)-1) = x
                    zr(nvale+zi(jddly+in)-1) = y
                    if (zi(jddlz+in) .ne. 0) zr(nvale+zi(jddlz+in)-1) = z
                    if (zi(jdlrx+in) .ne. 0) zr(nvale+zi(jdlrx+in)-1) = zr(ivale+zi(jdlrx+in)-1)
                    if (zi(jdlry+in) .ne. 0) zr(nvale+zi(jdlry+in)-1) = zr(ivale+zi(jdlry+in)-1)
                    if (zi(jdlrz+in) .ne. 0) zr(nvale+zi(jdlrz+in)-1) = zr(ivale+zi(jdlrz+in)-1)
                end if
            end do
!
        end do
        call jedetr('&&CHMIMA.VALE_P.NT')
!
    end if
!
58  continue
    if (valeur(1:4) .eq. 'INST') then
        if (typma .eq. 'NORM_TRA') then
            if (nbordr .ne. 1) then
                do in = 0, nbnoe-1
                    call rsadpa(nomsd, 'L', 1, 'INST', zi(inumer+zi(jddlx+in)-1), &
                                0, sjv=iad, styp=ctyp)
                    zr(nvale+zi(jddlx+in)-1) = zr(iad)
                    zr(nvale+zi(jddly+in)-1) = zr(iad)
                    if (zi(jddlz+in) .ne. 0) zr(nvale+zi(jddlz+in)-1) = zr(iad)
                    if (zi(jdlrx+in) .ne. 0) zr(nvale+zi(jdlrx+in)-1) = zr(iad)
                    if (zi(jdlry+in) .ne. 0) zr(nvale+zi(jdlry+in)-1) = zr(iad)
                    if (zi(jdlrz+in) .ne. 0) zr(nvale+zi(jdlrz+in)-1) = zr(iad)
                end do
                call jedetr('&&CHMIMA.DDL.DX')
                call jedetr('&&CHMIMA.DDL.DY')
                call jedetr('&&CHMIMA.DDL.DZ')
                call jedetr('&&CHMIMA.DDL.DRX')
                call jedetr('&&CHMIMA.DDL.DRY')
                call jedetr('&&CHMIMA.DDL.DRZ')
            else
                do j = 0, neq-1
                    call rsadpa(nomsd, 'L', 1, 'INST', zi(inumer+j), &
                                0, sjv=iad, styp=ctyp)
                    zr(nvale+j) = zr(iad)
                end do
            end if
        else
            do j = 0, neq-1
                call rsadpa(nomsd, 'L', 1, 'INST', zi(inumer+j), &
                            0, sjv=iad, styp=ctyp)
                zr(nvale+j) = zr(iad)
            end do
        end if
    else
        if (typma .eq. 'NORM_TRA') then
            call jedetr('&&CHMIMA.DDL.DX')
            call jedetr('&&CHMIMA.DDL.DY')
            call jedetr('&&CHMIMA.DDL.DZ')
            call jedetr('&&CHMIMA.DDL.DRX')
            call jedetr('&&CHMIMA.DDL.DRY')
            call jedetr('&&CHMIMA.DDL.DRZ')
        end if
    end if
!
    call jedetr('&&CHMIMA.INST')
    call jedetr(knum)

    if (typcha(1:4) .eq. 'NOEU' .and. typma .ne. 'NORM_TRA') then
        call cnscno(chs1, ' ', 'NON', 'G', nocha2, &
                    'F', iret)
        ASSERT(iret .eq. 0)
        call detrsd('CHAM_NO_S', chs1)
    end if
!
    call jedema()
end subroutine
