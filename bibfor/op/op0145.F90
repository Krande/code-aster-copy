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

subroutine op0145()
    implicit none
!-----------------------------------------------------------------------
!
!     OPERATEUR "DEFI_SPEC_TURB"
!
!-----------------------------------------------------------------------
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/reliem.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: ibid, dim, mxval, igrno, nbno, jnomno, jgrno
    character(len=1) :: typspe
    character(len=8) :: intspe, caelem, modele, nomzon, maillage
    character(len=16) :: concep, cmd, nommcf, mcfac(9)
    character(len=16) :: motcle(1), typmcl(1), nomno
    character(len=19) :: nomu
    character(len=24) :: vain, vare, vate, nnoe, chnumi
!
!-----------------------------------------------------------------------
    integer(kind=8) :: iangl, ifo, ifonct, iinter, imc, imcf, imci
    integer(kind=8) :: inat, inatur, inoeud, iocc, ispect
    integer(kind=8) :: lfon, lnat, lnom, long, lvain, lvare
    integer(kind=8) :: lvate, nbmcl
!-----------------------------------------------------------------------
    data mcfac/'SPEC_LONG_COR_1', 'SPEC_LONG_COR_2',&
     &            'SPEC_LONG_COR_3', 'SPEC_LONG_COR_4',&
     &            'SPEC_CORR_CONV_1', 'SPEC_CORR_CONV_2',&
     &            'SPEC_CORR_CONV_3',&
     &            'SPEC_FONC_FORME', 'SPEC_EXCI_POINT'/
    data motcle/'GROUP_NO'/
    data typmcl/'GROUP_NO'/
! ----------------------------------------------------------------------
    call jemarq()
    call infmaj()
!
    call getres(nomu, concep, cmd)
!
    do imcf = 1, 9
        call getfac(mcfac(imcf), iocc)
        if (iocc .eq. 1) exit
    end do
!   NBMCL EST AFFECTE A 12 , IL DOIT ETRE SUPERIEUR AU MAX DES NOMBRE
!   DE MC SIMPLES DES 9 MC FACTEURS
    nbmcl = 12
!
    nommcf = mcfac(imcf)
!
    read (nommcf(6:6), '(A1)') typspe
    if (typspe .eq. 'L') then
        read (nommcf(15:15), '(I1)') ispect
    else if (typspe .eq. 'F') then
        ispect = 11
    else if (typspe .eq. 'C') then
        read (nommcf(16:16), '(I1)') ispect
    else
        ispect = 21
    end if
!
! ----VERIFICATIONS AVANT EXECUTION----
!     =============================
!
    if (ispect .eq. 21) then
        call getvid(nommcf, 'INTE_SPEC', iocc=1, nbval=0, nbret=iinter)
        if (iinter .ne. 0) then
            call getvtx(nommcf, 'NATURE', iocc=1, nbval=0, nbret=inatur)
            call getvr8(nommcf, 'ANGLE', iocc=1, nbval=0, nbret=iangl)
            call getvtx(nommcf, 'NOEUD', iocc=1, nbval=0, nbret=inoeud)
            call getvtx(nommcf, 'GROUP_NO', iocc=1, nbval=0, nbret=igrno)
            if (abs(igrno) .gt. abs(inoeud)) inoeud = igrno
            if (inatur .ne. iangl .or. inatur .ne. inoeud .or. inoeud .ne. iangl) then
                call utmess('F', 'MODELISA5_66')
            end if
        else
            call getvtx(nommcf, 'NOEUD', iocc=1, nbval=0, nbret=inoeud)
            call getvtx(nommcf, 'GROUP_NO', iocc=1, nbval=0, nbret=igrno)
            if (abs(igrno) .gt. abs(inoeud)) inoeud = igrno
            if (abs(inoeud) .ne. 1) then
                call utmess('F', 'MODELISA5_67', sk=nommcf, si=inoeud)
            end if
        end if
    end if
!
! ----FIN DES VERIFICATIONS AVANT EXECUTION----
!     =====================================
!
!
! ----CREATION DES OBJETS ET REMPLISSAGE EN FONCTION DES----
! ----          DIFFERENTS TYPES DE SPECTRE             ----
!     ==================================================
!
! ----0.DENOMINATIONS DES OBJETS A CREER EN GLOBALE
!       -------------------------------------------
!
    vain = nomu//'.VAIN'
    vare = nomu//'.VARE'
    vate = nomu//'.VATE'
    nnoe = nomu//'.NNOE'
!
! ----VERIFICATIONS A L'EXECUTION----
!     ===========================
!
    nomno = '&&OP0145.GR_NO'
    if (ispect .eq. 11 .or. ispect .eq. 21) then
        call getvid(nommcf, 'INTE_SPEC', iocc=1, nbval=0, nbret=iinter)
        if (iinter .ne. 0) then
            call getvid(nommcf, 'INTE_SPEC', iocc=1, scal=intspe, nbret=ibid)
            chnumi = intspe//'.NUMI'
            call jelira(chnumi, 'LONMAX', mxval)
            if (ispect .eq. 11) then

                call getvid(nommcf, 'FONCTION', iocc=1, nbval=0, nbret=ifonct)
                dim = abs(ifonct)
                dim = dim*(dim+1)/2
                if (dim .ne. mxval) then
                    call utmess('F', 'MODELISA5_68')
                end if
            else
                call getvid(nommcf, 'MODELE', iocc=1, scal=modele, nbret=ibid)
                call dismoi('NOM_MAILLA', modele, 'MODELE', repk=maillage)
                call reliem(modele, maillage, 'NO_NOEUD', nommcf, iocc, &
                            1, motcle, typmcl, nomno, nbno)

                dim = nbno*(nbno+1)/2
                if (dim .ne. mxval) then
                    call utmess('F', 'MODELISA5_69')
                end if
            end if
        end if
    end if
!
! ----FIN DES VERIFICATIONS A L'EXECUTION----
!     ===================================

!
! ----1.MODELES "LONGUEUR DE CORRELATION"
!       ---------------------------------
!
    if (ispect .lt. 10 .or. nommcf(1:14) .eq. 'SPEC_CORR_CONV') then
!
! ------1.1.CREATION DES OBJETS SUR LA BASE GLOBALE
!
        call wkvect(vain, 'G V I', 1, lvain)
        call wkvect(vare, 'G V R', nbmcl, lvare)
        long = nbmcl+1
        call wkvect(vate, 'G V K16', long, lvate)
!
! ------1.2.CREATION D'OBJETS SUR LA BASE VOLATILE
!
        call wkvect('OP0145.TEMP.NOM', 'V V K16', nbmcl, lnom)
!
! ------1.3.REMPLISSAGE DES OBJETS
!
        if (nommcf .eq. 'SPEC_LONG_COR_1') then
            zk16(lnom) = 'LONG_COR        '
            zk16(lnom+1) = 'PROF_VITE_FLUI  '
            zk16(lnom+2) = 'VISC_CINE       '
        else if (nommcf .eq. 'SPEC_LONG_COR_2') then
            zk16(lnom) = 'LONG_COR        '
            zk16(lnom+1) = 'PROF_VITE_FLUI  '
            zk16(lnom+2) = 'FREQ_COUP       '
            zk16(lnom+3) = 'PHI0            '
            zk16(lnom+4) = 'BETA            '
        else if (nommcf .eq. 'SPEC_LONG_COR_3') then
            zk16(lnom) = 'LONG_COR        '
            zk16(lnom+1) = 'PROF_VITE_FLUI  '
            zk16(lnom+2) = 'FREQ_COUP       '
            zk16(lnom+3) = 'PHI0_1          '
            zk16(lnom+4) = 'BETA_1          '
            zk16(lnom+5) = 'PHI0_2          '
            zk16(lnom+6) = 'BETA_2          '
        else if (nommcf .eq. 'SPEC_LONG_COR_4') then
            zk16(lnom) = 'LONG_COR        '
            zk16(lnom+1) = 'PROF_VITE_FLUI  '
            zk16(lnom+2) = 'TAUX_VIDE       '
            zk16(lnom+3) = 'BETA            '
            zk16(lnom+4) = 'GAMMA           '
        else if (nommcf .eq. 'SPEC_CORR_CONV_1') then
            zk16(lnom) = 'LONG_COR_1      '
            zk16(lnom+1) = 'LONG_COR_2      '
            zk16(lnom+2) = 'VITE_FLUI       '
            zk16(lnom+3) = 'RHO_FLUI        '
            zk16(lnom+4) = 'FREQ_COUP       '
            zk16(lnom+5) = 'K               '
            zk16(lnom+6) = 'D_FLUI          '
            zk16(lnom+7) = 'COEF_VITE_FLUI_A'
            zk16(lnom+8) = 'COEF_VITE_FLUI_O'
            zk16(lnom+9) = 'METHODE         '
        else if (nommcf .eq. 'SPEC_CORR_CONV_2') then
            zk16(lnom) = 'FONCTION        '
            zk16(lnom+1) = 'LONG1_F         '
            zk16(lnom+2) = 'LONG2_F         '
            zk16(lnom+3) = 'METHODE         '
            zk16(lnom+4) = 'VITE_FLUI       '
            zk16(lnom+5) = 'COEF_VITE_FLUI_O'
            zk16(lnom+6) = 'FREQ_COUP       '
            zk16(lnom+7) = 'COEF_VITE_FLUI_A'
        else if (nommcf .eq. 'SPEC_CORR_CONV_3') then
            zk16(lnom) = 'TABLE_FONCTION  '
        else if (nommcf .eq. 'SPEC_FONC_FORME') then
            zk16(lnom) = 'INTE_SPEC       '
            zk16(lnom+1) = 'FONCTION        '
            zk16(lnom+2) = 'GRAPPE_1        '
            zk16(lnom+3) = 'NOEUD           '
            zk16(lnom+4) = 'CARA_ELEM       '
            zk16(lnom+5) = 'MODELE          '
        else if (nommcf .eq. 'SPEC_EXCI_POINT') then
            zk16(lnom) = 'INTE_SPEC      '
            zk16(lnom+1) = 'NATURE         '
            zk16(lnom+2) = 'ANGL           '
            zk16(lnom+3) = 'GRAPPE_2       '
            zk16(lnom+4) = 'RHO_FLUI       '
            zk16(lnom+5) = 'NOEUD          '
            zk16(lnom+6) = 'CARA_ELEM      '
            zk16(lnom+7) = 'MODELE         '
        end if

        zk16(lvate) = nommcf
        imci = 0
        do imc = 1, nbmcl
            if (zk16(lnom+imc-1) .eq. 'PROF_VITE_FLUI  ') then
                call getvid(nommcf, 'PROF_VITE_FLUI', iocc=1, scal=nomzon, nbret=ibid)
                zk16(lvate+imc) = nomzon
            else if (zk16(lnom+imc-1) .eq. 'METHODE         ') then
                call getvtx(nommcf, 'METHODE', iocc=1, scal=nomzon, nbret=ibid)
                zk16(lvate+imc) = nomzon
            else if (zk16(lnom+imc-1) .eq. 'FONCTION        ') then
                call getvid(nommcf, 'FONCTION', iocc=1, scal=nomzon, nbret=ibid)
                zk16(lvate+imc) = nomzon
            else if (zk16(lnom+imc-1) .eq. 'LONG1_F         ') then
                call getvid(nommcf, 'LONG1_F', iocc=1, scal=nomzon, nbret=ibid)
                zk16(lvate+imc) = nomzon
            else if (zk16(lnom+imc-1) .eq. 'LONG2_F         ') then
                call getvid(nommcf, 'LONG2_F', iocc=1, scal=nomzon, nbret=ibid)
                zk16(lvate+imc) = nomzon
            else if (zk16(lnom+imc-1) .eq. 'TABLE_FONCTION  ') then
                call getvid(nommcf, 'TABLE_FONCTION', iocc=1, scal=nomzon, nbret=ibid)
                zk16(lvate+imc) = nomzon
            else if (zk16(lnom+imc-1) .ne. '                ') then
                imci = imci+1
                zk16(lvate+imc) = zk16(lnom+imc-1)
                call getvr8(nommcf, zk16(lnom+imc-1), iocc=1, scal=zr(lvare+imci-1), nbret=ibid)
            end if
        end do
        zi(lvain) = ispect
!
! ----2.MODELES "FONCTIONS DE FORME" ET "EXCITATIONS PONCTUELLES"
!       ---------------------------------------------------------
!
    else
!
! ------2.1.CREATION ET REMPLISSAGE D'OBJETS COMMUNS
!
! ------2.1.1.OBJET .VAIN
!
        call wkvect(vain, 'G V I', 3, lvain)
        call getvid(nommcf, 'MODELE', iocc=1, scal=modele, nbret=ibid)
        call dismoi('NOM_MAILLA', modele, 'MODELE', repk=maillage)
        call reliem(modele, maillage, 'NO_NOEUD', nommcf, iocc, &
                    1, motcle, typmcl, nomno, nbno)
        call jeveuo(nomno, 'L', jnomno)
!
        if (nbno .ne. 1) then
            call utmess('F', 'MODELISA5_67', sk=nommcf, si=nbno)
        end if
!
        zi(lvain) = ispect
        if (iinter .eq. 0) then
            zi(lvain+1) = 1
        else
            zi(lvain+1) = 0
        end if
        zi(lvain+2) = nbno
!
! ------2.1.2.OBJET .NNOE
!
        call wkvect(nnoe, 'G V K8', nbno, jgrno)
!
        call jeveuo(nnoe, 'E', jgrno)
        zk8(jgrno) = zk8(jnomno)

!
! ------2.2.OBJETS .VATE ET .VARE
!
        long = 5
        if (iinter .ne. 0) then
            if (ispect .eq. 11) then
                ifonct = abs(ifonct)
                long = 4+ifonct
            else
                inoeud = abs(inoeud)
                long = 4+inoeud
            end if
        end if
!
        call wkvect(vate, 'G V K16', long, lvate)
        call getvid(nommcf, 'CARA_ELEM', iocc=1, scal=caelem, nbret=ibid)
        call getvid(nommcf, 'MODELE', iocc=1, scal=modele, nbret=ibid)
        zk16(lvate) = nommcf
        zk16(lvate+1) = caelem
        zk16(lvate+2) = modele

!
! ------2.2.1.MODELE "FONCTIONS DE FORME"
!
        if (ispect .eq. 11) then
            if (iinter .eq. 0) then
                zk16(lvate+3) = 'GRAPPE_1'
                call getvtx(nommcf, 'GRAPPE_1', iocc=1, scal=zk16(lvate+4), nbret=ibid)
            else
                call wkvect('OP0145.TEMP.FON', 'V V K8', ifonct, lfon)
                call getvid(nommcf, 'FONCTION', iocc=1, nbval=ifonct, vect=zk8(lfon), &
                            nbret=ibid)
                zk16(lvate+3) = intspe
                do ifo = 1, ifonct
                    zk16(lvate+3+ifo) = zk8(lfon+ifo-1)
                end do
            end if
!
! ------2.2.2.MODELE "EXCITATIONS PONCTUELLES"
!
        else
            if (iinter .eq. 0) then
                zk16(lvate+3) = 'GRAPPE_2'
                call getvtx(nommcf, 'GRAPPE_2', iocc=1, scal=zk16(lvate+4), nbret=ibid)
!
                call wkvect(vare, 'G V R', 1, lvare)
                call getvr8(nommcf, 'RHO_FLUI', iocc=1, scal=zr(lvare), nbret=ibid)
!
            else
                call wkvect('OP0145.TEMP.NAT', 'V V K8', inoeud, lnat)
                call getvtx(nommcf, 'NATURE', iocc=1, nbval=inoeud, vect=zk8(lnat), &
                            nbret=ibid)
                zk16(lvate+3) = intspe
                do inat = 1, inoeud
                    zk16(lvate+3+inat) = zk8(lnat+inat-1)
                end do
!
                call wkvect(vare, 'G V R', inoeud, lvare)
                call getvr8(nommcf, 'ANGLE', iocc=1, nbval=inoeud, vect=zr(lvare), &
                            nbret=ibid)
!
            end if
!
        end if
!
    end if
!
    call titre()
!
    call jedetr(nomno)
    call jedema()
end subroutine
