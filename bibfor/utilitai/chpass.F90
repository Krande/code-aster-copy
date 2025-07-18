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

subroutine chpass(tychr, ma, celmod, nomgd, prol0, &
                  chou)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/assert.h"
#include "asterfort/carces.h"
#include "asterfort/celces.h"
#include "asterfort/cescar.h"
#include "asterfort/cescel.h"
#include "asterfort/cesred.h"
#include "asterfort/chsfus.h"
#include "asterfort/chsut1.h"
#include "asterfort/cnocns.h"
#include "asterfort/cnscno.h"
#include "asterfort/cnsred.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/getvc8.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/reliem.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
!
    character(len=8) :: ma, chou
    character(len=*) :: celmod
    character(len=4) :: tychr, tych2
! person_in_charge: jacques.pellet at edf.fr
!     BUT : TRAITER :
!          - OPTION 'ASSE' DE LA COMMANDE CREA_CHAMP
!     -----------------------------------------------------------------
!
!
    integer(kind=8) :: n1, ib, nbocc, iocc, nbtrou, jnutro, nbmocl, lnom, ibid
    aster_logical :: chgcmp, cumul, lcumul(2), compOK, messConseil
    integer(kind=8) :: ncmp, jlicmp, gd, jcmpgd, iret, nncp, nchg
    integer(kind=8) :: jcesd, jcesc, i, ncmpdisp, j
    real(kind=8) :: coefr, lcoefr(2)
    complex(kind=8) :: coefc, lcoefc(2)
    character(len=8) :: kbid, modele, nomcmp
    character(len=8) :: nomgd, nomgd2, ma2
    character(len=3) :: prol0, tsca
    character(len=16) :: limocl(5), tymocl(5), typem
    character(len=19) :: chs1, chs2, nutrou, lichs(2), cesmod, option, champ
    character(len=19) :: chs3, ligrel
    character(len=24) :: cnom, valk(3)
!
    aster_logical :: lcoc, bool(1), iprem
    character(len=8), pointer :: licmp2(:) => null()
    character(len=8), pointer :: licmpdisp(:) => null()
!     -----------------------------------------------------------------
!
    call jemarq()
!
    chs1 = '&&CHPASS.CHS1'
    chs2 = '&&CHPASS.CHS2'
    chs3 = '&&CHPASS.CHS3'

    messConseil = ASTER_FALSE
!
!
! 1- CALCUL DE:
!      NOMGD : NOM DE LA GRANDEUR
!      GD    : NUMERO DE LA GRANDEUR
!      JCMPGD: ADRESSE DES CMPS DE LA GRANDEUR
!      PROL0 :/'OUI' POUR PROLONGER PAR ZERO LES CHAM_ELEM
!                    NON DEFINIS PARTOUT
!             /'NON' POUR ARRETER <F> DANS LE CAS PRECEDENT
!      LIGREL: NOM DU LIGREL ASSOCIE A CELMOD
!      OPTION: OPTION ASSOCIEE A CELMOD (SI CHAM_ELEM)
!      CESMOD: CHAM_ELEM_S EQUIVALENT A CELMOD (SI CHAM_ELEM)
!
!     ------------------------------------------------------------------
!
    if (ma .eq. ' ') then
        call utmess('F', 'UTILITAI_27')
    end if
!
    call jenonu(jexnom('&CATA.GD.NOMGD', nomgd), gd)
    if (gd .eq. 0) then
        call utmess('F', 'CALCULEL_67', sk=nomgd)
    end if
    call jeveuo(jexnum('&CATA.GD.NOMCMP', gd), 'L', jcmpgd)
!
    if (tychr(1:2) .eq. 'EL') then
        call dismoi('NOM_LIGREL', celmod, 'CHAM_ELEM', repk=ligrel)
        call dismoi('NOM_OPTION', celmod, 'CHAM_ELEM', repk=option)
        cesmod = '&&CHPASS.CESMOD'
        call celces(celmod, 'V', cesmod)
!
        modele = ligrel(1:8)
        call exisd('MODELE', modele, iret)
        if (iret .ne. 1) modele = ' '
!
        call jeveuo(cesmod//'.CESD', 'L', jcesd)
        call jeveuo(cesmod//'.CESC', 'L', jcesc)
        ncmpdisp = zi(jcesd+1)
        AS_ALLOCATE(vk8=licmpdisp, size=ncmpdisp)
        licmpdisp(1:ncmpdisp) = zk8(jcesc:jcesc-1+ncmpdisp)
!
    else
        ligrel = ' '
        option = ' '
        cesmod = ' '
        modele = ' '
        ncmpdisp = 0
    end if
!
!
!     2- BOUCLE DE VERIF SUR LES OCCURENCES DU MOT CLE "ASSE" :
!     ---------------------------------------------------------
    call getfac('ASSE', nbocc)
    cnom = '&&CHPASS.CHAM_GD_LISTE'
    call wkvect(cnom, 'V V K24', nbocc, lnom)
!
    do iocc = 1, nbocc
!
!       2.1 VERIFICATION DES CARACTERISTIQUES DU CHAMP :
!       ------------------------------------------------
        call getvid('ASSE', 'CHAM_GD', iocc=iocc, scal=champ, nbret=ib)
        zk24(lnom+iocc-1) = champ
!
        call dismoi('TYPE_CHAMP', champ, 'CHAMP', repk=tych2)
        call dismoi('NOM_MAILLA', champ, 'CHAMP', repk=ma2)
        if (ma .ne. ma2) then
            valk(1) = champ
            valk(2) = ma2
            valk(3) = ma
            call utmess('F', 'CALCULEL2_17', nk=3, valk=valk)
        end if
!
!
        valk(1) = champ
        valk(2) = tychr
        if (tychr .eq. 'NOEU') then
            if (tych2 .ne. 'NOEU') then
                call utmess('F', 'UTILITAI_28', nk=2, valk=valk)
            end if
!
        else if (tychr .eq. 'ELGA') then
            if ((tych2 .ne. 'CART') .and. (tych2 .ne. 'ELEM') .and. (tych2 .ne. 'ELGA')) then
                call utmess('F', 'UTILITAI_28', nk=2, valk=valk)
            end if
!
        else if (tychr .eq. 'ELNO') then
            if ((tych2 .ne. 'CART') .and. (tych2 .ne. 'ELNO')) then
                call utmess('F', 'UTILITAI_28', nk=2, valk=valk)
            end if
!
        else if (tychr .eq. 'ELEM') then
            if ((tych2 .ne. 'CART') .and. (tych2 .ne. 'ELEM')) then
                call utmess('F', 'UTILITAI_28', nk=2, valk=valk)
            end if
!
        else if (tychr .eq. 'CART') then
            if ((tych2 .ne. 'CART') .and. (tych2 .ne. 'ELEM')) then
                call utmess('F', 'UTILITAI_28', nk=2, valk=valk)
            end if
!
        else
            ASSERT(.false.)
        end if
!
    end do
!
!
!
!     3- ARGUMENTS POUR APPELER RELIEM :
!     ---------------------------------
    limocl(1) = 'TOUT'
    tymocl(1) = 'TOUT'
    limocl(2) = 'MAILLE'
    tymocl(2) = 'MAILLE'
    limocl(3) = 'GROUP_MA'
    tymocl(3) = 'GROUP_MA'
    limocl(4) = 'NOEUD'
    tymocl(4) = 'NOEUD'
    limocl(5) = 'GROUP_NO'
    tymocl(5) = 'GROUP_NO'
    nutrou = '&&CHPASS.NU_TROUVES'
    if (tychr .eq. 'NOEU') then
        typem = 'NU_NOEUD'
        nbmocl = 5
!
    else
        typem = 'NU_MAILLE'
        nbmocl = 3
    end if
!
!
!
!     4- BOUCLE SUR LES OCCURENCES DU MOT CLE "ASSE" :
!     -----------------------------------------------------
    nchg = 0
    iprem = .true.
    do iocc = 1, nbocc
        call getvid('ASSE', 'CHAM_GD', iocc=iocc, scal=champ, nbret=ib)
        call dismoi('TYPE_CHAMP', champ, 'CHAMP', repk=tych2)
!
        cumul = .false.
        call getvtx('ASSE', 'CUMUL', iocc=iocc, scal=kbid, nbret=ib)
        if (kbid .eq. 'OUI') cumul = .true.
!
!       4.0 CALCUL DE LA LISTE DES CMPS ET DU BOOLEEN CHGCMP
!        QUI INDIQUE QUE L'ON DOIT MODIFIER LES CMPS ET/OU LA GRANDEUR.
!       ---------------------------------------------------------------
        call getvtx('ASSE', 'NOM_CMP', iocc=iocc, nbval=0, nbret=n1)
        chgcmp = .false.
        if (n1 .lt. 0) then
            ncmp = -n1
            call wkvect('&&CHPASS.LICMP', 'V V K8', ncmp, jlicmp)
            call getvtx('ASSE', 'NOM_CMP', iocc=iocc, nbval=ncmp, vect=zk8(jlicmp), &
                        nbret=ib)
            AS_ALLOCATE(vk8=licmp2, size=ncmp)
            call getvtx('ASSE', 'NOM_CMP_RESU', iocc=iocc, nbval=0, nbret=n1)
            if (n1 .lt. 0) then
                chgcmp = .true.
                nchg = nchg+1
                if (n1 .ne. -ncmp) then
                    call utmess('F', 'UTILITAI_31')
                end if
                call getvtx('ASSE', 'NOM_CMP_RESU', iocc=iocc, nbval=ncmp, vect=licmp2, &
                            nbret=ib)
            else
                if (ncmpdisp > 0) licmp2(1:ncmp) = zk8(jlicmp:jlicmp-1+ncmp)
            end if

            if (ncmpdisp > 0 .and. nomgd .ne. 'VARI_R') then
                do i = 1, ncmp
                    compOK = ASTER_FALSE
                    do j = 1, ncmpdisp
                        if (licmp2(i) .eq. licmpdisp(j)) then
                            compOK = ASTER_TRUE
                            exit
                        end if
                    end do
                    if (.not. compOK) then
                        call utmess('A', 'UTILITAI_36', sk=licmp2(i))
                        messConseil = ASTER_TRUE
                    end if
                end do
            elseif (ncmpdisp > 0) then
                do i = 1, ncmp
                    do j = 1, ncmpdisp
                        nomcmp = licmp2(i)
                        if (nomcmp(1:1) .ne. 'V') then
                            call utmess('A', 'UTILITAI_36', sk=licmp2(i))
                        end if
                    end do
                end do
            end if
!
        else
            ncmp = 0
            jlicmp = 0
        end if
!
!       4.1 VERIFICATION DE LA GRANDEUR ASSOCIEE AU CHAMP
!       ------------------------------------------------------
        call dismoi('NOM_GD', champ, 'CHAMP', repk=nomgd2)
        if ((.not. chgcmp) .and. (nomgd2 .ne. nomgd)) then
            valk(1) = champ
            valk(2) = nomgd
            valk(3) = nomgd2
            call utmess('F', 'UTILITAI_32', nk=3, valk=valk)
        end if
!
!
        call dismoi('TYPE_SCA', nomgd2, 'GRANDEUR', repk=tsca)
        call getvc8('ASSE', 'COEF_C', iocc=iocc, scal=coefc, nbret=iret)
        if (iret .ne. 0) then
            if (tsca .ne. 'C') then
                call utmess('F', 'UTILITAI_33')
            end if
            lcoc = .true.
!
        else
            lcoc = .false.
            call getvr8('ASSE', 'COEF_R', iocc=iocc, scal=coefr, nbret=ib)
        end if
!
!       4.2 RECUPERATION DE LA LISTE DES NOEUDS OU MAILLES :
!       ----------------------------------------------------
        call reliem(modele, ma, typem, 'ASSE', iocc, &
                    nbmocl, limocl, tymocl, nutrou, nbtrou)
        if (nbtrou .eq. 0) cycle
        call jeveuo(nutrou, 'L', jnutro)
!
!
!       4.3 TRANSFORMATION ET REDUCTION DU CHAMP :
!       ------------------------------------------
        if (tych2 .eq. 'NOEU') then
            call cnocns(champ, 'V', chs1)
            call cnsred(chs1, nbtrou, zi(jnutro), ncmp, zk8(jlicmp), &
                        'V', chs2)
!
        else if (tych2(1:2) .eq. 'EL') then
            call celces(champ, 'V', chs1)

            call jeveuo(chs1//'.CESD', 'L', jcesd)
            call jeveuo(chs1//'.CESC', 'L', jcesc)

            if (nomgd2 .ne. 'VARI_R') then
                do i = 1, ncmp
                    compOK = ASTER_FALSE
                    do j = 1, zi(jcesd+1)
                        if (zk8(jlicmp-1+i) .eq. zk8(jcesc-1+j)) then
                            compOK = ASTER_TRUE
                            exit
                        end if
                    end do
                    if (.not. compOK) then
                        call utmess('F', 'UTILITAI_43', nk=2, valk=[zk8(jlicmp-1+i), champ])
                    end if
                end do
            end if

            call cesred(chs1, nbtrou, zi(jnutro), ncmp, zk8(jlicmp), &
                        'V', chs2)
!
        else if (tych2 .eq. 'CART') then
            if (tychr .eq. 'CART') then
                call carces(champ, 'ELEM', cesmod, 'V', chs1, &
                            'A', ib)
            else
                call carces(champ, tychr, cesmod, 'V', chs1, &
                            'A', ib)
            end if

            call jeveuo(chs1//'.CESD', 'L', jcesd)
            call jeveuo(chs1//'.CESC', 'L', jcesc)

            if (nomgd2 .ne. 'VARI_R') then
                do i = 1, ncmp
                    compOK = ASTER_FALSE
                    do j = 1, zi(jcesd+1)
                        if (zk8(jlicmp-1+i) .eq. zk8(jcesc-1+j)) then
                            compOK = ASTER_TRUE
                            exit
                        end if
                    end do
                    if (.not. compOK) then
                        call utmess('F', 'UTILITAI_43', nk=2, valk=[zk8(jlicmp-1+i), champ])
                    end if
                end do
            end if

            call cesred(chs1, nbtrou, zi(jnutro), ncmp, zk8(jlicmp), &
                        'V', chs2)
!
        else
            ASSERT(.false.)
        end if
!
!
!       4.4 SI ON DOIT CHANGER LES CMPS ET/OU LA GRANDEUR :
!       ----------------------------------------------------
        if (chgcmp) then
            call chsut1(chs2, nomgd, ncmp, zk8(jlicmp), licmp2, &
                        'V', chs2)
        end if
!
!
!       4.4 FUSION DU CHAMP REDUIT AVEC LE CHAMP RESULTAT :
!       ----------------------------------------------------
        if (iprem) then
            bool(1) = .false.
            call chsfus(1, chs2, bool(1), [coefr], [coefc], &
                        lcoc, 'V', chs3)
        else
            lichs(1) = chs3
            lichs(2) = chs2
            lcumul(1) = .false.
            lcumul(2) = cumul
            lcoefr(1) = 1.d0
            lcoefr(2) = coefr
            lcoefc(1) = 1.d0
            lcoefc(2) = coefc
            call chsfus(2, lichs, lcumul, lcoefr, lcoefc, &
                        lcoc, 'V', chs3)
        end if
        iprem = .false.
!
        call jedetr('&&CHPASS.LICMP')
        AS_DEALLOCATE(vk8=licmp2)
    end do
!
!     5 TRANSFORMATION DU CHAMP_S EN CHAMP :
!     ----------------------------------------------------
    if (tychr .eq. 'NOEU') then
        call cnscno(chs3, ' ', 'NON', 'G', chou, &
                    'F', ibid)
    else if (tychr(1:2) .eq. 'EL') then
        call cescel(chs3, ligrel, option, ' ', prol0, &
                    nncp, 'G', chou, 'F', ibid)
    else if (tychr .eq. 'CART') then
        call cescar(chs3, chou, 'G')
    else
        ASSERT(.false.)
    end if

    if (messConseil) call utmess('A', 'UTILITAI_42')
!
!
!     6- MENAGE :
!     -----------------------------------------------------
    if (ncmpdisp .gt. 0) AS_DEALLOCATE(vk8=licmpdisp)
    call detrsd('CHAM_NO_S', chs1)
    call detrsd('CHAM_NO_S', chs2)
    call detrsd('CHAM_NO_S', chs3)
    call detrsd('CHAM_ELEM_S', chs1)
    call detrsd('CHAM_ELEM_S', chs2)
    call detrsd('CHAM_ELEM_S', chs3)
!
    call detrsd('CHAM_ELEM_S', cesmod)
    call jedetr(nutrou)
!
    call jedema()
end subroutine
