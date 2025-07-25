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

subroutine speph0(nomu, table)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/gettco.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvem.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/posddl.h"
#include "asterfort/reliem.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsorac.h"
#include "asterfort/speph1.h"
#include "asterfort/speph2.h"
#include "asterfort/titre.h"
#include "asterfort/utchdl.h"
#include "asterfort/utmess.h"
#include "asterfort/utnono.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
!
    character(len=8) :: nomu, table
!   RESTITUTION D'UN INTERSPECTRE DE REPONSE MODALE DANS LA BASE
!   PHYSIQUE  OPERATEUR REST_SPEC_PHYS
!-----------------------------------------------------------------------
!
    integer(kind=8) :: ibid, nbmod1, nbtrou, lnumor, nbmode, ilmode, im, imod1, iad
    integer(kind=8) :: napexc, ilnoex, ncmpex, iret, ilcpex, idim1, idim0, nbn
    integer(kind=8) :: nbmail, i, inoeud, iddl, nupo, ivari, napex1, jvnomai
    integer(kind=8) :: nbmr, idim, imr, numod, in, nbpf, nbfo1, if1, ifor, ifoi, icham1
    integer(kind=8) :: isip, icham, tmod(1), nbmcl, ncmp
    integer(kind=8) :: i1, lnumi, lnumj, mxval, lfreq, lrefes, lfreqs, nbno, nbgrno, nbgrma
    real(kind=8) :: r8b, bande(2), freq1, epsi
    complex(kind=8) :: c16b
    aster_logical :: intmod, intphy, intdon
    character(len=8) :: k8b, modmec, modsta, noeud, noma, cmp, mail, noeud0
    character(len=8) :: maillage, limocl(4), nomcmp
    character(len=16) :: movrep, optcal, optcha, typcha, acces, typmec, nocham, acces0
    character(len=16) :: optch1, maille, typem
    character(len=16) :: option
    character(len=24) :: cham19, typba, group
    character(len=24) :: valk(3)
    character(len=24) :: chnumi, chnumj, chfreq
!
    character(len=3) :: toutor
    character(len=8), pointer :: maille_rep(:) => null()
    character(len=8), pointer :: nocmp_rep(:) => null()
    character(len=8), pointer :: noeud_rep(:) => null()
    character(len=24), pointer :: group_noeud(:) => null()
    integer(kind=8), pointer :: nume_ddl(:) => null()
    character(len=16), pointer :: refe(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
!
    call getvid(' ', 'MODE_MECA', scal=modmec, nbret=ibid)
    call gettco(modmec, typmec)
!
    epsi = 0.d0
    call rsorac(modmec, 'LONUTI', ibid, r8b, k8b, &
                c16b, epsi, k8b, tmod, 1, &
                nbtrou)
    nbmod1 = tmod(1)
    call wkvect('&&SPEPH0.NUMERO.ORDRE', 'V V I', nbmod1, lnumor)
    call rsorac(modmec, 'TOUT_ORDRE', ibid, r8b, k8b, &
                c16b, epsi, k8b, zi(lnumor), nbmod1, &
                nbtrou)
    call getvis(' ', 'NUME_ORDRE', nbval=0, nbret=nbmode)
    nbmode = -nbmode
    if (nbmode .eq. 0) then
        call getvtx(' ', 'TOUT_ORDRE', scal=toutor, nbret=ibid)
        if (toutor .eq. 'OUI') then
            nbmode = nbmod1
            call wkvect('&&SPEPH0.LISTEMODES', 'V V I', nbmod1, ilmode)
            do im = 1, nbmod1
                zi(ilmode-1+im) = im
                zi(ilmode-1+im) = zi(lnumor+im-1)
            end do
        else
            call getvr8(' ', 'BANDE', nbval=2, vect=bande, nbret=ibid)
            if (ibid .eq. 0) then
                call utmess('F', 'ALGORITH10_61')
            end if
            call wkvect('&&SPEPH0.LISTEMODES', 'V V I', nbmod1, ilmode)
            do im = 1, nbmod1
                imod1 = zi(lnumor+im-1)
                call rsadpa(modmec, 'L', 1, 'FREQ', imod1, &
                            0, sjv=iad, styp=k8b)
                freq1 = zr(iad)
                if ((freq1-bande(1))*(freq1-bande(2)) .le. 0.d0) then
                    nbmode = nbmode+1
                    zi(ilmode-1+nbmode) = imod1
                end if
            end do
            if (nbmode .eq. 0) then
                call utmess('F', 'ALGORITH10_31')
            end if
        end if
    else
        call getvis(' ', 'NUME_ORDRE', nbval=0, nbret=ibid)
        if (ibid .eq. 0) then
            call utmess('F', 'ALGORITH10_62')
        end if
        call wkvect('&&SPEPH0.LISTEMODES', 'V V I', nbmode, ilmode)
        call getvis(' ', 'NUME_ORDRE', nbval=nbmode, vect=zi(ilmode), nbret=ibid)
        do im = 1, nbmode
            if (zi(ilmode-1+im) .gt. nbmod1) then
                call utmess('F', 'ALGORITH10_32')
            end if
        end do
    end if
!

    call dismoi('NOM_MAILLA', modmec, 'RESULTAT', repk=maillage)

    napexc = 0
    movrep = 'RELATIF'
    call getvid(' ', 'MODE_STAT', scal=modsta, nbret=ibid)

    if (ibid .ne. 0) then
        call getvtx('EXCIT', 'NOM_CMP', iocc=1, nbval=0, nbret=ncmpex)
        ncmpex = -ncmpex
        if (ncmpex .ne. 0) then
            call wkvect('&&SPEPH0.LISTECMPEXC', 'V V K8', ncmpex, ilcpex)
            call getvtx('EXCIT', 'NOM_CMP', iocc=1, nbval=ncmpex, vect=zk8(ilcpex))
        else
            ASSERT(.false.)
        end if

        if (ncmpex .eq. 1) then
            typem = 'NO_NOEUD'
            nbmcl = 2
            limocl(1) = 'GROUP_NO'
            limocl(2) = 'NOEUD'
            call reliem(' ', maillage, typem, 'EXCIT', 1, &
                        nbmcl, limocl, limocl, '&&SPEPH0.LISTENOEEXC', napexc)
            call jeveuo('&&SPEPH0.LISTENOEEXC', 'L', ilnoex)
            if (napexc .gt. 1) then
                nomcmp = zk8(ilcpex)
                ncmpex = napexc
                call jedetr('&&SPEPH0.LISTECMPEXC')
                call wkvect('&&SPEPH0.LISTECMPEXC', 'V V K8', ncmpex, ilcpex)
                do i = 1, ncmpex
                    zk8(ilcpex) = nomcmp
                end do
            end if
        else
            call getvtx('EXCIT', 'NOEUD', iocc=1, nbval=0, nbret=napexc)
            if (napexc .eq. 0) then
                call getvtx('EXCIT', 'GROUP_NO', iocc=1, nbval=0, nbret=nbgrno)
                nbgrno = -nbgrno
                call wkvect('&&SPEPH0.LISTENOEEXC', 'V V K8', nbgrno, ilnoex)
                AS_ALLOCATE(nbgrno, vk24=group_noeud)
                call getvem(maillage, 'GROUP_NO', 'EXCIT', 'GROUP_NO', 1, &
                            nbgrno, group_noeud, iret)
                do i = 1, nbgrno
                    call utnono(' ', maillage, 'NOEUD', group_noeud(i), noeud, &
                                iret)
                    if (iret .eq. 10) then
                        call utmess('F', 'ELEMENTS_67', sk=group_noeud(i))
                    else if (iret .eq. 1) then
                        ASSERT(.false.)
                    else
                        zk8(ilnoex+i-1) = noeud
                    end if
                end do
                napexc = nbgrno
                AS_DEALLOCATE(vk24=group_noeud)
            else
                napexc = -napexc
                call wkvect('&&SPEPH0.LISTENOEEXC', 'V V K8', napexc, ilnoex)
                call getvtx('EXCIT', 'NOEUD', iocc=1, nbval=napexc, vect=zk8(ilnoex))
            end if

            if (napexc .ne. ncmpex) then
                ASSERT(.false.)
            end if
        end if

        call getvtx(' ', 'MOUVEMENT', scal=movrep, nbret=ibid)
    end if
!
    idim1 = nbmode+napexc
!
    chnumi = table//'.NUMI'
    chnumj = table//'.NUMJ'
    chfreq = table//'.DISC'
    call jeveuo(chnumi, 'L', lnumi)
    call jeveuo(chnumj, 'L', lnumj)
    call jeveuo(chfreq, 'L', lfreq)
    call jelira(chnumi, 'LONMAX', mxval)
    call jeveuo(table//'.REFE', 'L', vk16=refe)
    idim0 = 0
    do i1 = 1, mxval
        if (zi(lnumi-1+i1) .ge. idim0) then
            idim0 = zi(lnumi-1+i1)
        end if
    end do
!
    if (idim1 .ne. idim0) then
        call utmess('F', 'ALGORITH10_63')
    end if
!
!     --- OPTION DE RECOMBINAISON ---
!
    call getvtx(' ', 'NOM_CHAM', scal=optcha, nbret=ibid)
!
!     --- VERIFICATION DES DONNEES INTERSPECTRE ---
!
    nocham = refe(1)
!
    if (nocham .eq. 'ACCE_GENE') then
        if (optcha(1:4) .eq. 'ACCE') then
        else
            call utmess('F', 'ALGORITH10_64')
        end if
    else if (nocham .eq. 'VITE_GENE') then
        if (optcha(1:4) .eq. 'VITE') then
        else
            call utmess('F', 'ALGORITH10_65')
        end if
    else if (nocham .eq. 'DEPL_GENE') then
        if (optcha(1:4) .eq. 'ACCE') then
            call utmess('F', 'ALGORITH10_66')
        else if (optcha(1:4) .eq. 'VITE') then
            call utmess('F', 'ALGORITH10_67')
        end if
    end if
    optch1 = optcha
    if (optcha(1:4) .eq. 'VITE') optch1 = 'DEPL'
    if (optcha(1:4) .eq. 'ACCE') optch1 = 'DEPL'
!
!     --- RECUPERATION DES NOEUDS, NOM_CMP ET MAILLE ---
!
    call getvtx(' ', 'NOM_CMP', iocc=1, nbval=0, nbret=ncmp)
    ncmp = -ncmp
    if (ncmp .ne. 0) then
        AS_ALLOCATE(ncmp, vk8=nocmp_rep)
        call getvtx(' ', 'NOM_CMP', iocc=1, nbval=ncmp, vect=nocmp_rep)
    else
        ASSERT(.false.)
    end if

    jvnomai = 0

    if (ncmp .eq. 1) then
        typem = 'NO_NOEUD'
        nbmcl = 2
        limocl(1) = 'GROUP_NO'
        limocl(2) = 'NOEUD'
        call reliem(' ', maillage, typem, ' ', 1, &
                    nbmcl, limocl, limocl, '&&SPEPH0.LISTENOEREP', nbno)
        call jeveuo('&&SPEPH0.LISTENOEREP', 'L', vk8=noeud_rep)
        if (nbno .gt. 1) then
            nomcmp = nocmp_rep(1)
            ncmp = nbno
            AS_DEALLOCATE(vk8=nocmp_rep)
            AS_ALLOCATE(ncmp, vk8=nocmp_rep)
            do i = 1, ncmp
                nocmp_rep(i) = nomcmp
            end do
        end if

        typem = 'NO_MAILLE'
        nbmcl = 2
        limocl(1) = 'GROUP_MA'
        limocl(2) = 'MAILLE'
        call reliem(' ', maillage, typem, ' ', 1, &
                    nbmcl, limocl, limocl, '&&SPEPH0.LISTEMAEREP', nbmail)
        if (nbmail .ne. 0) then
            call jeveuo('&&SPEPH0.LISTEMAEREP', 'L', vk8=maille_rep)
            if (nbno .ne. nbmail) then
                call utmess('F', 'ALGORITH10_69')
            end if
        end if

        jvnomai = 1
    else
        call getvtx(' ', 'NOEUD', iocc=1, nbval=0, nbret=nbno)
        if (nbno .eq. 0) then
            call getvtx(' ', 'GROUP_NO', iocc=1, nbval=0, nbret=nbgrno)
            nbgrno = -nbgrno
            AS_ALLOCATE(nbgrno, vk8=noeud_rep)
            AS_ALLOCATE(nbgrno, vk24=group_noeud)
            call getvem(maillage, 'GROUP_NO', '', 'GROUP_NO', 1, &
                        nbgrno, group_noeud, iret)
            do i = 1, nbgrno
                call utnono(' ', maillage, 'NOEUD', group_noeud(i), noeud, &
                            iret)
                !print*, "noeud" , noeud, i
                if (iret .eq. 10) then
                    call utmess('F', 'ELEMENTS_67', sk=group_noeud(i))
                else if (iret .eq. 1) then
                    ASSERT(.false.)
                else
                    noeud_rep(i) = noeud
                end if
            end do
            nbno = nbgrno
            AS_DEALLOCATE(vk24=group_noeud)
        else
            nbno = -nbno
            AS_ALLOCATE(nbno, vk8=noeud_rep)
            call getvtx(' ', 'NOEUD', iocc=1, nbval=nbno, vect=noeud_rep)
        end if

        if (nbno .ne. ncmp) then
            call utmess('F', 'ALGORITH10_68', sk=group)
        end if

        call getvtx(' ', 'MAILLE', iocc=1, nbval=0, nbret=nbmail)
        if (nbmail .eq. 0) then
            call getvtx(' ', 'GROUP_MA', iocc=1, nbval=0, nbret=nbgrma)
            if (nbgrma .ne. 0) then
                nbgrma = -nbgrma
                AS_ALLOCATE(nbgrma, vk8=maille_rep)
                do i = 1, nbgrma
                    call getvem(maillage, 'GROUP_MA', 'EXCIT', 'GROUP_MA', i, &
                                0, group, iret)
                    call utnono(' ', maillage, 'MAILLE', group, mail, &
                                iret)
                    if (iret .eq. 10) then
                        call utmess('F', 'ELEMENTS_63', sk=group)
                    else if (iret .eq. 1) then
                        ASSERT(.false.)
                    else
                        maille_rep(i) = mail
                    end if
                end do
                nbmail = nbgrma
                if (nbno .ne. nbmail) then
                    call utmess('F', 'ALGORITH10_69')
                end if
            end if
        else
            nbmail = -nbmail
            AS_ALLOCATE(nbmail, vk8=maille_rep)
            call getvtx(' ', 'MAILLE', iocc=1, nbval=nbmail, vect=maille_rep)

            if (nbno .ne. nbmail) then
                call utmess('F', 'ALGORITH10_69')
            end if
        end if

    end if

    nbn = nbno
!
!
!     --- RECUPERATION DU NUMERO DU DDL ---
!
    call rsexch('F', modmec, optch1, zi(ilmode), cham19, &
                iret)
    AS_ALLOCATE(vi=nume_ddl, size=nbn)
    call dismoi('TYPE_SUPERVIS', cham19, 'CHAMP', repk=typcha)
!
    if (typcha(1:7) .eq. 'CHAM_NO') then
        do i = 1, nbn
            noeud = noeud_rep(i)
            cmp = nocmp_rep(i)
            call posddl('CHAM_NO', cham19, noeud, cmp, inoeud, &
                        iddl)
            if (inoeud .eq. 0) then
                call utmess('F', 'UTILITAI_92', sk=noeud)
            else if (iddl .eq. 0) then
                valk(1) = cmp
                valk(2) = noeud
                call utmess('F', 'UTILITAI_93', nk=2, valk=valk)
            end if
            nume_ddl(i) = iddl
        end do
!
    else if (typcha(1:9) .eq. 'CHAM_ELEM') then
        if (nbmail .eq. 0) then
            call utmess('F', 'ALGORITH10_72')
        end if
        call dismoi('NOM_MAILLA', cham19, 'CHAM_ELEM', repk=noma)
        nupo = 0
        ivari = 1
        do i = 1, nbn
            maille = maille_rep(i)
            noeud = noeud_rep(i)
            cmp = nocmp_rep(i)
            call utchdl(cham19, noma, maille, noeud, nupo, &
                        0, ivari, cmp, iddl)
            if (iddl .eq. 0) then
                valk(1) = cmp
                valk(2) = noeud
                valk(3) = maille
                call utmess('F', 'ALGORITH10_73', nk=3, valk=valk)
            end if
            nume_ddl(i) = iddl
        end do
    else
        call utmess('F', 'CALCULEL_17')
    end if
!
    call getvtx(' ', 'OPTION', scal=optcal, nbret=ibid)
    intphy = .false.
    intmod = .false.
    if (optcal(1:4) .eq. 'TOUT') intphy = .true.
    if (optcal(6:9) .eq. 'TOUT') intmod = .true.
!
! --- CARACTERISATION DU CONTENU DE LA TABLE   ---
! --- INTERSPECTRES OU AUTOSPECTRES UNIQUEMENT ---
!
    option = refe(2)
!
    intdon = .true.
    if (option(1:4) .eq. 'DIAG') intdon = .false.
    if (intmod .and. .not. intdon) then
        call utmess('F', 'MODELISA5_81')
    end if
!
!     --- ON NE PREND EN COMPTE QUE LES MODES DYNAMIQUES ---
!
    nbmod1 = nbmode
    napex1 = napexc
    if (movrep .eq. 'DIFFERENTIEL') nbmod1 = 0
    if (movrep .eq. 'RELATIF') napex1 = 0
    nbmr = napex1+nbmod1
    idim = nbmr*nbn
    call wkvect('&&SPEPH0_CHAM', 'V V R', idim, icham)
!
    do imr = 1, napex1
        noeud = zk8(ilnoex+imr-1)
        cmp = zk8(ilcpex+imr-1)
        acces = noeud//cmp
        call rsorac(modsta, 'NOEUD_CMP', ibid, r8b, acces, &
                    c16b, r8b, k8b, tmod, 1, &
                    nbtrou)
        numod = tmod(1)
        if (nbtrou .ne. 1 .and. noeud(1:1) .eq. 'N') then
            noeud0(1:7) = noeud(2:8)
            noeud0(8:8) = ' '
            acces0 = noeud0//cmp
            call rsorac(modsta, 'NOEUD_CMP', ibid, r8b, acces0, &
                        c16b, r8b, k8b, tmod, 1, &
                        nbtrou)
            numod = tmod(1)
        end if
        if (nbtrou .ne. 1) then
            valk(1) = modsta
            valk(2) = acces
            call utmess('F', 'ALGORITH14_63', nk=2, valk=valk)
        end if
        call rsexch('F', modsta, optch1, numod, cham19, &
                    iret)

        call dismoi('TYPE_SUPERVIS', cham19, 'CHAMP', repk=typcha)

        if (typcha(1:7) .eq. 'CHAM_NO') then
            call jeveuo(cham19(1:19)//'.VALE', 'L', isip)
            do in = 1, nbn
                noeud = noeud_rep(in)
                cmp = nocmp_rep(in)
                call posddl('CHAM_NO', cham19, noeud, cmp, inoeud, &
                            iddl)
                icham1 = icham+nbn*(imr-1)+in-1
                zr(icham1) = zr(isip+iddl-1)
            end do
        else if (typcha(1:9) .eq. 'CHAM_ELEM') then
            call jeveuo(cham19(1:19)//'.CELV', 'L', isip)
            call dismoi('NOM_MAILLA', cham19, 'CHAM_ELEM', repk=noma)
            nupo = 0
            ivari = 1
            do i = 1, nbn
                maille = maille_rep(i)
                noeud = noeud_rep(i)
                cmp = nocmp_rep(i)
                call utchdl(cham19, noma, maille, noeud, nupo, &
                            0, ivari, cmp, iddl)
                icham1 = icham+nbn*(imr-1)+in-1
                zr(icham1) = zr(isip+iddl-1)
            end do
        end if

    end do
!
    call dismoi('TYPE_BASE', modmec, 'RESU_DYNA', repk=typba, arret='C', &
                ier=iret)
!
    do imr = 1, nbmod1
        numod = zi(ilmode+imr-1)
        call rsexch('F', modmec, optch1, numod, cham19, &
                    iret)

        if (typcha(1:7) .eq. 'CHAM_NO') then
            call jeveuo(cham19(1:19)//'.VALE', 'L', isip)
        else if (typcha(1:9) .eq. 'CHAM_ELEM') then
            call jeveuo(cham19(1:19)//'.CELV', 'L', isip)
        end if

        if (typmec .eq. 'MODE_MECA_C') then
            do in = 1, nbn
                icham1 = icham+napex1*nbn+nbn*(imr-1)+in-1
                zr(icham1) = dble(zc(isip+nume_ddl(in)-1))
            end do
!  -------------------------------
!  si base modale, alors les nume_ddl des differents modes peuvent
!              etres differents
!
        else if (typba(1:1) .ne. ' ') then
            call dismoi('TYPE_SUPERVIS', cham19, 'CHAMP', repk=typcha)
            if (typcha(1:7) .eq. 'CHAM_NO') then
                do in = 1, nbn
                    noeud = noeud_rep(in)
                    cmp = nocmp_rep(in)
                    call posddl('CHAM_NO', cham19, noeud, cmp, inoeud, &
                                iddl)
                    icham1 = icham+nbn*(imr-1)+in-1
                    zr(icham1) = zr(isip+iddl-1)
                end do
            else if (typcha(1:9) .eq. 'CHAM_ELEM') then
                call dismoi('NOM_MAILLA', cham19, 'CHAM_ELEM', repk=noma)
                nupo = 0
                ivari = 1
                do i = 1, nbn
                    maille = maille_rep(i)
                    noeud = noeud_rep(i)
                    cmp = nocmp_rep(i)
                    call utchdl(cham19, noma, maille, noeud, nupo, &
                                0, ivari, cmp, iddl)
                    icham1 = icham+nbn*(imr-1)+in-1
                    zr(icham1) = zr(isip+iddl-1)
                end do
            end if
!  -------------------------------
        else
            do in = 1, nbn
                icham1 = icham+napex1*nbn+nbn*(imr-1)+in-1
                zr(icham1) = zr(isip+nume_ddl(in)-1)
            end do
        end if
    end do
!
!     --- CREATION DE LA TABLE DE SORTIE ---
!
    nbfo1 = (nbmr*(nbmr+1))/2
!
    call wkvect(nomu//'.REFE', 'G V K16', 3, lrefes)
    zk16(lrefes) = optcha
    zk16(lrefes+1) = optcal
    zk16(lrefes+2) = 'FREQ'
!
    call jelira(chfreq, 'LONMAX', nbpf)
    call wkvect(nomu//'.DISC', 'G V R', nbpf, lfreqs)
    call wkvect('&&SPEPH0.TEMP.FONR', 'V V R', nbpf*nbfo1, ifor)
    call wkvect('&&SPEPH0.TEMP.FONI', 'V V R', nbpf*nbfo1, ifoi)
!
    do if1 = 1, nbpf
        zr(lfreqs+if1-1) = zr(lfreq+if1-1)
    end do
!
    call speph2(movrep, napexc, nbmode, nbpf, intmod, &
                table, zr(ifor), zr(ifoi))
!
    call speph1(intphy, intmod, nomu, zr(icham), zr(ifor), &
                zr(ifoi), noeud_rep, nocmp_rep, nbmr, nbn, &
                nbpf)
!
    call titre()
!
    call jedetr('&&SPEPH0.NUMERO.ORDRE')
    call jedetr('&&SPEPH0.LISTEMODES')
    AS_DEALLOCATE(vi=nume_ddl)
    call jedetr('&&SPEPH0_CHAM')
    call jedetr('&&SPEPH0.TEMP.FONR')
    call jedetr('&&SPEPH0.TEMP.FONI')
    if (jvnomai .eq. 0) then
        AS_DEALLOCATE(vk8=noeud_rep)
        AS_DEALLOCATE(vk8=maille_rep)
    end if
    AS_DEALLOCATE(vk8=nocmp_rep)
    call jeexin('&&SPEPH0.LISTENOEEXC', iret)
    if (iret .ne. 0) call jedetr('&&SPEPH0.LISTENOEEXC')
    call jeexin('&&SPEPH0.LISTECMPEXC', iret)
    if (iret .ne. 0) call jedetr('&&SPEPH0.LISTECMPEXC')
!
    call jedema()
end subroutine
