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
subroutine rvpara(nomtab, mcf, nbpost)
! IN  NOMTAB  : NOM DE LA TABLE PRINCIPALE PRODUITE PAR LA COMMANDE
! IN  MCF     : MOT-CLE FACTEUR
! IN  NBPOST  : NOMBRE DE POST-TRAITEMENT A CONSIDERER
! ----------------------------------------------------------------------
!     INITIALISE LA TABLE DE POST_RELEVE_T ASSOCIEE A LA TABLE DE
!     REFERENCE NOMTAB
!     ------------------------------------------------------------------
!
    implicit none
!
! 0.1. ==> ARGUMENTS
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/gettco.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/juveca.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsorac.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
#include "asterfort/titrea.h"
#include "asterfort/utmess.h"
#include "asterfort/utncmp.h"
#include "asterfort/wkvect.h"
    character(len=6) :: mcf
    character(len=8) :: nomtab
    integer(kind=8) :: nbpost
!
!
! 0.3. ==> VARIABLES LOCALES
!
    character(len=6) :: nompro
    parameter(nompro='RVPARA')
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: iocc, ibid, tord(1)
    integer(kind=8) :: jcham, jresu, jncmp, ncmp, i
    integer(kind=8) :: jinva, jprin, jmome, jmail, jmoye, j, jtrad
    integer(kind=8) :: jtran, n1, n2, n3, jcmp1, jcmp2, jcmp3, nbc, nume
    integer(kind=8) :: iret, nbp, jinst, jordr, jmode, jabsc, jfreq
    integer(kind=8) :: jnoeu, n11, n12, n13, n14, n15, n16, n17, n18
    integer(kind=8) :: jncas, jangl, jnocp, numcmp, jnucp, nbordr, jnume
    real(kind=8) :: r8b
    aster_logical :: lmima, lmoye, lextr, lmoygr
    complex(kind=8) :: c16b
    character(len=8) :: k8b, resu, typara(100), nomcmp
    character(len=16) :: k16b, nomsy, tysd
    character(len=24) :: nomobj, chextr, nopara(100), knume
    character(len=24) :: valk(3)
    character(len=24) :: k24bid
!
    character(len=24) :: nocmp
    integer(kind=8) :: jnocmp, ncmpmx
!     ------------------------------------------------------------------
!
!====
! 1. PREALABLES
!====
!
    call jemarq()
!
    call infmaj()
    call infniv(ifm, niv)
!
    if (niv .ge. 2) then
        call utmess('I', 'POSTRELE_8', sk=nomtab)
    end if
!
    nocmp = '&&'//nompro//'_NOM_CMP_TABLE  '
    ncmpmx = 100
    call wkvect(nocmp, 'V V K8', ncmpmx, jnocmp)
!
    jabsc = 0
    jcham = 0
    jresu = 0
    jordr = 0
    jmode = 0
    jinst = 0
    jfreq = 0
    jncas = 0
    jangl = 0
    jnocp = 0
    jncmp = 0
    jinva = 0
    jprin = 0
    jmome = 0
    jnoeu = 0
    jmail = 0
    jmoye = 0
    jtran = 0
    jtrad = 0
    ncmp = 0
!
!====
! 2. ON PARCOURT TOUTES LES ACTIONS DEMANDEES
!====
!
    do iocc = 1, nbpost
!
! 2.1. ==> ON CHERCHE SI C'EST LA BONNE TABLE
!
!
        call getvid(mcf, 'CHAM_GD', iocc=iocc, nbval=0, nbret=n2)
        if (n2 .ne. 0) jcham = jcham+1
!
        call getvid(mcf, 'RESULTAT', iocc=iocc, nbval=0, nbret=n3)
        if (n3 .ne. 0) then
            jresu = jresu+1
            call getvid(mcf, 'RESULTAT', iocc=iocc, scal=k8b, nbret=n3)
            call gettco(k8b, tysd)
            if (tysd .eq. 'EVOL_ELAS' .or. tysd .eq. 'EVOL_THER' .or. tysd .eq. 'EVOL_NOLI' &
                .or. tysd .eq. 'EVOL_CHAR' .or. tysd .eq. 'DYNA_TRANS') then
                jinst = jinst+1
            elseif (tysd .eq. 'DYNA_HARMO' .or. tysd .eq. &
                    'HARM_GENE' .or. tysd .eq. 'ACOU_HARMO') then
                jfreq = jfreq+1
            elseif (tysd .eq. 'MODE_MECA' .or. tysd .eq. 'MODE_GENE' &
                    .or. tysd .eq. 'MODE_ACOU') then
                jfreq = jfreq+1
                jmode = jmode+1
                jnocp = jnocp+1
            else if (tysd .eq. 'MULT_ELAS') then
                jncas = jncas+1
            else if (tysd(1:8) .eq. 'FOURIER_') then
                jmode = jmode+1
            else if (tysd .eq. 'COMB_FOURIER') then
                jangl = jangl+1
            end if
        end if
!
        call getvid(mcf, 'LIST_ORDRE', iocc=iocc, nbval=0, nbret=n11)
        if (n11 .ne. 0) jordr = jordr+1
!
        call getvis(mcf, 'NUME_ORDRE', iocc=iocc, nbval=0, nbret=n12)
        if (n12 .ne. 0) jordr = jordr+1
!
        call getvid(mcf, 'LIST_MODE', iocc=iocc, nbval=0, nbret=n13)
        if (n13 .ne. 0) jmode = jmode+1
!
        call getvis(mcf, 'NUME_MODE', iocc=iocc, nbval=0, nbret=n14)
        if (n14 .ne. 0) jmode = jmode+1
!
        call getvid(mcf, 'LIST_INST', iocc=iocc, nbval=0, nbret=n15)
        if (n15 .ne. 0) jinst = jinst+1
!
        call getvr8(mcf, 'INST', iocc=iocc, nbval=0, nbret=n16)
        if (n16 .ne. 0) jinst = jinst+1
!
        call getvid(mcf, 'LIST_FREQ', iocc=iocc, nbval=0, nbret=n17)
        if (n17 .ne. 0) jfreq = jfreq+1
!
        call getvr8(mcf, 'FREQ', iocc=iocc, nbval=0, nbret=n18)
        if (n18 .ne. 0) jfreq = jfreq+1
!
        if ((n2+n11+n12+n13+n14+n15+n16+n17+n18) .eq. 0) jordr = jordr+1
!
        call getvtx(mcf, 'TOUT_CMP', iocc=iocc, nbval=0, nbret=n1)
        if (n1 .ne. 0) then
            jncmp = jncmp+1
            nomobj = '&&'//nompro//'.NCMP'
            if (n2 .ne. 0) then
                call getvid(mcf, 'CHAM_GD', iocc=iocc, scal=nomsy, nbret=n2)
                call utncmp(nomsy, nbc, nomobj)
            else
                call getvid(mcf, 'RESULTAT', iocc=iocc, scal=resu, nbret=n3)
                call getvtx(mcf, 'NOM_CHAM', iocc=iocc, scal=nomsy, nbret=n1)
!
                call rsorac(resu, 'LONUTI', 0, r8b, k8b, &
                            c16b, r8b, k8b, tord, 1, &
                            ibid)
                nbordr = tord(1)
                knume = '&&'//nompro//'.NUME_ORDRE'
                call wkvect(knume, 'V V I', nbordr, jnume)
                call rsorac(resu, 'TOUT_ORDRE', 0, r8b, k8b, &
                            c16b, r8b, k8b, zi(jnume), nbordr, &
                            ibid)
                do i = 1, nbordr
                    nume = zi(jnume+i-1)
                    call rsexch(' ', resu, nomsy, nume, chextr, &
                                iret)
                    if (iret .eq. 0) goto 16
                end do
                call utmess('F', 'POSTRELE_9', sk=nomsy)
16              continue
                call jedetr(knume)
                call utncmp(chextr, nbc, nomobj)
            end if
            if (nbc .eq. 0) then
                call utmess('F', 'POSTRELE_59')
            end if
            call jeveuo(nomobj, 'L', jcmp1)
            do i = 1, nbc
                do j = 1, ncmp
                    if (zk8(jnocmp-1+j) .eq. zk8(jcmp1+i-1)) goto 10
                end do
                ncmp = ncmp+1
                if (ncmp .gt. ncmpmx) then
                    ncmpmx = 2*ncmpmx
                    call juveca(nocmp, ncmpmx)
                    call jeveuo(nocmp, 'E', jnocmp)
                end if
                zk8(jnocmp-1+ncmp) = zk8(jcmp1+i-1)
10              continue
            end do
            call jedetr(nomobj)
        end if
!
        call getvtx(mcf, 'NOM_CMP', iocc=iocc, nbval=0, nbret=n1)
        if (n1 .ne. 0) then
!
            call getvtx(mcf, 'TRAC_NOR', iocc=iocc, nbval=0, nbret=n12)
            if (n12 .ne. 0) jtran = jtran+1
!
            call getvtx(mcf, 'TRAC_DIR', iocc=iocc, nbval=0, nbret=n14)
            if (n14 .ne. 0) jtrad = jtrad+1
!
            if ((n12+n14) .ne. 0) goto 24
            jncmp = jncmp+1
            nbc = -n1
            call wkvect('&&'//nompro//'.NCMP', 'V V K8', nbc, jcmp2)
            call getvtx(mcf, 'NOM_CMP', iocc=iocc, nbval=nbc, vect=zk8(jcmp2), &
                        nbret=n1)
!           CALL GETVIS ( MCF, 'NUME_CMP', IOCC,1,0, IBID,N11)
            n11 = 0
            if (n11 .ne. 0) then
                numcmp = -n11
                call wkvect('&&'//nompro//'.NU_CMP', 'V V I', numcmp, jnucp)
!           CALL GETVIS(MCF,'NUME_CMP',IOCC,IARG,NUMCMP,ZI(JNUCP),N11)
                n11 = 0
                if (zk8(jcmp2) (1:4) .eq. 'VARI') then
                    ASSERT(nbc .eq. 1)
                    do i = 1, numcmp
                        call codent(zi(jnucp+i-1), 'G', k8b)
                        nomcmp = 'VARI_'//k8b(1:3)
                        do j = 1, ncmp
                            if (zk8(jnocmp-1+j) .eq. nomcmp) goto 120
                        end do
                        ncmp = ncmp+1
                        if (ncmp .gt. ncmpmx) then
                            ncmpmx = 2*ncmpmx
                            call juveca(nocmp, ncmpmx)
                            call jeveuo(nocmp, 'E', jnocmp)
                        end if
                        zk8(jnocmp-1+ncmp) = nomcmp
120                     continue
                    end do
                else
                    do i = 1, nbc
                        do j = 1, ncmp
                            if (zk8(jnocmp-1+j) .eq. zk8(jcmp2+i-1)) goto 124
                        end do
                        ncmp = ncmp+1
                        if (ncmp .gt. ncmpmx) then
                            ncmpmx = 2*ncmpmx
                            call juveca(nocmp, ncmpmx)
                            call jeveuo(nocmp, 'E', jnocmp)
                        end if
                        zk8(jnocmp-1+ncmp) = zk8(jcmp2+i-1)
124                     continue
                    end do
                end if
                call jedetr('&&'//nompro//'.NU_CMP')
            else
                do i = 1, nbc
                    do j = 1, ncmp
                        if (zk8(jnocmp-1+j) .eq. zk8(jcmp2+i-1)) goto 20
                    end do
                    ncmp = ncmp+1
                    if (ncmp .gt. ncmpmx) then
                        ncmpmx = 2*ncmpmx
                        call juveca(nocmp, ncmpmx)
                        call jeveuo(nocmp, 'E', jnocmp)
                    end if
                    zk8(jnocmp-1+ncmp) = zk8(jcmp2+i-1)
20                  continue
                end do
            end if
            call jedetr('&&'//nompro//'.NCMP')
        end if
24      continue
!
        call getvtx(mcf, 'ELEM_PRINCIPAUX', iocc=iocc, nbval=0, nbret=n1)
        if (n1 .ne. 0) jprin = jprin+1
!
        call getvtx(mcf, 'RESULTANTE', iocc=iocc, nbval=0, nbret=n1)
        call getvtx(mcf, 'MOMENT', iocc=iocc, nbval=0, nbret=n2)
        if ((n1 .ne. 0) .and. (n2 .ne. 0)) jmome = jmome+1
        if ((n1 .ne. 0) .and. (n2 .eq. 0)) then
            jncmp = jncmp+1
            nbc = -n1
            call wkvect('&&'//nompro//'.NCMP', 'V V K8', nbc, jcmp3)
            call getvtx(mcf, 'RESULTANTE', iocc=iocc, nbval=nbc, vect=zk8(jcmp3), &
                        nbret=n1)
            do i = 1, nbc
                do j = 1, ncmp
                    if (zk8(jnocmp-1+j) .eq. zk8(jcmp3+i-1)) goto 30
                end do
                ncmp = ncmp+1
                if (ncmp .gt. ncmpmx) then
                    ncmpmx = 2*ncmpmx
                    call juveca(nocmp, ncmpmx)
                    call jeveuo(nocmp, 'E', jnocmp)
                end if
                zk8(jnocmp-1+ncmp) = zk8(jcmp3+i-1)
30              continue
            end do
            call jedetr('&&'//nompro//'.NCMP')
        end if
!
        lmima = .false.
        lmoygr = .false.
        lmoye = .false.
        lextr = .false.
        call getvtx(mcf, 'OPERATION', iocc=iocc, scal=k16b, nbret=n3)
        if (k16b .eq. 'EXTREMA') lmima = .true.
        if (k16b .eq. 'MOYENNE_ARITH') lmoygr = .true.
        if (k16b .eq. 'MOYENNE') then
            jmoye = jmoye+1
            lmoye = .true.
        end if
        if (k16b .eq. 'EXTRACTION') then
            lextr = .true.
!
            call getvtx(mcf, 'INVARIANT', iocc=iocc, nbval=0, nbret=n2)
            if (n2 .ne. 0) jinva = jinva+1
!
            if (n1 .eq. 0) jabsc = jabsc+1
!
            call getvtx(mcf, 'NOEUD', iocc=iocc, nbval=0, nbret=n2)
            if ((n1 .eq. 0) .and. (n2 .ne. 0)) jnoeu = jnoeu+1
!
            call getvtx(mcf, 'GROUP_NO', iocc=iocc, nbval=0, nbret=n2)
            if ((n1 .eq. 0) .and. (n2 .ne. 0)) jnoeu = jnoeu+1
        end if
!
        call getvtx(mcf, 'MOYE_NOEUD', iocc=iocc, nbval=0, nbret=n1)
        if (n1 .ne. 0) then
            call getvtx(mcf, 'MOYE_NOEUD', iocc=iocc, scal=k8b, nbret=n1)
            if (k8b(1:3) .eq. 'NON') jmail = jmail+1
        end if
!
    end do
!
!====
! 3. CONNAISSANT LES CARACTERISTIQUES DE LA TABLE, ON INITIALISE
!====
!
! 3.1. ==> MISE EN PLACE DES PARAMETRES
!
    nbp = 1
    nopara(nbp) = 'INTITULE'
    typara(nbp) = 'K16'
    if ((lextr .or. lmoye) .and. jnoeu .ne. 0) then
        nbp = nbp+1
        nopara(nbp) = 'NOEUD'
        typara(nbp) = 'K8'
    end if
    if (jcham .ne. 0) then
        nbp = nbp+1
        nopara(nbp) = 'CHAM_GD'
        typara(nbp) = 'K8'
    end if
    if (jresu .ne. 0) then
        nbp = nbp+1
        nopara(nbp) = 'RESU'
        typara(nbp) = 'K8'
        nbp = nbp+1
        nopara(nbp) = 'NOM_CHAM'
        typara(nbp) = 'K16'
    end if
!
    if ((jresu .ne. 0 .and. jordr .eq. 0) .or. (jordr .ne. 0)) then
        nbp = nbp+1
        nopara(nbp) = 'NUME_ORDRE'
        typara(nbp) = 'I'
    end if
    if (jmode .ne. 0) then
        nbp = nbp+1
        nopara(nbp) = 'NUME_MODE'
        typara(nbp) = 'I'
    end if
    if (jinst .ne. 0) then
        nbp = nbp+1
        nopara(nbp) = 'INST'
        typara(nbp) = 'R'
    end if
    if (jfreq .ne. 0) then
        nbp = nbp+1
        nopara(nbp) = 'FREQ'
        typara(nbp) = 'R'
    end if
    if ((lextr .or. lmoye) .and. jncas .ne. 0) then
        nbp = nbp+1
        nopara(nbp) = 'NOM_CAS'
        typara(nbp) = 'K16'
    end if
    if ((lextr .or. lmoye) .and. jnocp .ne. 0) then
        nbp = nbp+1
        nopara(nbp) = 'NOEUD_CMP'
        typara(nbp) = 'K16'
    end if
    if ((lextr .or. lmoye) .and. jangl .ne. 0) then
        nbp = nbp+1
        nopara(nbp) = 'ANGL'
        typara(nbp) = 'R'
    end if
    if ((lextr .or. lmoye) .and. jmail .ne. 0) then
        nbp = nbp+1
        nopara(nbp) = 'MAILLE'
        typara(nbp) = 'K8'
    end if
    if ((lextr .or. lmoye) .and. jabsc .ne. 0) then
        nbp = nbp+1
        nopara(nbp) = 'ABSC_CURV'
        typara(nbp) = 'R'
        nbp = nbp+1
        nopara(nbp) = 'COOR_X'
        typara(nbp) = 'R'
        nbp = nbp+1
        nopara(nbp) = 'COOR_Y'
        typara(nbp) = 'R'
        nbp = nbp+1
        nopara(nbp) = 'COOR_Z'
        typara(nbp) = 'R'
    end if
    if ((lextr .or. lmoye) .and. jncmp .ne. 0) then
        do i = 1, ncmp
            nbp = nbp+1
            nopara(nbp) = zk8(jnocmp-1+i)
            typara(nbp) = 'R'
        end do
    end if
    if ((lextr .or. lmoye) .and. jinva .ne. 0) then
        nbp = nbp+1
        nopara(nbp) = 'VMIS'
        typara(nbp) = 'R'
        nbp = nbp+1
        nopara(nbp) = 'TRESCA'
        typara(nbp) = 'R'
        nbp = nbp+1
        nopara(nbp) = 'TRACE'
        typara(nbp) = 'R'
        nbp = nbp+1
        nopara(nbp) = 'DETER'
        typara(nbp) = 'R'
    end if
    if ((lextr .or. lmoye) .and. jprin .ne. 0) then
        nbp = nbp+1
        nopara(nbp) = 'PRIN_1'
        typara(nbp) = 'R'
        nbp = nbp+1
        nopara(nbp) = 'PRIN_2'
        typara(nbp) = 'R'
        nbp = nbp+1
        nopara(nbp) = 'PRIN_3'
        typara(nbp) = 'R'
        nbp = nbp+1
        nopara(nbp) = 'VECT_1_X'
        typara(nbp) = 'R'
        nbp = nbp+1
        nopara(nbp) = 'VECT_1_Y'
        typara(nbp) = 'R'
        nbp = nbp+1
        nopara(nbp) = 'VECT_1_Z'
        typara(nbp) = 'R'
        nbp = nbp+1
        nopara(nbp) = 'VECT_2_X'
        typara(nbp) = 'R'
        nbp = nbp+1
        nopara(nbp) = 'VECT_2_Y'
        typara(nbp) = 'R'
        nbp = nbp+1
        nopara(nbp) = 'VECT_2_Z'
        typara(nbp) = 'R'
        nbp = nbp+1
        nopara(nbp) = 'VECT_3_X'
        typara(nbp) = 'R'
        nbp = nbp+1
        nopara(nbp) = 'VECT_3_Y'
        typara(nbp) = 'R'
        nbp = nbp+1
        nopara(nbp) = 'VECT_3_Z'
        typara(nbp) = 'R'
    end if
    if ((lextr .or. lmoye) .and. jmome .ne. 0) then
        nbp = nbp+1
        nopara(nbp) = 'RESULT_X'
        typara(nbp) = 'R'
        nbp = nbp+1
        nopara(nbp) = 'RESULT_Y'
        typara(nbp) = 'R'
        nbp = nbp+1
        nopara(nbp) = 'RESULT_Z'
        typara(nbp) = 'R'
        nbp = nbp+1
        nopara(nbp) = 'MOMENT_X'
        typara(nbp) = 'R'
        nbp = nbp+1
        nopara(nbp) = 'MOMENT_Y'
        typara(nbp) = 'R'
        nbp = nbp+1
        nopara(nbp) = 'MOMENT_Z'
        typara(nbp) = 'R'
    end if
    if ((lextr .or. lmoye) .and. jtran .ne. 0) then
        nbp = nbp+1
        nopara(nbp) = 'TRAC_NOR'
        typara(nbp) = 'R'
        nbp = nbp+1
        nopara(nbp) = 'TR_NOR_1'
        typara(nbp) = 'R'
        nbp = nbp+1
        nopara(nbp) = 'TR_NOR_2'
        typara(nbp) = 'R'
        nbp = nbp+1
        nopara(nbp) = 'TR_NOR_3'
        typara(nbp) = 'R'
    end if
    if ((lextr .or. lmoye) .and. jtrad .ne. 0) then
        nbp = nbp+1
        nopara(nbp) = 'TRAC_DIR'
        typara(nbp) = 'R'
        nbp = nbp+1
        nopara(nbp) = 'TR_DIR_1'
        typara(nbp) = 'R'
        nbp = nbp+1
        nopara(nbp) = 'TR_DIR_2'
        typara(nbp) = 'R'
        nbp = nbp+1
        nopara(nbp) = 'TR_DIR_3'
        typara(nbp) = 'R'
    end if
    if ((lextr .or. lmoye) .and. jmoye .ne. 0) then
        nbp = nbp+1
        nopara(nbp) = 'QUANTITE'
        typara(nbp) = 'K16'
    end if
    if (lmima) then
        nbp = nbp+1
        nopara(nbp) = 'EXTREMA'
        typara(nbp) = 'K8'
        nbp = nbp+1
        nopara(nbp) = 'MAILLE'
        typara(nbp) = 'K8'
        nbp = nbp+1
        nopara(nbp) = 'NOEUD'
        typara(nbp) = 'K8'
        nbp = nbp+1
        nopara(nbp) = 'CMP'
        typara(nbp) = 'K8'
        nbp = nbp+1
        nopara(nbp) = 'VALE'
        typara(nbp) = 'R'
    end if
!
    if (lmoygr) then
        nbp = nbp+1
        nopara(nbp) = 'CMP'
        typara(nbp) = 'K8'
        nbp = nbp+1
        nopara(nbp) = 'MOYENNE'
        typara(nbp) = 'R'
    end if
!
    if (niv .ge. 2) then
        do 1789, n1 = 1, nbp
            valk(1) = nopara(n1)
            valk(2) = typara(n1)
            call utmess('I', 'POSTRELE_10', nk=2, valk=valk)
1789        continue
            end if
!
! 3.2. ==> CREATION/INITIALISATION DE LA TABLE
!
            call tbcrsd(nomtab, 'G')
            call tbajpa(nomtab, nbp, nopara, typara)
!
            k24bid = nomtab(1:8)//'           .TITR'
            call titrea('T', nomtab, nomtab, k24bid, 'C', &
                        ' ', 0, 'G', '(1PE12.5)')
!
            call jedetr(nocmp)
!
            call jedema()
!
            end subroutine
