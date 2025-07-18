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

subroutine carbe3(charge)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterfort/aflrch.h"
#include "asterfort/afrela.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvem.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/lxlgut.h"
#include "asterfort/mgauss.h"
#include "asterfort/pmppr.h"
#include "asterfort/utbtab.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/int_to_char8.h"
#include "asterfort/char8_to_int.h"
!
    character(len=8) :: charge
!
!     TRAITER LE MOT CLE LIAISON_RBE3 DE AFFE_CHAR_MECA
!     ET ENRICHIR LA CHARGE (CHARGE) AVEC LES RELATIONS LINEAIRES
! IN/JXVAR : CHARGE : NOM D'UNE SD CHARGE
! ----------------------------------------------------------------------
    integer(kind=8) :: vali(2)
!
    character(len=4) :: typcoe, typval
    character(len=7) :: typcha
    character(len=8) :: k8bid, mode, noma, nomres, noemai, nomnoe, ddltrr(6)
    character(len=8) :: ddlcod, ddlmac(6), betaf
    character(len=15) :: coordo
    character(len=16) :: motfac, concep, nomcmd
    character(len=19) :: lisrel
    character(len=24) :: ddlstr, grouno, gromai
    complex(kind=8) :: betac
    integer(kind=8) :: ifm, niv, iret
    integer(kind=8) :: idxrbe, idxlig, idxcol, idxvec, idxnoe, idxgro, idxter
    integer(kind=8) :: idxddl
    integer(kind=8) :: posesc, posmai, cntlig, cntddl, cntnoe, inilig
    integer(kind=8) :: jlises, jcofes, jddles, jcescl, jcoore
    integer(kind=8) :: jnorel, jddl, jcmur, jcmuc, jcmuf, jdirec, jdime
    integer(kind=8) :: jnogro, jnoesc, jnzddl, jnznor
    integer(kind=8) :: nbrbe3, nbdles, nbcfes, nbddl, nblign, nbcol, nbgrou, nbent
    integer(kind=8) :: nbnoeu, nbdlma, maxesc, maxles, maxddl, dime
    aster_logical :: fincod, ddlesc(6), ddlmai(6), frstco, dime2d
    real(kind=8) :: rbid, coomai(3), cooesc(3), lc, norme, lcsqua, stws(6, 6)
    real(kind=8) :: cofesc, beta, x(6, 6)
    real(kind=8), pointer :: b(:) => null()
    real(kind=8), pointer :: s(:) => null()
    real(kind=8), pointer :: w(:) => null()
    real(kind=8), pointer :: xab(:) => null()
    real(kind=8), pointer :: vale(:) => null()
! ----------------------------------------------------------------------
!
    motfac = 'LIAISON_RBE3    '
    call getfac(motfac, nbrbe3)
!
    if (nbrbe3 .eq. 0) goto 999
!
    beta = 0.0d0
    betac = (1.0d0, 0.0d0)
!
    typcoe = 'REEL'
    typval = 'REEL'
    ddltrr(1) = 'DX'
    ddltrr(2) = 'DY'
    ddltrr(3) = 'DZ'
    ddltrr(4) = 'DRX'
    ddltrr(5) = 'DRY'
    ddltrr(6) = 'DRZ'
!
    call infniv(ifm, niv)
    call jemarq()
!
    if (niv .eq. 2) then
        write (ifm, *) 'NOMBRE RELATIONS RBE3 : ', nbrbe3
!       CALL JXVERI(' ')
    end if
!
    call getres(nomres, concep, nomcmd)
    call dismoi('TYPE_CHARGE', charge, 'CHARGE', repk=typcha)
    call dismoi('NOM_MODELE', charge, 'CHARGE', repk=mode)
    call dismoi('DIM_GEOM', mode, 'MODELE', repi=dime)
    dime2d = (dime .eq. 2)
    if (niv .eq. 2) then
        write (ifm, *) 'MODELE 2D : ', dime2d
    end if
!
    call dismoi('NOM_MAILLA', charge, 'CHARGE', repk=noma)
!
    grouno = noma//'.GROUPENO'
    coordo = noma//'.COORDO'
    call jeveuo(coordo//'    .VALE', 'L', vr=vale)
!
!     -- CALCUL DE MAXLES : NBRE DE TERMES MAXI D'UNE LISTE
!        DE GROUP_NO_ESCL OU DE NOEUD_ESCL
!        --------------------------------------------------
    maxles = 0
    do idxrbe = 1, nbrbe3
        call getvem(noma, 'GROUP_NO', motfac, 'GROUP_NO_ESCL', idxrbe, &
                    0, k8bid, nbgrou)
        maxles = max(maxles, -nbgrou)
        call getvem(noma, 'NOEUD', motfac, 'NOEUD_ESCL', idxrbe, &
                    0, k8bid, nbnoeu)
        maxles = max(maxles, -nbnoeu)
    end do
!
    if (maxles .eq. 0) then
        call utmess('F', 'MODELISA10_7')
    end if
    call wkvect('&&CARBE3.LISTESCL', 'V V K24', maxles, jlises)
!
!     -- CALCUL DE MAXDDL ET VERIFICATION DES NOEUDS ET GROUP_NO
!        MAXDDL EST LE NOMBRE MAXI DE NOEUDS IMPLIQUES DANS UNE
!        RELATION LINEAIRE
!        -------------------------------------------------------
    maxddl = 0
    maxesc = 0
    do idxrbe = 1, nbrbe3
!
        call getvtx(motfac, 'GROUP_NO_ESCL', iocc=idxrbe, nbval=maxles, vect=zk24(jlises), &
                    nbret=nbgrou)
        if (nbgrou .ne. 0) then
            do idxgro = 1, nbgrou
                call jelira(jexnom(grouno, zk24(jlises-1+idxgro)), 'LONUTI', nbnoeu)
                if (nbnoeu .eq. 0) then
                    call utmess('F', 'MODELISA10_8', sk=zk24(jlises-1+idxgro))
                end if
                maxesc = maxesc+nbnoeu
            end do
        else
            call getvtx(motfac, 'NOEUD_ESCL', iocc=idxrbe, nbval=0, nbret=nbnoeu)
            nbnoeu = -nbnoeu
            maxesc = maxesc+nbnoeu
        end if
!
        cntddl = 1
        cntddl = cntddl+maxesc*6
        maxddl = max(maxddl, cntddl)
    end do
!
    if (niv .eq. 2) then
        write (ifm, *) 'MAXESC : ', maxesc
        write (ifm, *) 'MAXDDL : ', maxddl
    end if
!
!     -- ALLOCATION DES TABLEAUX DE TRAVAIL
!     -------------------------------------
    lisrel = '&&CARBE3.RLLISTE'
    call wkvect('&&CARBE3.LISNOREL', 'V V K8', maxddl, jnorel)
    call wkvect('&&CARBE3.LISNZNOR', 'V V K8', maxddl, jnznor)
    call wkvect('&&CARBE3.LISNOESC', 'V V K8', maxesc, jnoesc)
    call wkvect('&&CARBE3.DDL  ', 'V V K8', maxddl, jddl)
    call wkvect('&&CARBE3.NZDDL', 'V V K8', maxddl, jnzddl)
    call wkvect('&&CARBE3.COEMUR', 'V V R', maxddl, jcmur)
    call wkvect('&&CARBE3.COEMUC', 'V V C', maxddl, jcmuc)
    call wkvect('&&CARBE3.COEMUF', 'V V K8', maxddl, jcmuf)
    call wkvect('&&CARBE3.DIRECT', 'V V R', 3*maxddl, jdirec)
    call wkvect('&&CARBE3.DIMENSION', 'V V I', maxddl, jdime)
    call wkvect('&&CARBE3.CESCL', 'V V L', maxesc*6, jcescl)
    call wkvect('&&CARBE3.COEFES', 'V V R', maxesc, jcofes)
    call wkvect('&&CARBE3.COOREL', 'V V R', maxesc*3, jcoore)
    call wkvect('&&CARBE3.DDLESCL', 'V V K24', maxesc, jddles)
!
!     BOUCLE SUR LES RELATIONS RBE3
!     -----------------------------------
    do idxrbe = 1, nbrbe3
        if (niv .eq. 2) then
            write (ifm, *) 'INDEX RELATION RBE3 : ', idxrbe
        end if
!
        call getvtx(motfac, 'GROUP_NO_MAIT', iocc=idxrbe, nbval=0, nbret=nbent)
        nbent = -nbent
        if (nbent .ne. 0) then
            call getvem(noma, 'GROUP_NO', motfac, 'GROUP_NO_MAIT', idxrbe, &
                        1, gromai, nbent)
            call jeveuo(jexnom(grouno, gromai), 'L', jnogro)
            call jelira(jexnom(grouno, gromai), 'LONUTI', nbent)
            if (nbent .ne. 1) then
                call utmess('F', 'MODELISA10_9', sk=gromai, si=nbent)
            end if
            noemai = int_to_char8(zi(jnogro-1+1))
        end if
!
        call getvtx(motfac, 'NOEUD_MAIT', iocc=idxrbe, nbval=0, nbret=nbent)
        if (nbent .ne. 0) then
            call getvem(noma, 'NOEUD', motfac, 'NOEUD_MAIT', idxrbe, &
                        1, noemai, nbent)
        end if
!
        posmai = char8_to_int(noemai)
        coomai(1) = vale(3*(posmai-1)+1)
        coomai(2) = vale(3*(posmai-1)+2)
        coomai(3) = vale(3*(posmai-1)+3)
        if (niv .eq. 2) then
            write (ifm, *) 'COORDS : ', coomai, ' DU NOEUD MAITRE : ', &
                noemai
        end if
!
        call getvtx(motfac, 'DDL_MAIT', iocc=idxrbe, nbval=6, vect=ddlmac, &
                    nbret=nbdlma)
!
        do idxlig = 1, 6
            ddlmai(idxlig) = .false.
        end do
!
        do idxlig = 1, nbdlma
            ddlcod = ddlmac(idxlig) (1:lxlgut(ddlmac(idxlig)))
            if (ddltrr(1) .eq. ddlcod) then
                ddlmai(1) = .true.
            else if (ddltrr(2) .eq. ddlcod) then
                ddlmai(2) = .true.
            else if (ddltrr(3) .eq. ddlcod) then
                ddlmai(3) = .true.
            else if (ddltrr(4) .eq. ddlcod) then
                ddlmai(4) = .true.
            else if (ddltrr(5) .eq. ddlcod) then
                ddlmai(5) = .true.
            else if (ddltrr(6) .eq. ddlcod) then
                ddlmai(6) = .true.
            end if
        end do
!
        call getvtx(motfac, 'GROUP_NO_ESCL', iocc=idxrbe, nbval=0, nbret=nbgrou)
        if (nbgrou .ne. 0) then
            nbgrou = -nbgrou
            nbnoeu = 0
            call getvtx(motfac, 'GROUP_NO_ESCL', iocc=idxrbe, nbval=nbgrou, vect=zk24(jlises), &
                        nbret=nbent)
            cntnoe = 0
            do idxgro = 1, nbgrou
                call jeveuo(jexnom(grouno, zk24(jlises-1+idxgro)), 'L', jnogro)
                call jelira(jexnom(grouno, zk24(jlises-1+idxgro)), 'LONUTI', nbent)
                nbnoeu = nbnoeu+nbent
                do idxnoe = 1, nbent
                    cntnoe = cntnoe+1
                    nomnoe = int_to_char8(zi(jnogro-1+idxnoe))
                    zk8(jnoesc+cntnoe-1) = nomnoe
                end do
            end do
        end if
!
        call getvtx(motfac, 'NOEUD_ESCL', iocc=idxrbe, nbval=0, nbret=nbent)
        if (nbent .ne. 0) then
            nbnoeu = -nbent
            call getvtx(motfac, 'NOEUD_ESCL', iocc=idxrbe, nbval=nbnoeu, vect=zk8(jnoesc), &
                        nbret=nbent)
        end if
!
        if (niv .eq. 2) then
            write (ifm, *) 'LISTE DES ', nbnoeu, ' NOEUDS ESCLAVES'
            write (ifm, *) (zk8(jnoesc+idxlig-1), idxlig=1, nbnoeu)
        end if
!
        call getvtx(motfac, 'DDL_ESCL', iocc=idxrbe, nbval=nbnoeu, vect=zk24(jddles), &
                    nbret=nbddl)
!
        if (nbddl .ne. 1 .and. nbddl .ne. nbnoeu) then
            vali(1) = nbddl
            vali(2) = nbnoeu
            call utmess('F', 'MODELISA10_10', ni=2, vali=vali)
        end if
!
!       BOUCLE SUR LES NOEUDS ESCLAVES POUR EXTRAIRE LES DDLS
!       -----------------------------------------------------
        nbdles = 0
        do idxnoe = 1, nbnoeu
            if (nbddl .gt. 1 .or. idxnoe .eq. 1) then
                if (nbddl .eq. 1) then
                    ddlstr = zk24(jddles-1+1)
                else
                    ddlstr = zk24(jddles-1+idxnoe)
                end if
!
!           EXTRACTION DDL_ESCL
!           -------------------------------------------------------
                do idxlig = 1, 6
                    ddlesc(idxlig) = .false.
                end do
!
                idxcol = 1
                do idxlig = 1, lxlgut(ddlstr)
                    if (ddlstr(idxlig:idxlig) .eq. '-') then
                        ddlcod = ddlstr(idxcol:idxlig-1)
                        idxcol = idxlig+1
                        fincod = .true.
                    else if (idxlig .eq. lxlgut(ddlstr)) then
                        ddlcod = ddlstr(idxcol:idxlig)
                        fincod = .true.
                    else
                        fincod = .false.
                    end if
                    if (fincod) then
                        if (ddltrr(1) .eq. ddlcod) then
                            ddlesc(1) = .true.
                        else if (ddltrr(2) .eq. ddlcod) then
                            ddlesc(2) = .true.
                        else if (ddltrr(3) .eq. ddlcod) then
                            ddlesc(3) = .true.
                        else if (ddltrr(4) .eq. ddlcod) then
                            ddlesc(4) = .true.
                        else if (ddltrr(5) .eq. ddlcod) then
                            ddlesc(5) = .true.
                        else if (ddltrr(6) .eq. ddlcod) then
                            ddlesc(6) = .true.
                        else
                            call utmess('F', 'MODELISA10_11', sk=ddlcod)
                        end if
                    end if
                end do
            end if
            do idxlig = 1, 6
                if (ddlesc(idxlig)) then
                    nbdles = nbdles+1
                    zk8(jnorel-1+nbdles) = zk8(jnoesc-1+idxnoe)
                    zk8(jddl-1+nbdles) = ddltrr(idxlig)
                end if
                zl(jcescl-1+(idxnoe-1)*6+idxlig) = ddlesc(idxlig)
            end do
        end do
!
        if (niv .eq. 2) then
            write (ifm, *) 'NOMBRE DDL NOEUDS ESCLAVES : ', nbdles
        end if
!
!       BOUCLE SUR LES NOEUDS ESCLAVES POUR CALCULER Lc
!       -----------------------------------------------
        lc = 0
        do idxnoe = 1, nbnoeu
            posesc = char8_to_int(zk8(jnoesc-1+idxnoe))
            cooesc(1) = vale(3*(posesc-1)+1)
            cooesc(2) = vale(3*(posesc-1)+2)
            cooesc(3) = vale(3*(posesc-1)+3)
            zr(jcoore-1+3*(idxnoe-1)+1) = cooesc(1)-coomai(1)
            zr(jcoore-1+3*(idxnoe-1)+2) = cooesc(2)-coomai(2)
            zr(jcoore-1+3*(idxnoe-1)+3) = cooesc(3)-coomai(3)
            norme = zr(jcoore-1+3*(idxnoe-1)+1)*zr(jcoore-1+3*(idxnoe-1) &
                                                   +1)+zr(jcoore-1+3*(idxnoe-1)+2)* &
                    zr(jcoore-1+3*(idxnoe-1)+ &
                       2)+zr(jcoore-1+3*(idxnoe-1)+3)* &
                    zr(jcoore-1+3*(idxnoe-1)+3)
            if (norme .ne. 0.0d0) then
                norme = sqrt(norme)
            end if
            lc = lc+norme
        end do
        lc = lc/nbnoeu
        if (niv .eq. 2) then
            write (ifm, *) 'LC : ', lc
        end if
        lcsqua = lc*lc
!
!       BOUCLE SUR LES NOEUDS ESCLAVES POUR CALCULER W
!       -------------------------------------------------------
        AS_ALLOCATE(vr=w, size=nbdles*nbdles)
        call getvr8(motfac, 'COEF_ESCL', iocc=idxrbe, nbval=nbnoeu, vect=zr(jcofes), &
                    nbret=nbcfes)
        if (nbcfes .lt. 0) then
            nbcfes = -nbcfes
        end if
!
        if (nbddl .ne. 1 .and. nbcfes .ne. 1 .and. nbddl .ne. nbcfes) then
            vali(1) = nbddl
            vali(2) = nbcfes
            call utmess('F', 'MODELISA10_12', ni=2, vali=vali)
        end if
!
        if (nbcfes .ne. 1 .and. nbcfes .ne. nbnoeu) then
            vali(1) = nbcfes
            vali(2) = nbnoeu
            call utmess('F', 'MODELISA10_13', ni=2, vali=vali)
        end if
!
        nbcol = 0
        inilig = 0
        do idxnoe = 1, nbnoeu
            if (nbcfes .eq. 1) then
                cofesc = zr(jcofes-1+1)
            else
                cofesc = zr(jcofes-1+idxnoe)
            end if
!
            frstco = .true.
            cntlig = 0
            do idxcol = 1, 6
                if (zl(jcescl-1+(idxnoe-1)*6+idxcol)) then
                    nbcol = nbcol+1
                    nblign = inilig
                    do idxlig = 1, 6
                        if (zl(jcescl-1+(idxnoe-1)*6+idxlig)) then
                            nblign = nblign+1
                            if (frstco) then
                                cntlig = cntlig+1
                            end if
                            idxvec = nbdles*(nbcol-1)+nblign
                            if (idxlig .ne. idxcol) then
                                w(idxvec) = 0.0d0
                            else if (idxcol .le. 3) then
                                w(idxvec) = cofesc
                            else
                                w(idxvec) = cofesc*lcsqua
                            end if
                        end if
                    end do
                    frstco = .false.
                end if
            end do
            inilig = inilig+cntlig
        end do
!
        if (niv .eq. 2) then
            write (ifm, *) 'IMPRESSION W'
            do idxlig = 1, nbdles
                write (ifm, *) 'LIGNE : ', idxlig, (w(idxlig+ &
                                                      nbdles*(idxcol-1)), idxcol=1, nbdles)
            end do
        end if
!
!       BOUCLE SUR LES NOEUDS ESCLAVES POUR CALCULER S
!       -------------------------------------------------------
        nblign = 0
        inilig = 0
        AS_ALLOCATE(vr=s, size=nbdles*6)
        do idxnoe = 1, nbnoeu
            frstco = .true.
            cntlig = 0
            do idxcol = 1, 6
                nblign = inilig
                do idxlig = 1, 6
                    if (zl(jcescl-1+(idxnoe-1)*6+idxlig)) then
                        nblign = nblign+1
                        if (frstco) then
                            cntlig = cntlig+1
                        end if
                        idxvec = nbdles*(idxcol-1)+nblign
                        if ((idxlig .le. 3 .and. idxcol .le. 3) .or. &
                            (idxlig .ge. 4 .and. idxcol .ge. 4)) then
                            if (idxlig .eq. idxcol) then
                                s(idxvec) = 1.0d0
                            else
                                s(idxvec) = 0.0d0
                            end if
                        else if (idxlig .ge. 4 .and. idxcol .le. 3) then
                            s(idxvec) = 0.0d0
                        else
                            if (idxlig .eq. 1 .and. idxcol .eq. 5) then
                                s(idxvec) = zr(jcoore-1+3*(idxnoe-1)+3)
                            else if (idxlig .eq. 1 .and. idxcol .eq. 6) &
                                then
                                s(idxvec) = -zr(jcoore-1+3*(idxnoe-1)+2)
                            else if (idxlig .eq. 2 .and. idxcol .eq. 4) &
                                then
                                s(idxvec) = -zr(jcoore-1+3*(idxnoe-1)+3)
                            else if (idxlig .eq. 2 .and. idxcol .eq. 6) &
                                then
                                s(idxvec) = zr(jcoore-1+3*(idxnoe-1)+1)
                            else if (idxlig .eq. 3 .and. idxcol .eq. 4) &
                                then
                                s(idxvec) = zr(jcoore-1+3*(idxnoe-1)+2)
                            else if (idxlig .eq. 3 .and. idxcol .eq. 5) &
                                then
                                s(idxvec) = -zr(jcoore-1+3*(idxnoe-1)+1)
                            else
                                s(idxvec) = 0.0d0
                            end if
                        end if
                    end if
                end do
                frstco = .false.
            end do
            inilig = inilig+cntlig
        end do
!
        if (niv .eq. 2) then
            write (ifm, *) 'IMPRESSION S'
            do idxlig = 1, nbdles
                write (ifm, *) 'LIGNE : ', idxlig, (s(idxlig+ &
                                                      nbdles*(idxcol-1)), idxcol=1, 6)
            end do
        end if
!
        AS_ALLOCATE(vr=xab, size=nbdles*6)
        call utbtab('ZERO', nbdles, 6, w, s, &
                    xab, stws)
!
        if (niv .eq. 2) then
            write (ifm, *) 'IMPRESSION MATRICE MGAUSS'
            do idxlig = 1, 6
                write (ifm, *) 'LIGNE : ', idxlig, (stws(idxcol, &
                                                         idxlig), idxcol=1, 6)
            end do
        end if
!
        do idxlig = 1, 6
            do idxcol = 1, 6
                if (idxcol .eq. idxlig) then
                    x(idxcol, idxlig) = 1.0d0
                else
                    x(idxcol, idxlig) = 0.0d0
                end if
            end do
        end do
!
        call mgauss('NFSP', stws, x, 6, 6, &
                    6, rbid, iret)
!
        if (iret .ne. 0) then
            ASSERT(.false.)
        end if
!
        if (niv .eq. 2) then
            write (ifm, *) 'IMPRESSION MATRICE X'
            do idxlig = 1, 6
                write (ifm, *) 'LIGNE : ', idxlig, (x(idxcol, idxlig), &
                                                    idxcol=1, 6)
            end do
        end if
!
        call pmppr(s, nbdles, 6, -1, w, &
                   nbdles, nbdles, -1, xab, 6, &
                   nbdles)
!
        AS_DEALLOCATE(vr=w)
        AS_DEALLOCATE(vr=s)
!
        AS_ALLOCATE(vr=b, size=nbdles*6)
!
        call pmppr(x, 6, 6, -1, xab, &
                   6, nbdles, 1, b, 6, &
                   nbdles)
!
        AS_DEALLOCATE(vr=xab)
!
        do idxlig = 1, 6
            if (ddlmai(idxlig)) then
                idxter = 1
                zk8(jnznor-1+idxter) = noemai
                zr(jcmur-1+idxter) = -1.0d0
                zk8(jnzddl-1+idxter) = ddltrr(idxlig)
                idxddl = 1
                do idxcol = 1, nbdles
                    idxvec = 6*(idxcol-1)+idxlig
                    if (b(idxvec) .ne. 0) then
                        idxter = idxter+1
                        zk8(jnzddl-1+idxter) = zk8(jddl-1+idxddl)
                        zk8(jnznor-1+idxter) = zk8(jnorel-1+idxddl)
                        zr(jcmur-1+idxter) = b(idxvec)
                    end if
                    idxddl = idxddl+1
                end do
                call afrela(zr(jcmur), zc(jcmuc), zk8(jnzddl), zk8(jnznor), zi(jdime), &
                            zr(jdirec), idxter, beta, betac, betaf, &
                            typcoe, typval, 0.d0, lisrel)
            end if
        end do
!
        AS_DEALLOCATE(vr=b)
!
    end do
!
!     -- AFFECTATION DE LA LISTE_RELA A LA CHARGE :
!     ---------------------------------------------
    call aflrch(lisrel, charge, 'NLIN')
!
    call jedema()
!
!     IF (NIV.EQ.2) THEN
!       CALL JXVERI(' ')
!     ENDIF
!
999 continue
!
end subroutine
