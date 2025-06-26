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
subroutine ircerl(ifi, nbel, ligrel, nbgrel, longr, &
                  ncmpmx, vale, nomcmp, nomel, loc, &
                  celd, connex, point, nomnos, nbcmpt, &
                  nucmpu, nbnot, numnoe, nbmat, nummai, &
                  lsup, borsup, linf, borinf, lmax, &
                  lmin, lcor, ndim, coor, nolili, &
                  formr, ncmpv, nucmp)
! aslint: disable=W1501,W1504
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8vide.h"
#include "asterfort/codent.h"
#include "asterfort/dgmode.h"
#include "asterfort/digdel.h"
#include "asterfort/exisdg.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/lxlgut.h"
#include "asterfort/nbec.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: ifi, nbel, ligrel(*), nbgrel, longr(*), ncmpmx, nbnot, nbcmpt
    integer(kind=8) :: nucmpu(*), celd(*), connex(*), point(*), numnoe(*), nbmat, ndim
    integer(kind=8) :: nummai(*), ncmpv, nucmp(*)
    real(kind=8) :: borsup, borinf, coor(*), vale(*)
    character(len=*) :: nomcmp(*), nomel(*), loc, nomnos(*), formr
    character(len=19) :: nolili
    aster_logical :: lsup, linf, lmax, lmin, lcor
!        ECRITURE D'UN CHAMELEM SUR LISTING
!        A VALEURS REELLES
!  ENTREE:
!     IFI   : UNITE LOGIQUE DU FICHIER
!     NBEL  : NOMBRE D'ELEMENTS DU LIGREL ( DU MAILLAGE)
!     LIGREL: LIGREL COMPLET
!     NBGREL: NOMBRE DE GRELS
!     LONGR : POINTEUR DE LONGUEUR DE LIGREL
!     NCMPMX: NOMBRE MAXI DE CMP DE LA GRANDEUR NOMGD
!     VALE  : VALEURS DU CHAM_NO
!     NOMCMP: NOMS DES CMP
!     NOMEL : NOMS DES MAILLES SUPPORTS DES ELEMENTS
!     LOC   : LOCALISATION DES VALEURS (ELNO OU ELGA OU ELEM)
!     CELD  : DESCRIPTEUR DU CHAM_ELEM (MODES LOCAUX,ADRESSES->.CELV)
!     CONNEX: CONNECTIVITES DES MAILLES
!     POINT : POINTEUR DANS LES CONNECTIVITES
!     NOMNOS: NOMS DES NOEUDS
!     NBCMPT: NOMBRE DE COMPOSANTES A IMPRIMER
!     NUCMPU: NUMEROS DES COMPOSANTES A IMPRIMER
!     NBMAT : NOMBRE DE MAILLES OU ON DESIRE IMPRIMER LE CHAMELEM
!     NUMMAI: NUMEROS DES MAILLES OU ON DESIRE IMPRIMER LE CHAMELEM
!     LSUP  : =.TRUE.  INDIQUE PRESENCE D'UNE BORNE SUPERIEURE
!     BORSUP: VALEUR DE LA BORNE SUPERIEURE
!     LINF  : =.TRUE.  INDIQUE PRESENCE D'UNE BORNE INFERIEURE
!     BORINF: VALEUR DE LA BORNE INFERIEURE
!     LMAX  : =.TRUE.  INDIQUE IMPRESSION VALEUR MAXIMALE
!     LMIN  : =.TRUE.  INDIQUE IMPRESSION VALEUR MINIMALE
!     LCOR  : =.TRUE.  IMPRESSION DES COORDONNEES DE NOEUDS DEMANDEE
!     NDIM  : DIMENSION DU PROBLEME
!     COOR  : TABLEAU DES COORDONNEES DE NOEUDS
!     NOLILI: NOM DU LIGREL
!     FORMR : FORMAT D'ECRITURE DES REELS SUR "RESULTAT"
!     ------------------------------------------------------------------
    integer(kind=8) :: imodel, ilong
    real(kind=8) :: rundf
    character(len=7) :: cbid
    character(len=8) :: nomno, nomcp, forcmp, nomcor(3)
    character(len=10) :: format
    character(len=24) :: nrepe
    character(len=50) :: fmt, fmv, fmt1, fmt2, form1
    aster_logical :: limpr
!     ------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, i2, iachml, iad, iadr, iaec, icm
    integer(kind=8) :: icmp, icmp2, icoe, icoef, icoef2, icomax, icomp2
    integer(kind=8) :: icou, id, iel, ielg, if, igre, igrel
    integer(kind=8) :: iino, ilig, imai, imail, imax, imin, in
    integer(kind=8) :: inmax, inmin, inom, inop, inot, inu, ipca
    integer(kind=8) :: ipo2, ipoin, ipoin1, iposg, iposv, irepe, ires
    integer(kind=8) :: iva, ival, ivmax, ivmin, j, jco, jmod
    integer(kind=8) :: lgr, mode, modsau, nbcpt, nbno, ncmp
    integer(kind=8) :: ncmp2, ncmpp, ncou, nec, npcalc, nsca, nscal
    integer(kind=8) :: nuno
!-----------------------------------------------------------------------
    call jemarq()
    rundf = r8vide()
    nomcor(1) = 'X'
    nomcor(2) = 'Y'
    nomcor(3) = 'Z'
    format = formr
    lgr = lxlgut(format)
    id = 0
    if = 0
    call jeveuo('&CATA.TE.MODELOC', 'L', imodel)
    call jeveuo(jexatr('&CATA.TE.MODELOC', 'LONCUM'), 'L', ilong)
    do i = 1, lgr-1
        if (format(i:i) .eq. 'D' .or. format(i:i) .eq. 'E' .or. format(i:i) .eq. 'F' .or. &
            format(i:i) .eq. 'G') then
            id = i+1
            goto 2
        end if
        if (format(i:i) .eq. '.') then
            if = i-1
            goto 2
        end if
2       continue
    end do
    if (id .ne. 0 .and. if .ge. id) then
        forcmp = 'A'//format(id:if)
    else
        forcmp = 'A12'
    end if
!
!  -- DETERMINATION DU NOMBRE MAXIMUM DE SOUS_POINTS ---
    icomax = 0
    do igre = 1, nbgrel
        icoef = max(1, celd(4))
        if (icoef .gt. icomax) icomax = icoef
    end do
    ncmp = ncmpv
    if (ncmp .gt. 0) then
        ncmp = 0
        do i = 1, ncmpv
            if (nucmp(i) .le. icomax) then
                ncmp = ncmp+1
            else
                call codent(nucmp(i), 'G', cbid)
                nomcp = 'V'//cbid
                call utmess('A', 'PREPOST_74', sk=nomcp)
            end if
        end do
        if (ncmp .eq. 0) then
            call utmess('A', 'PREPOST_75')
            goto 999
        end if
        icomax = ncmp
    end if
!
    if (lmax .or. lmin) then
        call jedetr('&&IRCERL.NCMT')
        call wkvect('&&IRCERL.NCMT', 'V V K16', ncmpmx*icomax, inot)
        do i = 1, ncmpmx
            if (icomax .gt. 1 .or. ncmp .ge. 1) then
                do jco = 1, icomax
                    if (ncmp .gt. 0) then
                        call codent(nucmp(jco), 'G', cbid)
                    else
                        call codent(jco, 'G', cbid)
                    end if
                    nomcp = nomcmp(i)
                    zk16(inot-1+(i-1)*icomax+jco) = 'V'//cbid
                end do
            else
                zk16(inot-1+i) = nomcmp(i)
            end if
        end do
    end if
    if (lmax) then
        call jedetr('&&IRCERL.MAX')
        call wkvect('&&IRCERL.MAX', 'V V R', ncmpmx*icomax, imax)
        call jedetr('&&IRCERL.MAIMAX')
        call wkvect('&&IRCERL.MAIMAX', 'V V K8', ncmpmx*icomax, inmax)
        call jedetr('&&IRCERL.NBVMAX')
        call wkvect('&&IRCERL.NBVMAX', 'V V I', ncmpmx*icomax, ivmax)
        do i = 1, ncmpmx*icomax
            zr(imax-1+i) = rundf
        end do
    end if
    if (lmin) then
        call jedetr('&&IRCERL.MIN')
        call wkvect('&&IRCERL.MIN', 'V V R', ncmpmx*icomax, imin)
        call jedetr('&&IRCERL.MAIMIN')
        call wkvect('&&IRCERL.MAIMIN', 'V V K8', ncmpmx*icomax, inmin)
        call jedetr('&&IRCERL.NBVMIN')
        call wkvect('&&IRCERL.NBVMIN', 'V V I', ncmpmx*icomax, ivmin)
        do i = 1, ncmpmx*icomax
            zr(imin-1+i) = rundf
        end do
    end if
    if (loc .eq. 'ELGA' .or. loc .eq. 'ELEM' .or. .not. lcor) ndim = 0
    nrepe = nolili//'.REPE'
    call jeveuo(nrepe, 'L', irepe)
    if (nbmat .ne. 0) nbel = nbmat
    modsau = 0
    do imai = 1, nbel
        if (nbmat .ne. 0) then
            imail = nummai(imai)
        else
            imail = imai
        end if
        igrel = zi(irepe+2*(imail-1)+1-1)
        if (igrel .eq. 0) goto 12
        ielg = zi(irepe+2*(imail-1)+2-1)
        mode = celd(celd(4+igrel)+2)
        if (mode .eq. 0) goto 12
        if (mode .ne. modsau) then
            ipoin1 = longr(igrel)
            jmod = imodel+zi(ilong-1+mode)-1
            nec = nbec(zi(jmod-1+2))
            call jedetr('&&IRCERL.ENT_COD')
            call wkvect('&&IRCERL.ENT_COD', 'V V I', nec, iaec)
            call dgmode(mode, imodel, ilong, nec, zi(iaec))
            iad = celd(celd(4+igrel)+8)
            nscal = digdel(mode)
            icoef = max(1, celd(4))
            nsca = nscal*icoef
            icoef2 = icoef
            if (ncmp .gt. 0) icoef2 = ncmp
            ncmpp = 0
            ncmp2 = 0
!
! -- IPOSG : POSITION DE LA COMPOSANTE DANS LA GRANDEUR
! -- IPOSV : POSITION DE LA COMPOSANTE DANS LE .VALE
!
            call jedetr('&&IRCERL.POSG')
            call wkvect('&&IRCERL.POSG', 'V V I', ncmpmx*icoef2, iposg)
            call jedetr('&&IRCERL.POSV')
            call wkvect('&&IRCERL.POSV', 'V V I', ncmpmx, iposv)
            call jedetr('&&IRCERL.COEF')
            call wkvect('&&IRCERL.COEF', 'V V I', ncmpmx*icoef2, icoe)
            call jedetr('&&IRCERL.NCMP')
            call wkvect('&&IRCERL.NCMP', 'V V K16', ncmpmx*icoef2, inom)
            if (lsup .or. linf) then
                call jedetr('&&IRCERL.NCPP')
                call wkvect('&&IRCERL.NCPP', 'V V K16', ncmpmx*icoef2, inop)
                call jedetr('&&IRCERL.PO2')
                call wkvect('&&IRCERL.PO2', 'V V I', ncmpmx*icoef2, ipo2)
            end if
            call jedetr('&&IRCERL.VAL')
            call wkvect('&&IRCERL.VAL', 'V V R', ncmpmx*icoef2, ival)
            do i = 1, ncmpmx*icoef2
                zi(iposg-1+i) = 0
            end do
            do i = 1, ncmpmx
                zi(iposv-1+i) = 0
            end do
            do i = 1, ncmpmx
                if (exisdg(zi(iaec), i)) then
                    ncmpp = ncmpp+1
                    if (nbcmpt .ne. 0) then
                        do icm = 1, nbcmpt
                            icmp2 = nucmpu(icm)
                            if (i .eq. icmp2) then
                                ncmp2 = ncmp2+1
                                do jco = 1, icoef2
                                    zi(iposg-1+(icm-1)*icoef2+jco) = i
                                end do
                                zi(iposv-1+icm) = ncmpp
                            end if
                        end do
                    else
                        do jco = 1, icoef2
                            zi(iposg-1+(ncmpp-1)*icoef2+jco) = i
                        end do
                    end if
                end if
            end do
            if (nbcmpt .eq. 0) ncmp2 = ncmpp
            npcalc = nscal/ncmpp
!
! --- RETASSAGE DU TABLEAU DES POSITIONS DES COMPOSANTES DANS GRANDEUR-
!
            if (nbcmpt .ne. 0) then
                i2 = 0
                do i = 1, nbcmpt*icoef2
                    if (zi(iposg-1+i) .ne. 0) then
                        i2 = i2+1
                        zi(iposg-1+i2) = zi(iposg-1+i)
                    end if
                end do
            end if
!
! --- STOCKAGE DES NOMS DE COMPOSANTES ---
            do i = 1, ncmp2
                if (icoef2 .gt. 1 .or. ncmp .ge. 1) then
                    do jco = 1, icoef2
                        if (ncmp .gt. 0) then
                            call codent(nucmp(jco), 'G', cbid)
                        else
                            call codent(jco, 'G', cbid)
                        end if
                        nomcp = nomcmp(zi(iposg-1+i))
                        zk16(inom-1+(i-1)*icoef2+jco) = 'V'//cbid
                    end do
                else
                    zk16(inom-1+i) = nomcmp(zi(iposg-1+i))
                end if
            end do
!
! --- CREATION DES FORMATS D'ECRITURE ---
!
            if (.not. lmax .and. .not. lmin) then
                ilig = (ncmp2*icoef2+ndim)/6
                ires = (ncmp2*icoef2+ndim)-ilig*6
                fmt = ' '
                fmv = ' '
                if (ires .ne. 0) then
                    fmt = '(1X,A8,6(1X,'//forcmp//'),30(/,9X,6(1X,'//forcmp//')))'
                    if (loc .eq. 'ELNO') then
                        fmv = '(1X,A8,6(1X,'//format//'),30(/,9X,6(1X,'//format//')))'
                    else if (loc .eq. 'ELGA') then
                        fmv = '(2X,I7,6(1X,'//format//'),30(/,9X,6(1X,'//format//')))'
                    else if (loc .eq. 'ELEM') then
                        fmv = '(9X,6(1X,'//format//'),30(/,9X,6(1X,'//format//')))'
                    end if
                else if (ires .eq. 0 .and. ilig .eq. 1) then
                    fmt = '(1X,A8,6(1X,'//forcmp//'))'
                    if (loc .eq. 'ELNO') then
                        fmv = '(1X,A8,6(1X,'//format//'))'
                    else if (loc .eq. 'ELGA') then
                        fmv = '(2X,I7,6(1X,'//format//'))'
                    else if (loc .eq. 'ELEM') then
                        fmv = '(9X,6(1X,'//format//'))'
                    end if
                else
                    write (fmt, '(A,A8,A,I2,A,A8,A)') '(1X,A8,6(1X,', forcmp,&
     &                   '),', (ilig-1), '(/,9X,6(1X,', forcmp, ')))'
                    if (loc .eq. 'ELNO') then
                        write (fmv, '(A,A10,A,I2,A,A10,A)') '(1X,A8,6(1X,',&
     &               format, '),', (ilig-1), '(/,9X,6(1X,', format, ')))'
                    else if (loc .eq. 'ELGA') then
                        write (fmv, '(A,A10,A,I2,A,A10,A)') '(2X,I7,6(1X,',&
     &               format, '),', (ilig-1), '(/,9X,6(1X,', format, ')))'
                    else if (loc .eq. 'ELEM') then
                        write (fmv, '(A,A10,A,I2,A,A10,A)') '(9X,6(1X,', &
                            format, '),', (ilig-1), '(/,9X,6(1X,', format, &
                            ')))'
                    end if
                end if
            end if
        end if
!
! --- BOUCLE SUR LES ELEMENTS ---
!
        iel = ligrel(ipoin1+ielg-1)
        limpr = .true.
        if (.not. lsup .and. .not. linf .and. .not. lmax .and. .not. lmin) then
            if (ndim .eq. 0) then
                write (ifi, fmt) nomel(iel), (zk16(inom-1+i) (1:11), &
                                              i=1, icoef2*ncmp2)
            else
                write (ifi, fmt) nomel(iel), (nomcor(i), i=1, ndim), &
                    (zk16(inom-1+i) (1:11), i=1, icoef2*ncmp2)
            end if
        end if
        iachml = iad+nsca*(ielg-1)
        if (loc .eq. 'ELGA' .or. loc .eq. 'ELEM') then
            do ipca = 1, npcalc
                j = iachml-1+ncmpp*icoef*(ipca-1)
                if (nbcmpt .eq. 0) then
                    do i = 1, ncmp2
                        if (ncmp .gt. 0) then
                            do jco = 1, icoef2
                                zr(ival-1+(i-1)*icoef2+jco) = vale(j+i+ &
                                                                   (nucmp(jco)-1)*ncmpp)
                                zi(icoe-1+(i-1)*icoef2+jco) = jco
                            end do
                        else
                            do jco = 1, icoef2
                                zr(ival-1+(i-1)*icoef2+jco) = vale(j+i+ &
                                                                   (jco-1)*ncmpp)
                                zi(icoe-1+(i-1)*icoef2+jco) = jco
                            end do
                        end if
                    end do
                else
                    do i = 1, ncmp2
                        inu = zi(iposv-1+i)
                        if (ncmp .gt. 0) then
                            do jco = 1, icoef2
                                zr(ival-1+(i-1)*icoef2+jco) = vale(j+ &
                                                                   inu+(nucmp(jco)-1)*ncmpp)
                                zi(icoe-1+(i-1)*icoef2+jco) = jco
                            end do
                        else
                            do jco = 1, icoef2
                                zr(ival-1+(i-1)*icoef2+jco) = vale(j+ &
                                                                   inu+(jco-1)*ncmpp)
                                zi(icoe-1+(i-1)*icoef2+jco) = jco
                            end do
                        end if
                    end do
                end if
!
! --  TRI DES COMPOSANTES DANS L'INTERVALLE BORINF,BORSUP
!
                if (lsup .or. linf) then
                    do iva = 1, icoef2*ncmp2
                        if (lsup) then
                            if ((zr(ival-1+iva)-borsup) .gt. 0.d0) zi(icoe-1+iva) = 0
                        end if
                        if (linf) then
                            if ((zr(ival-1+iva)-borinf) .lt. 0.d0) zi(icoe-1+iva) = 0
                        end if
                    end do
!
! --- RETASSAGE POUR IMPRIMER COMPOSANTES PRESENTES DANS L'INTERVALLE --
!
                    icomp2 = 0
                    do i = 1, icoef2*ncmp2
                        if (zi(icoe-1+i) .ne. 0) then
                            icomp2 = icomp2+1
                            zi(icoe-1+icomp2) = zi(icoe-1+i)
                            zi(ipo2-1+icomp2) = zi(iposg-1+i)
                            zr(ival-1+icomp2) = zr(ival-1+i)
                            zk16(inop-1+icomp2) = zk16(inom-1+i)
                        end if
                    end do
                    if (icomp2 .eq. 0) goto 16
!
! -- IMPRESSION ----
!
                    if (.not. lmax .and. .not. lmin) then
                        ilig = (icomp2)/6
                        ires = (icomp2)-ilig*6
                        fmt1 = ' '
                        fmt2 = ' '
                        if (loc .eq. 'ELGA') then
                            if (ires .ne. 0) then
                                fmt1 = '(9X,6(1X,'//forcmp//'),30(/,9X,6(1X,'//forcmp//')))'
                                fmt2 = '(2X,I7,6(1X,'//format//'),30(/,9X,6(1X,'//format//')))'
                            else if (ires .eq. 0 .and. ilig .eq. 1) then
                                fmt1 = '(9X,6(1X,'//forcmp//'))'
                                fmt2 = '(2X,I7,6(1X,'//format//'))'
                            else
                                write (fmt1, '(A,A8,A,I2,A,A8,A)') '(1X,A8,6(1X,',&
     &                   forcmp, '),', (ilig-1), '(/,9X,6(1X,', forcmp, ')))'
                                write (fmt2, '(A,A10,A,I2,A,A10,A)') '(2X,I7,6(1X,',&
     &                   format, '),', (ilig-1), '(/,9X,6(1X,', format, ')))'
                            end if
                        else
                            if (ires .ne. 0) then
                                fmt1 = '(9X,6(1X,'//forcmp//'),30(/,9X,6(1X,'//forcmp//')))'
                                fmt2 = '(9X,6(1X,'//format//'),30(/,9X,6(1X,'//format//')))'
                            else if (ires .eq. 0 .and. ilig .eq. 1) then
                                fmt1 = '(9X,6(1X,'//forcmp//'))'
                                fmt2 = '(9X,6(1X,'//format//'))'
                            else
                                write (fmt1, '(A,A8,A,I2,A,A8,A)') '(1X,A8,6(1X,',&
     &                   forcmp, '),', (ilig-1), '(/,9X,6(1X,', forcmp, ')))'
                                write (fmt2, '(A,A10,A,I2,A,A10,A)') '(9X,6(1X,',&
     &                   format, '),', (ilig-1), '(/,9X,6(1X,', format, ')))'
                            end if
                        end if
                        if (lsup .or. linf) then
                            if (limpr) then
                                write (ifi, '(A,I2,A)') nomel(iel)
                                limpr = .false.
                            end if
                        end if
                        if (loc .eq. 'ELGA') then
                            write (ifi, fmt1) (zk16(inop-1+i) (1:11), i=1, &
                                               icomp2)
                            write (ifi, fmt2) ipca, (zr(ival-1+icmp), &
                                                     icmp=1, icomp2)
                        else
                            write (ifi, fmt1) (zk16(inop-1+i) (1:11), i=1, &
                                               icomp2)
                            write (ifi, fmt2) (zr(ival-1+icmp), icmp=1, &
                                               icomp2)
                        end if
                    end if
                    nbcpt = icomp2
                else
                    if (.not. lmax .and. .not. lmin) then
                        if (loc .eq. 'ELGA') then
                            write (ifi, fmv) ipca, (zr(ival-1+icmp), &
                                                    icmp=1, icoef2*ncmp2)
                        else
                            write (ifi, fmv) (zr(ival-1+icmp), icmp=1, &
                                              icoef2*ncmp2)
                        end if
                    end if
                    nbcpt = icoef2*ncmp2
                end if
!
! -- RECHERCHE DE LA VALEUR MAXIMALE ---
!
                if (lmax) then
                    do i = 1, nbcpt
                        if (lsup .or. linf) then
                            iadr = (zi(ipo2-1+i)-1)*icoef2+zi(icoe-1+i)
                        else
                            iadr = (zi(iposg-1+i)-1)*icoef2+zi(icoe-1+i)
                        end if
                        if (zr(imax-1+iadr) .eq. rundf) then
                            zr(imax-1+iadr) = zr(ival-1+i)
                            zk8(inmax-1+iadr) = nomel(iel)
                            zi(ivmax-1+iadr) = 1
                        else if (zr(ival-1+i) .gt. zr(imax-1+iadr)) then
                            zr(imax-1+iadr) = zr(ival-1+i)
                            zk8(inmax-1+iadr) = nomel(iel)
                            zi(ivmax-1+iadr) = 1
                        else if (zr(ival-1+i) .eq. zr(imax-1+iadr)) then
                            zi(ivmax-1+iadr) = zi(ivmax-1+iadr)+1
                        end if
                    end do
                end if
!
! -- RECHERCHE DE LA VALEURE MINIMALE ---
!
                if (lmin) then
                    do i = 1, nbcpt
                        if (lsup .or. linf) then
                            iadr = (zi(ipo2-1+i)-1)*icoef2+zi(icoe-1+i)
                        else
                            iadr = (zi(iposg-1+i)-1)*icoef2+zi(icoe-1+i)
                        end if
                        if (zr(imin-1+iadr) .eq. rundf) then
                            zr(imin-1+iadr) = zr(ival-1+i)
                            zk8(inmin-1+iadr) = nomel(iel)
                            zi(ivmin-1+iadr) = 1
                        else if (zr(ival-1+i) .lt. zr(imin-1+iadr)) then
                            zr(imin-1+iadr) = zr(ival-1+i)
                            zk8(inmin-1+iadr) = nomel(iel)
                            zi(ivmin-1+iadr) = 1
                        else if (zr(ival-1+i) .eq. zr(imin-1+iadr)) then
                            zi(ivmin-1+iadr) = zi(ivmin-1+iadr)+1
                        end if
                    end do
                end if
16              continue
            end do
!CCCCC
        else if (loc .eq. 'ELNO') then
            ipoin = point(iel)
            nbno = point(iel+1)-ipoin
            ncou = npcalc/nbno
            do icou = 1, ncou
                if (ncou .gt. 1) then
                    if (.not. lmax .and. .not. lmin) then
                        if (ncou .eq. 2) then
                            if (icou .eq. 1) write (ifi, '(A)') ' PEAU INTERNE'
                            if (icou .eq. 2) write (ifi, '(A)') ' PEAU EXTERNE'
                        else
                            write (ifi, '(A,I3)') ' COUCHE NUMERO:', &
                                icou
                        end if
                    end if
                end if
                do in = 1, nbno
                    nuno = connex(ipoin-1+in)
                    if (nbnot .ne. 0) then
                        do iino = 1, nbnot
                            if (nuno .eq. numnoe(iino)) goto 189
                        end do
                        goto 18
189                     continue
                    end if
                    nomno = nomnos(nuno)
                    j = iachml-1+ncmpp*icoef*(in-1)+(icou-1)*ncmpp* &
                        icoef*nbno
                    if (nbcmpt .eq. 0) then
                        do i = 1, ncmp2
                            if (ncmp .gt. 0) then
                                do jco = 1, icoef2
                                    zr(ival-1+(i-1)*icoef2+jco) = &
                                        vale(j+i+(nucmp(jco)-1)*ncmpp)
                                    zi(icoe-1+(i-1)*icoef2+jco) = jco
                                end do
                            else
                                do jco = 1, icoef2
                                    zr(ival-1+(i-1)*icoef2+jco) = &
                                        vale(j+i+(jco-1)*ncmpp)
                                    zi(icoe-1+(i-1)*icoef2+jco) = jco
                                end do
                            end if
                        end do
                    else
                        do i = 1, ncmp2
                            inu = zi(iposv-1+i)
                            if (ncmp .gt. 0) then
                                do jco = 1, icoef2
                                    zr(ival-1+(i-1)*icoef2+jco) = &
                                        vale(j+inu+(nucmp(jco)-1)*ncmpp)
                                    zi(icoe-1+(i-1)*icoef2+jco) = jco
                                end do
                            else
                                do jco = 1, icoef2
                                    zr(ival-1+(i-1)*icoef2+jco) = &
                                        vale(j+inu+(jco-1)*ncmpp)
                                    zi(icoe-1+(i-1)*icoef2+jco) = jco
                                end do
                            end if
                        end do
                    end if
!
! --  TRI DES COMPOSANTES DANS L'INTERVALLE BORINF,BORSUP
!
                    if (lsup .or. linf) then
                        do iva = 1, icoef2*ncmp2
                            if (lsup) then
                                if ((zr(ival-1+iva)-borsup) .gt. 0.d0) zi(icoe-1+iva) = 0
                            end if
                            if (linf) then
                                if ((zr(ival-1+iva)-borinf) .lt. 0.d0) zi(icoe-1+iva) = 0
                            end if
                        end do
!
! --- RETASSAGE POUR IMPRIMER COMPOSANTES PRESENTES DANS L'INTERVALLE --
!
                        icomp2 = 0
                        do i = 1, icoef2*ncmp2
                            if (zi(icoe-1+i) .ne. 0) then
                                icomp2 = icomp2+1
                                zi(icoe-1+icomp2) = zi(icoe-1+i)
                                zi(ipo2-1+icomp2) = zi(iposg-1+i)
                                zr(ival-1+icomp2) = zr(ival-1+i)
                                zk16(inop-1+icomp2) = zk16(inom-1+i)
                            end if
                        end do
                        if (icomp2 .eq. 0) goto 18
!
! -- IMPRESSION  --
!
                        if (.not. lmax .and. .not. lmin) then
                            ilig = (icomp2+ndim)/6
                            ires = (icomp2+ndim)-ilig*6
                            fmt1 = ' '
                            fmt2 = ' '
                            if (ires .ne. 0) then
                                fmt1 = '(9X,6(1X,'//forcmp//'),30(/,9X,6(1X,'//forcmp//')))'
                                fmt2 = '(1X,A8,6(1X,'//format//'),30(/,9X,6(1X,'//format//')))'
                            else if (ires .eq. 0 .and. ilig .eq. 1) then
                                fmt1 = '(9X,6(1X,'//forcmp//'))'
                                fmt2 = '(1X,A8,6(1X,'//format//'))'
                            else
                                write (fmt1, '(A,A8,A,I2,A,A8,A)') '(9X,6(1X,',&
     &                   forcmp, '),', (ilig-1), '(/,9X,6(1X,', forcmp, ')))'
                                write (fmt2, '(A,A10,A,I2,A,A10,A)') &
                                    '(1X,A8,6(1X,', format, '),', (ilig-1) &
                                    , '(/,9X,6(1X,', format, ')))'
                            end if
                            if (lsup .or. linf) then
                                if (limpr) then
                                    write (ifi, '(A,I2,A)') nomel(iel)
                                    limpr = .false.
                                end if
                            end if
                            if (ndim .eq. 0) then
                                write (ifi, fmt1) (zk16(inop-1+i) (1:11), &
                                                   i=1, icomp2)
                                write (ifi, fmt2) nomno, (zr(ival-1+icmp) &
                                                          , icmp=1, icomp2)
                            else
                                write (ifi, fmt1) (nomcor(i), i=1, ndim), &
                                    (zk16(inop-1+i) (1:11), i=1, icomp2)
                                write (ifi, fmt2) nomno, (coor((nuno-1)* &
                                                         3+i), i=1, ndim), (zr(ival-1+icmp), icmp= &
                                                                                  1, icomp2)
                            end if
                        end if
                        nbcpt = icomp2
                    else
                        if (.not. lmax .and. .not. lmin) then
                            if (ndim .eq. 0) then
                                write (ifi, fmv) nomno, (zr(ival-1+icmp), &
                                                         icmp=1, icoef2*ncmp2)
                            else
                                write (ifi, fmv) nomno, (coor((nuno-1)*3+ &
                                                              i), i=1, ndim), (zr(ival-1+icmp), &
                                                                               icmp=1, icoef2*ncmp2)
                            end if
                        end if
                        nbcpt = icoef2*ncmp2
                    end if
!
! -- RECHERCHE DE LA VALEUR MAXIMALE ---
!
                    if (lmax) then
                        do i = 1, nbcpt
                            if (lsup .or. linf) then
                                iadr = (zi(ipo2-1+i)-1)*icoef2+zi( &
                                       icoe-1+i)
                            else
                                iadr = (zi(iposg-1+i)-1)*icoef2+zi( &
                                       icoe-1+i)
                            end if
                            if (zr(imax-1+iadr) .eq. rundf) then
                                zr(imax-1+iadr) = zr(ival-1+i)
                                zk8(inmax-1+iadr) = nomel(iel)
                                zi(ivmax-1+iadr) = 1
                            elseif (zr(ival-1+i) .gt. zr(imax-1+iadr)) &
                                then
                                zr(imax-1+iadr) = zr(ival-1+i)
                                zk8(inmax-1+iadr) = nomel(iel)
                                zi(ivmax-1+iadr) = 1
                            elseif (zr(ival-1+i) .eq. zr(imax-1+iadr)) &
                                then
                                zi(ivmax-1+iadr) = zi(ivmax-1+iadr)+1
                            end if
                        end do
                    end if
!
! -- RECHERCHE DE LA VALEURE MINIMALE ---
!
                    if (lmin) then
                        do i = 1, nbcpt
                            if (lsup .or. linf) then
                                iadr = (zi(ipo2-1+i)-1)*icoef2+zi( &
                                       icoe-1+i)
                            else
                                iadr = (zi(iposg-1+i)-1)*icoef2+zi( &
                                       icoe-1+i)
                            end if
                            if (zr(imin-1+iadr) .eq. rundf) then
                                zr(imin-1+iadr) = zr(ival-1+i)
                                zk8(inmin-1+iadr) = nomel(iel)
                                zi(ivmin-1+iadr) = 1
                            elseif (zr(ival-1+i) .lt. zr(imin-1+iadr)) &
                                then
                                zr(imin-1+iadr) = zr(ival-1+i)
                                zk8(inmin-1+iadr) = nomel(iel)
                                zi(ivmin-1+iadr) = 1
                            elseif (zr(ival-1+i) .eq. zr(imin-1+iadr)) &
                                then
                                zi(ivmin-1+iadr) = zi(ivmin-1+iadr)+1
                            end if
                        end do
                    end if
18                  continue
                end do
            end do
        end if
12      continue
    end do
    write (ifi, *) ' '
!
! --- IMPRESSION DE LA VALEUR MAXIMALE ---
!
    if (lmax) then
        do i = 1, ncmpmx*icoef2
            if (zr(imax-1+i) .ne. rundf) then
                form1 = '(1X,3A,1X,'//format//',A,I4,2A)'
                write (ifi, form1) 'LA VALEUR MAXIMALE DE ', zk16(inot-1+i),&
     &       ' EST ', zr(imax-1+i),&
     &       ' EN ', zi(ivmax-1+i), ' MAILLE(S) : ', zk8(inmax-1+i)
            end if
        end do
    end if
!
! --- IMPRESSION DE LA VALEUR MINIMALE ---
!
    if (lmin) then
        do i = 1, ncmpmx*icoef2
            if (zr(imin-1+i) .ne. rundf) then
                form1 = '(1X,3A,1X,'//format//',A,I4,2A)'
                write (ifi, form1) 'LA VALEUR MINIMALE DE ', zk16(inot-1+i),&
     &       ' EST ', zr(imin-1+i),&
     &       ' EN ', zi(ivmin-1+i), ' MAILLE(S) : ', zk8(inmin-1+i)
            end if
        end do
    end if
!
    call jedetr('&&IRCERL.NCMT')
    call jedetr('&&IRCERL.MAX')
    call jedetr('&&IRCERL.MAIMAX')
    call jedetr('&&IRCERL.NBVMAX')
    call jedetr('&&IRCERL.MIN')
    call jedetr('&&IRCERL.MAIMIN')
    call jedetr('&&IRCERL.NBVMIN')
    call jedetr('&&IRCERL.ENT_COD')
    call jedetr('&&IRCERL.POSG')
    call jedetr('&&IRCERL.POSV')
    call jedetr('&&IRCERL.COEF')
    call jedetr('&&IRCERL.NCMP')
    call jedetr('&&IRCERL.NCPP')
    call jedetr('&&IRCERL.PO2')
    call jedetr('&&IRCERL.VAL')
!
999 continue
    call jedema()
end subroutine
