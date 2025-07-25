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
subroutine ircnrl(ifi, nbno, prno, nueq, nec, &
                  dg, ncmpmx, vale, nomcmp, nomnoe, &
                  lcor, ndim, coor, numnoe, nbcmpt, &
                  nucmpu, lsup, borsup, linf, borinf, &
                  lmax, lmin, formr)
! aslint: disable=W1504
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8vide.h"
#include "asterfort/exisdg.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/lxlgut.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: ifi, nbno, prno(*), nueq(*), nec, dg(*), ncmpmx
    integer(kind=8) :: ndim, numnoe(*), nbcmpt, nucmpu(*)
    real(kind=8) :: borsup, borinf, coor(*), vale(*)
    character(len=*) :: nomcmp(*), nomnoe(*), formr
    aster_logical :: lcor, lsup, linf, lmax, lmin
!
!        ECRITURE D'UN CHAM_NO SUR FICHIER IFI AU FORMAT 'RESULTAT'
!        A VALEURS REELLES
!      ENTREE:
!         IFI   : UNITE LOGIQUE DU FICHIER UNIVERSEL
!         NBNO  : NOMBRE DE NOEUDS A IMPRIMER
!         PRNO  : OBJET .PRNO(ILIGREL) D'UN NUME_EQUA
!         NUEQ  : OBJET .NUEQ D'UN NUME_EQUA
!         NEC   : NOMBRE D'ENTIERS-CODES
!         DG    : TABLEAU DES ENTIERS CODES
!         NCMPMX: NOMBRE MAXI DE CMP DE LA GRANDEUR NOMGD
!         VALE  : VALEURS DU CHAM_NO
!         NOMCMP: NOMS DES CMP
!         NOMNOE: NOMS DES NOEUDS
!         LCOR  : IMPRESSION DES COORDONNEES .TRUE. IMPRESSION
!         NDIM  : DIMENSION DU MAILLAGE
!         COOR  : COORDONNEES D'UN MAILLAGE
!         NUMNOE: NUMEROS DES NOEUDS A IMPRIMER
!         NBCMPT: NOMBRE DE COMPOSANTES A IMPRIMER
!         NOCMPU: NUMEROS DES COMPOSANTES A IMPRIMER
!         LSUP  : =.TRUE. INDIQUE PRESENCE D'UNE BORNE SUPERIEURE
!         BORSUP: VALEUR DE LA BORNE SUPERIEURE
!         LINF  : =.TRUE. INDIQUE PRESENCE D'UNE BORNE INFERIEURE
!         BORINF: VALEUR DE LA BORNE INFERIEURE
!         LMAX  : =.TRUE. INDIQUE IMPRESSION VALEUR MAXIMALE
!         LMIN  : =.TRUE. INDIQUE IMPRESSION VALEUR MINIMALE
!         FORMR : FORMAT D'ECRITURE DES REELS SUR "RESULTAT"
!     ------------------------------------------------------------------
!     ATTENTION EN CAS DE MODIFICATION DE CE SS-PGME, PENSER A IRCNC8
!     ------------------------------------------------------------------
!
    real(kind=8) :: rundf
    integer(kind=8) :: impre
    character(len=8) :: nomcor(3), forcmp
    character(len=10) :: format
    character(len=50) :: fmt, form1
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, icm, icmp, icmp2, icomp2, icompt, id
    integer(kind=8) :: iec, ieq, if, ilign, imax, imin, inec
    integer(kind=8) :: inmax, inmin, inno, ino, ipos, ipres, irest
    integer(kind=8) :: irval, iva, ival, ivmax, ivmin, lgr
    integer(kind=8) :: ncmp
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
! -- ALLOCATION DES TABLEAUX DE TRAVAIL ---
!
    call jedetr('&&IRCNRL.VAL')
    call wkvect('&&IRCNRL.VAL', 'V V R', ncmpmx, irval)
    call jedetr('&&IRCNRL.POS')
    call wkvect('&&IRCNRL.POS', 'V V I', ncmpmx, ipos)
    if (nec .gt. 0) then
        call jedetr('&&IRCNRL.ENT')
        call wkvect('&&IRCNRL.ENT', 'V V I', nec, inec)
        do iec = 1, nec
            zi(inec-1+iec) = 0
        end do
    end if
    if (lmax) then
        call jedetr('&&IRCNRL.MAX')
        call wkvect('&&IRCNRL.MAX', 'V V R', ncmpmx, imax)
        call jedetr('&&IRCNRL.NOEMAX')
        call wkvect('&&IRCNRL.NOEMAX', 'V V K8', ncmpmx, inmax)
        call jedetr('&&IRCNRL.NBVMAX')
        call wkvect('&&IRCNRL.NBVMAX', 'V V I', ncmpmx, ivmax)
        do i = 1, ncmpmx
            zr(imax-1+i) = rundf
        end do
    end if
    if (lmin) then
        call jedetr('&&IRCNRL.MIN')
        call wkvect('&&IRCNRL.MIN', 'V V R', ncmpmx, imin)
        call jedetr('&&IRCNRL.NOEMIN')
        call wkvect('&&IRCNRL.NOEMIN', 'V V K8', ncmpmx, inmin)
        call jedetr('&&IRCNRL.NBVMIN')
        call wkvect('&&IRCNRL.NBVMIN', 'V V I', ncmpmx, ivmin)
        do i = 1, ncmpmx
            zr(imin-1+i) = rundf
        end do
    end if
!
    do inno = 1, nbno
        ino = numnoe(inno)
        do iec = 1, nec
            dg(iec) = prno((ino-1)*(nec+2)+2+iec)
        end do
!
!        NCMP : NOMBRE DE CMPS SUR LE NOEUD INO
!        IVAL : ADRESSE DU DEBUT DU NOEUD INO DANS .NUEQ
!
        ival = prno((ino-1)*(nec+2)+1)
        ncmp = prno((ino-1)*(nec+2)+2)
        if (ncmp .eq. 0) goto 11
!
        do i = 1, ncmpmx
            zi(ipos-1+i) = 0
        end do
        icompt = 0
        impre = 0
        ipres = 0
        do icmp = 1, ncmpmx
            if (exisdg(dg, icmp)) then
                ipres = ipres+1
                ieq = nueq(ival-1+ipres)
                if (nbcmpt .ne. 0) then
                    do icm = 1, nbcmpt
                        icmp2 = nucmpu(icm)
                        if (icmp .eq. icmp2) then
                            zr(irval-1+icm) = vale(ieq)
                            zi(ipos-1+icm) = icmp
                            goto 12
                        end if
                    end do
                else
                    icompt = ipres
                    zr(irval-1+icompt) = vale(ieq)
                    zi(ipos-1+icompt) = icmp
                end if
            end if
12          continue
        end do
!
! --- RETASSAGE POUR IMPRIMER COMPOSANTES ORDRE UTILISATEUR---
!
        if (nbcmpt .ne. 0) then
            icompt = 0
            do i = 1, nbcmpt
                if (zi(ipos-1+i) .ne. 0) then
                    icompt = icompt+1
                    zi(ipos-1+icompt) = zi(ipos-1+i)
                    zr(irval-1+icompt) = zr(irval-1+i)
                end if
            end do
        end if
        do iec = 1, nec
            if (dg(iec) .ne. zi(inec-1+iec)) then
                impre = 1
                zi(inec-1+iec) = dg(iec)
            end if
        end do
!
! --  TRI DES COMPOSANTES DANS L'INTERVALLE BORINF,BORSUP
!
        if (lsup .or. linf) then
            do iva = 1, icompt
                if (lsup) then
                    if ((zr(irval-1+iva)-borsup) .gt. 0.d0) zi(ipos-1+iva) = 0
                end if
                if (linf) then
                    if ((zr(irval-1+iva)-borinf) .lt. 0.d0) zi(ipos-1+iva) = 0
                end if
            end do
!
! --- RETASSAGE POUR IMPRIMER COMPOSANTES PRESENTES DANS L'INTERVALLE --
!
            icomp2 = 0
            do i = 1, icompt
                if (zi(ipos-1+i) .ne. 0) then
                    icomp2 = icomp2+1
                    zi(ipos-1+icomp2) = zi(ipos-1+i)
                    zr(irval-1+icomp2) = zr(irval-1+i)
                end if
            end do
            icompt = icomp2
        end if
        if (icompt .eq. 0) then
            goto 11
        end if
!
! -- RECHERCHE DE LA VALEURE MAXIMALE ---
!
        if (lmax) then
            do i = 1, icompt
                if (zr(imax-1+zi(ipos-1+i)) .eq. rundf) then
                    zr(imax-1+zi(ipos-1+i)) = zr(irval-1+i)
                    zk8(inmax-1+zi(ipos-1+i)) = nomnoe(inno)
                    zi(ivmax-1+zi(ipos-1+i)) = 1
                else if (zr(irval-1+i) .gt. zr(imax-1+zi(ipos-1+i))) then
                    zr(imax-1+zi(ipos-1+i)) = zr(irval-1+i)
                    zk8(inmax-1+zi(ipos-1+i)) = nomnoe(inno)
                    zi(ivmax-1+zi(ipos-1+i)) = 1
                else if (zr(irval-1+i) .eq. zr(imax-1+zi(ipos-1+i))) then
                    zi(ivmax-1+zi(ipos-1+i)) = zi(ivmax-1+zi(ipos-1+i))+ &
                                               1
                end if
            end do
        end if
!
! -- RECHERCHE DE LA VALEURE MINIMALE ---
!
        if (lmin) then
            do i = 1, icompt
                if (zr(imin-1+zi(ipos-1+i)) .eq. rundf) then
                    zr(imin-1+zi(ipos-1+i)) = zr(irval-1+i)
                    zk8(inmin-1+zi(ipos-1+i)) = nomnoe(inno)
                    zi(ivmin-1+zi(ipos-1+i)) = 1
                else if (zr(irval-1+i) .lt. zr(imin-1+zi(ipos-1+i))) then
                    zr(imin-1+zi(ipos-1+i)) = zr(irval-1+i)
                    zk8(inmin-1+zi(ipos-1+i)) = nomnoe(inno)
                    zi(ivmin-1+zi(ipos-1+i)) = 1
                else if (zr(irval-1+i) .eq. zr(imin-1+zi(ipos-1+i))) then
                    zi(ivmin-1+zi(ipos-1+i)) = zi(ivmin-1+zi(ipos-1+i))+ &
                                               1
                end if
            end do
        end if
!
! - IMPRESSION DES VALEURS ---
!
        if (.not. lmax .and. .not. lmin .and. lcor) then
            ilign = (icompt+ndim)/6
            irest = (icompt+ndim)-ilign*6
            if (impre .eq. 1 .or. lsup .or. linf) then
                fmt = ' '
                if (irest .ne. 0) then
                    fmt = '( 1X,A,6(1X,'//forcmp//'),30(/,9X,6(1X,'//forcmp//')) )'
                else if (irest .eq. 0 .and. ilign .eq. 1) then
                    fmt = '(1X,A,6(1X,'//forcmp//'))'
                else
                    write (fmt, '(A,A8,A,I2,A,A8,A)') '(1X,A,6(1X,', forcmp,&
     &                     '),', (ilign-1), '(/,9X,6(1X,', forcmp, ')))'
                end if
                write (ifi, fmt) 'NOEUD   ', (nomcor(i), i=1, ndim),&
     &                        (nomcmp(zi(ipos-1+i)), i=1, icompt)
            end if
            fmt = ' '
            if (irest .ne. 0) then
                fmt = '(1X,A,6(1X,'//format//'),30(/,9X,6(1X,'//format//')))'
            else if (irest .eq. 0 .and. ilign .eq. 1) then
                fmt = '(1X,A,6(1X,'//format//'))'
            else
                write (fmt, '(A,A10,A,I2,A,A10,A)') '(1X,A,6(1X,', &
                    format, '),', (ilign-1), '(/,9X,6(1X,', format, ')))'
            end if
            write (ifi, fmt) nomnoe(inno), (coor((ino-1)*3+i), i=1, ndim) &
                , (zr(irval-1+i), i=1, icompt)
        else if (.not. lmax .and. .not. lmin) then
            ilign = (icompt)/6
            irest = (icompt)-ilign*6
            if (impre .eq. 1 .or. lsup .or. linf) then
                fmt = ' '
                if (irest .ne. 0) then
                    fmt = '( 1X,A,6(1X,'//forcmp//'),30(/,9X,6(1X,'//forcmp//')) )'
                else if (irest .eq. 0 .and. ilign .eq. 1) then
                    fmt = '(1X,A,6(1X,'//forcmp//'))'
                else
                    write (fmt, '(A,A8,A,I2,A,A8,A)') '(1X,A,6(1X,', forcmp,&
     &                    '),', (ilign-1), '(/,9X,6(1X,', forcmp, ')))'
                end if
                write (ifi, fmt) 'NOEUD   ',&
     &                        (nomcmp(zi(ipos-1+i)), i=1, icompt)
            end if
            fmt = ' '
            if (irest .ne. 0) then
                fmt = '(1X,A,6(1X,'//format//'),30(/,9X,6(1X,'//format//')))'
            else if (irest .eq. 0 .and. ilign .eq. 1) then
                fmt = '(1X,A,6(1X,'//format//'))'
            else
                write (fmt, '(A,A10,A,I2,A,A10,A)') '(1X,A,6(1X,', &
                    format, '),', (ilign-1), '(/,9X,6(1X,', format, ')))'
            end if
            write (ifi, fmt) nomnoe(inno), (zr(irval-1+i), i=1, icompt)
        end if
11      continue
    end do
    write (ifi, '(A)') ' '
!
! --- IMPRESSION DE LA VALEUR MAXIMALE ---
!
    if (lmax) then
        do i = 1, ncmpmx
            if (zr(imax-1+i) .ne. rundf) then
                form1 = '(1X,3A,1X,'//format//',A,I4,A,A8)'
                write (ifi, form1) 'LA VALEUR MAXIMALE DE ', nomcmp(i),&
     &       ' EST', zr(imax-1+i),&
     &       ' EN ', zi(ivmax-1+i), ' NOEUD(S) : ', zk8(inmax-1+i)
            end if
        end do
    end if
!
! --- IMPRESSION DE LA VALEUR MINIMALE ---
!
    if (lmin) then
        do i = 1, ncmpmx
            if (zr(imin-1+i) .ne. rundf) then
                form1 = '(1X,3A,1X,'//format//',A,I4,A,A8)'
                write (ifi, form1) 'LA VALEUR MINIMALE DE ', nomcmp(i),&
     &       ' EST', zr(imin-1+i),&
     &       ' EN ', zi(ivmin-1+i), ' NOEUD(S) : ', zk8(inmin-1+i)
            end if
        end do
    end if
!
    call jedetr('&&IRCNRL.VAL')
    call jedetr('&&IRCNRL.POS')
    call jedetr('&&IRCNRL.ENT')
    call jedetr('&&IRCNRL.MAX')
    call jedetr('&&IRCNRL.NOEMAX')
    call jedetr('&&IRCNRL.NBVMAX')
    call jedetr('&&IRCNRL.MIN')
    call jedetr('&&IRCNRL.NOEMIN')
    call jedetr('&&IRCNRL.NBVMIN')
    call jedema()
end subroutine
