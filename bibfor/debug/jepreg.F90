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
subroutine jepreg(cunit, clas, numerg, cmess, info)
! person_in_charge: j-pierre.lefebvre at edf.fr
    implicit none
#include "asterf_types.h"
#include "jeveux_private.h"
#include "asterfort/iunifi.h"
#include "asterfort/jjalls.h"
#include "asterfort/jjlidy.h"
#include "asterfort/jxdeps.h"
#include "asterfort/jxliro.h"
    character(len=*) :: cunit, clas, cmess
    integer(kind=8) :: numerg, info
! ----------------------------------------------------------------------
! ROUTINE UTILISATEUR D'IMPRESSION DU CONTENU D'UN ENREGISTREMENT
! DU FICHIER D'ACCES DIRECT ASSOCIE A UNE BASE
!
! IN  CUNIT  : NOM LOCAL DU FICHIER DES IMPRESSIONS
! IN  CLAS   : CLASSE ASSOCIEE A LA BASE ( ' ' : TOUTES LES CLASSES )
! IN  NUMERG : NUMERO DE L'ENREGISTREMENT
! IN  CMESS  : MESSAGE D'INFORMATION
! IN  INFO   : NIVEAU DES IMPRESSIONS
!              SI > 1 ALORS ON IMPRIME LE CONTENU DE L'ENREGISTREMENT
!              SINON ON IMPRIME UNIQUEMENT LE CHAINAGE
! ----------------------------------------------------------------------
    integer(kind=8) :: lbis, lois, lols, lor8, loc8
    common/ienvje/lbis, lois, lols, lor8, loc8
    integer(kind=8) :: lk1zon, jk1zon, liszon, jiszon
    common/izonje/lk1zon, jk1zon, liszon, jiszon
!-----------------------------------------------------------------------
    integer(kind=8) :: iadmo, jdocu, jgenr, jorig, jrnom
    integer(kind=8) :: jtype, julist, jusadi, k, l, lgbl, n
    integer(kind=8) :: ncla1, ncla2
!-----------------------------------------------------------------------
    parameter(n=5)
!
    common/jkatje/jgenr(n), jtype(n), jdocu(n), jorig(n), jrnom(n)
!
    integer(kind=8) :: nblmax, nbluti, longbl, kitlec, kitecr, kiadm, iitlec, iitecr
    integer(kind=8) :: nitecr, kmarq
    common/ificje/nblmax(n), nbluti(n), longbl(n),&
     &                 kitlec(n), kitecr(n), kiadm(n),&
     &                 iitlec(n), iitecr(n), nitecr(n), kmarq(n)
    common/jusadi/jusadi(n)
    character(len=2) :: dn2
    character(len=5) :: classe
    character(len=8) :: nomfic, kstout, kstini
    common/kficje/classe, nomfic(n), kstout(n), kstini(n),&
     &                 dn2(n)
    character(len=8) :: nombas
    common/kbasje/nombas(n)
    integer(kind=8) :: istat
    common/istaje/istat(4)
    real(kind=8) :: svuse, smxuse
    common/statje/svuse, smxuse
! ----------------------------------------------------------------------
    character(len=1) :: kclas
    integer(kind=8) :: itp(1), jitp, iaditp, iadyn, idco, idos, idec, icomp
    integer(kind=8) :: ji, nl, nd, iaddi(2), idosl, idcol, lgl, ic
! DEB ------------------------------------------------------------------
    kclas = clas(1:min(1, len(clas)))
    julist = iunifi(cunit)
    if (julist .eq. 0) goto 999
    if (kclas .eq. ' ') then
        ncla1 = 1
        ncla2 = index(classe, '$')-1
        if (ncla2 .lt. 0) ncla2 = n
    else
        ncla1 = index(classe, kclas)
        ncla2 = ncla1
    end if
    do ic = ncla1, ncla2
        lgbl = 1024*longbl(ic)*lois
        write (julist, '(    1X,4A)') ('--------------------', k=1, 4)
        write (julist, *) ' '
        write (julist, '(1X,A)') cmess
        write (julist, *) ' NOM DE LA BASE                    : ',&
     &        nombas(ic)
        write (julist, *) ' NB D''ENREGISTREMENTS MAXIMUM      : ',&
     &        nblmax(ic)
        write (julist, *) ' NB D''ENREGISTREMENTS UTILISES     : ',&
     &        nbluti(ic)
        write (julist, *) ' LONGUEUR D''ENREGISTREMENT (OCTETS): ', &
            lgbl
        write (julist, *) '                                  '
        write (julist, '(    1X,4A)') ('--------------------', k=1, 4)
!
        if (numerg .le. nbluti(ic)) then
            lgbl = 1024*longbl(ic)*lois
            call jjalls(lgbl, 0, 'V', 'I', lois, &
                        'INIT', itp, jitp, iaditp, iadyn)
            iszon(jiszon+iaditp-1) = istat(2)
            iszon(jiszon+iszon(jiszon+iaditp-4)-4) = istat(4)
            svuse = svuse+(iszon(jiszon+iaditp-4)-iaditp+4)
            if (iadyn .ne. 0) svuse = svuse+1
            iaddi(1) = numerg
            iaddi(2) = 0
            iadmo = (iaditp-1)*lois+iszon(jiszon+iaditp-3)+1
            if (numerg .eq. iitlec(ic)) then
                call jxdeps(kitlec(ic)+1, iadmo, lgbl)
            else if (numerg .eq. iitecr(ic)) then
                call jxdeps(kitecr(ic)+1, iadmo, lgbl)
            else
                call jxliro(ic, iaditp, iaddi, lgbl)
            end if
            if (info .gt. 1) then
                ji = jiszon+iaditp
                nl = lgbl/(10*lois)
                nd = mod(lgbl, (10*lois))/lois
                write (julist, *) 'CONTENU BRUT DE L''ENREGISTREMENT ', &
                    numerg
                write (julist, 1001) (10*(l-1)+1, (iszon(ji+10*(l- &
                                                                1)+k-1), k=1, 10), l=1, nl)
                if (nd .ne. 0) then
                    write (julist, 1001) 10*nl+1, (iszon(ji+10*nl+ &
                                                         k-1), k=1, nd)
                end if
            end if
            idco = iusadi(jusadi(ic)+3*numerg-2)
            idos = iusadi(jusadi(ic)+3*numerg-1)
            if (idos .gt. 0 .or. idco .gt. 0) then
!
! ------- L'ENREGISTREMENT CONTIENT UN OU UNE PARTIE D'UN "GROS" OBJET
!
                write (julist, *) ' '
                write (julist, *) 'ENREGISTREMENT NUMERO : ', numerg
                write (julist, *) 'ENREGISTREMENT AFFECTE A UN UNIQUE OBJET'
                if (idco .eq. 0) then
                    write (julist, *) 'OBJET SIMPLE DE NOM :  ',&
     &                        rnom(jrnom(ic)+idos), idos
                else
                    write (julist, *) 'OBJET DE COLLECTION DE NOM : ', &
                        rnom(jrnom(ic)+idco), idco, ' NUMERO ', idos
                end if
            else if (idco .lt. 0 .or. idos .lt. 0) then
!
! ------- L'ENREGISTREMENT CORRESPOND A UN OBJET DETRUIT
!
                write (julist, *) ' '
                write (julist, *) 'ENREGISTREMENT NUMERO : ', numerg
                write (julist, *) 'ENREGISTREMENT LIBERE', idco, idos
            else if (idco .eq. 0 .and. idos .eq. 0) then
!
! ------- L'ENREGISTREMENT CONTIENT DES PETITS OBJETS
!
                write (julist, *) ' '
                write (julist, *) 'CHAINAGE DE L''ENREGISTREMENT NUMERO :',&
     &                     numerg
                idec = 0
                icomp = 0
300             continue
                icomp = icomp+1
                idcol = iszon(jiszon+iaditp+idec)
                idosl = iszon(jiszon+iaditp+idec+1)
                lgl = iszon(jiszon+iaditp+idec+2)
                if (idcol .eq. 0 .and. idosl .eq. 0) then
                    write (julist, 1002) icomp, idcol, idosl, lgl,&
     &                        ' FIN DE L''ENREGISTREMENT ATTEINT'
                    goto 350
                else if (idcol .lt. 0 .or. idosl .lt. 0) then
                    write (julist, 1002) icomp, idcol, idosl, lgl,&
     &                        ' OBJET DETRUIT'
                    goto 320
                end if
                if (idcol .eq. 0) then
                    write (julist, 1002) icomp, idcol, idosl, lgl,&
     &                   ' OBJET SIMPLE DE NOM : ', rnom(jrnom(ic)+idosl)
                else
                    write (julist, 1003) icomp, idcol, idosl, lgl,&
     &                    ' OBJET DE COLLECTION DE NOM : ',&
     &                   rnom(jrnom(ic)+idcol), ' ET DE NUMERO : ', idosl
                end if
320             continue
                idec = idec+lgl+3
                goto 300
350             continue
            end if
        else
            write (julist, *) 'NUMERO D''ENREGISTREMENT INEXISTANT ', &
                numerg
            goto 999
        end if
        write (julist, '(    1X,4A)') ('--------------------', k=1, 4)
!
        if (iadyn .ne. 0) call jjlidy(iadyn, iaditp)
    end do
999 continue
!
1001 format((i7, ' - ', 10(i12, 1x)))
1002 format(i8, 1x, 'ID COLLECTION:', i8, ' ID OBJET:', i8, ' LONGUEUR:',&
   &       i8, a, a, i6)
1003 format(i8, 1x, 'ID COLLECTION:', i8, ' ID OBJET:', i8, ' LONGUEUR:',&
   &       i8, a, a, a, i8)
! FIN -----------------------------------------------------------------
end subroutine
