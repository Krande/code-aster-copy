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
subroutine jxcopy(clsinz, nominz, clsouz, nmoutz, nbext)
! person_in_charge: j-pierre.lefebvre at edf.fr
! aslint: disable=W1303
! for the path name
    implicit none
#include "jeveux_private.h"
#include "asterc/cpfile.h"
#include "asterc/readdr.h"
#include "asterc/rmfile.h"
#include "asterc/writdr.h"
#include "asterfort/get_jvbasename.h"
#include "asterfort/jeinif.h"
#include "asterfort/jjalls.h"
#include "asterfort/jjlidy.h"
#include "asterfort/jxferm.h"
#include "asterfort/jxouvr.h"
#include "asterfort/utmess.h"
    character(len=*) :: clsinz, nominz, clsouz, nmoutz
    character(len=1) :: clasin, clasou
    character(len=8) :: nomin, nomout
!     ------------------------------------------------------------------
!     RECOPIE D'UNE BASE DE DONNEES APRES ELIMINATION DES
!     ENREGISTREMENTS DEVENUS INACCESSIBLES
!     ------------------------------------------------------------------
! IN  CLSINZ : NOM DE CLASSE ASSOCIEE EN ENTREE
! IN  NOMINZ : NOM DE LA BASE EN ENTREE
! IN  CLSOUZ : NOM DE CLASSE ASSOCIEE EN SORTIE
! IN  NMOUTZ : NOM DE LA BASE EN SORTIE
! OUT NBEXT  : NOMBRE D'"EXTENDS" UTILISES APRES RETASSAGE
!     ------------------------------------------------------------------
    integer(kind=8) :: lk1zon, jk1zon, liszon, jiszon
    common/izonje/lk1zon, jk1zon, liszon, jiszon
!     ------------------------------------------------------------------
    integer(kind=8) :: lbis, lois, lols, lor8, loc8
    common/ienvje/lbis, lois, lols, lor8, loc8
    integer(kind=8) :: istat
    common/istaje/istat(4)
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: iadloc, iadyn, ici, ico, ierr, k
    integer(kind=8) :: lbloc, n, nbext, nbloc, nrep, numext
!-----------------------------------------------------------------------
    parameter(n=5)
!
    integer(kind=8) :: nblmax, nbluti, longbl, kitlec, kitecr, kiadm, iitlec, iitecr
    integer(kind=8) :: nitecr, kmarq
    common/ificje/nblmax(n), nbluti(n), longbl(n),&
     &                 kitlec(n), kitecr(n), kiadm(n),&
     &                 iitlec(n), iitecr(n), nitecr(n), kmarq(n)
    integer(kind=8) :: idn, iext, nbenrg
    common/iextje/idn(n), iext(n), nbenrg(n)
    character(len=2) :: dn2
    character(len=5) :: classe
    character(len=8) :: nomfic, kstout, kstini
    common/kficje/classe, nomfic(n), kstout(n), kstini(n),&
     &                 dn2(n)
    integer(kind=8) :: nrhcod, nremax, nreuti
    common/icodje/nrhcod(n), nremax(n), nreuti(n)
    real(kind=8) :: svuse, smxuse
    common/statje/svuse, smxuse
    character(len=128) :: repglo, repvol
    common/banvje/repglo, repvol
    integer(kind=8) :: lrepgl, lrepvo
    common/balvje/lrepgl, lrepvo
!     ------------------------------------------------------------------
    character(len=1) :: kclas
    character(len=8) :: nomba1, nomba2
    character(len=512) :: noml1, noml2
    integer(kind=8) :: itp(1), jitp, iaditp, lgbl1, lgbl2, info, iret
! DEB ------------------------------------------------------------------
    nomin = nominz
    clasin = clsinz
    nomout = nmoutz
    clasou = clsouz
!
    kclas = clasin
    call jeinif('POURSUITE', 'DETRUIT', nomin, kclas, 1, &
                1, 1)
    ici = index(classe, kclas)
    kclas = clasou
    nrep = nremax(ici)
    nbloc = nbenrg(ici)
    lbloc = longbl(ici)
    call get_jvbasename(nomout(1:4), -1, noml1)
    info = 1
    call rmfile(noml1, info, iret)
    call jeinif('DEBUT', 'SAUVE', nomout, kclas, nrep, &
                nbloc, lbloc)
    ico = index(classe, kclas)
    nomba1 = nomfic(ici) (1:4)//'.   '
    nomba2 = nomfic(ico) (1:4)//'.   '
!
    lgbl1 = 1024*longbl(ici)*lois
    lgbl2 = 1024*longbl(ico)*lois
    call jjalls(lgbl1, 0, ' ', 'I', lois, &
                'INIT', itp, jitp, iaditp, iadyn)
    iszon(jiszon+iaditp-1) = istat(1)
    iszon(jiszon+iszon(jiszon+iaditp-4)-4) = istat(4)
    svuse = svuse+(iszon(jiszon+iaditp-4)-iaditp+4)
    smxuse = max(smxuse, svuse)
    do k = 1, (nbluti(ici)-1)/nbenrg(ici)
        call jxouvr(ico, k+1, mode=2)
        iext(ico) = iext(ico)+1
    end do
!
    do k = 1, nbluti(ici)
        numext = (k-1)/nbenrg(ici)
        iadloc = k-(numext*nbenrg(ici))
        call get_jvbasename(nomba1, numext+1, noml1)
        call readdr(noml1, iszon(jiszon+iaditp), lgbl1, iadloc, ierr)
        if (ierr .ne. 0) then
            call utmess('F', 'JEVEUX_47')
        end if
        call get_jvbasename(nomba2, numext+1, noml2)
        call writdr(noml2, iszon(jiszon+iaditp), lgbl2, iadloc, ierr)
        if (ierr .ne. 0) then
            call utmess('F', 'JEVEUX_48')
        end if
    end do
    nbext = numext+1
    call jxferm(ici)
    call jxferm(ico)
    call jjlidy(iadyn, iaditp)
    classe(ico:ico) = ' '
    classe(ici:ici) = ' '
!
! Destruction de tous les receptacles avant recopie.
!
    call get_jvbasename(nomba1, -2, noml1)
    call rmfile(noml1, info, iret)
!
    do k = 1, nbext
        call get_jvbasename(nomba1, k, noml1)
        call get_jvbasename(nomba2, k, noml2)
        call cpfile('M', noml2, noml1)
    end do
! FIN ------------------------------------------------------------------
end subroutine
