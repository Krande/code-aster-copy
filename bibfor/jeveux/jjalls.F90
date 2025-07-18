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
subroutine jjalls(lonoi, ic, genri, typei, lty, &
                  ci, itab, jitab, iadmi, iadyn)
! person_in_charge: j-pierre.lefebvre at edf.fr
! aslint: disable=C1002
    implicit none
#include "asterf_types.h"
#include "jeveux_private.h"
#include "asterc/hpalloc.h"
#include "asterfort/assert.h"
#include "asterfort/jjldyn.h"
#include "asterfort/jxlocs.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: lonoi, lty, itab(*), jitab, iadmi, iadyn
    character(len=*) :: genri, typei, ci
! ----------------------------------------------------------------------
! ALLOUE UN SEGMENT DE VALEUR EN MEMOIRE
!
! IN  LONOI  : LONGUEUR EN OCTETS DU SEGMENT DE VALEUR
! IN  IC     : CLASSE DE L'OBJET
! IN  GENRI  : GENRE DE L'OBJET JEVEUX
! IN  TYPEI  : TYPE DE L'OBJET JEVEUX
! IN  LTY    : LONGUEUR DU TYPE DE L'OBJET JEVEUX
! IN  CI     : = 'INIT' POUR INITIALISER LE SEGMENT DE VALEUR
! IN  ITAB   : TABLEAU PAR RAPPORT AUQUEL ON DETERMINE JITAB
! OUT JITAB  : ADRESSE DANS ITAB DU SEGMENT DE VALEUR
! OUT IADMI  : ADRESSE DU PREMIER MOT DU SEGMENT DE VALEUR
! OUT IADYN  : ADRESSE DU TABLEAU ALLOUE DYNAMIQUEMENT
! ----------------------------------------------------------------------
    integer(kind=8) :: lk1zon, jk1zon, liszon, jiszon
    common/izonje/lk1zon, jk1zon, liszon, jiszon
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iada, idm, ierr, iessai, ildyna, jcara
    integer(kind=8) :: jdate, jhcod, jiadd, jiadm, jiszo2, jlong, jlono
    integer(kind=8) :: jltyp, jluti, jmarq, lgbl, lsi, lso
    integer(kind=8) :: ltot, n, nde
!-----------------------------------------------------------------------
    parameter(n=5)
    common/jiatje/jltyp(n), jlong(n), jdate(n), jiadd(n), jiadm(n),&
     &                 jlono(n), jhcod(n), jcara(n), jluti(n), jmarq(n)
    integer(kind=8) :: nblmax, nbluti, longbl, kitlec, kitecr, kiadm, iitlec, iitecr
    integer(kind=8) :: nitecr, kmarq
    common/ificje/nblmax(n), nbluti(n), longbl(n),&
     &                 kitlec(n), kitecr(n), kiadm(n),&
     &                 iitlec(n), iitecr(n), nitecr(n), kmarq(n)
! ----------------------------------------------------------------------
    integer(kind=8) :: istat
    common/istaje/istat(4)
    integer(kind=8) :: lbis, lois, lols, lor8, loc8
    common/ienvje/lbis, lois, lols, lor8, loc8
    integer(kind=8) :: ldyn, lgdyn, nbdyn, nbfree
    common/idynje/ldyn, lgdyn, nbdyn, nbfree
    real(kind=8) :: mxdyn, mcdyn, mldyn, vmxdyn, vmet, lgio, cuvtrav
    common/r8dyje/mxdyn, mcdyn, mldyn, vmxdyn, vmet, lgio(2), cuvtrav
! ----------------------------------------------------------------------
    integer(kind=8) :: init, iblanc, valloc, lsic
    integer(kind=8) :: ic, ival(4), unmega
    aster_logical :: linit, ldeps
    character(len=8) :: cblanc
    equivalence(cblanc, iblanc)
    parameter(nde=6)
! ----------------------------------------------------------------------
! REMARQUE : LE PARAMETER NDE EST AUSSI DEFINI DANS JXLIRO JXECRO
! ----------------------------------------------------------------------
    data cblanc/'        '/
! DEB ------------------------------------------------------------------
    ltot = 0
    jitab = 0
    iadmi = 0
    iadyn = 0
    linit = (ci(1:4) .eq. 'INIT')
    lso = lonoi
!
!     LA TAILLE DU SEGMENT DE VALEURS EST AJUSTEE POUR S'ALLIGNER
!     SUIVANT LA LONGUEUR DU TYPE (SI SUPERIEUR A L'ENTIER)
!
!
    if (lty .ne. lois) then
        lso = lso+lty
        if (mod(lso, lois) .ne. 0) lso = (1+lso/lois)*lois
    end if
!
!     LA TAILLE DU SEGMENT DE VALEURS EST AJUSTEE A LA LONGUEUR DE BLOC
!     SI ON EST COMPRIS ENTRE LGBL-(NDE*LOIS) ET LGBL POUR DISPOSER DE
!     LA PLACE MINIMUM NECESSAIRE POUR LES GROS OBJETS
!
    if (ic .ne. 0) then
        if (longbl(ic) .gt. 1) then
            lgbl = 1024*longbl(ic)*lois
            if (lso .ge. lgbl-nde*lois .and. lso .lt. lgbl) then
                lso = lgbl
            end if
        end if
    end if
    ASSERT(lso .ne. 0)
    lsi = lso/lois
!
!     LE SEGMENT DE VALEURS EST ALLOUE DYNAMIQUEMENT
!
    iessai = 0
    ildyna = 0
!
    lsic = lsi+8
50  continue
    ildyna = ildyna+1
!
!   ON TESTE SI LE CUMUL DES ALLOCATIONS RESTE INFERIEUR A LA LIMITE
!   UNE EXCEPTION EST FAITE LORSQUE L'APPELANT EST jjagod  MAIS DANS
!   CE CAS PRECIS ON EVITE DE FAIRE APPEL A jjldyn
!
    if (mcdyn+lsic .gt. vmxdyn) then
        if (ildyna .gt. 1) then
            unmega = 1048576
            ival(1) = (lsic*lois)/unmega
            ival(2) = nint(vmxdyn*lois)/unmega
            ival(3) = nint(mcdyn*lois)/unmega
            ival(4) = (ltot*lois)/unmega
            call utmess('F', 'JEVEUX_62', ni=4, vali=ival)
        else
!
!  ON APPELLE JJLDYN AVEC L'ARGUMENT -2 POUR EVITER DE REACTUALISER
!  LA LIMITE MEMOIRE VMXDYN
!
            call jjldyn(2, -2, ltot)
            if (mcdyn+lsic .gt. vmxdyn) then
                call jjldyn(0, -1, ltot)
            end if
            goto 50
        end if
    end if
    iessai = iessai+1
    call hpalloc(iada, lsic, ierr, 0)
    if (ierr .eq. 0) then
        valloc = loc(iszon)
        jiszo2 = (iada-valloc)/lois
        iadmi = jiszo2+5-jiszon
        idm = jiszo2+1
        iadyn = iada
        mcdyn = mcdyn+lsic
        mxdyn = max(mxdyn, mcdyn*lois)
        nbdyn = nbdyn+1
    else
        if (iessai .gt. 1) then
            ival(1) = lsic*lois
            ival(2) = ltot*lois
            call utmess('F', 'JEVEUX_60', ni=2, vali=ival)
        else
            call jjldyn(2, -2, ltot)
            if (mcdyn+lsic .gt. vmxdyn) then
                call jjldyn(0, -1, ltot)
            end if
            goto 50
        end if
    end if
!
    iszon(idm) = idm+lsi+8-jiszon
    iszon(idm+1) = 0
    iszon(idm+2) = 0
    iszon(idm+3) = istat(1)
    iszon(idm+lsi+4) = istat(1)
    iszon(idm+lsi+5) = 0
    iszon(idm+lsi+6) = 0
    iszon(idm+lsi+7) = 0
!
    ldeps = .true.
    call jxlocs(itab, genri, lty, lonoi, iadmi, &
                ldeps, jitab)
!
    if (linit) then
        init = 0
        if (typei(1:1) .eq. 'K') init = iblanc
        do i = 1, lsi
            iszon(jiszon+iadmi+i-1) = init
        end do
    end if
! FIN ------------------------------------------------------------------
end subroutine
