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
subroutine te0443(option, nomte)
!
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/vdefro.h"
#include "asterfort/vdrep2.h"
#include "asterfort/vdrepe.h"
#include "asterfort/vdsiro.h"
#include "asterfort/vectan.h"
#include "asterfort/vectgt.h"
    character(len=16) :: option, nomte
!......................................................................
!
!    - FONCTION REALISEE: CHANGEMENT DE REPERE POUR LES COQUE_3D
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                        'REPE_TENS'  :  TENSEURS
!                        'REPE_GENE'  :  QUANTITES GENERALISEES
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
!
    integer(kind=8) :: ndim, nno, nnos, npg, ipoids, ivf, idfdx, jgano, iret(4)
    integer(kind=8) :: jgeom, jin, jout, jang, np, i, itab(7), iret1, iret2, nbsp
    integer(kind=8) :: vali(2)
!-----------------------------------------------------------------------
    integer(kind=8) :: intsn, j, jcara, k, lzi, lzr, nb1
    integer(kind=8) :: nb2, ncmp, ncpmax, npgsn, nptmax, nspmax
    real(kind=8) :: epais, zero
!-----------------------------------------------------------------------
    parameter(nptmax=9, ncpmax=8, nspmax=162)
    real(kind=8) :: alpha, beta, s
    real(kind=8) :: matvn1(2, 2, 10), matvg1(2, 2, 10)
    real(kind=8) :: matvn2(2, 2, 10), matvg2(2, 2, 10)
    real(kind=8) :: vecta(9, 2, 3), vectn(9, 3), vectpt(9, 2, 3)
    real(kind=8) :: vectg(2, 3), vectt(3, 3), conin(nptmax*ncpmax*nspmax)
    integer(kind=8)      :: repere_nature
    character(len=8) :: pain, paout
    character(len=24) :: messk(2)
! --------------------------------------------------------------------------------------------------
!
    zero = 0.0d0
!
    call jevech('PCACOQU', 'L', jcara)
    epais = zr(jcara)
!
    if (option .ne. 'REPE_TENS' .and. option .ne. 'REPE_GENE') then
!       OPTION DE CALCUL INVALIDE
        ASSERT(.false.)
    end if
!
    if (option .eq. 'REPE_TENS') then
        ncmp = 6
        call tecach('ONO', 'PCOGAIN', 'L', iret(1), nval=7, itab=itab)
        call tecach('ONO', 'PCONOIN', 'L', iret(2), nval=7, itab=itab)
        call tecach('ONO', 'PDEGAIN', 'L', iret(3), nval=7, itab=itab)
        call tecach('ONO', 'PDENOIN', 'L', iret(4), nval=7, itab=itab)
        iret1 = iret(1)+iret(2)+iret(3)+iret(4)
        ASSERT(iret1 .eq. 6)
!
        if (iret(1) .eq. 0) then
            pain = 'PCOGAIN'
            paout = 'PCOGAOUT'
        else if (iret(2) .eq. 0) then
            pain = 'PCONOIN'
            paout = 'PCONOOUT'
        else if (iret(3) .eq. 0) then
            pain = 'PDEGAIN'
            paout = 'PDEGAOUT'
        else if (iret(4) .eq. 0) then
            pain = 'PDENOIN'
            paout = 'PDENOOUT'
        end if
!
    else if (option .eq. 'REPE_GENE') then
        ncmp = 8
        call tecach('ONO', 'PEFGAIN', 'L', iret(1), nval=7, itab=itab)
        call tecach('ONO', 'PEFNOIN', 'L', iret(2), nval=7, itab=itab)
        call tecach('ONO', 'PDGGAIN', 'L', iret(3), nval=7, itab=itab)
        call tecach('ONO', 'PDGNOIN', 'L', iret(4), nval=7, itab=itab)
        iret1 = iret(1)+iret(2)+iret(3)+iret(4)
        ASSERT(iret1 .eq. 6)
!
        if (iret(1) .eq. 0) then
            pain = 'PEFGAIN'
            paout = 'PEFGAOUT'
        else if (iret(2) .eq. 0) then
            pain = 'PEFNOIN'
            paout = 'PEFNOOUT'
        else if (iret(3) .eq. 0) then
            pain = 'PDGGAIN'
            paout = 'PDGGAOUT'
        else if (iret(4) .eq. 0) then
            pain = 'PDGNOIN'
            paout = 'PDGNOOUT'
        end if
    end if
!
!  appel a elrefe_info pour recuperer nno et npg
    call elrefe_info(fami='MASS', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfdx, jgano=jgano)
!
    if (pain(4:5) .eq. 'NO') then
        np = nno
    else if (pain(4:5) .eq. 'GA') then
        np = npg
    end if
!
    call jevech('PGEOMER', 'L', jgeom)
    call jevech('PANGREP', 'L', jang)
    call jevech(pain, 'L', jin)
    call jevech(paout, 'E', jout)
    call tecach('OOO', pain, 'L', iret2, nval=7, itab=itab)
    nbsp = itab(7)
    if ((nbsp .ne. 1) .and. (mod(nbsp, 3) .ne. 0)) then
        call utmess('F', 'ELEMENTS5_54', si=nbsp)
    end if
!
    alpha = zr(jang)
    beta = zr(jang+1)
    repere_nature = nint(zr(jang+2))
    call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
    nb1 = zi(lzi-1+1)
    nb2 = zi(lzi-1+2)
    npgsn = zi(lzi-1+4)
    call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)
!
!   CALCUL DES MATRICES DE PASSAGE GLOBAL-INTRINSEQUE
    call vectan(nb1, nb2, zr(jgeom), zr(lzr), vecta, vectn, vectpt)
!
!   DETERMINATION DES REPERES  LOCAUX DE L'ELEMENT AUX POINTS
!   D'INTEGRATION ET STOCKAGE DE CES REPERES DANS LE VECTEUR .DESR
    k = 0
    do intsn = 1, npgsn
        call vectgt(1, nb1, zr(jgeom), zero, intsn, &
                    zr(lzr), epais, vectn, vectg, vectt)
        do j = 1, 3
            do i = 1, 3
                k = k+1
                zr(lzr+2000+k-1) = vectt(i, j)
            end do
        end do
    end do
!
    ASSERT(ncmp .le. ncpmax)
    ASSERT(np .le. nptmax)
    vali(1) = nspmax
    vali(2) = nbsp
    if (nbsp .gt. nspmax) then
        call utmess('F', 'ELEMENTS5_4', ni=2, vali=vali)
    end if
!
!   LE TABLEAU CONIN A ETE ALLOUE DE FACON STATIQUE POUR
!   OPTIMISER LE CPU CAR LES APPELS A WKVECT DANS LES TE SONT COUTEUX.
    call vdrepe(nomte, matvn1, matvg1)
!
!   ON PREND L INVERSE DES MATRICES (CAR ON REVIENT EN REPERE INTRINSEQUE)
    if ((repere_nature .eq. 0) .or. (repere_nature .eq. 2)) then
        if (pain(4:5) .eq. 'NO') then
            do i = 1, np
                s = matvn1(1, 2, i)
                matvn1(2, 1, i) = s
                matvn1(1, 2, i) = -s
            end do
        else if (pain(4:5) .eq. 'GA') then
            do i = 1, np
                s = matvg1(1, 2, i)
                matvg1(2, 1, i) = s
                matvg1(1, 2, i) = -s
            end do
        end if
    end if
!
    if (repere_nature .eq. 0) then
!       PASSAGE DES CONTRAINTES DU REPERE LOCAL 1
!       A L'ELEMENT AU REPERE INTRINSEQUE DE LA COQUE
        if (option .eq. 'REPE_TENS') then
            if (pain(4:5) .eq. 'NO') then
                call vdsiro(np, nbsp, matvn1, 'IU', 'N', zr(jin), conin)
            else if (pain(4:5) .eq. 'GA') then
                call vdsiro(np, nbsp, matvg1, 'IU', 'G', zr(jin), conin)
            end if
        else if (option .eq. 'REPE_GENE') then
            if (pain(4:5) .eq. 'NO') then
                call vdefro(np, matvn1, zr(jin), conin)
            else if (pain(4:5) .eq. 'GA') then
                call vdefro(np, matvg1, zr(jin), conin)
            end if
        end if
!       CALCUL DES MATRICES DE PASSAGE DU CHGT DE REPERE
        call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
        call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)
        call vdrep2(alpha, beta, zi(lzi), zr(lzr), matvn2, matvg2)
!
!       PASSAGE DES QUANTITES DU REPERE INTRINSEQUE
!       A L'ELEMENT AU REPERE LOCAL DE LA COQUE
        if (option .eq. 'REPE_TENS') then
            if (pain(4:5) .eq. 'NO') then
                call vdsiro(np, nbsp, matvn2, 'IU', 'N', conin, zr(jout))
            else if (pain(4:5) .eq. 'GA') then
                call vdsiro(np, nbsp, matvg2, 'IU', 'G', conin, zr(jout))
            end if
        else if (option .eq. 'REPE_GENE') then
            if (pain(4:5) .eq. 'NO') then
                call vdefro(np, matvn2, conin, zr(jout))
            else if (pain(4:5) .eq. 'GA') then
                call vdefro(np, matvg2, conin, zr(jout))
            end if
        end if
!
    else if (repere_nature .eq. 1) then
!       PASSAGE DES CONTRAINTES DU REPERE INTRINSEQUE
!       A L'ELEMENT AU REPERE LOCAL 1 DE LA COQUE REPERE = 'COQUE_INTR_UTIL'
        if (option .eq. 'REPE_TENS') then
            if (pain(4:5) .eq. 'NO') then
                call vdsiro(np, nbsp, matvn1, 'IU', 'N', zr(jin), zr(jout))
            else if (pain(4:5) .eq. 'GA') then
                call vdsiro(np, nbsp, matvg1, 'IU', 'G', zr(jin), zr(jout))
            end if
        else if (option .eq. 'REPE_GENE') then
            if (pain(4:5) .eq. 'NO') then
                call vdefro(np, matvn1, zr(jin), zr(jout))
            else if (pain(4:5) .eq. 'GA') then
                call vdefro(np, matvg1, zr(jin), zr(jout))
            end if
        end if
!
    else if (repere_nature .eq. 2) then
!       PASSAGE DES CONTRAINTES DU REPERE LOCAL 1
!       A L'ELEMENT AU REPERE INTRINSEQUE DE LA COQUE REPERE = 'COQUE_UTIL_INTR'
        if (option .eq. 'REPE_TENS') then
            if (pain(4:5) .eq. 'NO') then
                call vdsiro(np, nbsp, matvn1, 'IU', 'N', zr(jin), zr(jout))
            else if (pain(4:5) .eq. 'GA') then
                call vdsiro(np, nbsp, matvg1, 'IU', 'G', zr(jin), zr(jout))
            end if
        else if (option .eq. 'REPE_GENE') then
            if (pain(4:5) .eq. 'NO') then
                call vdefro(np, matvn1, zr(jin), zr(jout))
            else if (pain(4:5) .eq. 'GA') then
                call vdefro(np, matvg1, zr(jin), zr(jout))
            end if
        end if
    else
        messk(1) = nomte
        messk(2) = 'COQUE_UTIL_CYL'
        call utmess('F', 'ALGORITH12_41', nk=2, valk=messk)
    end if
!
end subroutine
