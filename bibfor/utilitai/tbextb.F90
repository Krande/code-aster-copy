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
subroutine tbextb(tabin, basout, tabout, npacri, lipacr, &
                  lcrpa, vi, vr, vc, vk, &
                  lprec, lcrit, iret)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/ismaem.h"
#include "asterc/r8maem.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
#include "asterfort/utmess.h"
!
    integer(kind=8) :: npacri, vi(*), iret
    real(kind=8) :: vr(*), lprec(*)
    complex(kind=8) :: vc(*)
    character(len=*) :: tabin, basout, tabout, lipacr(*), lcrpa(*), vk(*)
    character(len=*) :: lcrit(*)
!     FILTRAGE ET EXTRACTION D'UNE NOUVELLE TABLE.
! ----------------------------------------------------------------------
! IN  : TABIN  : NOM DE LA TABLE DONT ON VEUT EXTRAIRE DES LIGNES
! IN  : BASOUT : BASE DE CREATION DE "TABOUT"
! IN  : TABOUT : NOM DE LA TABLE QUI CONTIENDRA LES LIGNES EXTRAITES
! IN  : NPACRI : NOMBRE DE PARAMETRES IMPLIQUES DANS LES CRITERES
! IN  : LIPACR : LISTE DES PARAMETRES CRITERES
! IN  : LCRPA  : LISTE DES CRITERES DE COMPARAISON
! IN  : VI     : LISTE DES CRITERES POUR LES PARAMETRES "I"
! IN  : VR     : LISTE DES CRITERES POUR LES PARAMETRES "R"
! IN  : VC     : LISTE DES CRITERES POUR LES PARAMETRES "C"
! IN  : VK     : LISTE DES CRITERES POUR LES PARAMETRES "K"
! IN  : LPREC  : PRECISION POUR LES PARAMETRES "R"
! IN  : LCRIT  : CRITERE POUR LES PARAMETRES "R"
! OUT : IRET   : =  0 , OK
!                = 10 , LE PARAMETRE N'EXISTE PAS
!                = 20 , PAS DE LIGNES POUR LE PARAMETRE DONNE
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
    integer(kind=8) :: irt, nbpara, nblign, nbpu
    integer(kind=8) :: i, j, k, n, jvale, itrouv, itrou2
    integer(kind=8) :: ki, kr, kc, kk, jvall, nbp
    integer(kind=8) :: imax, imin
    real(kind=8) :: prec, refr, rmax, rmin
    complex(kind=8) :: cmin, cmax
    character(len=1) :: base
    character(len=4) :: type, crit
    character(len=8) :: rela
    character(len=19) :: nomtab
    character(len=24) :: nomjv, nomjvl, inpar, jnpar
    aster_logical :: lok
    integer(kind=8), pointer :: numero(:) => null()
    character(len=24), pointer :: para_r(:) => null()
    character(len=8), pointer :: type_r(:) => null()
    complex(kind=8), pointer :: vale_c(:) => null()
    integer(kind=8), pointer :: vale_i(:) => null()
    character(len=80), pointer :: vale_k(:) => null()
    real(kind=8), pointer :: vale_r(:) => null()
    integer(kind=8), pointer :: tbnp(:) => null()
    character(len=24), pointer :: tblp(:) => null()
! ----------------------------------------------------------------------
!
    call jemarq()
!
    nomtab = tabin
    base = basout(1:1)
    iret = 0
!
!     --- VERIFICATION DE LA BASE ---
!
    ASSERT(base .eq. 'V' .or. base .eq. 'G')
!
!     --- VERIFICATION DE LA TABLE ---
!
    call jeexin(nomtab//'.TBBA', irt)
    if (irt .eq. 0) then
        call utmess('F', 'UTILITAI4_64')
    end if
!
    call jeveuo(nomtab//'.TBNP', 'E', vi=tbnp)
    nbpara = tbnp(1)
    nblign = tbnp(2)
    if (nbpara .eq. 0) then
        call utmess('F', 'UTILITAI4_65')
    end if
    if (nblign .eq. 0) then
        call utmess('F', 'UTILITAI4_66')
    end if
!
!     --- VERIFICATION QUE LES PARAMETRES EXISTENT DANS LA TABLE ---
!
    call jeveuo(nomtab//'.TBLP', 'L', vk24=tblp)
    do i = 1, npacri
        inpar = lipacr(i)
        do j = 1, nbpara
            jnpar = tblp(1+4*(j-1))
            if (inpar .eq. jnpar) goto 10
        end do
        iret = 10
        goto 999
10      continue
    end do
!
    nbpu = nblign
    AS_ALLOCATE(vi=numero, size=nbpu)
    do i = 1, nbpu
        numero(i) = i
    end do
!
    ki = 0
    kr = 0
    kc = 0
    kk = 0
    do i = 1, npacri
        itrouv = 0
        inpar = lipacr(i)
        rela = lcrpa(i)
        do j = 1, nbpara
            jnpar = tblp(1+4*(j-1))
            if (inpar .eq. jnpar) then
                type = tblp(1+4*(j-1)+1)
                nomjv = tblp(1+4*(j-1)+2)
                nomjvl = tblp(1+4*(j-1)+3)
                call jeveuo(nomjv, 'L', jvale)
                call jeveuo(nomjvl, 'L', jvall)
                if (rela .eq. 'VIDE') then
                    do k = 1, nbpu
                        n = numero(k)
                        numero(k) = 0
                        if (zi(jvall+n-1) .eq. 0) then
                            itrouv = itrouv+1
                            numero(itrouv) = n
                        end if
                    end do
                    goto 24
                else if (rela .eq. 'NON_VIDE') then
                    do k = 1, nbpu
                        n = numero(k)
                        numero(k) = 0
                        if (zi(jvall+n-1) .eq. 1) then
                            itrouv = itrouv+1
                            numero(itrouv) = n
                        end if
                    end do
                    goto 24
                else if (rela .eq. 'MAXI') then
                    itrou2 = 0
                    if (type(1:1) .eq. 'I') then
                        imax = -ismaem()
                        do k = 1, nbpu
                            n = numero(k)
                            numero(k) = 0
                            if (zi(jvall+n-1) .eq. 0) goto 52
                            if (zi(jvale+n-1) .gt. imax) then
                                imax = zi(jvale+n-1)
                                itrou2 = n
                            end if
52                          continue
                        end do
                    else if (type(1:1) .eq. 'R') then
                        rmax = -r8maem()
                        do k = 1, nbpu
                            n = numero(k)
                            numero(k) = 0
                            if (zi(jvall+n-1) .eq. 0) goto 53
                            if (zr(jvale+n-1) .gt. rmax) then
                                rmax = zr(jvale+n-1)
                                itrou2 = n
                            end if
53                          continue
                        end do
                    end if
                    if (itrou2 .ne. 0) then
                        itrouv = itrouv+1
                        numero(itrouv) = itrou2
                    end if
                    goto 24
                else if (rela .eq. 'MAXI_ABS') then
                    itrou2 = 0
                    if (type(1:1) .eq. 'I') then
                        imax = -ismaem()
                        do k = 1, nbpu
                            n = numero(k)
                            numero(k) = 0
                            if (zi(jvall+n-1) .eq. 0) goto 54
                            if (abs(zi(jvale+n-1)) .gt. imax) then
                                imax = abs(zi(jvale+n-1))
                                itrou2 = n
                            end if
54                          continue
                        end do
                    else if (type(1:1) .eq. 'R') then
                        rmax = -r8maem()
                        do k = 1, nbpu
                            n = numero(k)
                            numero(k) = 0
                            if (zi(jvall+n-1) .eq. 0) goto 55
                            if (abs(zr(jvale+n-1)) .gt. rmax) then
                                rmax = abs(zr(jvale+n-1))
                                itrou2 = n
                            end if
55                          continue
                        end do
                    else if (type(1:1) .eq. 'C') then
                        cmax = dcmplx(-1.d-50, -1.d-50)
                        do k = 1, nbpu
                            n = numero(k)
                            numero(k) = 0
                            if (zi(jvall+n-1) .eq. 0) goto 60
                            if (abs(zc(jvale+n-1)) .gt. abs(cmax)) then
                                cmax = zc(jvale+n-1)
                                itrou2 = n
                            end if
60                          continue
                        end do
                    end if
                    if (itrou2 .ne. 0) then
                        itrouv = itrouv+1
                        numero(itrouv) = itrou2
                    end if
                    goto 24
                else if (rela .eq. 'MINI') then
                    itrou2 = 0
                    if (type(1:1) .eq. 'I') then
                        imin = ismaem()
                        do k = 1, nbpu
                            n = numero(k)
                            numero(k) = 0
                            if (zi(jvall+n-1) .eq. 0) goto 56
                            if (zi(jvale+n-1) .lt. imin) then
                                imin = zi(jvale+n-1)
                                itrou2 = n
                            end if
56                          continue
                        end do
                    else if (type(1:1) .eq. 'R') then
                        rmin = r8maem()
                        do k = 1, nbpu
                            n = numero(k)
                            numero(k) = 0
                            if (zi(jvall+n-1) .eq. 0) goto 57
                            if (zr(jvale+n-1) .lt. rmin) then
                                rmin = zr(jvale+n-1)
                                itrou2 = n
                            end if
57                          continue
                        end do
                    end if
                    if (itrou2 .ne. 0) then
                        itrouv = itrouv+1
                        numero(itrouv) = itrou2
                    end if
                    goto 24
                else if (rela .eq. 'MINI_ABS') then
                    itrou2 = 0
                    if (type(1:1) .eq. 'I') then
                        imin = ismaem()
                        do k = 1, nbpu
                            n = numero(k)
                            numero(k) = 0
                            if (zi(jvall+n-1) .eq. 0) goto 58
                            if (abs(zi(jvale+n-1)) .lt. imin) then
                                imin = abs(zi(jvale+n-1))
                                itrou2 = n
                            end if
58                          continue
                        end do
                    else if (type(1:1) .eq. 'R') then
                        rmin = r8maem()
                        do k = 1, nbpu
                            n = numero(k)
                            numero(k) = 0
                            if (zi(jvall+n-1) .eq. 0) goto 59
                            if (abs(zr(jvale+n-1)) .lt. rmin) then
                                rmin = abs(zr(jvale+n-1))
                                itrou2 = n
                            end if
59                          continue
                        end do
                    else if (type(1:1) .eq. 'C') then
                        cmin = dcmplx(+1.d+50, +1.d+50)
                        do k = 1, nbpu
                            n = numero(k)
                            numero(k) = 0
                            if (zi(jvall+n-1) .eq. 0) goto 61
                            if (abs(zc(jvale+n-1)) .lt. abs(cmin)) then
                                cmin = zc(jvale+n-1)
                                itrou2 = n
                            end if
61                          continue
                        end do
                    end if
                    if (itrou2 .ne. 0) then
                        itrouv = itrouv+1
                        numero(itrouv) = itrou2
                    end if
                    goto 24
                end if
                if (type(1:1) .eq. 'I') then
                    ki = ki+1
                    do k = 1, nbpu
                        n = numero(k)
                        numero(k) = 0
                        if (zi(jvall+n-1) .eq. 0) goto 30
                        if (rela .eq. 'EQ') then
                            if (zi(jvale+n-1) .eq. vi(ki)) then
                                itrouv = itrouv+1
                                numero(itrouv) = n
                            end if
                        else if (rela .eq. 'LT') then
                            if (zi(jvale+n-1) .lt. vi(ki)) then
                                itrouv = itrouv+1
                                numero(itrouv) = n
                            end if
                        else if (rela .eq. 'GT') then
                            if (zi(jvale+n-1) .gt. vi(ki)) then
                                itrouv = itrouv+1
                                numero(itrouv) = n
                            end if
                        else if (rela .eq. 'NE') then
                            if (zi(jvale+n-1) .ne. vi(ki)) then
                                itrouv = itrouv+1
                                numero(itrouv) = n
                            end if
                        else if (rela .eq. 'LE') then
                            if (zi(jvale+n-1) .le. vi(ki)) then
                                itrouv = itrouv+1
                                numero(itrouv) = n
                            end if
                        else if (rela .eq. 'GE') then
                            if (zi(jvale+n-1) .ge. vi(ki)) then
                                itrouv = itrouv+1
                                numero(itrouv) = n
                            end if
                        end if
30                      continue
                    end do
                    goto 24
                else if (type(1:1) .eq. 'R') then
                    kr = kr+1
                    prec = lprec(kr)
                    crit = lcrit(kr)
                    do k = 1, nbpu
                        n = numero(k)
                        numero(k) = 0
                        if (zi(jvall+n-1) .eq. 0) goto 31
                        if (rela .eq. 'EQ') then
                            refr = zr(jvale+n-1)
                            if (crit .eq. 'RELA') then
                                lok = (abs(vr(kr)-refr) .le. prec*abs(refr))
                            else if (crit .eq. 'EGAL') then
                                lok = (vr(kr) .eq. refr)
                            else
                                lok = (abs(vr(kr)-refr) .le. prec)
                            end if
                            if (lok) then
                                itrouv = itrouv+1
                                numero(itrouv) = n
                            end if
                        else if (rela .eq. 'LT') then
                            if (zr(jvale+n-1) .lt. vr(kr)) then
                                itrouv = itrouv+1
                                numero(itrouv) = n
                            end if
                        else if (rela .eq. 'GT') then
                            if (zr(jvale+n-1) .gt. vr(kr)) then
                                itrouv = itrouv+1
                                numero(itrouv) = n
                            end if
                        else if (rela .eq. 'NE') then
                            if (zr(jvale+n-1) .ne. vr(kr)) then
                                itrouv = itrouv+1
                                numero(itrouv) = n
                            end if
                        else if (rela .eq. 'LE') then
                            if (zr(jvale+n-1) .le. vr(kr)) then
                                itrouv = itrouv+1
                                numero(itrouv) = n
                            end if
                        else if (rela .eq. 'GE') then
                            if (zr(jvale+n-1) .ge. vr(kr)) then
                                itrouv = itrouv+1
                                numero(itrouv) = n
                            end if
                        end if
31                      continue
                    end do
                    goto 24
                else if (type(1:1) .eq. 'C') then
                    kc = kc+1
                    do k = 1, nbpu
                        n = numero(k)
                        numero(k) = 0
                        if (zi(jvall+n-1) .eq. 0) goto 32
                        if (rela .eq. 'EQ') then
                            if (zc(jvale+n-1) .eq. vc(kc)) then
                                itrouv = itrouv+1
                                numero(itrouv) = n
                            end if
                        else if (rela .eq. 'NE') then
                            if (zc(jvale+n-1) .ne. vc(kc)) then
                                itrouv = itrouv+1
                                numero(itrouv) = n
                            end if
                        end if
32                      continue
                    end do
                    goto 24
                else if (type(1:3) .eq. 'K80') then
                    kk = kk+1
                    do k = 1, nbpu
                        n = numero(k)
                        numero(k) = 0
                        if (zi(jvall+n-1) .eq. 0) goto 33
                        if (rela .eq. 'EQ') then
                            if (zk80(jvale+n-1) .eq. vk(kk)) then
                                itrouv = itrouv+1
                                numero(itrouv) = n
                            end if
                        else if (rela .eq. 'NE') then
                            if (zk80(jvale+n-1) .ne. vk(kk)) then
                                itrouv = itrouv+1
                                numero(itrouv) = n
                            end if
                        end if
33                      continue
                    end do
                    goto 24
                else if (type(1:3) .eq. 'K32') then
                    kk = kk+1
                    do k = 1, nbpu
                        n = numero(k)
                        numero(k) = 0
                        if (zi(jvall+n-1) .eq. 0) goto 34
                        if (rela .eq. 'EQ') then
                            if (zk32(jvale+n-1) .eq. vk(kk)) then
                                itrouv = itrouv+1
                                numero(itrouv) = n
                            end if
                        else if (rela .eq. 'NE') then
                            if (zk32(jvale+n-1) .ne. vk(kk)) then
                                itrouv = itrouv+1
                                numero(itrouv) = n
                            end if
                        end if
34                      continue
                    end do
                    goto 24
                else if (type(1:3) .eq. 'K24') then
                    kk = kk+1
                    do k = 1, nbpu
                        n = numero(k)
                        numero(k) = 0
                        if (zi(jvall+n-1) .eq. 0) goto 35
                        if (rela .eq. 'EQ') then
                            if (zk24(jvale+n-1) .eq. vk(kk)) then
                                itrouv = itrouv+1
                                numero(itrouv) = n
                            end if
                        else if (rela .eq. 'NE') then
                            if (zk24(jvale+n-1) .ne. vk(kk)) then
                                itrouv = itrouv+1
                                numero(itrouv) = n
                            end if
                        end if
35                      continue
                    end do
                    goto 24
                else if (type(1:3) .eq. 'K16') then
                    kk = kk+1
                    do k = 1, nbpu
                        n = numero(k)
                        numero(k) = 0
                        if (zi(jvall+n-1) .eq. 0) goto 36
                        if (rela .eq. 'EQ') then
                            if (zk16(jvale+n-1) .eq. vk(kk)) then
                                itrouv = itrouv+1
                                numero(itrouv) = n
                            end if
                        else if (rela .eq. 'NE') then
                            if (zk16(jvale+n-1) .ne. vk(kk)) then
                                itrouv = itrouv+1
                                numero(itrouv) = n
                            end if
                        end if
36                      continue
                    end do
                    goto 24
                else if (type(1:2) .eq. 'K8') then
                    kk = kk+1
                    do k = 1, nbpu
                        n = numero(k)
                        numero(k) = 0
                        if (zi(jvall+n-1) .eq. 0) goto 37
                        if (rela .eq. 'EQ') then
                            if (zk8(jvale+n-1) .eq. vk(kk)) then
                                itrouv = itrouv+1
                                numero(itrouv) = n
                            end if
                        else if (rela .eq. 'NE') then
                            if (zk8(jvale+n-1) .ne. vk(kk)) then
                                itrouv = itrouv+1
                                numero(itrouv) = n
                            end if
                        end if
37                      continue
                    end do
                    goto 24
                end if
            end if
        end do
24      continue
        if (itrouv .eq. 0) then
            iret = 20
            goto 999
        end if
        nbpu = itrouv
    end do
!
!     --- ON RECUPERE LES PARAMETRES ET LEURS TYPES  ---
!
    AS_ALLOCATE(vk8=type_r, size=nbpara)
    AS_ALLOCATE(vk24=para_r, size=nbpara)
    do i = 1, nbpara
        para_r(i) = tblp(1+4*(i-1))
        type_r(i) = tblp(1+4*(i-1)+1)
    end do
!
!     --- CREATION DE LA TABLE ---
!
    call tbcrsd(tabout, basout)
    call tbajpa(tabout, nbpara, para_r, type_r)
    AS_ALLOCATE(vi=vale_i, size=nbpara)
    AS_ALLOCATE(vr=vale_r, size=nbpara)
    AS_ALLOCATE(vc=vale_c, size=nbpara)
    AS_ALLOCATE(vk80=vale_k, size=nbpara)
    do k = 1, nbpu
        ki = 0
        kr = 0
        kc = 0
        kk = 0
        n = numero(k)
        nbp = 0
        do j = 1, nbpara
            jnpar = tblp(1+4*(j-1))
            type = tblp(1+4*(j-1)+1)
            nomjv = tblp(1+4*(j-1)+2)
            nomjvl = tblp(1+4*(j-1)+3)
            call jeveuo(nomjv, 'L', jvale)
            call jeveuo(nomjvl, 'L', jvall)
            if (zi(jvall+n-1) .eq. 0) goto 42
            nbp = nbp+1
            para_r(nbp) = jnpar
            if (type(1:1) .eq. 'I') then
                ki = ki+1
                vale_i(ki) = zi(jvale+n-1)
            else if (type(1:1) .eq. 'R') then
                kr = kr+1
                vale_r(kr) = zr(jvale+n-1)
            else if (type(1:1) .eq. 'C') then
                kc = kc+1
                vale_c(kc) = zc(jvale+n-1)
            else if (type(1:3) .eq. 'K80') then
                kk = kk+1
                vale_k(kk) = zk80(jvale+n-1)
            else if (type(1:3) .eq. 'K32') then
                kk = kk+1
                vale_k(kk) = zk32(jvale+n-1)
            else if (type(1:3) .eq. 'K24') then
                kk = kk+1
                vale_k(kk) = zk24(jvale+n-1)
            else if (type(1:3) .eq. 'K16') then
                kk = kk+1
                vale_k(kk) = zk16(jvale+n-1)
            else if (type(1:2) .eq. 'K8') then
                kk = kk+1
                vale_k(kk) = zk8(jvale+n-1)
            end if
42          continue
        end do
        if (nbp .eq. 0) goto 40
        call tbajli(tabout, nbp, para_r, vale_i, vale_r, &
                    vale_c, vale_k, 0)
40      continue
    end do
!
999 continue
    AS_DEALLOCATE(vi=numero)
    AS_DEALLOCATE(vk8=type_r)
    AS_DEALLOCATE(vk24=para_r)
    AS_DEALLOCATE(vi=vale_i)
    AS_DEALLOCATE(vr=vale_r)
    AS_DEALLOCATE(vc=vale_c)
    AS_DEALLOCATE(vk80=vale_k)
!
    call jedema()
end subroutine
