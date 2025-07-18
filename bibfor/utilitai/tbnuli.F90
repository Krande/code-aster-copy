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
subroutine tbnuli(tabin, npacri, lipacr, vi, vr, &
                  vc, vk, lprec, lcrit, nume)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
!
    integer(kind=8) :: npacri, vi(*), nume
    real(kind=8) :: vr(*), lprec(*)
    complex(kind=8) :: vc(*)
    character(len=*) :: tabin, lipacr(*), vk(*), lcrit(*)
!     RECUPERATION D'UN NUMERO DE LIGNE
! ----------------------------------------------------------------------
! IN  : TABIN  : NOM DE LA TABLE DONT ON VEUT RECUPERER UNE LIGNE
! IN  : NPACRI : NOMBRE DE PARAMETRES IMPLIQUES DANS LES CRITERES
! IN  : LIPACR : LISTE DES PARAMETRES CRITERES
! IN  : VI     : LISTE DES CRITERES POUR LES PARAMETRES "I"
! IN  : VR     : LISTE DES CRITERES POUR LES PARAMETRES "R"
! IN  : VC     : LISTE DES CRITERES POUR LES PARAMETRES "C"
! IN  : VK     : LISTE DES CRITERES POUR LES PARAMETRES "K"
! IN  : LPREC  : PRECISION POUR LES PARAMETRES "R"
! IN  : LCRIT  : CRITERE POUR LES PARAMETRES "R"
! OUT : NUME   : = 0 , LA LIGNE N'A PAS PU ETRE RECUPERE
!                = I , ON A RECUPERE LA LIGNE
!                < 0 , PLUSIEURS LIGNES RECUPEREES
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
    integer(kind=8) :: iret, nbpara, nblign, nbpu, i, j, k, n
    integer(kind=8) :: jvale, itrouv, ki, kr, kc, kk, jvall
    real(kind=8) :: prec, refr
    character(len=4) :: type, crit
    character(len=19) :: nomtab
    character(len=24) :: nomjv, nomjvl, inpar, jnpar
    character(len=24) :: valk
    aster_logical :: lok
    integer(kind=8), pointer :: numero(:) => null()
    character(len=24), pointer :: tblp(:) => null()
    integer(kind=8), pointer :: tbnp(:) => null()
! ----------------------------------------------------------------------
!
    call jemarq()
!
    nume = 0
    nomtab = tabin
!
!     --- VERIFICATION DE LA TABLE ---
!
    call jeexin(nomtab//'.TBBA', iret)
    if (iret .eq. 0) then
        call utmess('F', 'UTILITAI4_64')
    end if
!
    call jeveuo(nomtab//'.TBNP', 'E', vi=tbnp)
    nbpara = tbnp(1)
    nblign = tbnp(2)
    if (nbpara .eq. 0) then
        call utmess('F', 'UTILITAI4_65')
    end if
    if (nblign .eq. 0) goto 999
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
        valk = inpar
        call utmess('F', 'UTILITAI6_89', sk=valk)
10      continue
    end do
!
    nomjv = tblp(3)
    call jelira(nomjv, 'LONUTI', nbpu)
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
        do j = 1, nbpara
            jnpar = tblp(1+4*(j-1))
            if (inpar .eq. jnpar) then
                type = tblp(1+4*(j-1)+1)
                nomjv = tblp(1+4*(j-1)+2)
                nomjvl = tblp(1+4*(j-1)+3)
                call jeveuo(nomjv, 'L', jvale)
                call jeveuo(nomjvl, 'L', jvall)
                if (type(1:1) .eq. 'I') then
                    ki = ki+1
                    do k = 1, nbpu
                        n = numero(k)
                        numero(k) = 0
                        if (zi(jvall+n-1) .eq. 0) goto 30
                        if (zi(jvale+n-1) .eq. vi(ki)) then
                            itrouv = itrouv+1
                            numero(itrouv) = n
                        end if
30                      continue
                    end do
                else if (type(1:1) .eq. 'R') then
                    kr = kr+1
                    prec = lprec(kr)
                    crit = lcrit(kr)
                    do k = 1, nbpu
                        n = numero(k)
                        numero(k) = 0
                        if (zi(jvall+n-1) .eq. 0) goto 31
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
31                      continue
                    end do
                else if (type(1:1) .eq. 'C') then
                    kc = kc+1
                    do k = 1, nbpu
                        n = numero(k)
                        numero(k) = 0
                        if (zi(jvall+n-1) .eq. 0) goto 32
                        if (zc(jvale+n-1) .eq. vc(kc)) then
                            itrouv = itrouv+1
                            numero(itrouv) = n
                        end if
32                      continue
                    end do
                else if (type(1:3) .eq. 'K80') then
                    kk = kk+1
                    do k = 1, nbpu
                        n = numero(k)
                        numero(k) = 0
                        if (zi(jvall+n-1) .eq. 0) goto 33
                        if (zk80(jvale+n-1) .eq. vk(kk)) then
                            itrouv = itrouv+1
                            numero(itrouv) = n
                        end if
33                      continue
                    end do
                else if (type(1:3) .eq. 'K32') then
                    kk = kk+1
                    do k = 1, nbpu
                        n = numero(k)
                        numero(k) = 0
                        if (zi(jvall+n-1) .eq. 0) goto 34
                        if (zk32(jvale+n-1) .eq. vk(kk)) then
                            itrouv = itrouv+1
                            numero(itrouv) = n
                        end if
34                      continue
                    end do
                else if (type(1:3) .eq. 'K24') then
                    kk = kk+1
                    do k = 1, nbpu
                        n = numero(k)
                        numero(k) = 0
                        if (zi(jvall+n-1) .eq. 0) goto 35
                        if (zk24(jvale+n-1) .eq. vk(kk)) then
                            itrouv = itrouv+1
                            numero(itrouv) = n
                        end if
35                      continue
                    end do
                else if (type(1:3) .eq. 'K16') then
                    kk = kk+1
                    do k = 1, nbpu
                        n = numero(k)
                        numero(k) = 0
                        if (zi(jvall+n-1) .eq. 0) goto 36
                        if (zk16(jvale+n-1) .eq. vk(kk)) then
                            itrouv = itrouv+1
                            numero(itrouv) = n
                        end if
36                      continue
                    end do
                else if (type(1:2) .eq. 'K8') then
                    kk = kk+1
                    do k = 1, nbpu
                        n = numero(k)
                        numero(k) = 0
                        if (zi(jvall+n-1) .eq. 0) goto 37
                        if (zk8(jvale+n-1) .eq. vk(kk)) then
                            itrouv = itrouv+1
                            numero(itrouv) = n
                        end if
37                      continue
                    end do
                end if
            end if
        end do
        nbpu = itrouv
    end do
!
    if (nbpu .eq. 1) then
        nume = numero(1)
    else if (nbpu .gt. 1) then
        nume = -nbpu
    end if
!
    AS_DEALLOCATE(vi=numero)
!
999 continue
    call jedema()
end subroutine
