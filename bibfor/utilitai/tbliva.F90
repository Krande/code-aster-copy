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
subroutine tbliva(nomta, npacri, lipacr, vi, vr, &
                  vc, vk, crit, prec, para, &
                  ctype, vali, valr, valc, valk, &
                  ier)
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
    integer(kind=8) :: npacri, vi(*), vali, ier
    real(kind=8) :: vr(*), valr, prec(*)
    complex(kind=8) :: vc(*), valc
    character(len=*) :: nomta, lipacr(*), vk(*), valk, crit(*), ctype, para
!      LECTURE D'UNE VALEUR D'UNE CELLULE DE LA TABLE.
! ----------------------------------------------------------------------
! IN  : NOMTA  : NOM DE LA STRUCTURE "TABLE".
! IN  : NPACRI : NOMBRE DE PARAMETRES IMPLIQUES DANS LES CRITERES
! IN  : LIPACR : LISTE DES PARAMETRES CRITERES
! IN  : VI     : LISTE DES CRITERES POUR LES PARAMETRES "I"
! IN  : VR     : LISTE DES CRITERES POUR LES PARAMETRES "R"
! IN  : VC     : LISTE DES CRITERES POUR LES PARAMETRES "C"
! IN  : VK     : LISTE DES CRITERES POUR LES PARAMETRES "K"
! IN  : CRIT   : CRITERE POUR LES PARAMETRES REELS
! IN  : PREC   : PRECISION POUR LES PARAMETRES REELS
! IN  : PARA   : PARAMETRE A TROUVER
! OUT : CTYPE  : TYPE DE LA VALEUR TROUVEE
! OUT : VALI   : VALEUR TROUVEE SI PARAMETRES "I"
! OUT : VALR   : VALEUR TROUVEE SI PARAMETRES "R"
! OUT : VALC   : VALEUR TROUVEE SI PARAMETRES "C"
! OUT : VALK   : VALEUR TROUVEE SI PARAMETRES "K"
! OUT : IER    : CODE RETOUR 0 : OK
!                            1 : PARA N'EXISTE PAS
!                            2 : PAS DE LIGNE TROUVEE
!                            3 : PLUSIEURS LIGNES TROUVEES
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
    integer(kind=8) :: iret, nbpara, nblign, nbpu
    integer(kind=8) :: i, j, k, n, jvale, itrouv
    integer(kind=8) :: ki, kr, kc, k8, jvall
    real(kind=8) :: refr, xr, epsi
    complex(kind=8) :: refc, xc
    character(len=4) :: rela, type
    character(len=19) :: nomtab
    character(len=24) :: nomjv, nomjvl, inpar, jnpar
    aster_logical :: lok
    integer(kind=8), pointer :: numero(:) => null()
    integer(kind=8), pointer :: tbnp(:) => null()
    character(len=24), pointer :: tblp(:) => null()
! ----------------------------------------------------------------------
!
    call jemarq()
!
    ier = 0
    ctype = '?'
!
    vali = 0
    valr = 0.d0
    valc = dcmplx(0.d0, 0.d0)
    valk = ' '
!
    nomtab = nomta
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
    if (nblign .eq. 0) then
        call utmess('F', 'UTILITAI4_66')
    end if
!
    call jeveuo(nomtab//'.TBLP', 'L', vk24=tblp)
!
!     --- VERIFICATION QUE LES PARAMETRES EXISTENT DANS LA TABLE ---
!
    do i = 1, npacri
        inpar = lipacr(i)
        do j = 1, nbpara
            jnpar = tblp(1+4*(j-1))
            if (inpar .eq. jnpar) goto 10
        end do
        ier = 1
        goto 999
10      continue
    end do
    inpar = para
    do j = 1, nbpara
        jnpar = tblp(1+4*(j-1))
        if (inpar .eq. jnpar) goto 16
    end do
    ier = 1
    goto 999
16  continue
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
    k8 = 0
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
                        if (zi(jvall+n-1) .eq. 0) goto 30
                        if (zi(jvale+n-1) .eq. vi(ki)) then
                            itrouv = itrouv+1
                            numero(itrouv) = n
                        end if
30                      continue
                    end do
                else if (type(1:1) .eq. 'R') then
                    kr = kr+1
                    rela = crit(kr)
                    epsi = prec(kr)
                    xr = vr(kr)
                    do k = 1, nbpu
                        n = numero(k)
                        if (zi(jvall+n-1) .eq. 0) goto 31
                        refr = zr(jvale+n-1)
                        if (rela .eq. 'RELA') then
                            lok = (abs(xr-refr) .le. epsi*abs(refr))
                        else if (rela .eq. 'EGAL') then
                            lok = (refr .eq. xr)
                        else
                            lok = (abs(xr-refr) .le. epsi)
                        end if
                        if (lok) then
                            itrouv = itrouv+1
                            numero(itrouv) = n
                        end if
31                      continue
                    end do
                else if (type(1:1) .eq. 'C') then
                    kc = kc+1
                    rela = crit(kc)
                    epsi = prec(kc)
                    xc = vc(kc)
                    do k = 1, nbpu
                        n = numero(k)
                        if (zi(jvall+n-1) .eq. 0) goto 32
                        refc = zc(jvale+n-1)
                        if (rela .eq. 'RELA') then
                            lok = (abs(xc-refc) .le. epsi*abs(refc))
                        else
                            lok = (abs(xc-refc) .le. epsi)
                        end if
                        if (lok) then
                            itrouv = itrouv+1
                            numero(itrouv) = n
                        end if
32                      continue
                    end do
                else if (type(1:3) .eq. 'K80') then
                    k8 = k8+1
                    do k = 1, nbpu
                        n = numero(k)
                        if (zi(jvall+n-1) .eq. 0) goto 33
                        if (zk80(jvale+n-1) .eq. vk(k8)) then
                            itrouv = itrouv+1
                            numero(itrouv) = n
                        end if
33                      continue
                    end do
                else if (type(1:3) .eq. 'K32') then
                    k8 = k8+1
                    do k = 1, nbpu
                        n = numero(k)
                        if (zi(jvall+n-1) .eq. 0) goto 34
                        if (zk32(jvale+n-1) .eq. vk(k8)) then
                            itrouv = itrouv+1
                            numero(itrouv) = n
                        end if
34                      continue
                    end do
                else if (type(1:3) .eq. 'K24') then
                    k8 = k8+1
                    do k = 1, nbpu
                        n = numero(k)
                        if (zi(jvall+n-1) .eq. 0) goto 35
                        if (zk24(jvale+n-1) .eq. vk(k8)) then
                            itrouv = itrouv+1
                            numero(itrouv) = n
                        end if
35                      continue
                    end do
                else if (type(1:3) .eq. 'K16') then
                    k8 = k8+1
                    do k = 1, nbpu
                        n = numero(k)
                        if (zi(jvall+n-1) .eq. 0) goto 36
                        if (zk16(jvale+n-1) .eq. vk(k8)) then
                            itrouv = itrouv+1
                            numero(itrouv) = n
                        end if
36                      continue
                    end do
                else if (type(1:2) .eq. 'K8') then
                    k8 = k8+1
                    do k = 1, nbpu
                        n = numero(k)
                        if (zi(jvall+n-1) .eq. 0) goto 37
                        if (zk8(jvale+n-1) .eq. vk(k8)) then
                            itrouv = itrouv+1
                            numero(itrouv) = n
                        end if
37                      continue
                    end do
                end if
            end if
        end do
        if (itrouv .eq. 0) then
            ier = 2
            goto 999
        end if
        nbpu = itrouv
    end do
!
    itrouv = 0
    inpar = para
    do j = 1, nbpara
        jnpar = tblp(1+4*(j-1))
        if (inpar .eq. jnpar) then
            type = tblp(1+4*(j-1)+1)
            nomjv = tblp(1+4*(j-1)+2)
            nomjvl = tblp(1+4*(j-1)+3)
            call jeveuo(nomjv, 'L', jvale)
            call jeveuo(nomjvl, 'L', jvall)
            if (type(1:1) .eq. 'I') then
                do k = 1, nbpu
                    n = numero(k)
                    if (zi(jvall+n-1) .eq. 0) goto 50
                    itrouv = itrouv+1
                    vali = zi(jvale+n-1)
                    ctype = 'I'
50                  continue
                end do
            else if (type(1:1) .eq. 'R') then
                do k = 1, nbpu
                    n = numero(k)
                    if (zi(jvall+n-1) .eq. 0) goto 51
                    itrouv = itrouv+1
                    valr = zr(jvale+n-1)
                    ctype = 'R'
51                  continue
                end do
            else if (type(1:1) .eq. 'C') then
                do k = 1, nbpu
                    n = numero(k)
                    if (zi(jvall+n-1) .eq. 0) goto 52
                    itrouv = itrouv+1
                    valc = zc(jvale+n-1)
                    ctype = 'C'
52                  continue
                end do
            else if (type(1:3) .eq. 'K80') then
                k8 = k8+1
                do k = 1, nbpu
                    n = numero(k)
                    if (zi(jvall+n-1) .eq. 0) goto 53
                    itrouv = itrouv+1
                    valk = zk80(jvale+n-1)
                    ctype = 'K'
53                  continue
                end do
            else if (type(1:3) .eq. 'K32') then
                do k = 1, nbpu
                    n = numero(k)
                    if (zi(jvall+n-1) .eq. 0) goto 54
                    itrouv = itrouv+1
                    valk = zk32(jvale+n-1)
                    ctype = 'K'
54                  continue
                end do
            else if (type(1:3) .eq. 'K24') then
                do k = 1, nbpu
                    n = numero(k)
                    if (zi(jvall+n-1) .eq. 0) goto 55
                    itrouv = itrouv+1
                    valk = zk24(jvale+n-1)
                    ctype = 'K'
55                  continue
                end do
            else if (type(1:3) .eq. 'K16') then
                do k = 1, nbpu
                    n = numero(k)
                    if (zi(jvall+n-1) .eq. 0) goto 56
                    itrouv = itrouv+1
                    valk = zk16(jvale+n-1)
                    ctype = 'K'
56                  continue
                end do
            else if (type(1:2) .eq. 'K8') then
                do k = 1, nbpu
                    n = numero(k)
                    if (zi(jvall+n-1) .eq. 0) goto 57
                    itrouv = itrouv+1
                    valk = zk8(jvale+n-1)
                    ctype = 'K'
57                  continue
                end do
            end if
        end if
    end do
!
    if (itrouv .eq. 0) ier = 2
    if (itrouv .gt. 1) ier = 3
!
999 continue
    AS_DEALLOCATE(vi=numero)
!
    call jedema()
end subroutine
