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
subroutine tbajli(nomta, nbpar, nompar, vi, vr, &
                  vc, vk, nume)
    implicit none
#include "jeveux.h"
#include "asterc/ismaem.h"
#include "asterc/r8vide.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/juveca.h"
#include "asterfort/tbadap.h"
#include "asterfort/utmess.h"
!
    integer(kind=8) :: nbpar, vi(*), nume
    real(kind=8) :: vr(*)
    complex(kind=8) :: vc(*)
    character(len=*) :: nomta, nompar(*), vk(*)
!      AJOUTER UNE LIGNE A UNE TABLE.
! ----------------------------------------------------------------------
! IN  : NOMTA  : NOM DE LA STRUCTURE "TABLE".
! IN  : NBPAR  : NOMBRE DE PARAMETRES DE NOMPAR
! IN  : NOMPAR : NOMS DES PARAMETRES DE LA LIGNE
! IN  : VI     : LISTE DES VALEURS POUR LES PARAMETRES "I"
! IN  : VR     : LISTE DES VALEURS POUR LES PARAMETRES "R"
! IN  : VC     : LISTE DES VALEURS POUR LES PARAMETRES "C"
! IN  : VK     : LISTE DES VALEURS POUR LES PARAMETRES "K"
! IN  : NUME   : NUMERO DE LIGNE
!                = 0 : ON AJOUTE UNE LIGNE
!                > 0 : ON SURCHARGE UNE LIGNE
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
    integer(kind=8) :: iret, nbpara, nblign, nbpm, nbpu
    integer(kind=8) :: ndim, i, j, jvale, jlogq, ki, kr, kc, kk
    character(len=3) :: type
    character(len=19) :: nomtab
    character(len=24) :: nomjv, nomjvl, inpar, jnpar
    character(len=24), pointer :: tblp(:) => null()
    integer(kind=8), pointer :: tbnp(:) => null()
! ----------------------------------------------------------------------
!
    call jemarq()
!
    call tbadap(nomta, nbpar, nompar, vi, vr, &
                vc, vk)
!
    nomtab = ' '
    nomtab = nomta
    call jeexin(nomtab//'.TBBA', iret)
    if (iret .eq. 0) then
        call utmess('F', 'UTILITAI4_64')
    end if
    if (nomtab(18:19) .ne. '  ') then
        call utmess('F', 'UTILITAI4_68')
    end if
!
    call jeveuo(nomtab//'.TBNP', 'E', vi=tbnp)
    nbpara = tbnp(1)
    nblign = tbnp(2)
    if (nbpara .eq. 0) then
        call utmess('F', 'UTILITAI4_65')
    end if
    if (nume .lt. 0) then
        call utmess('F', 'UTILITAI4_70')
    end if
    if (nume .gt. nblign) then
        call utmess('F', 'UTILITAI4_74')
    end if
!
    call jeveuo(nomtab//'.TBLP', 'L', vk24=tblp)
    nomjv = tblp(3)
    call jelira(nomjv, 'LONMAX', nbpm)
    call jelira(nomjv, 'LONUTI', nbpu)
!
    ndim = nbpu+1
    if (ndim .gt. nbpm) then
        ndim = nint((ndim*3.d0)/2.d0)
        do i = 1, nbpara
            nomjv = tblp(1+4*(i-1)+2)
            call juveca(nomjv, ndim)
            nomjv = tblp(1+4*(i-1)+3)
            call juveca(nomjv, ndim)
        end do
    end if
!
    if (nume .eq. 0) then
        nblign = nblign+1
        tbnp(2) = nblign
!
        do i = 1, nbpara
            nomjv = tblp(1+4*(i-1)+2)
            call jeecra(nomjv, 'LONUTI', nblign)
            nomjv = tblp(1+4*(i-1)+3)
            call jeecra(nomjv, 'LONUTI', nblign)
        end do
!
        ki = 0
        kr = 0
        kc = 0
        kk = 0
        do j = 1, nbpar
            jnpar = nompar(j)
            do i = 1, nbpara
                inpar = tblp(1+4*(i-1))
                if (jnpar .eq. inpar) then
                    type = tblp(1+4*(i-1)+1) (1:3)
                    nomjv = tblp(1+4*(i-1)+2)
                    nomjvl = tblp(1+4*(i-1)+3)
                    call jeveuo(nomjv, 'E', jvale)
                    call jeveuo(nomjvl, 'E', jlogq)
                    if (type(1:1) .eq. 'I') then
                        ki = ki+1
                        if (vi(ki) .eq. ismaem()) then
                            zi(jlogq+nblign-1) = 0
                        else
                            zi(jvale+nblign-1) = vi(ki)
                            zi(jlogq+nblign-1) = 1
                        end if
                    else if (type(1:1) .eq. 'R') then
                        kr = kr+1
                        if (vr(kr) .eq. r8vide()) then
                            zi(jlogq+nblign-1) = 0
                        else
                            zr(jvale+nblign-1) = vr(kr)
                            zi(jlogq+nblign-1) = 1
                        end if
                    else if (type(1:1) .eq. 'C') then
                        kc = kc+1
                        if (dble(vc(kc)) .eq. r8vide() .and. dimag(vc(kc)) .eq. r8vide()) then
                            zi(jlogq+nblign-1) = 0
                        else
                            zc(jvale+nblign-1) = vc(kc)
                            zi(jlogq+nblign-1) = 1
                        end if
                    else if (type(1:3) .eq. 'K80') then
                        kk = kk+1
                        if (vk(kk) (1:7) .eq. '???????') then
                            zi(jlogq+nblign-1) = 0
                        else
                            zk80(jvale+nblign-1) = vk(kk)
                            zi(jlogq+nblign-1) = 1
                        end if
                    else if (type(1:3) .eq. 'K32') then
                        kk = kk+1
                        if (vk(kk) (1:7) .eq. '???????') then
                            zi(jlogq+nblign-1) = 0
                        else
                            zk32(jvale+nblign-1) = vk(kk)
                            zi(jlogq+nblign-1) = 1
                        end if
                    else if (type(1:3) .eq. 'K24') then
                        kk = kk+1
                        if (vk(kk) (1:7) .eq. '???????') then
                            zi(jlogq+nblign-1) = 0
                        else
                            zk24(jvale+nblign-1) = vk(kk)
                            zi(jlogq+nblign-1) = 1
                        end if
                    else if (type(1:3) .eq. 'K16') then
                        kk = kk+1
                        if (vk(kk) (1:7) .eq. '???????') then
                            zi(jlogq+nblign-1) = 0
                        else
                            zk16(jvale+nblign-1) = vk(kk)
                            zi(jlogq+nblign-1) = 1
                        end if
                    else if (type(1:2) .eq. 'K8') then
                        kk = kk+1
                        if (vk(kk) (1:7) .eq. '???????') then
                            zi(jlogq+nblign-1) = 0
                        else
                            zk8(jvale+nblign-1) = vk(kk)
                            zi(jlogq+nblign-1) = 1
                        end if
                    end if
                    goto 34
                end if
            end do
            call utmess('F', 'TABLE0_1', sk=jnpar)
34          continue
        end do
!
    else
        ki = 0
        kr = 0
        kc = 0
        kk = 0
        do j = 1, nbpar
            jnpar = nompar(j)
            do i = 1, nbpara
                inpar = tblp(1+4*(i-1))
                if (jnpar .eq. inpar) then
                    type = tblp(1+4*(i-1)+1) (1:3)
                    nomjv = tblp(1+4*(i-1)+2)
                    nomjvl = tblp(1+4*(i-1)+3)
                    call jeveuo(nomjv, 'E', jvale)
                    call jeveuo(nomjvl, 'E', jlogq)
                    if (type(1:1) .eq. 'I') then
                        ki = ki+1
                        if (vi(ki) .eq. ismaem()) then
                            zi(jlogq+nblign-1) = 0
                        else
                            zi(jvale+nume-1) = vi(ki)
                            zi(jlogq+nume-1) = 1
                        end if
                    else if (type(1:1) .eq. 'R') then
                        kr = kr+1
                        if (vr(kr) .eq. r8vide()) then
                            zi(jlogq+nume-1) = 0
                        else
                            zr(jvale+nume-1) = vr(kr)
                            zi(jlogq+nume-1) = 1
                        end if
                    else if (type(1:1) .eq. 'C') then
                        kc = kc+1
                        if (dble(vc(kc)) .eq. r8vide() .and. dimag(vc(kc)) .eq. r8vide()) then
                            zi(jlogq+nume-1) = 0
                        else
                            zc(jvale+nume-1) = vc(kc)
                            zi(jlogq+nume-1) = 1
                        end if
                    else if (type(1:3) .eq. 'K80') then
                        kk = kk+1
                        if (vk(kk) (1:7) .eq. '???????') then
                            zi(jlogq+nblign-1) = 0
                        else
                            zk80(jvale+nume-1) = vk(kk)
                            zi(jlogq+nume-1) = 1
                        end if
                    else if (type(1:3) .eq. 'K32') then
                        kk = kk+1
                        if (vk(kk) (1:7) .eq. '???????') then
                            zi(jlogq+nblign-1) = 0
                        else
                            zk32(jvale+nume-1) = vk(kk)
                            zi(jlogq+nume-1) = 1
                        end if
                    else if (type(1:3) .eq. 'K24') then
                        kk = kk+1
                        if (vk(kk) (1:7) .eq. '???????') then
                            zi(jlogq+nblign-1) = 0
                        else
                            zk24(jvale+nume-1) = vk(kk)
                            zi(jlogq+nume-1) = 1
                        end if
                    else if (type(1:3) .eq. 'K16') then
                        kk = kk+1
                        if (vk(kk) (1:7) .eq. '???????') then
                            zi(jlogq+nblign-1) = 0
                        else
                            zk16(jvale+nume-1) = vk(kk)
                            zi(jlogq+nume-1) = 1
                        end if
                    else if (type(1:2) .eq. 'K8') then
                        kk = kk+1
                        if (vk(kk) (1:7) .eq. '???????') then
                            zi(jlogq+nblign-1) = 0
                        else
                            zk8(jvale+nume-1) = vk(kk)
                            zi(jlogq+nume-1) = 1
                        end if
                    end if
                    goto 44
                end if
            end do
            call utmess('F', 'TABLE0_1', sk=jnpar)
44          continue
        end do
!
    end if
!
    call jedema()
end subroutine
