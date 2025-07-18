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
subroutine tbfutb(tabout, basout, ntab, ltabin, para, &
                  typpar, vi, vr, vc, vk)
    implicit none
#include "jeveux.h"
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
    integer(kind=8) :: ntab, vi(*)
    real(kind=8) :: vr(*)
    complex(kind=8) :: vc(*)
    character(len=*) :: tabout, basout, ltabin(*), para, typpar, vk(*)
!     FUSIONNER PLUSIEURS TABLES EN UNE SEULE TABLE.
! ----------------------------------------------------------------------
! IN  : TABOUT : NOM DE LA TABLE QUE L'ON VEUT OBTENIR
! IN  : BASOUT : BASE DE CREATION DE "TABOUT"
! IN  : NTAB   : NOMBRE DE TABLES QUE L'ON VEUT FUSIONNER
! IN  : LTABIN : NOMS DES TABLES QUE L'ON VEUT FUSIONNER
! IN  : PARA   : NOUVEAU PARAMETRE NECESSAIRE
! IN  : TYPPAR : TYPE DU NOUVEAU PARAMETRE
! IN  : VI     : LISTE DES CRITERES POUR LES PARAMETRES "I"
! IN  : VR     : LISTE DES CRITERES POUR LES PARAMETRES "R"
! IN  : VC     : LISTE DES CRITERES POUR LES PARAMETRES "C"
! IN  : VK     : LISTE DES CRITERES POUR LES PARAMETRES "K"
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
    integer(kind=8) :: iret, nbpara, nblign, jtbnp, nbpu, nbpart, ipar
    integer(kind=8) :: jtblp, i, j, k, jvale
    integer(kind=8) :: ki, kr, kc, kk, jvall
!
    character(len=1) :: base
    character(len=4) :: type, ktype
    character(len=19) :: nomtab
    character(len=24) :: nomjv, nomjvl, inpar, jnpar, knpar
    character(len=24) :: valk(3)
    character(len=24), pointer :: para_r(:) => null()
    character(len=8), pointer :: type_r(:) => null()
    complex(kind=8), pointer :: vale_c(:) => null()
    integer(kind=8), pointer :: vale_i(:) => null()
    character(len=80), pointer :: vale_k(:) => null()
    real(kind=8), pointer :: vale_r(:) => null()
! ----------------------------------------------------------------------
!
    call jemarq()
!
    base = basout(1:1)
!
!     --- VERIFICATION DE LA BASE ---
!
    ASSERT(base .eq. 'V' .or. base .eq. 'G')
!
!     --- VERIFICATION DES TABLES ---
!
    inpar = para
    nbpart = 0
    nbpu = 0
    do i = 1, ntab
        nomtab = ltabin(i)
        call jeexin(nomtab//'.TBBA', iret)
        if (iret .eq. 0) then
            call utmess('F', 'UTILITAI4_64')
        end if
!
        call jeveuo(nomtab//'.TBNP', 'L', jtbnp)
        nbpara = zi(jtbnp)
        nblign = zi(jtbnp+1)
        nbpart = nbpart+nbpara
        nbpu = max(nbpu, nblign)
        if (nbpara .eq. 0) then
            call utmess('F', 'UTILITAI4_65')
        end if
        if (nblign .eq. 0) then
            call utmess('F', 'UTILITAI4_66')
        end if
!
        call jeveuo(nomtab//'.TBLP', 'L', jtblp)
        do j = 1, nbpara
            jnpar = zk24(jtblp+4*(j-1))
            if (inpar .eq. jnpar) then
                valk(1) = jnpar
                valk(2) = nomtab
                call utmess('F', 'UTILITAI8_20', nk=2, valk=valk)
            end if
        end do
!
    end do
!
!     --- ON ELIMINE LES PARAMETRES DOUBLONS ---
!
    nbpart = nbpart+1
    AS_ALLOCATE(vk8=type_r, size=nbpart)
    AS_ALLOCATE(vk24=para_r, size=nbpart)
    ipar = 1
    if (para(1:1) .ne. ' ') then
        para_r(1) = para
        type_r(1) = typpar
    else
        para_r(1) = zk24(jtblp)
        type_r(1) = zk24(jtblp+1)
    end if
    do i = 1, ntab
        nomtab = ltabin(i)
        call jeveuo(nomtab//'.TBNP', 'L', jtbnp)
        call jeveuo(nomtab//'.TBLP', 'L', jtblp)
        nbpara = zi(jtbnp)
        do j = 1, nbpara
            jnpar = zk24(jtblp+4*(j-1))
            type = zk24(jtblp+4*(j-1)+1)
            do k = 1, ipar
                knpar = para_r(k)
                ktype = type_r(k)
                if (knpar .eq. jnpar) then
                    if (type .ne. ktype) then
                        valk(1) = jnpar
                        valk(2) = jnpar
                        valk(3) = knpar
                        call utmess('F', 'UTILITAI8_21', nk=3, valk=valk)
                    end if
                    goto 22
                end if
            end do
            ipar = ipar+1
            para_r(ipar) = jnpar
            type_r(ipar) = type
22          continue
        end do
    end do
    nbpart = ipar
!
!     --- CREATION DE LA TABLE ---
!
    call tbcrsd(tabout, basout)
    call tbajpa(tabout, nbpart, para_r, type_r)
    AS_ALLOCATE(vi=vale_i, size=nbpu)
    AS_ALLOCATE(vr=vale_r, size=nbpu)
    AS_ALLOCATE(vc=vale_c, size=nbpu)
    AS_ALLOCATE(vk80=vale_k, size=nbpu)
    do i = 1, ntab
        nomtab = ltabin(i)
        call jeveuo(nomtab//'.TBNP', 'L', jtbnp)
        call jeveuo(nomtab//'.TBLP', 'L', jtblp)
        nbpara = zi(jtbnp)
        nblign = zi(jtbnp+1)
        do k = 1, nblign
            ki = 0
            kr = 0
            kc = 0
            kk = 0
            if (para .ne. ' ') then
                ipar = 1
                para_r(1) = para
            else
                ipar = 0
            end if
            if (typpar(1:1) .eq. 'I') then
                ki = ki+1
                vale_i(ki) = vi(i)
            else if (typpar(1:1) .eq. 'R') then
                kr = kr+1
                vale_r(kr) = vr(i)
            else if (typpar(1:1) .eq. 'C') then
                kc = kc+1
                vale_c(kc) = vc(i)
            else if (typpar(1:3) .eq. 'K80') then
                kk = kk+1
                vale_k(kk) = vk(i)
            else if (typpar(1:3) .eq. 'K32') then
                kk = kk+1
                vale_k(kk) = vk(i)
            else if (typpar(1:3) .eq. 'K24') then
                kk = kk+1
                vale_k(kk) = vk(i)
            else if (typpar(1:3) .eq. 'K16') then
                kk = kk+1
                vale_k(kk) = vk(i)
            else if (typpar(1:2) .eq. 'K8') then
                kk = kk+1
                vale_k(kk) = vk(i)
            end if
            do j = 1, nbpara
                jnpar = zk24(jtblp+4*(j-1))
                type = zk24(jtblp+4*(j-1)+1)
                nomjv = zk24(jtblp+4*(j-1)+2)
                nomjvl = zk24(jtblp+4*(j-1)+3)
                call jeveuo(nomjv, 'L', jvale)
                call jeveuo(nomjvl, 'L', jvall)
                if (zi(jvall+k-1) .eq. 0) goto 42
                ipar = ipar+1
                para_r(ipar) = jnpar
                if (type(1:1) .eq. 'I') then
                    ki = ki+1
                    vale_i(ki) = zi(jvale+k-1)
                else if (type(1:1) .eq. 'R') then
                    kr = kr+1
                    vale_r(kr) = zr(jvale+k-1)
                else if (type(1:1) .eq. 'C') then
                    kc = kc+1
                    vale_c(kc) = zc(jvale+k-1)
                else if (type(1:3) .eq. 'K80') then
                    kk = kk+1
                    vale_k(kk) = zk80(jvale+k-1)
                else if (type(1:3) .eq. 'K32') then
                    kk = kk+1
                    vale_k(kk) = zk32(jvale+k-1)
                else if (type(1:3) .eq. 'K24') then
                    kk = kk+1
                    vale_k(kk) = zk24(jvale+k-1)
                else if (type(1:3) .eq. 'K16') then
                    kk = kk+1
                    vale_k(kk) = zk16(jvale+k-1)
                else if (type(1:2) .eq. 'K8') then
                    kk = kk+1
                    vale_k(kk) = zk8(jvale+k-1)
                end if
42              continue
            end do
            call tbajli(tabout, ipar, para_r, vale_i, vale_r, &
                        vale_c, vale_k, 0)
        end do
    end do
!
!
!
    AS_DEALLOCATE(vk8=type_r)
    AS_DEALLOCATE(vk24=para_r)
    AS_DEALLOCATE(vi=vale_i)
    AS_DEALLOCATE(vr=vale_r)
    AS_DEALLOCATE(vc=vale_c)
    AS_DEALLOCATE(vk80=vale_k)
!
    call jedema()
end subroutine
