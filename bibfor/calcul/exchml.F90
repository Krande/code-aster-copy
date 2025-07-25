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

subroutine exchml(imodat, iparg)

    use calcul_module, only: ca_iachii_, ca_iachik_, ca_iachin_, ca_iachlo_, &
                             ca_iamloc_, ca_iawlo2_, ca_igr_, ca_iichin_, &
                             ca_ilchlo_, ca_ilmloc_, ca_nbelgr_, ca_nbgr_, &
                             ca_ncmpmx_, ca_nec_, ca_typegd_, &
                             ca_lparal_, ca_paral_, ca_iel_, ca_iachid_

    implicit none

! person_in_charge: jacques.pellet at edf.fr

#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/chloet.h"
#include "asterfort/exisdg.h"
#include "asterfort/jacopo.h"
#include "asterfort/utmess.h"

    integer(kind=8) :: iparg, imodat
!----------------------------------------------------------------------
!     entrees:
!        imodat : mode local attendu
!        iparg  : numero du parametre dans l'option
!----------------------------------------------------------------------
    integer(kind=8) :: jceld, mode, debgr2, lggre2, iaux1
    integer(kind=8) :: itypl1, modlo1, nbpoi1, lgcata
    integer(kind=8) :: itypl2, modlo2, nbpoi2
    integer(kind=8) :: ncmp1, ncmp2
    integer(kind=8) ::  debugr, lggrel
    integer(kind=8) :: jec, ncmp, jad1, jad2, ipt2, k, ipt1
    integer(kind=8) :: nbpoi, icmp1, icmp2, kcmp, ipt
    aster_logical :: etendu, lverec
    character(len=8) :: tych, cas
! --------------------------------------------------------------------------------------------------
!
    tych = zk8(ca_iachik_-1+2*(ca_iichin_-1)+1)
    ASSERT(tych(1:4) .eq. 'CHML')
!
    jceld = zi(ca_iachii_-1+ca_iachid_*(ca_iichin_-1)+4)
    lggre2 = zi(jceld-1+zi(jceld-1+4+ca_igr_)+4)
    debgr2 = zi(jceld-1+zi(jceld-1+4+ca_igr_)+8)
!
    mode = zi(jceld-1+zi(jceld-1+4+ca_igr_)+2)
!
    lgcata = zi(ca_iawlo2_-1+5*(ca_nbgr_*(iparg-1)+ca_igr_-1)+2)
    lggrel = zi(ca_iawlo2_-1+5*(ca_nbgr_*(iparg-1)+ca_igr_-1)+4)
    debugr = zi(ca_iawlo2_-1+5*(ca_nbgr_*(iparg-1)+ca_igr_-1)+5)
!
!   si mode=0 : il faut mettre champ_loc.exis a .false.
    if (mode .eq. 0) then
        do k = 1, lggrel
            zl(ca_ilchlo_-1+debugr-1+k) = .false.
        end do
        goto 999
    end if
!
!   si le champ a le mode attendu : on recopie
    if (mode .eq. imodat) then
        call jacopo(lggrel, ca_typegd_, ca_iachin_-1+debgr2, ca_iachlo_+debugr-1)
        goto 998
    end if
!   si le champ n'a pas le mode attendu ...
    call chloet(iparg, etendu, jceld)
    if (etendu) then
        call utmess('F', 'CALCUL_8')
    end if
!
    modlo1 = ca_iamloc_-1+zi(ca_ilmloc_-1+mode)
    modlo2 = ca_iamloc_-1+zi(ca_ilmloc_-1+imodat)
    itypl1 = zi(modlo1-1+1)
    itypl2 = zi(modlo2-1+1)
    ASSERT(itypl1 .le. 3)
    ASSERT(itypl2 .le. 3)
    nbpoi1 = zi(modlo1-1+4)
    nbpoi2 = zi(modlo2-1+4)
!
    ncmp1 = lggre2/(nbpoi1*ca_nbelgr_)
    ncmp2 = lgcata/nbpoi2
!   on verifie que les points ne sont pas "diff__" :
    ASSERT(nbpoi1 .lt. 10000)
    ASSERT(nbpoi2 .lt. 10000)
!   Dans quel cas de figure se trouve-t-on ?
    lverec = .true.
    if (nbpoi1 .eq. nbpoi2) then
        if (ncmp1 .eq. ncmp2) then
!           le cas "copie" est bizarre : il s'agit de 2 modes locaux
!           de meme contenu mais de noms differents. Faut-il l'interdire ?
            cas = 'COPIE'
            ncmp = ncmp1
!           quelques verifications
            ASSERT(itypl1 .eq. itypl2)
!           pour les champs elga, on verifie que c'est la meme famille
            if (itypl1 .eq. 3) then
                ASSERT(zi(modlo1+4+ca_nec_) .eq. zi(modlo2+4+ca_nec_))
            end if
        else
            cas = 'TRICMP'
            lverec = .false.
        end if
    else
        ASSERT(ncmp1 .eq. ncmp2)
        ncmp = ncmp1
        if (nbpoi1 .eq. 1) then
            cas = 'EXPAND'
        else if (nbpoi2 .eq. 1) then
            cas = 'MOYENN'
        else
            ASSERT(.false.)
        end if
    end if
    if (lverec) then
!       on verifie que les cmps sont les memes :
!          (sinon il faudrait trier ... => a faire (trigd) )
        do jec = 1, ca_nec_
            ASSERT(zi(modlo1-1+4+jec) .eq. zi(modlo2-1+4+jec))
        end do
    end if
!
    if (cas .eq. 'EXPAND' .or. cas .eq. 'COPIE') then
!       cas "expand" ou "copie"
        do ca_iel_ = 1, ca_nbelgr_
            if (ca_lparal_) then
                if (.not. ca_paral_(ca_iel_)) cycle
            end if
            if (cas .eq. 'EXPAND') then
                jad1 = ca_iachin_-1+debgr2+(ca_iel_-1)*ncmp
                do ipt2 = 1, nbpoi2
                    jad2 = ca_iachlo_+debugr-1+((ca_iel_-1)*nbpoi2+ipt2-1)*ncmp
                    call jacopo(ncmp, ca_typegd_, jad1, jad2)
                end do
            else if (cas .eq. 'COPIE') then
                ASSERT(nbpoi1 .eq. nbpoi2)
                jad1 = ca_iachin_-1+debgr2+(ca_iel_-1)*ncmp*nbpoi1
                jad2 = ca_iachlo_-1+debugr+(ca_iel_-1)*ncmp*nbpoi1
                call jacopo(ncmp*nbpoi1, ca_typegd_, jad1, jad2)
            end if
        end do
    else if (cas .eq. 'TRICMP') then
!       cas "tricmp"
        nbpoi = nbpoi1
        icmp1 = 0
        icmp2 = 0
        do kcmp = 1, ca_ncmpmx_
            if (exisdg(zi(modlo2-1+5), kcmp)) then
                icmp2 = icmp2+1
                if (exisdg(zi(modlo1-1+5), kcmp)) then
                    icmp1 = icmp1+1
                else
!                   -- a faire ... (gestion de zl)
                    ASSERT(.false.)
                end if
            else
                cycle
            end if
            ASSERT(icmp1 .ge. 1 .and. icmp1 .le. ncmp1)
            ASSERT(icmp2 .ge. 1 .and. icmp2 .le. ncmp2)
            do ca_iel_ = 1, ca_nbelgr_
                if (ca_lparal_) then
                    if (.not. ca_paral_(ca_iel_)) cycle
                end if

                do ipt = 1, nbpoi
                    jad1 = ca_iachin_+debgr2-1+((ca_iel_-1)*nbpoi+ipt-1)*ncmp1
                    jad2 = ca_iachlo_+debugr-1+((ca_iel_-1)*nbpoi+ipt-1)*ncmp2
                    jad1 = jad1-1+icmp1
                    jad2 = jad2-1+icmp2
                    call jacopo(1, ca_typegd_, jad1, jad2)
                end do
            end do
        end do
    else if (nbpoi2 .eq. 1) then
!       cas "moyenn"
        if (ca_typegd_ .eq. 'R') then
            if (ca_lparal_) then
                do ca_iel_ = 1, ca_nbelgr_
                    if (ca_paral_(ca_iel_)) then
                        iaux1 = ca_iachlo_+debugr-1+(ca_iel_-1)*ncmp
                        do k = 1, ncmp
                            zr(iaux1-1+k) = 0.d0
                        end do
                    end if
                end do
            else
                do k = 1, ca_nbelgr_*ncmp
                    zr(ca_iachlo_+debugr-1-1+k) = 0.d0
                end do
            end if
        else if (ca_typegd_ .eq. 'C') then
            if (ca_lparal_) then
                do ca_iel_ = 1, ca_nbelgr_
                    if (ca_paral_(ca_iel_)) then
                        iaux1 = ca_iachlo_+debugr-1+(ca_iel_-1)*ncmp
                        do k = 1, ncmp
                            zc(iaux1-1+k) = (0.d0, 0.d0)
                        end do
                    end if
                end do
            else
                do k = 1, ca_nbelgr_*ncmp
                    zc(ca_iachlo_+debugr-1-1+k) = (0.d0, 0.d0)
                end do
            end if
        else
            ASSERT(.false.)
        end if
!
        do ca_iel_ = 1, ca_nbelgr_
            if (ca_lparal_) then
                if (.not. ca_paral_(ca_iel_)) cycle
            end if
            jad2 = ca_iachlo_+debugr-1+(ca_iel_-1)*ncmp
            do ipt1 = 1, nbpoi1
                jad1 = ca_iachin_-1+debgr2+((ca_iel_-1)*nbpoi1+ipt1-1)*ncmp
                do k = 0, ncmp-1
                    if (ca_typegd_ .eq. 'R') then
                        zr(jad2+k) = zr(jad2+k)+zr(jad1+k)/dble(nbpoi1)
                    else if (ca_typegd_ .eq. 'C') then
                        zc(jad2+k) = zc(jad2+k)+zc(jad1+k)/dble(nbpoi1)
                    end if
                end do
            end do
        end do
    else
!       Autres cas pas encore programmes
        ASSERT(.false.)
    end if
!
998 continue
    do k = 1, lggrel
        zl(ca_ilchlo_-1+debugr-1+k) = .true.
    end do
!
999 continue
end subroutine
