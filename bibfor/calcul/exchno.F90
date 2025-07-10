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
subroutine exchno(imodat, iparg)
!
    use calcul_module, only: ca_iachii_, ca_iachlo_, ca_ialiel_, ca_iamaco_, &
                             ca_iamloc_, ca_iamsco_, ca_iawlo2_, ca_igr_, ca_iichin_, &
                             ca_illiel_, ca_ilmaco_, ca_ilmloc_, ca_ilmsco_, ca_nbelgr_, &
                             ca_nbgr_, ca_nec_, ca_typegd_, ca_lparal_, &
                             ca_paral_, ca_iel_, ca_iachid_
!
    implicit none
!
! person_in_charge: jacques.pellet at edf.fr
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/trigd.h"
!
    integer(kind=8) :: imodat, iparg
!----------------------------------------------------------------------
!     entrees:
!        imodat  : indice dans la collection modeloc
!        igr    : numero du grel a traiter.
!     sorties:
!        ecriture dans le champ local
!----------------------------------------------------------------------
!
    integer(kind=8) :: ima, ino, nno, nugl
    integer(kind=8) :: prno1, prno2, modloc, ityplo
    integer(kind=8) :: deb1, deb2, idg1, idg2, nbpt, nbpt2, lgcata, ncmp
    integer(kind=8) :: iaux1, k, debugr
    aster_logical :: diff, moyenn
!
! --------------------------------------------------------------------------------------------------
!
!   fonctions formules :
!       numail(igr,iel)=numero de la maille associee a l'element iel
#define numail(ca_igr_,ca_iel_) zi(ca_ialiel_-1+zi(ca_illiel_+ca_igr_-1)+ca_iel_-1)
!       numglm(ima,ino)=numero global du noeud ino de la maille ima
!                     ima etant une maille du maillage.
#define numglm(ima,ino) zi(ca_iamaco_-1+zi(ca_ilmaco_+ima-1)+ino-1)
!       numgls(ima,ino)=numero global du noeud ino de la maille ima
!                     ima etant une maille supplementaire du ligrel
#define numgls(ima,ino) zi(ca_iamsco_-1+zi(ca_ilmsco_+ima-1)+ino-1)
!
! --------------------------------------------------------------------------------------------------
!
    modloc = ca_iamloc_-1+zi(ca_ilmloc_-1+imodat)
    ityplo = zi(modloc-1+1)
    debugr = zi(ca_iawlo2_-1+5*(ca_nbgr_*(iparg-1)+ca_igr_-1)+5)
    lgcata = zi(ca_iawlo2_-1+5*(ca_nbgr_*(iparg-1)+ca_igr_-1)+2)
!
    ASSERT(ityplo .lt. 4)
!
!   1-  cas: chno -> elga :
!   -----------------------
!   le cas ityplo=3 n est pas prevu : developpement a faire ...
    ASSERT(ityplo .ne. 3)
!
!   2-  cas: chno -> elno :
!       cas: chno -> elem (moyenne)
!   --------------------------------
    if ((ityplo .eq. 2) .or. (ityplo .eq. 1)) then
        if (ityplo .eq. 2) then
            moyenn = .false.
        else
            moyenn = .true.
        end if
!       2.1 on cherche nno sur le 1er element :
!       ---------------------------------------
        ima = numail(ca_igr_, 1)
        ASSERT(ima .ne. 0)
        if (ima .gt. 0) then
            nno = zi(ca_ilmaco_-1+ima+1)-zi(ca_ilmaco_-1+ima)
        else
            nno = zi(ca_ilmsco_-1-ima+1)-zi(ca_ilmsco_-1-ima)-1
        end if
!       2.2 on recupere le debut du descripteur grandeur :
!       --------------------------------------------------
        nbpt = zi(modloc-1+4)
        nbpt2 = mod(nbpt, 10000)
        if (nbpt .ne. nbpt2) then
            diff = .true.
        else
            diff = .false.
            idg2 = 5
        end if
!       moyenn => (nbpt2=1)
        ASSERT((.not. moyenn) .or. (nbpt2 .eq. 1))
!       .not.moyenn => (nbpt2=nno)
        ASSERT(moyenn .or. (nbpt2 .eq. nno))
!       2.3 si moyenn, il faut mettre a zero le champ local
!           (pour pouvoir cumuler)
!       --------------------------------------------------
        if (moyenn) then
            ncmp = lgcata
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
        end if
!
!           -- c'est 1 champ avec profil_noeud:
!           ------------------------------------
        prno1 = zi(ca_iachii_-1+ca_iachid_*(ca_iichin_-1)+8)
        prno2 = zi(ca_iachii_-1+ca_iachid_*(ca_iichin_-1)+9)
        deb2 = debugr
        do ca_iel_ = 1, ca_nbelgr_
            if (ca_lparal_) then
                if (.not. ca_paral_(ca_iel_)) then
                    deb2 = deb2+lgcata
                    goto 110
                end if
            end if
            ima = numail(ca_igr_, ca_iel_)
            ASSERT(ima .ne. 0)
            do ino = 1, nno
                if (diff) idg2 = 5+ca_nec_*(ino-1)
                if (ima .gt. 0) then
                    nugl = numglm(ima, ino)
                else
                    nugl = numgls((-ima), ino)
                end if
                deb1 = (abs(nugl)-1)*(ca_nec_+2)+1
                idg1 = (abs(nugl)-1)*(ca_nec_+2)+3
!
                if (nugl .gt. 0) then
                    call trigd(zi(prno1-1+idg1), zi(prno1-1+deb1), zi(modloc-1+idg2), deb2, &
                               moyenn, ino, nno)
                else
                    call trigd(zi(prno2-1+idg1), zi(prno2-1+deb1), zi(modloc-1+idg2), deb2, &
                               moyenn, ino, nno)
                end if
            end do
110         continue
        end do
!
        if (moyenn) then
            ncmp = lgcata
            if (ca_typegd_ .eq. 'R') then
                if (ca_lparal_) then
                    do ca_iel_ = 1, ca_nbelgr_
                        if (ca_paral_(ca_iel_)) then
                            iaux1 = ca_iachlo_+debugr-1+(ca_iel_-1)*ncmp
                            do k = 1, ncmp
                                zr(iaux1-1+k) = zr(iaux1-1+k)/dble(nno)
                            end do
                        end if
                    end do
                else
                    do k = 1, ca_nbelgr_*ncmp
                        zr(ca_iachlo_-1+k) = zr(ca_iachlo_+debugr-1-1+k)/dble(nno)
                    end do
                end if
            else if (ca_typegd_ .eq. 'C') then
                if (ca_lparal_) then
                    do ca_iel_ = 1, ca_nbelgr_
                        if (ca_paral_(ca_iel_)) then
                            iaux1 = ca_iachlo_+debugr-1+(ca_iel_-1)*ncmp
                            do k = 1, ncmp
                                zc(iaux1-1+k) = zc(iaux1-1+k)/dble(nno)
                            end do
                        end if
                    end do
                else
                    do k = 1, ca_nbelgr_*ncmp
                        zc(ca_iachlo_-1+k) = zc(ca_iachlo_+debugr-1-1+k)/dble(nno)
                    end do
                end if
            else
                ASSERT(.false.)
            end if
        end if
!
    end if
!
!
end subroutine
