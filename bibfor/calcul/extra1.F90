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
subroutine extra1(nin, lchin, lpain)
!
    use calcul_module, only: ca_iachii_, ca_iachik_, ca_iachin_, ca_iachix_, ca_iachlo_, &
                             ca_ianueq_, ca_iaoppa_, ca_iawlo2_, ca_iawloc_, ca_igd_, ca_igr_, &
                             ca_iichin_, ca_ilchlo_, ca_itypgd_, ca_lprno_, ca_nbgr_, ca_ncmpmx_, &
                             ca_nec_, ca_npario_, ca_typegd_, ca_nute_, ca_nuop_, ca_iachid_
!
    implicit none
!
! person_in_charge: jacques.pellet at edf.fr
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/assert.h"
#include "asterfort/excart.h"
#include "asterfort/exchml.h"
#include "asterfort/exchgo.h"
#include "asterfort/exchno.h"
#include "asterfort/exresl.h"
#include "asterfort/nbpara.h"
#include "asterfort/nopara.h"
#include "asterfort/utmess.h"
!
    integer(kind=8) :: nin
    character(len=*) :: lchin(*)
    character(len=8) :: lpain(*)
!-----------------------------------------------------------------------
!     but: preparer les champs locaux "in"
!-----------------------------------------------------------------------
    integer(kind=8) :: debugr, lggrel
    character(len=19) :: chin
    character(len=4) :: type
    character(len=8) :: nompar
    integer(kind=8) :: k, iparg, imodat
    integer(kind=8) :: ipar, npin, iparin
    aster_logical :: exich
!-------------------------------------------------------------------
!
    npin = nbpara(ca_nuop_, ca_nute_, 'IN ')
    do ipar = 1, npin
        nompar = nopara(ca_nuop_, ca_nute_, 'IN ', ipar)
        iparg = indik8(zk8(ca_iaoppa_), nompar, 1, ca_npario_)
        iparin = indik8(lpain, nompar, 1, nin)
        exich = ((iparin .gt. 0) .and. zl(ca_iachix_-1+iparin))
        if (.not. exich) then
            zi(ca_iawloc_-1+3*(iparg-1)+1) = -1
            zi(ca_iawlo2_-1+5*(ca_nbgr_*(iparg-1)+ca_igr_-1)+2) = 0
            goto 90
        end if
!
        ASSERT(iparin .ne. 0)
        chin = lchin(iparin)
        if (chin(1:1) .eq. ' ') then
            call utmess('E', 'CALCUL_13', sk=nompar)
        end if
!
!
        ca_iichin_ = iparin
        ca_igd_ = zi(ca_iachii_-1+ca_iachid_*(iparin-1)+1)
        ca_nec_ = zi(ca_iachii_-1+ca_iachid_*(iparin-1)+2)
        ca_ncmpmx_ = zi(ca_iachii_-1+ca_iachid_*(iparin-1)+3)
        ca_iachin_ = zi(ca_iachii_-1+ca_iachid_*(iparin-1)+5)
        ca_ianueq_ = zi(ca_iachii_-1+ca_iachid_*(iparin-1)+10)
        ca_lprno_ = zi(ca_iachii_-1+ca_iachid_*(iparin-1)+11)
        iparg = indik8(zk8(ca_iaoppa_), nompar, 1, ca_npario_)
        ca_iachlo_ = zi(ca_iawloc_-1+3*(iparg-1)+1)
        ca_ilchlo_ = zi(ca_iawloc_-1+3*(iparg-1)+2)
        imodat = zi(ca_iawlo2_-1+5*(ca_nbgr_*(iparg-1)+ca_igr_-1)+1)
        ASSERT((ca_iachlo_ .lt. -2) .or. (ca_iachlo_ .gt. 0))
        ASSERT(ca_ilchlo_ .ne. -1)
        type = zk8(ca_iachik_-1+2*(iparin-1)+1) (1:4)
        ca_typegd_ = zk8(ca_iachik_-1+2*(iparin-1)+2)
        if (ca_typegd_ .eq. 'R') then
            ca_itypgd_ = 1
        else if (ca_typegd_ .eq. 'C') then
            ca_itypgd_ = 2
        else if (ca_typegd_ .eq. 'I') then
            ca_itypgd_ = 3
        else if (ca_typegd_ .eq. 'K8') then
            ca_itypgd_ = 4
        else if (ca_typegd_ .eq. 'K16') then
            ca_itypgd_ = 5
        else if (ca_typegd_ .eq. 'K24') then
            ca_itypgd_ = 6
        else
            ASSERT(.false.)
        end if
!
!
!       1- mise a .false. du champ_loc.EXIS
!       -----------------------------------------------------
        lggrel = zi(ca_iawlo2_-1+5*(ca_nbgr_*(iparg-1)+ca_igr_-1)+4)
        debugr = zi(ca_iawlo2_-1+5*(ca_nbgr_*(iparg-1)+ca_igr_-1)+5)
        do k = 1, lggrel
            zl(ca_ilchlo_-1+debugr-1+k) = .false.
        end do
!
!
!       2- on fait l'extraction:
!       -------------------------------------------
        if (type .eq. 'CART') call excart(imodat, iparg)
        if (type .eq. 'CHML') call exchml(imodat, iparg)
        if (type .eq. 'CHNO') call exchno(imodat, iparg)
        if (type .eq. 'CHGO') call exchgo(imodat, iparg)
        if (type .eq. 'RESL') call exresl(imodat, iparg, chin)
90      continue
    end do
!
!
!
end subroutine
