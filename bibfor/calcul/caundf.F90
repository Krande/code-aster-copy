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

subroutine caundf(code, opt, te)
    use calcul_module, only: ca_iaoppa_, ca_iawlo2_, ca_iawloc_, &
                             ca_iawtyp_, ca_igr_, ca_nbgr_, ca_npario_
    implicit none
! person_in_charge: jacques.pellet at edf.fr

#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterc/isnnem.h"
#include "asterc/r8nnem.h"
#include "asterfort/assert.h"
#include "asterfort/jenuno.h"
#include "asterfort/jexnum.h"
#include "asterfort/nbpara.h"
#include "asterfort/nopara.h"
#include "asterfort/utmess.h"

    integer(kind=8) :: opt, te
    character(len=5) :: code
!-----------------------------------------------------------------------
!     entrees:
!        code :  / 'ECRIT' : on ecrit une valeur undef au bout des chloc
!                / 'VERIF' : on verifie la valeur undef au bout des chloc
!        opt  : option
!        te   : type_element
!-----------------------------------------------------------------------
    integer(kind=8) :: innem
    integer(kind=8) :: np, ipar
    integer(kind=8) ::  iparg, lggrel, iachlo
    character(len=3) :: typsca
    character(len=8) :: nompar
    aster_logical :: arret, ecras
    character(len=16) :: nomte, nomopt
    integer(kind=8) :: ich, debugr, lgcata
    real(kind=8) :: rnnem
    character(len=8) :: knnem
    character(len=24) :: valk(3)
!-------------------------------------------------------------------

    innem = isnnem()
    rnnem = r8nnem()
    knnem = '????????'

    ASSERT((code .eq. 'ECRIT') .or. (code .eq. 'VERIF'))

    if (code .eq. 'ECRIT') then
!   ---------------------------

!        -- champs "in" et "out" :
        do iparg = 1, ca_npario_
            lgcata = zi(ca_iawlo2_-1+5*(ca_nbgr_*(iparg-1)+ca_igr_-1)+2)
            if (lgcata .le. 0) cycle
            iachlo = zi(ca_iawloc_-1+3*(iparg-1)+1)
            if ((iachlo .eq. -1) .or. (iachlo .eq. -2)) cycle

            typsca = zk8(ca_iawtyp_-1+iparg) (1:3)
            lggrel = zi(ca_iawlo2_-1+5*(ca_nbgr_*(iparg-1)+ca_igr_-1)+4)
            debugr = zi(ca_iawlo2_-1+5*(ca_nbgr_*(iparg-1)+ca_igr_-1)+5)

            if (typsca .eq. 'R') then
                zr(iachlo-1+debugr-1+lggrel+1) = rnnem
            else if (typsca .eq. 'C') then
                zc(iachlo-1+debugr-1+lggrel+1) = dcmplx(rnnem, rnnem)
            else if (typsca .eq. 'I') then
                zi(iachlo-1+debugr-1+lggrel+1) = innem
            else if (typsca .eq. 'K8') then
                zk8(iachlo-1+debugr-1+lggrel+1) = knnem
            else if (typsca .eq. 'K16') then
                zk16(iachlo-1+debugr-1+lggrel+1) = knnem
            else if (typsca .eq. 'K24') then
                zk24(iachlo-1+debugr-1+lggrel+1) = knnem
            else
                ASSERT(.false.)
            end if
        end do

    else if (code .eq. 'VERIF') then
!   ------------------------------

!        -- champs "out" :
        arret = .false.
        np = nbpara(opt, te, 'OUT')
        do ipar = 1, np
            ecras = .false.
            nompar = nopara(opt, te, 'OUT', ipar)
            iparg = indik8(zk8(ca_iaoppa_), nompar, 1, ca_npario_)
            lgcata = zi(ca_iawlo2_-1+5*(ca_nbgr_*(iparg-1)+ca_igr_-1)+2)
            if (lgcata .le. 0) cycle
            ich = zi(ca_iawloc_-1+3*(iparg-1)+3)
            if (ich .eq. 0) cycle
            iachlo = zi(ca_iawloc_-1+3*(iparg-1)+1)
            if ((iachlo .eq. -1) .or. (iachlo .eq. -2)) cycle

            typsca = zk8(ca_iawtyp_-1+iparg) (1:3)
            lggrel = zi(ca_iawlo2_-1+5*(ca_nbgr_*(iparg-1)+ca_igr_-1)+4)
            debugr = zi(ca_iawlo2_-1+5*(ca_nbgr_*(iparg-1)+ca_igr_-1)+5)

            if (typsca .eq. 'R') then
                if (.not. isnan(zr(iachlo-1+debugr-1+lggrel+1))) ecras = .true.
            else if (typsca .eq. 'C') then
                if (.not. isnan(dble(zc(iachlo-1+debugr-1+lggrel+1)))) ecras = .true.
                if (.not. isnan(dimag(zc(iachlo-1+debugr-1+lggrel+1)))) ecras = .true.
            else if (typsca .eq. 'I') then
                if (zi(iachlo-1+debugr-1+lggrel+1) .ne. innem) ecras = .true.
            else
                ASSERT(.false.)
            end if

            if (ecras) then
                arret = .true.
                call jenuno(jexnum('&CATA.TE.NOMTE', te), nomte)
                call jenuno(jexnum('&CATA.OP.NOMOPT', opt), nomopt)
                valk(1) = nomte
                valk(2) = nomopt
                valk(3) = nompar
                call utmess('E', 'CALCUL_40', nk=3, valk=valk)
            end if

        end do

        ASSERT(.not. arret)

    end if

end subroutine
