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
subroutine monte1(te2, nout, lchout, lpaout, igr2)
!
    use calcul_module, only: ca_iachoi_, ca_iadsgd_, ca_iaoppa_, ca_iawlo2_, ca_iawloc_, &
                             ca_iawtyp_, ca_nbelgr_, ca_nbgr_, ca_npario_, &
                             ca_paral_, ca_lparal_, ca_nuop_, ca_iel_
!
    implicit none
! person_in_charge: jacques.pellet at edf.fr
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/digde2.h"
#include "asterfort/grdeur.h"
#include "asterfort/jacopo.h"
#include "asterfort/jedema.h"
#include "asterfort/jelibe.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/modatt.h"
#include "asterfort/nbpara.h"
#include "asterfort/nopara.h"
!
    integer(kind=8) :: nout, te2, igr2
    character(len=19) :: ch19
    character(len=*) :: lchout(*)
    character(len=8) :: lpaout(*)
!-----------------------------------------------------------------------
!     entrees:
!     igr2   : numero du grel dont on sauve les champs locaux
!
!     sorties:
!     met a jour les champs globaux de sortie de l ca_nuop_ion ca_nuop_
!-----------------------------------------------------------------------
    integer(kind=8) :: ipar, np, mod1, jpar, gd, iaux1, iaux2, iaux0
    integer(kind=8) :: iparg, iachlo, lggrel, jcelv, jresl
    integer(kind=8) :: descgd, jceld, code, debugr, ncmpel, debgr2
    character(len=8) :: nompar, typsca
!-----------------------------------------------------------------------
!
    call jemarq()
!
!
    np = nbpara(ca_nuop_, te2, 'OUT')
    do ipar = 1, np
        nompar = nopara(ca_nuop_, te2, 'OUT', ipar)
        iparg = indik8(zk8(ca_iaoppa_), nompar, 1, ca_npario_)
        iachlo = zi(ca_iawloc_-1+3*(iparg-1)+1)
        if (iachlo .eq. -1) cycle
!
        gd = grdeur(nompar)
        descgd = ca_iadsgd_+7*(gd-1)
        code = zi(descgd-1+1)
!
        mod1 = modatt(ca_nuop_, te2, 'OUT', ipar)
        jpar = indik8(lpaout, nompar, 1, nout)
        ch19 = lchout(jpar)
!
!
        typsca = zk8(ca_iawtyp_-1+iparg)
        lggrel = zi(ca_iawlo2_-1+5*(ca_nbgr_*(iparg-1)+igr2-1)+4)
        debugr = zi(ca_iawlo2_-1+5*(ca_nbgr_*(iparg-1)+igr2-1)+5)
        if (lggrel .eq. 0) cycle
!
!
        if (code .eq. 1) then
!           -- cas : cham_elem
            jceld = zi(ca_iachoi_-1+2*(jpar-1)+1)
            debgr2 = zi(jceld-1+zi(jceld-1+4+igr2)+8)
            jcelv = zi(ca_iachoi_-1+2*(jpar-1)+2)
            call jacopo(lggrel, typsca, iachlo+debugr-1, jcelv-1+debgr2)
!
!
        else
!           -- cas : resuelem
            call jeveuo(jexnum(ch19//'.RESL', igr2), 'E', jresl)
!
            if (ca_lparal_) then
                ncmpel = digde2(mod1)
                do ca_iel_ = 1, ca_nbelgr_
                    if (ca_paral_(ca_iel_)) then
                        iaux0 = (ca_iel_-1)*ncmpel
                        iaux1 = iachlo+debugr-1+iaux0
                        iaux2 = jresl+iaux0
                        call jacopo(ncmpel, typsca, iaux1, iaux2)
                    end if
                end do
            else
                call jacopo(lggrel, typsca, iachlo+debugr-1, jresl)
            end if
!
            call jelibe(jexnum(ch19//'.RESL', igr2))
        end if
    end do
!
    call jedema()
end subroutine
