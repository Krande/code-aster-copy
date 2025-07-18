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

subroutine vpini2(eigsol, lcomod, nbvecg, nfreqg, nbpark, &
                  nbpari, nbparr, vecrer, vecrei, vecrek, &
                  vecvp, mxresf)
!
! CREATION ET INITIALISATION DES SDS RESULTATS DE OP0045.
! RQ. ILS SONT DETRUITS DANS VPPOST.
! -------------------------------------------------------------------------------------------------
! person_in_charge: olivier.boiteau at edf.fr
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/isnnem.h"
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/vecini.h"
#include "asterfort/vecink.h"
#include "asterfort/vecint.h"
#include "asterfort/vplecs.h"
#include "asterfort/wkvect.h"
!
!
! --- INPUT
!
    integer(kind=8), intent(in) :: nbvecg, nfreqg, nbpark, nbpari, nbparr
    aster_logical, intent(in) :: lcomod
    character(len=19), intent(in) :: eigsol
    character(len=24), intent(in) :: vecrer, vecrei, vecrek, vecvp
!
! --- OUTPUT
!
    integer(kind=8), intent(out) :: mxresf
!
! --- INPUT/OUTPUT
! None
!
! --- VARIABLES LOCALES
!
    integer(kind=8) :: nbvect, nfreq, iauxr, iauxi, iauxk, lresui, lresur, lresuk
    integer(kind=8) :: lraide, neq, indf, lvec
    real(kind=8) :: undf
    character(len=19) :: raide
    character(len=24) :: kzero
    aster_logical :: lc, lkr, lns
!
! -----------------------
! --- CORPS DE LA ROUTINE
! -----------------------
!
!
! --  INITS.
    call jemarq()
    undf = r8vide()
    indf = isnnem()
    kzero = ' '
!
! --  LECTURE DES PARAMETRES MODAUX
    call vplecs(eigsol, nbvect_=nbvect, nfreq_=nfreq, &
                raide_=raide, lc_=lc, lkr_=lkr, lns_=lns)
    call jeveuo(raide//'.&INT', 'E', lraide)
    neq = zi(lraide+2)
!
! --  CREATION DES SDS
    if (lcomod) then
        if (lc .or. lns .or. (.not. lkr)) then
            ASSERT(.false.)
        end if
        mxresf = nfreqg
        iauxr = nbparr*nbvecg
        iauxi = nbpari*nbvecg
        iauxk = nbpark*nbvecg
    else
        mxresf = nfreq
        iauxr = nbparr*nbvect
        iauxi = nbpari*nbvect
        iauxk = nbpark*nbvect
    end if
!
    call wkvect(vecrei, 'V V I', iauxi, lresui)
    call vecint(iauxi, indf, zi(lresui))
    call wkvect(vecrer, 'V V R', iauxr, lresur)
    call vecini(iauxr, undf, zr(lresur))
    call wkvect(vecrek, 'V V K24', iauxk, lresuk)
    call vecink(iauxk, kzero, zk24(lresuk))
!
!
    if (lkr .and. (.not. lc) .and. (.not. lns)) then
        if (lcomod) then
            call wkvect(vecvp, 'V V R', neq*nbvecg, lvec)
        else
            call wkvect(vecvp, 'V V R', neq*nbvect, lvec)
        end if
    else
! --  CAS GENERALISE COMPLEXE OU QUADRATIQUE REEL ET COMPLEXE ---
        if (lcomod) then
            ASSERT(.false.)
        end if
        call wkvect(vecvp, 'V V C', neq*nbvect, lvec)
    end if
!
    call jedema()
!
!     FIN DE VPINI2
!
end subroutine
