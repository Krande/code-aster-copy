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

subroutine xfem_store_pc(matass, base, nonu, neq, deeq, &
                         nbnoxfem, nbnomax, ino_xfem, ieq_loc, neq_mloc, &
                         maxi_ddl, iglob_ddl, deca, tab_mloc, pc, kstruct)
!
!-----------------------------------------------------------------------
! BUT : CREATION D UNE MATR_ASSE A PARTIR D UNE MATRICE BLOC <DENSE>
!-----------------------------------------------------------------------
!
! ARGUMENTS :
!------------
!
!
!  SORTIE :
!     - PC : NOM DE MATR_ASSE
!-----------------------------------------------------------------------
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jexnum.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mtdefs.h"
#include "asterfort/wkvect.h"
!-----------------------------------------------------------------------
!
    character(len=19) :: matass, pc
    character(len=14) :: nonu
    character(len=1) :: base
    character(len=5) :: kstruct
    integer(kind=8) :: neq, nbnoxfem, maxi_ddl, deca, nbnomax
    integer(kind=8) :: ino_xfem(nbnomax), deeq(*)
    integer(kind=8) :: ieq_loc(neq), neq_mloc(nbnoxfem), iglob_ddl(maxi_ddl*nbnoxfem)
    real(kind=8) :: tab_mloc(deca*nbnoxfem)
!
!-----------------------------------------------------------------------
    character(len=1) :: kbid
    character(len=14) :: nu_pc
    character(len=24), pointer :: refa_pc(:) => null()
    integer(kind=8) :: cumul_ilig, jcoll, jvale_sup, jvale_inf, nunoj
    integer(kind=8) :: ipos, decaj, jsmhc_pc, jsmdi_pc, jadr, iexi, nvale
!-----------------------------------------------------------------------
!
    call jemarq()
!
    ASSERT(pc(1:11) .eq. pc)
    nu_pc = pc(1:11)//'_NU'
!
!    - PREMIERE PASSE : CALCUL DE L ESPACE <DENSE> NECESSAIRE
!                          POUR STOCKER LES MATRICES LOCALES TRIANFULAIRES SUPERIEURES
!
    if (kstruct .eq. 'D_P_B') then
        cumul_ilig = 0
        do jcoll = 1, neq
            if (ieq_loc(jcoll) .eq. 0) then
                cumul_ilig = cumul_ilig+1
            else
                cumul_ilig = cumul_ilig+ieq_loc(jcoll)
            end if
        end do
        nvale = 2
    elseif (kstruct .eq. 'DIAGO') then
        cumul_ilig = neq
        nvale = 1
    else
        ASSERT(.false.)
    end if
!
!    - DEUXIEME PASSE : ALLOCATION ET ECRITURE
!
    call jeexin(pc//'.REFA', iexi)
    if (iexi .gt. 0) call detrsd('MATR_ASSE', pc)
    call mtdefs(pc, matass, base, ' ')
    call jedetr(pc//'.VALM')
    call jecrec(pc//'.VALM', base//' V R', 'NU', 'DISPERSE', 'CONSTANT', nvale)
! TRIANGULAIRE SUPERIEURE
    call jecroc(jexnum(pc//'.VALM', 1))
    call jeecra(pc//'.VALM', 'LONMAX', cumul_ilig, kbid)
    call jeveuo(jexnum(pc//'.VALM', 1), 'E', jvale_sup)
! TRIANGULAIRE INFERIEURE
    if (nvale .eq. 2) then
        call jecroc(jexnum(pc//'.VALM', 2))
        call jeveuo(jexnum(pc//'.VALM', 2), 'E', jvale_inf)
    end if
!
    call jeveuo(pc//'.REFA', 'E', vk24=refa_pc)
!
    call jeexin(nu_pc, iexi)
    if (iexi .gt. 0) call detrsd('NUME_DDL', nu_pc)
    call copisd('NUME_DDL', base, nonu, nu_pc)
    refa_pc(2) = nu_pc
    call jedetr(nu_pc//'.SMOS.SMDI')
    call wkvect(nu_pc//'.SMOS.SMDI', base//' V I', neq, jsmdi_pc)
    jadr = 0
    do jcoll = 1, neq
        if (ieq_loc(jcoll) .eq. 0 .or. kstruct .eq. 'DIAGO') then
            jadr = jadr+1
        else
            jadr = jadr+ieq_loc(jcoll)
        end if
        zi(jsmdi_pc-1+jcoll) = jadr
    end do
    ASSERT(jadr .eq. cumul_ilig)
!   - ALLOCATION DU .SMOS.SMHC
    call jedetr(nu_pc//'.SMOS.SMHC')
    call wkvect(nu_pc//'.SMOS.SMHC', base//' V S', cumul_ilig, jsmhc_pc)
!
    decaj = 0
    do jcoll = 1, neq
        if (ieq_loc(jcoll) .eq. 0) then
            decaj = decaj+1
            zi4(jsmhc_pc-1+decaj) = int(jcoll, 4)
            zr(jvale_sup-1+decaj) = 1.d0
            if (nvale .eq. 2) zr(jvale_inf-1+decaj) = 1.d0
        else
!
            nunoj = ino_xfem(deeq(2*(jcoll-1)+1))
            if (kstruct .eq. 'D_P_B') then
                do ipos = 1, ieq_loc(jcoll)
                    decaj = decaj+1
                    zi4(jsmhc_pc-1+decaj) = int(iglob_ddl(maxi_ddl*(nunoj-1)+ipos), 4)
                    zr(jvale_sup-1+decaj) = tab_mloc(deca*(nunoj-1)+neq_mloc(nunoj)* &
                                                     (ipos-1)+ieq_loc(jcoll))
                    zr(jvale_inf-1+decaj) = tab_mloc(deca*(nunoj-1)+neq_mloc(nunoj)* &
                                                     (ieq_loc(jcoll)-1)+ipos)
                end do
            else
                decaj = decaj+1
                zi4(jsmhc_pc-1+decaj) = int(jcoll, 4)
                zr(jvale_sup-1+decaj) = tab_mloc(deca*(nunoj-1)+ieq_loc(jcoll))
            end if
!
        end if
    end do
!
    ASSERT(decaj .eq. cumul_ilig)
!
    call jedema()
!
end subroutine
