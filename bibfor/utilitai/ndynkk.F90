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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine ndynkk(sddyna, chaine, nomsd)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/ndynlo.h"
!
    character(len=19) :: sddyna
    character(len=*) :: chaine
    character(len=19) :: nomsd
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Dynamic management
!
! Get parameters (string) for dynamic datastructure
!
! --------------------------------------------------------------------------------------------------
!
! In  sddyna           : datastructure for dynamic
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24) :: nosd, tcha
    integer(kind=8) :: jnosd, jtcha
    character(len=24) :: vecent, vecabs
    integer(kind=8) :: jvecen, jvecab
    character(len=24) :: veol, vaol
    integer(kind=8) :: jveol, jvaol
    aster_logical :: ldyna
    character(len=24) :: cham24
    character(len=15) :: sdmuap, sdprmo, sdexso
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATIONS
!
    ldyna = ndynlo(sddyna, 'DYNAMIQUE')
    nomsd = ' '
!
! --- ACCES AUX OBJETS DE LA SD SDDYNA
!
    if (ldyna) then
        nosd = sddyna(1:15)//'.NOM_SD'
        veol = sddyna(1:15)//'.VEEL_OLD'
        vaol = sddyna(1:15)//'.VEAS_OLD'
        tcha = sddyna(1:15)//'.TYPE_CHA'
        vecent = sddyna(1:15)//'.VECENT'
        vecabs = sddyna(1:15)//'.VECABS'
        call jeveuo(vecent, 'L', jvecen)
        call jeveuo(vecabs, 'L', jvecab)
        call jeveuo(nosd, 'L', jnosd)
        call jeveuo(veol, 'E', jveol)
        call jeveuo(vaol, 'E', jvaol)
        call jeveuo(tcha, 'L', jtcha)
    else
        goto 999
    end if
!
    if (chaine(1:6) .eq. 'DEPENT') then
        cham24 = zk24(jvecen+1-1)
    else if (chaine(1:6) .eq. 'VITENT') then
        cham24 = zk24(jvecen+2-1)
    else if (chaine(1:6) .eq. 'ACCENT') then
        cham24 = zk24(jvecen+3-1)
    else if (chaine(1:6) .eq. 'DEPABS') then
        cham24 = zk24(jvecab+1-1)
    else if (chaine(1:6) .eq. 'VITABS') then
        cham24 = zk24(jvecab+2-1)
    else if (chaine(1:6) .eq. 'ACCABS') then
        cham24 = zk24(jvecab+3-1)
!
    else if (chaine(1:6) .eq. 'STADYN') then
        cham24 = zk24(jnosd-1+4)
!
    else if (chaine(1:11) .eq. 'OLDP_VEFEDO') then
        cham24 = zk24(jveol+1-1)
    else if (chaine(1:11) .eq. 'OLDP_VEFSDO') then
        cham24 = zk24(jveol+2-1)
    else if (chaine(1:11) .eq. 'OLDP_VEDIDO') then
        cham24 = zk24(jveol+3-1)
    else if (chaine(1:11) .eq. 'OLDP_VEDIDI') then
        cham24 = zk24(jveol+4-1)
    else if (chaine(1:11) .eq. 'OLDP_VEFINT') then
        cham24 = zk24(jveol+5-1)
    else if (chaine(1:11) .eq. 'OLDP_VEONDP') then
        cham24 = zk24(jveol+6-1)
    else if (chaine(1:11) .eq. 'OLDP_VESSTF') then
        cham24 = zk24(jveol+8-1)
!
    else if (chaine(1:11) .eq. 'OLDP_CNFEDO') then
        cham24 = zk24(jvaol+1-1)
    else if (chaine(1:11) .eq. 'OLDP_CNFSDO') then
        cham24 = zk24(jvaol+2-1)
    else if (chaine(1:11) .eq. 'OLDP_CNDIDO') then
        cham24 = zk24(jvaol+3-1)
    else if (chaine(1:11) .eq. 'OLDP_CNDIDI') then
        cham24 = zk24(jvaol+4-1)
    else if (chaine(1:11) .eq. 'OLDP_CNFINT') then
        cham24 = zk24(jvaol+5-1)
    else if (chaine(1:11) .eq. 'OLDP_CNONDP') then
        cham24 = zk24(jvaol+6-1)
    else if (chaine(1:11) .eq. 'OLDP_CNLAPL') then
        cham24 = zk24(jvaol+7-1)
    else if (chaine(1:11) .eq. 'OLDP_CNSSTF') then
        cham24 = zk24(jvaol+8-1)
    else if (chaine(1:11) .eq. 'OLDP_CNCINE') then
        cham24 = zk24(jvaol+9-1)
    else if (chaine(1:11) .eq. 'OLDP_CNVISS') then
        cham24 = zk24(jvaol+10-1)
    else if (chaine(1:11) .eq. 'OLDP_CNHYST') then
        cham24 = zk24(jvaol+11-1)
    else if (chaine(1:11) .eq. 'OLDP_CNSSTR') then
        cham24 = zk24(jvaol+12-1)
    else if (chaine(1:11) .eq. 'OLDP_CNELTC') then
        cham24 = zk24(jvaol+13-1)
    else if (chaine(1:11) .eq. 'OLDP_CNELTF') then
        cham24 = zk24(jvaol+14-1)
!
    else if (chaine(1:6) .eq. 'CHONDP') then
        cham24 = zk24(jtcha+1-1)
!
! --- SD EXCIT_SOL
!
    else if (chaine(1:6) .eq. 'SDEXSO') then
        sdexso = zk24(jnosd-1+5) (1:15)
        cham24 = sdexso
! - Static modes for multi-support
    else if (chaine .eq. 'multSuppMode') then
        cham24 = zk24(jnosd-1+6)
!
! --- SD PROJECTION MODALE
!
    else if (chaine(1:11) .eq. 'PRMO_DEPGEM') then
        sdprmo = zk24(jnosd-1+3) (1:15)
        cham24 = sdprmo(1:15)//'.DGM'
    else if (chaine(1:11) .eq. 'PRMO_DEPGEP') then
        sdprmo = zk24(jnosd-1+3) (1:15)
        cham24 = sdprmo(1:15)//'.DGP'
    else if (chaine(1:11) .eq. 'PRMO_VITGEM') then
        sdprmo = zk24(jnosd-1+3) (1:15)
        cham24 = sdprmo(1:15)//'.VGM'
    else if (chaine(1:11) .eq. 'PRMO_VITGEP') then
        sdprmo = zk24(jnosd-1+3) (1:15)
        cham24 = sdprmo(1:15)//'.VGP'
    else if (chaine(1:11) .eq. 'PRMO_ACCGEM') then
        sdprmo = zk24(jnosd-1+3) (1:15)
        cham24 = sdprmo(1:15)//'.AGM'
    else if (chaine(1:11) .eq. 'PRMO_ACCGEP') then
        sdprmo = zk24(jnosd-1+3) (1:15)
        cham24 = sdprmo(1:15)//'.AGP'
    else if (chaine(1:11) .eq. 'PRMO_BASMOD') then
        sdprmo = zk24(jnosd-1+3) (1:15)
        cham24 = sdprmo(1:15)//'.BAM'
    else if (chaine(1:11) .eq. 'PRMO_MASGEN') then
        sdprmo = zk24(jnosd-1+3) (1:15)
        cham24 = sdprmo(1:15)//'.MAG'
    else if (chaine(1:11) .eq. 'PRMO_RIGGEN') then
        sdprmo = zk24(jnosd-1+3) (1:15)
        cham24 = sdprmo(1:15)//'.RIG'
    else if (chaine(1:11) .eq. 'PRMO_AMOGEN') then
        sdprmo = zk24(jnosd-1+3) (1:15)
        cham24 = sdprmo(1:15)//'.AMOG'
    else if (chaine(1:11) .eq. 'PRMO_FONGEN') then
        sdprmo = zk24(jnosd-1+3) (1:15)
        cham24 = sdprmo(1:15)//'.FOG'
    else if (chaine(1:11) .eq. 'PRMO_FORGEN') then
        sdprmo = zk24(jnosd-1+3) (1:15)
        cham24 = sdprmo(1:15)//'.FRG'
    else if (chaine(1:11) .eq. 'PRMO_ACCGCN') then
        sdprmo = zk24(jnosd-1+3) (1:15)
        cham24 = sdprmo(1:15)//'.AGN'
    else if (chaine(1:11) .eq. 'PRMO_VALFON') then
        sdprmo = zk24(jnosd-1+3) (1:15)
        cham24 = sdprmo(1:15)//'.VAF'
    else if (chaine(1:11) .eq. 'PRMO_FMODAL') then
        sdprmo = zk24(jnosd-1+3) (1:15)
        cham24 = sdprmo(1:15)//'.FMD'
!
! --- SD MULTI_APPUi
!
    else if (chaine(1:11) .eq. 'MUAP_MAFDEP') then
        sdmuap = zk24(jnosd-1+1) (1:15)
        cham24 = sdmuap(1:15)//'.FDP'
    else if (chaine(1:11) .eq. 'MUAP_MAFVIT') then
        sdmuap = zk24(jnosd-1+1) (1:15)
        cham24 = sdmuap(1:15)//'.FVT'
    else if (chaine(1:11) .eq. 'MUAP_MAFACC') then
        sdmuap = zk24(jnosd-1+1) (1:15)
        cham24 = sdmuap(1:15)//'.FAC'
    else if (chaine(1:11) .eq. 'MUAP_MAMULA') then
        sdmuap = zk24(jnosd-1+1) (1:15)
        cham24 = sdmuap(1:15)//'.MUA'
    else if (chaine(1:11) .eq. 'MUAP_MAPSID') then
        sdmuap = zk24(jnosd-1+1) (1:15)
        cham24 = sdmuap(1:15)//'.PSD'
!
! --- Carte fluide (damping)
!
    else if (chaine(1:9) .eq. 'AMOR_FLUI') then
        cham24 = zk24(jnosd-1+8)
    else
        ASSERT(.false.)
    end if
!
    nomsd = cham24(1:19)
!
999 continue
!
    call jedema()
!
end subroutine
