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

subroutine vrcref(modele, chmat, carele, chvref, basez)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/alchml.h"
#include "asterfort/assert.h"
#include "asterfort/carces.h"
#include "asterfort/celces.h"
#include "asterfort/cescel.h"
#include "asterfort/cesexi.h"
#include "asterfort/cesvar.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/juvinn.h"
#include "asterfort/utmess.h"
    character(len=8) :: modele, chmat, carele
    character(len=19) :: chvref
    character(len=1), intent(in), optional :: basez
! ======================================================================
!   BUT : FABRIQUER LE CHAMP DE VARIABLES DE COMMANDE DE "REFERENCE"
!   ARGUMENTS :
!   MODELE (K8)  IN/JXIN : SD MODELE
!   CHMAT  (K8)  IN/JXIN : SD CHAM_MATER
!   CARELE (K8)  IN/JXIN : SD CARA_ELEM (SOUS-POINTS)
!   CHVREF (K19) IN/JXOUT: SD CHAM_ELEM/ELGA (VARC DE "REFERENCE")
! ----------------------------------------------------------------------
!
!
    character(len=8) :: models, chmats, carels
    character(len=19) :: chvres
    integer(kind=8) :: n1, iad, isp, ipt
    integer(kind=8) :: k, k2, nbma, ncmp, icmp, jcesl1, jcesd1
    integer(kind=8) :: jcesd, jcesl, ima, nbpt, nbsp, nbcvrc, ibid
    integer(kind=8) :: jdcld, jdcll, nncp, iret
    character(len=8) :: varc, noma1, noma2
    character(len=19) :: dceli, celmod, cart1, ces1, ligrmo, csvref
    character(len=24) :: valk(4)
    real(kind=8) :: valref
    aster_logical :: avrc
    real(kind=8), pointer :: cesv1(:) => null()
    real(kind=8), pointer :: cesv(:) => null()
    integer(kind=8), pointer :: dclv(:) => null()
    character(len=8), pointer :: cvrc(:) => null()
    character(len=8), pointer :: cvvar(:) => null()
    character(len=1) :: base
    save models, chmats, carels, chvres
    data models/' '/, chmats/' '/, carels/' '/, chvres/' '/
! ----------------------------------------------------------------------
!
    call jemarq()
!
    if (present(basez)) then
        base = basez
    else
        base = 'V'
    end if
!
!     -- SI LE CHAMP A PRODUIRE EXISTE DEJA ET QUE LES ARGUMENTS
!        SONT LES MEMES QUE LA FOIS PRECEDENTE, ON SORT RAPIDEMENT :
    call exisd('CHAM_ELEM', chvref, iret)
    if (iret .gt. 0) then
        if (modele .ne. models) goto 5
        if (chmat .ne. chmats) goto 5
        if (carele .ne. carels) goto 5
        if (chvref .ne. chvres) goto 5
        goto 999
    end if
5   continue
!     -- SAUVERGARDE DES ARGUMENTS POUR LE PROCHAIN APPEL
    models = modele
    chmats = chmat
    carels = carele
    chvres = chvref
!
!
!     -- ON SE PROTEGE DES MODELES QUI NE CONNAISSENT PAS LES VRC :
    celmod = '&&VRCREF.CELMOD'
    ligrmo = modele//'.MODELE'
    call jeexin(ligrmo//'.LIEL', iret)
    if (iret .eq. 0) goto 999
    call alchml(ligrmo, 'INIT_VARC', 'PVARCPR', 'V', celmod, &
                iret, ' ')
    call detrsd('CHAMP', celmod)
    if (iret .eq. 1) goto 999
!
!
    call jeexin(chmat//'.CVRCVARC', iret)
!     AVRC : .TRUE. SI AFFE_MATERIAU/AFFE_VARC EST UTILISE
    avrc = .false.
    if (iret .gt. 0) then
        call jeveuo(chmat//'.CVRCVARC', 'L', vk8=cvrc)
        call jelira(chmat//'.CVRCVARC', 'LONMAX', nbcvrc)
        do k = 1, nbcvrc
            if (cvrc(k) .ne. ' ') then
                avrc = .true.
                goto 11
            end if
        end do
!
11      continue
    end if
!
!
!
!     -- CAS : PAS DE AFFE_VARC :
!     ------------------------------------------
    if (.not. avrc) goto 999
!
!     -- CAS AFFE_VARC  :
!     ------------------------
    call jeveuo(chmat//'.CVRCVARC', 'L', vk8=cvvar)
    call jelira(chmat//'.CVRCVARC', 'LONMAX', nbcvrc)
!
!
!     0. VERIFICATION DE LA COHERENCE DE MODELE  ET CHMAT (FICHE 19507)
!     ------------------------------------------------------------------
    call dismoi('NOM_MAILLA', modele, 'MODELE', repk=noma1)
    call dismoi('NOM_MAILLA', chmat//'.CHAMP_MAT', 'CHAMP', repk=noma2)
    if (noma1 .ne. noma2) then
        valk(1) = noma1
        valk(2) = noma2
        call utmess('F', 'CALCULEL4_23', nk=2, valk=valk)
    end if
!
!
!     1. ALLOCATION DE CSVREF
!     ------------------------------------------
    dceli = '&&VRCREF.DCELI'
    csvref = '&&VRCREF.CSVREF'
    call cesvar(carele, ' ', ligrmo, dceli)
!
!     -- MODIFICATION DE DCELI : TOUTES LES MAILLES ONT
!        NBCVRC COMPOSANTES.
    call jeveuo(dceli//'.CESD', 'L', jdcld)
    call jeveuo(dceli//'.CESL', 'L', jdcll)
    call jeveuo(dceli//'.CESV', 'E', vi=dclv)
    nbma = zi(jdcld-1+1)
!
    do ima = 1, nbma
        nbpt = zi(jdcld-1+5+4*(ima-1)+1)
        nbsp = max(1, zi(jdcld-1+5+4*(ima-1)+2))
        ASSERT(nbpt .eq. 1)
        ASSERT(nbsp .eq. 1)
        call cesexi('C', jdcld, jdcll, ima, 1, &
                    1, 2, iad)
        if (iad .gt. 0) dclv(iad) = nbcvrc
    end do
!
    call alchml(ligrmo, 'INIT_VARC', 'PVARCPR', 'V', celmod, &
                iret, dceli)
    ASSERT(iret .eq. 0)
    call detrsd('CHAMP', dceli)
    call celces(celmod, 'V', csvref)
    call detrsd('CHAMP', celmod)
!
    call jelira(csvref//'.CESV', 'LONMAX', n1)
!
    call jeveuo(csvref//'.CESD', 'L', jcesd)
    call jeveuo(csvref//'.CESL', 'E', jcesl)
    call jeveuo(csvref//'.CESV', 'E', vr=cesv)
    call jelira(csvref//'.CESL', 'LONMAX', n1)
    do k = 1, n1
        zl(jcesl-1+k) = .false.
    end do
!
!
!
!     2. REMPLISSAGE DE CSVREF.CESV :
!     ------------------------------------------
    varc = ' '
    do k = 1, nbcvrc
        if (cvvar(k) .eq. varc) goto 1
        varc = cvvar(k)
        cart1 = chmat//'.'//varc//'.1'
        ces1 = '&&VRCREF.CES1'
        call carces(cart1, 'ELEM', ' ', 'V', ces1, &
                    'A', iret)
        ASSERT(iret .eq. 0)
!
        call jeveuo(ces1//'.CESD', 'L', jcesd1)
        call jeveuo(ces1//'.CESV', 'L', vr=cesv1)
        call jeveuo(ces1//'.CESL', 'L', jcesl1)
!
        nbma = zi(jcesd-1+1)
        ASSERT(nbma .eq. zi(jcesd1-1+1))
!
!       -- CALCUL DE NCMP
        ncmp = 0
        do k2 = k, nbcvrc
            if (cvvar(k2) .eq. varc) ncmp = ncmp+1
        end do
!
        do ima = 1, nbma
            nbpt = zi(jcesd-1+5+4*(ima-1)+1)
            nbsp = max(1, zi(jcesd-1+5+4*(ima-1)+2))
!
            call cesexi('C', jcesd1, jcesl1, ima, 1, &
                        1, 1, iad)
            if (iad .le. 0) goto 70
            valref = cesv1(iad)
!
            do ipt = 1, nbpt
                do isp = 1, nbsp
                    do icmp = 1, ncmp
                        call cesexi('C', jcesd, jcesl, ima, ipt, &
                                    isp, k-1+icmp, iad)
                        ASSERT(iad .le. 0)
                        if (iad .eq. 0) goto 51
                        iad = -iad
                        zl(jcesl-1+iad) = .true.
                        cesv(iad) = valref
51                      continue
                    end do
                end do
            end do
70          continue
        end do
        call detrsd('CHAMP', ces1)
1       continue
    end do
!
!
!     3. RECOPIE DU CHAMP SIMPLE DANS LE CHAMP CHVREF
!     -----------------------------------------------------
!
!     LE CHAMP DE TEMP_REF PEUT CONTENIR DES VALEURS "VIDES"
!     ON LES TRANSFORME EN "NAN"
    call juvinn(csvref//'.CESV')
    call cescel(csvref, ligrmo, 'INIT_VARC', 'PVARCPR', 'NAN', &
                nncp, base, chvref, 'F', ibid)
    call detrsd('CHAM_ELEM_S', csvref)
!
999 continue
    call jedema()
end subroutine
