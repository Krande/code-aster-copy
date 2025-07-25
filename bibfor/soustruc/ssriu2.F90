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

subroutine ssriu2(nomu)
    implicit none
!
!     ARGUMENTS:
!     ----------
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/crmeri.h"
#include "asterfort/dismoi.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/mtdsc2.h"
#include "asterfort/mtdscr.h"
#include "asterfort/rldlr8.h"
#include "asterfort/tldlg3.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=8) :: nomu
! ----------------------------------------------------------------------
!     BUT:
!        1)FACTORISER PARTIELLEMENT LA MATR_ASSE DE RIGIDITE
!             "K_II**(-1)"
!        2)CALCULER PHI_IE=  (K_II**(-1))*K_IE
!     ===>  ATTENTION LE PHI_IE CALCULE EST L'OPPOSE DE CELUI D'IMBERT
!        3)CALCULER KP_EE = K_EE - K_EI*PHI_IE
!
!
!     IN: NOMU   : NOM DU MACR_ELEM_STAT
!
!     OUT:  PHI_IE, KP_EE,
!           "K_II**(-1)" (DANS LE DEBUT DE .RIGIMECA.UALF)
!
!
! ----------------------------------------------------------------------
!
!
    integer(kind=8) :: i, scdi, schc, iblo
    character(len=8) :: promes
    aster_logical :: modif
!
    real(kind=8) :: rtbloc
    character(len=19) :: nu, matas, stock
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: iakpee, iaphi0, iaphie, iascbl, iascdi
    integer(kind=8) :: iblold, iblph, ier, ii, iiblph, isingu, j
    integer(kind=8) :: jualf, k, kk, lgblph
    integer(kind=8) :: lmat, nbbloc, nblph, nddle, nddli, ndeci, nlblph
    integer(kind=8) :: npvneg
    integer(kind=8), pointer :: scib(:) => null()
    character(len=24), pointer :: refa(:) => null()
    integer(kind=8), pointer :: desm(:) => null()
    real(kind=8), pointer :: varm(:) => null()
    integer(kind=8), pointer :: vschc(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
    nu = nomu
    nu = nu(1:14)//'.NUME'
    stock = nu(1:14)//'.SLCS'
    matas = nomu//'.RIGIMECA'
!
    modif = .true.
    call dismoi('NOM_PROJ_MESU', nomu, 'MACR_ELEM_STAT', repk=promes)
    if (promes .eq. ' ') modif = .false.
!
    call jeveuo(nomu//'.VARM', 'E', vr=varm)
    rtbloc = varm(1)
    call jeveuo(nomu//'.DESM', 'E', vi=desm)
    nddle = desm(4)
!     NDDLE = 50
    nddli = desm(5)
!
!     -- FACTORISATION PARTIELLE DE LA MATRICE DE RIGIDITE:
!     -----------------------------------------------------
    call mtdscr(matas)
    call jeveuo(matas(1:19)//'.&INT', 'E', lmat)
    call tldlg3('LDLT', ' ', 1, lmat, 1, nddli, 0, &
                ndeci, isingu, npvneg, ier, ' ')
    if (ier .gt. 0) then
        call utmess('F', 'ALGORITH5_19')
    end if
!
!
!     -- ALLOCATION DE PHI_IE ET INITIALISATION PAR K_IE
!     -- ALLOCATION DE KP_EE  ET INITIALISATION PAR K_EE:
!     -------------------------------------------------------
!
    call mtdsc2(zk24(zi(lmat+1)), 'SCDI', 'L', iascdi)
    call jeveuo(zk24(zi(lmat+1)) (1:19)//'.REFA', 'L', vk24=refa)
    call jeveuo(refa(2) (1:14)//'.SLCS.SCHC', 'L', vi=vschc)
    call mtdsc2(zk24(zi(lmat+1)), 'SCBL', 'L', iascbl)
    call jelira(matas//'.UALF', 'NMAXOC', nbbloc)
    call jeveuo(stock//'.SCIB', 'L', vi=scib)
!
!     NLBLPH : NOMBRE DE LIGNES DE PHI_IE QUE L'ON PEUT REGROUPER
!              DANS UN BLOC DE LONGUEUR 5*RTBLOC
    nlblph = max(1, min(int(5*rtbloc*1024)/nddli, nddle))
!
!     LGBLPH : LONGUEUR DES BLOCS DE PHI_IE :
    lgblph = nlblph*nddli
!
!     NBLPH : NOMBRE DE BLOCS DE PHI_IE :
    nblph = (nddle*nddli-1)/lgblph+1
!
!
    call jecrec(nomu//'.PHI_IE', 'G V R', 'NU', 'DISPERSE', 'CONSTANT', &
                nblph)
    call jeecra(nomu//'.PHI_IE', 'LONMAX', lgblph)
!
    call jecrec(nomu//'.MAEL_RAID_VALE', 'G V R', 'NU', 'DISPERSE', &
                'CONSTANT', 1)
    call jeecra(nomu//'.MAEL_RAID_VALE', 'LONMAX', (nddle*(nddle+1)/2))
    call jecroc(jexnum(nomu//'.MAEL_RAID_VALE', 1))
    call jeveuo(jexnum(nomu//'.MAEL_RAID_VALE', 1), 'E', iakpee)

!   call wkvect(nomu//'.MAEL_RAID_VALE', 'G V R', (nddle*(nddle+1)/2), iakpee)
!
    iblold = 0
    j = 0
    do iblph = 1, nblph
        call jecroc(jexnum(nomu//'.PHI_IE', iblph))
        call jeveuo(jexnum(nomu//'.PHI_IE', iblph), 'E', iaphi0)
        do iiblph = 1, nlblph
            j = j+1
            if (j .gt. nddle) goto 40
            iaphie = iaphi0+(iiblph-1)*nddli
            iblo = scib(nddli+j)
            scdi = zi(iascdi-1+nddli+j)
            schc = vschc(nddli+j)
            if (iblo .ne. iblold) then
                if (iblold .gt. 0) call jelibe(jexnum(matas//'.UALF', iblold))
                call jeveuo(jexnum(matas//'.UALF', iblo), 'L', jualf)
            end if
            iblold = iblo
            k = 0
!
            if (modif) then
                do i = nddli+j+1-schc, nddli
                    k = k+1
                    zr(iaphie-1+i) = 0.d0
                end do
            else
!
                do i = nddli+j+1-schc, nddli
                    k = k+1
                    zr(iaphie-1+i) = zr(jualf-1+scdi-schc+k)
                end do
!
                do i = max(1, j+1-schc), j
                    ii = ((j-1)*j)/2+i
                    zr(iakpee-1+ii) = zr(jualf-1+scdi+i-j)
                end do
            end if
!
        end do
40      continue
        call jelibe(jexnum(nomu//'.PHI_IE', iblph))
    end do
    if (iblold .gt. 0) call jelibe(jexnum(matas//'.UALF', iblold))
!
!
!     -- CALCUL DE PHI_IE = (K_II**(-1))*K_IE:
!     ----------------------------------------
    if (modif) then
    else
        do iblph = 1, nblph
            call jeveuo(jexnum(nomu//'.PHI_IE', iblph), 'E', iaphi0)
            call rldlr8(zk24(zi(lmat+1)), vschc, zi(iascdi), zi(iascbl), nddli, &
                        nbbloc, zr(iaphi0), nlblph)
            call jelibe(jexnum(nomu//'.PHI_IE', iblph))
        end do
    end if
!
!
!     -- CALCUL DE KP_EE:
!     -------------------
    if (modif) then
        call crmeri(promes, iakpee)
    else
        iblold = 0
        do j = 1, nddle
            iblo = scib(nddli+j)
            scdi = zi(iascdi-1+nddli+j)
            schc = vschc(nddli+j)
            if (iblo .ne. iblold) then
                if (iblold .gt. 0) call jelibe(jexnum(matas//'.UALF', iblold))
                call jeveuo(jexnum(matas//'.UALF', iblo), 'L', jualf)
            end if
            iblold = iblo
!
            i = 0
            do iblph = 1, nblph
                call jeveuo(jexnum(nomu//'.PHI_IE', iblph), 'L', iaphi0)
                do iiblph = 1, nlblph
                    i = i+1
                    if (i .gt. j) goto 90
                    iaphie = iaphi0+(iiblph-1)*nddli
                    ii = (j*(j-1)/2)+i
                    kk = 0
                    do k = nddli+j+1-schc, nddli
                        kk = kk+1
                        zr(iakpee-1+ii) = zr(iakpee-1+ii)-zr(iaphie-1+k)*zr(jualf-1+scdi-schc+&
                                          &kk)
                    end do
                end do
90              continue
                call jelibe(jexnum(nomu//'.PHI_IE', iblph))
            end do
!
        end do
        if (iblold .gt. 0) call jelibe(jexnum(matas//'.UALF', iblold))
!
! FIN TEST SUR MODIF
    end if
!
!
    call jelibe(refa(2) (1:14)//'.SLCS.SCHC')
!
    call jedema()
end subroutine
