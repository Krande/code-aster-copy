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
subroutine peenca(champ, long, vr, nbmail, nummai)
    implicit none
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/celver.h"
#include "asterfort/digdel.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nbelem.h"
#include "asterfort/nbgrel.h"
#include "asterfort/scalai.h"
#include "asterfort/utmess.h"
    character(len=*) :: champ
    integer(kind=8) :: long, nbmail, nummai(*)
    real(kind=8) :: vr(long)
!     FAIRE DES OPERATIONS SUR UN CHAM_ELEM DE TYPE ENERGIE
!            (NOTION D'INTEGRALE DU CHAMP SUR LE MODELE)
!     ------------------------------------------------------------------
! IN  : CHAMP  : NOM DU CHAM_ELEM
! IN  : LONG   : LONGUEUR DU VECTEUR VR
! OUT : VR     : VECTEUR CONTENANT LES RESULATTS GLOBAUX
! IN  : NBMAIL : = 0 , CALCUL SUR TOUT LE CHAM_ELEM
!                SINON CALCUL SUR UN NOMBRE DE MAILLES
! IN  : NUMMAI : NUMEROS DES MAILLES
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
    integer(kind=8) :: longt, mode
    real(kind=8) :: rzero, ztot
    character(len=4) :: docu
    character(len=8) :: scal
    character(len=19) :: champ2, ligrel
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ibid, icoef, idecgr, iel, im, inum
    integer(kind=8) :: j, k, nbgr
    integer(kind=8) :: nel
    character(len=24), pointer :: celk(:) => null()
    real(kind=8), pointer :: celv(:) => null()
    integer(kind=8), pointer :: liel(:) => null()
    integer(kind=8), pointer :: celd(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
    champ2 = champ
    rzero = 0.0d0
!
!     -- ON VERIFIE QUE LE CHAM_ELEM N'EST PAS TROP DYNAMIQUE :
    call celver(champ2, 'NBVARI_CST', 'STOP', ibid)
    call celver(champ2, 'NBSPT_1', 'STOP', ibid)
!
    call jelira(champ2//'.CELD', 'DOCU', cval=docu)
    if (docu .ne. 'CHML') then
        call utmess('F', 'CALCULEL3_52')
    end if
    call jeveuo(champ2//'.CELK', 'L', vk24=celk)
    ligrel = celk(1) (1:19)
!
    call jeveuo(champ2//'.CELD', 'L', vi=celd)
!
!     --- TYPE DE LA GRANDEUR ---
    scal = scalai(celd(1))
!
    nbgr = nbgrel(ligrel)
!
!     -- ON MET A ZERO LE VECTEUR "VSCAL":
    if (scal(1:1) .eq. 'R') then
        do i = 1, long
            vr(i) = rzero
        end do
    else
        call utmess('F', 'CALCULEL3_74', sk=scal)
    end if
!
    call jeveuo(champ2//'.CELV', 'L', vr=celv)
    if (nbmail .le. 0) then
        do j = 1, nbgr
            mode = celd(celd(4+j)+2)
            if (mode .eq. 0) goto 30
            longt = digdel(mode)
            icoef = max(1, celd(4))
            longt = longt*icoef
            nel = nbelem(ligrel, j)
            idecgr = celd(celd(4+j)+8)
            do k = 1, nel
!
!              --- TOTALE ---
                i = 1
                ztot = celv(idecgr+(k-1)*longt+i-1)
                vr(1) = vr(1)+ztot
            end do
30          continue
        end do
        vr(2) = 100.0d0
    else
        ztot = rzero
        do j = 1, nbgr
            mode = celd(celd(4+j)+2)
            if (mode .eq. 0) goto 34
            longt = digdel(mode)
            icoef = max(1, celd(4))
            longt = longt*icoef
            nel = nbelem(ligrel, j)
            idecgr = celd(celd(4+j)+8)
            do k = 1, nel
                ztot = ztot+celv(idecgr+(k-1)*longt)
            end do
34          continue
        end do
        call jeveuo(ligrel//'.LIEL', 'L', vi=liel)
        do im = 1, nbmail
            inum = 0
            do j = 1, nbgr
                mode = celd(celd(4+j)+2)
                nel = nbelem(ligrel, j)
!
                if (mode .eq. 0) then
                    inum = inum+nel+1
                    goto 42
                end if
                longt = digdel(mode)
                icoef = max(1, celd(4))
                longt = longt*icoef
!
                idecgr = celd(celd(4+j)+8)
                do k = 1, nel
                    iel = liel(1+inum+k-1)
                    if (iel .ne. nummai(im)) goto 44
!
!                 --- TOTALE ---
                    i = 1
                    vr(1) = vr(1)+celv(idecgr+(k-1)*longt+i-1)
                    goto 40
44                  continue
                end do
                inum = inum+nel+1
42              continue
            end do
40          continue
        end do
        if ((vr(1) .lt. r8prem()) .and. (ztot .lt. r8prem())) then
            vr(2) = 0.0d0
        else
            vr(2) = 100.0d0*vr(1)/ztot
        end if
    end if
!
    call jedema()
end subroutine
