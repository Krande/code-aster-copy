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
subroutine gchs2f(char1, char2, char3)
    implicit none
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/jedema.h"
#include "asterfort/jedupo.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/lxlgut.h"
#include "asterfort/nbec.h"
#include "asterfort/wkvect.h"
!
    character(len=19) :: char1, char2, char3
!
!     BUT : TRANSFORME :
!              CHARGE 'SCALAIRE' --> CHARGE 'FONCTION' (CONSTANTE)
!
!           (ROUTINE SPECIFIQUE A L'OPERATEUR CALC_G,
!            APPELEE PAR GCHARF, DONT LE BUT EST DE
!            FUSIONNER CHAR1 ET CHAR2)
!
!     IN       :    CHAR1  :  CHARGE 'SCALAIRE'
!              :    CHAR2  :  CHARGE 'FONCTION'
!     IN/OUT   :    CHAR3  :  CHARGE 'FONCTION'
!
! ======================================================================
! ----------------------------------------------------------------------
    integer(kind=8) :: ncmp1, ncmp2, i, jval3, k, izo, jfpro
    integer(kind=8) :: jfval, kk, iec, reste, code, jncmp1, jncmp2, j, nec, ior
    real(kind=8) :: epsi
    character(len=8) :: nocmp1, nomfon
    character(len=19) :: nomf19
    real(kind=8), pointer :: vale(:) => null()
    integer(kind=8), pointer :: des1(:) => null()
    integer(kind=8), pointer :: des2(:) => null()
    integer(kind=8), pointer :: des3(:) => null()
!
    call jemarq()
!
    nomfon = char3(1:8)
    epsi = r8prem()
!
!     DUPLICATION AVANT MISE A JOUR
    call jedupo(char1//'.DESC', 'V', char3//'.DESC', .false._1)
    call jedupo(char1//'.NOMA', 'V', char3//'.NOMA', .false._1)
    call jedupo(char1//'.NOLI', 'V', char3//'.NOLI', .false._1)
    call jedupo(char1//'.LIMA', 'V', char3//'.LIMA', .false._1)
!
!     DESC (MAJ 1/2)
    call jeveuo(char1//'.DESC', 'L', vi=des1)
    call jeveuo(char2//'.DESC', 'L', vi=des2)
    call jeveuo(char3//'.DESC', 'E', vi=des3)
    des3(1) = des2(1)
!
    call jelira(jexnum('&CATA.GD.NOMCMP', des1(1)), 'LONMAX', ncmp1)
    call jelira(jexnum('&CATA.GD.NOMCMP', des2(1)), 'LONMAX', ncmp2)
    call jeveuo(jexnum('&CATA.GD.NOMCMP', des1(1)), 'L', jncmp1)
    call jeveuo(jexnum('&CATA.GD.NOMCMP', des2(1)), 'L', jncmp2)
!
!     VALE
    call jeveuo(char1//'.VALE', 'L', vr=vale)
    call wkvect(char3//'.VALE', 'V V K8', ncmp2*des3(3), jval3)
!
    k = 0
    kk = 0
    do i = 1, ncmp1
        nocmp1 = zk8(jncmp1-1+i)
        j = indik8(zk8(jncmp2), nocmp1, 1, ncmp2)
        if (j .ne. 0) then
            k = k+1
            if (k .eq. 1) then
                iec = (i-1)/30+1
                reste = i-30*(iec-1)
                code = 2**reste
            end if
            do izo = 1, des1(3)
                if (abs(vale(1+(izo-1)*ncmp1+i-1)) .gt. epsi) then
                    kk = kk+1
                    call codent(kk, 'D0', nomfon(7:8))
                    nomf19 = nomfon
                    ASSERT(lxlgut(nomf19) .le. 24)
                    call wkvect(nomf19//'.PROL', 'V V K24', 6, jfpro)
                    zk24(jfpro) = 'CONSTANT'
                    zk24(jfpro+1) = 'LIN LIN'
                    zk24(jfpro+2) = 'TOUTPARA'
                    zk24(jfpro+3) = 'TOUTRESU'
                    zk24(jfpro+4) = 'CC'
                    zk24(jfpro+5) = nomf19
                    call wkvect(nomf19//'.VALE', 'V V R', 2, jfval)
                    zr(jfval) = 1.d0
                    zr(jfval+1) = vale(1+(izo-1)*ncmp1+i-1)
                    zk8(jval3+(izo-1)*ncmp2+j-1) = nomfon
                else
                    zk8(jval3+(izo-1)*ncmp2+j-1) = '&FOZERO'
                end if
            end do
        end if
    end do
!
!     DESC (MAJ 2/2)
    nec = nbec(des3(1))
    do i = 1, nec*des3(3)
        des3(1+3+2*des3(3)+i-1) = ior(des3(1+3+2*des3(1+3- &
                                                      1)+i-1), code)
    end do
!
    call jedema()
!
end subroutine
