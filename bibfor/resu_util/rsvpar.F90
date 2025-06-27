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
! aslint: disable=W0413
!
subroutine rsvpar(resultNameZ, numeStore, paraNameZ, &
                  paraValeI_, paraValeR_, paraValeKZ_, &
                  ier_)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsnopa.h"
!
    character(len=*), intent(in) :: resultNameZ
    integer(kind=8), intent(in) :: numeStore
    character(len=*), intent(in) :: paraNameZ
    integer(kind=8), optional, intent(in) :: paraValeI_
    real(kind=8), optional, intent(in) :: paraValeR_
    character(len=*), optional, intent(in) :: paraValeKZ_
    integer(kind=8), optional, intent(out) :: ier_
! --------------------------------------------------------------------------------------------------
!      VERIFICATION DE L'EXISTANCE D'UN NOM DE PARAMETRE ET DE
!      SA VALEUR DANS UN RESULTAT COMPOSE
! --------------------------------------------------------------------------------------------------
! IN  : NOMSD  : NOM DE LA STRUCTURE "RESULTAT".
! IN  : IORDR  : NUMERO D'ORDRE A TRAITER.
! IN  : NOMPAR : NOM DU PARAMETRE A VERIFIER.
! IN  : IPAR   : VALEUR DU PARAMETRE ( TYPE INTEGER )
! IN  : RPAR   : VALEUR DU PARAMETRE ( TYPE REAL )
! IN  : KPAR   : VALEUR DU PARAMETRE ( TYPE CHARACTER )
! OUT : IER    : = 0    CE N'EST PAS UN PARAMETRE DU "RESULTAT".
!              : = 110  LA VALEUR DU PARAMETRE N'EST PAS CORRECTE.
!              : = 100  LA VALEUR DU PARAMETRE EST CORRECTE.
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iPara, nbPara, nbAcce, ier
    character(len=3) :: ctype
    integer(kind=8) :: paraValeI
    real(kind=8) :: paraValeR
    character(len=80) :: paraValeK
    integer(kind=8) :: jadr
    character(len=16) :: paraName
    character(len=16), pointer :: listParaName(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    ier = 0

    paraName = paraNameZ
    paraValeI = 0
    paraValeR = 0.d0
    paraValeK = " "
    if (present(paraValeI_)) then
        paraValeI = paraValeI_
    end if
    if (present(paraValeR_)) then
        paraValeR = paraValeR_
    end if
    if (present(paraValeKZ_)) then
        paraValeK = paraValeKZ_
    end if
!
    call rsnopa(resultNameZ, 1, '&&RSVPAR.NOMS_PARA', nbAcce, nbPara)
    call jeveuo('&&RSVPAR.NOMS_PARA', 'L', vk16=listParaName)
    do iPara = 1, nbPara
        if (paraName .eq. listParaName(iPara)) then
            goto 12
        end if
    end do
    goto 999
!
12  continue
    ier = 110
    call rsadpa(resultNameZ, 'L', 1, paraName, numeStore, &
                1, sjv=jadr, styp=ctype, istop=0)
    if (ctype(1:1) .eq. 'I') then
        if (zi(jadr) .eq. paraValeI) ier = 100
    else if (ctype(1:1) .eq. 'R') then
        if (zr(jadr) .eq. paraValeR) ier = 100
    else if (ctype(1:3) .eq. 'K80') then
        if (zk80(jadr) .eq. paraValeK) ier = 100
    else if (ctype(1:3) .eq. 'K32') then
        if (zk32(jadr) .eq. paraValeK) ier = 100
    else if (ctype(1:3) .eq. 'K24') then
        if (zk24(jadr) .eq. paraValeK) ier = 100
    else if (ctype(1:3) .eq. 'K16') then
        if (zk16(jadr) .eq. paraValeK) ier = 100
    else if (ctype(1:2) .eq. 'K8') then
        if (zk8(jadr) .eq. paraValeK) ier = 100
    end if
!
999 continue
    if (present(ier_)) then
        ier_ = ier
    end if
    call jedetr('&&RSVPAR.NOMS_PARA')
    call jedema()
end subroutine
