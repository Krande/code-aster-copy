! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
subroutine asdir(monoap, muapde, id, neq, nbsup, &
                 nsupp, tcosup, recmod, repdir)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/wkvect.h"
    integer :: id, neq, nbsup, nsupp(*), tcosup(nbsup, *)
    real(kind=8) :: recmod(nbsup, neq, *), repdir(neq, *)
    aster_logical :: monoap, muapde
!     COMMANDE : COMB_SISM_MODAL
!        CALCUL DES REPONSES DIRECTIONNELLES DES SUPPORTS
!     ------------------------------------------------------------------
! IN  : MONOAP : =.TRUE.  , CAS DU MONO-SUPPORT
!                =.FALSE. , CAS DU MULTI-SUPPORT
! IN  : MUAPDE : =.TRUE.  , CAS DU MULTI-SUPPORTS DECORRELES
!                =.FALSE. , CAS DU MULTI-SUPPORTS CORRELES
! IN  : ID     : LA DIRECTION
! IN  : NEQ    : NOMBRE D'EQUATIONS
! IN  : NBSUP  : NOMBRE DE SUPPORTS
! IN  : NSUPP  : MAX DU NOMBRE DE SUPPORT PAR DIRECTION
! IN  : TCOSUP : VECTEUR DES TYPES DE RECOMBINAISON DES SUPPORTS
! IN  : RECMOD : VECTEUR DES RECOMBINAISONS MODALES PAR APPUIS
! OUT : REPDIR : VECTEUR DES RECOMBINAISONS PAR DIRECTIONS
!     ------------------------------------------------------------------
    integer :: in, is, jabs
    real(kind=8) :: xxx, xx1, xx2
    real(kind=8), pointer :: line(:) => null()
    real(kind=8), pointer :: quad(:) => null()
!     ------------------------------------------------------------------
!
    call jemarq()
!
    if (monoap .or. .not. muapde) then
        do in = 1, neq
            repdir(in, id) = recmod(1, in, id)
        end do
    else
        AS_ALLOCATE(vr=quad, size=neq)
        AS_ALLOCATE(vr=line, size=neq)
        call wkvect('&&ASDIR.ABS ', 'V V R', neq, jabs)
        do is = 1, nsupp(id)
            if (tcosup(is, id) .eq. 1) then
!              --- COMBINAISON QUADRATIQUE ---
                do in = 1, neq
                    xxx = recmod(is, in, id)
                    quad(in) = quad(in)+xxx
                end do
            else if (tcosup(is, id) .eq. 2) then
!              --- COMBINAISON LINEAIRE ---
                do in = 1, neq
                    if (recmod(is, in, id) .ge. 0.d0) then
                        xxx = sqrt(recmod(is, in, id))
                        line(in) = line(in)+xxx
                    end if
                end do
            else
!              --- COMBINAISON VALEUR ABSOLUE ---
                do in = 1, neq
                    xxx = sqrt(abs(recmod(is, in, id)))
                    zr(jabs+in-1) = zr(jabs+in-1)+xxx
                end do
            end if
        end do
        do in = 1, neq
            xx1 = line(in)*line(in)
            xx2 = zr(jabs+in-1)*zr(jabs+in-1)
            repdir(in, id) = quad(in)+xx1+xx2
        end do
        AS_DEALLOCATE(vr=quad)
        AS_DEALLOCATE(vr=line)
        call jedetr('&&ASDIR.ABS ')
    end if
    call jedema()
end subroutine
