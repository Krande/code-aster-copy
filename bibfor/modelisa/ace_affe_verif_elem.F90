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

subroutine ace_affe_verif_elem(noma, jdme, lesmailles, ng, grpmail, grpnbma, &
                               ace_nu, mclf, coderet, TFVCode)
!
    use cara_elem_parameter_module
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jenonu.h"
#include "asterfort/utmess.h"
#include "asterfort/int_to_char8.h"
!
    character(len=8), intent(in)          :: noma
    integer(kind=8), intent(in)           :: ng, jdme, ace_nu
    integer(kind=8), intent(inout)        :: lesmailles(*)
    character(len=24), intent(inout)      :: grpmail(*)
    integer(kind=8), intent(inout)        :: grpnbma(*)
    character(len=*), intent(in)          :: mclf
    integer(kind=8), intent(out)          :: coderet
    integer(kind=8), optional, intent(in) :: TFVCode
!
    aster_logical     :: l_parallel_mesh
!
    integer(kind=8)   :: ii, jj, kk, jdgm, nbmagr, nummai, nutyel
    integer(kind=8)   :: iret, code0, code9, TFVCode0
    character(len=8)  :: nommai
    character(len=24) :: mlggma, vmessk(5)
!
!-----------------------------------------------------------------------
!
    mlggma = noma//'.GROUPEMA'
!
    coderet = 0
!   Toutes les mailles possibles de la liste des groupes de mailles
    if (ng .le. 0) goto 999
!
!   Il y a un VERIMA de fait avant. On est sur un processeur donc si en // :
!       '.GROUPEMA' n'existe pas, ce n'est pas une erreur
!       grpmail(ii) n'existe pas, ce n'est pas une erreur
!
    l_parallel_mesh = isParallelMesh(noma)
    if (l_parallel_mesh) then
        call jeexin(mlggma, iret)
        if (iret .eq. 0) goto 999
    end if
!
!   Pour ne pas faire le test dans la boucle
    if (present(TFVCode)) then
!
        bii2: do ii = 1, ng
            grpnbma(ii) = 0
            if (l_parallel_mesh) then
                call jenonu(jexnom(mlggma, grpmail(ii)), iret)
                if (iret .eq. 0) cycle bii2
            end if
            call jeveuo(jexnom(mlggma, grpmail(ii)), 'L', jdgm)
            call jelira(jexnom(mlggma, grpmail(ii)), 'LONUTI', nbmagr)
            grpnbma(ii) = nbmagr
            coderet = max(coderet, nbmagr)
            bgroup2: do jj = 1, nbmagr
                nummai = zi(jdgm+jj-1)
                nutyel = zi(jdme+nummai-1)
                do kk = elem_supp%aceind(ace_nu, 1), elem_supp%aceind(ace_nu, 2)
                    if (nutyel .eq. elem_supp%catanum(kk)) then
!                       Si la maille est déjà traitée ==> lesmailles(nummai) < 0
!                       On doit vérifier que la surcharge est autorisée
!                       Si son code est différent, la surcharge est interdite
                        code0 = lesmailles(nummai)
                        code9 = nutyel+2*elem_supp%MaxCataNum+TFVCode
!                       Si on est dans la surcharge
                        if (code0 .lt. 0) then
                            if (abs(code0) .ne. code9) then
                                TFVCode0 = abs(code0)-nutyel-2*elem_supp%MaxCataNum
                                vmessk(1) = trim(ACE_NM_ELEMENT(elem_supp%acenum(kk)))
                                call ACE_DeCodeNomFormeVaria(TFVCode0, vmessk(2:3))
                                call ACE_DeCodeNomFormeVaria(TFVCode, vmessk(4:5))
                                call utmess('F', 'AFFECARAELEM_34', nk=5, valk=vmessk)
                            end if
                        end if
                        lesmailles(nummai) = -code9
                        cycle bgroup2
                    end if
                end do
                nommai = int_to_char8(nummai)
                vmessk(1) = mclf
                vmessk(2) = nommai
                call utmess('F', 'MODELISA_8', nk=2, valk=vmessk)
            end do bgroup2
        end do bii2
!
    else
!
        bii1: do ii = 1, ng
            grpnbma(ii) = 0
            if (l_parallel_mesh) then
                call jenonu(jexnom(mlggma, grpmail(ii)), iret)
                if (iret .eq. 0) cycle bii1
            end if
            call jeveuo(jexnom(mlggma, grpmail(ii)), 'L', jdgm)
            call jelira(jexnom(mlggma, grpmail(ii)), 'LONUTI', nbmagr)
            grpnbma(ii) = nbmagr
            coderet = max(coderet, nbmagr)
            bgroup1: do jj = 1, nbmagr
                nummai = zi(jdgm+jj-1)
                nutyel = zi(jdme+nummai-1)
                do kk = elem_supp%aceind(ace_nu, 1), elem_supp%aceind(ace_nu, 2)
                    if (nutyel .eq. elem_supp%catanum(kk)) then
!                       La maille va être traitée ou retraitée
!                       La surcharge est autorisée, on ne vérifie donc pas son "code"
                        code0 = lesmailles(nummai)
                        ASSERT(abs(code0) .eq. nutyel)
                        lesmailles(nummai) = -nutyel
                        cycle bgroup1
                    end if
                end do
                nommai = int_to_char8(nummai)
                vmessk(1) = mclf
                vmessk(2) = nommai
                call utmess('F', 'MODELISA_8', nk=2, valk=vmessk)
            end do bgroup1
        end do bii1
!
    end if
!
999 continue
!
end subroutine ace_affe_verif_elem
