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
! person_in_charge: matthieu.le-cren at edf.fr
!
subroutine cgHighOrder(cgField, cgTheta, cgStudy, cgStat, chscer, chseli, v_cer, v_eli, &
                       jcesd2, jcesl2, chcer, cheli)
!
    use calcG_type
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/calcG_type.h"
#include "asterfort/cescel.h"
#include "asterfort/cesexi.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/nbelem.h"
#include "asterfort/typele.h"
!
    type(CalcG_field), intent(in) :: cgField
    type(CalcG_theta), intent(in) :: cgTheta
    type(CalcG_Study), intent(in) :: cgStudy
    type(CalcG_stat), intent(inout) :: cgStat
    character(len=19), intent(in) :: chscer, chseli
    integer(kind=8), intent(in) :: jcesd2, jcesl2
    real(kind=8), pointer :: v_cer(:)
    real(kind=8), pointer :: v_eli(:)
    character(len=19), intent(inout) :: chcer
    character(len=19), intent(inout) :: cheli
!
! --------------------------------------------------------------------------------------------------
!
!     CALC_G --- Utilities
!
!    Define dicretization field
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: nomte
    character(len=19) :: ligrmo
    integer(kind=8) :: igr, nbgrel, nel, nute, iel, iad, ima, nncp, iret
    integer(kind=8), pointer :: v_liel(:) => null()
    real(kind=8) :: start, finish
!
    call cpu_time(start)
    call jemarq()
!
    if (cgField%ndim .eq. 3) then
!
        call dismoi('NOM_LIGREL', cgStudy%model, 'MODELE', repk=ligrmo)
        call jelira(ligrmo//'.LIEL', 'NMAXOC', nbgrel)
!
        do igr = 1, nbgrel
!
!           Récupération nombre d'éléments dans le groupe
            nel = nbelem(ligrmo, igr)
            call jeveuo(jexnum(ligrmo//'.LIEL', igr), 'L', vi=v_liel)
!
!           Récupération du nom du type d'éléments du groupe
            nute = typele(ligrmo, igr)
            call jenuno(jexnum('&CATA.TE.NOMTE', nute), nomte)
!
!           Boucle sur les éléments du groupe
            do iel = 1, nel
                ima = v_liel(iel)
                if (ima .lt. 0) cycle
!
                if (cgTheta%form_fiss .eq. 'CERCLE') then
!
                    call cesexi('C', jcesd2, jcesl2, ima, 1, 1, 1, iad)
                    iad = abs(iad)
                    zl(jcesl2-1+iad) = ASTER_TRUE
                    v_cer(iad) = cgTheta%rayon
!
                elseif (cgTheta%form_fiss .eq. 'ELLIPSE') then
!
                    call cesexi('C', jcesd2, jcesl2, ima, 1, 1, 1, iad)
                    iad = abs(iad)
                    zl(jcesl2-1+iad) = ASTER_TRUE
                    v_eli(iad) = cgTheta%demi_grand_axe
!
                    call cesexi('C', jcesd2, jcesl2, ima, 1, 1, 2, iad)
                    iad = abs(iad)
                    zl(jcesl2-1+iad) = ASTER_TRUE
                    v_eli(iad) = cgTheta%demi_petit_axe
!
                end if
            end do
        end do
!
        if (cgTheta%form_fiss .eq. 'CERCLE') then
!           Conversion du champ par élément simple en champ par éléments
            call cescel(chscer, ligrmo, 'CALC_G', 'PCER', 'NON', &
                        nncp, 'V', chcer, 'F', iret)
!
        elseif (cgTheta%form_fiss .eq. 'ELLIPSE') then
!           Conversion du champ par élément simple en champ par éléments
            call cescel(chseli, ligrmo, 'CALC_G', 'PELI', 'NON', &
                        nncp, 'V', cheli, 'F', iret)
!
        else
            ASSERT(ASTER_FALSE)
        end if

    end if
!
    call jedema()
    call cpu_time(finish)
    cgStat%cgCmpGtheta_disc = cgStat%cgCmpGtheta_disc+finish-start
    cgStat%nb_cgCmpGtheta_disc = cgStat%nb_cgCmpGtheta_disc+1
!
end subroutine
