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
! person_in_charge: nicolas.pignet at edf.fr
!
subroutine cgDiscrField(cgField, cgTheta, cgStudy, cgStat, chsdeg, chslag, v_absc, v_basf, v_cesv, &
                        jcesd, jcesl, i_theta, lpain, lchin, nchin)
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
    character(len=19), intent(in) :: chsdeg, chslag
    integer(kind=8), intent(in) :: jcesd, jcesl, i_theta
    real(kind=8), pointer :: v_basf(:)
    real(kind=8), pointer :: v_absc(:)
    integer(kind=8), pointer  :: v_cesv(:)
    character(len=24), intent(inout) :: lchin(*)
    character(len=8), intent(inout) :: lpain(*)
    integer(kind=8), intent(inout) :: nchin
!
! --------------------------------------------------------------------------------------------------
!
!     CALC_G --- Utilities
!
!    Define dicretization field
!
! --------------------------------------------------------------------------------------------------
!
    character(len=19), parameter :: chdeg = '&&cgtheta.CHDEG'
    character(len=19), parameter :: chlag = '&&cgtheta.CHLAG'
    character(len=16) :: nomte
    character(len=19) :: ligrmo
    integer(kind=8) :: igr, nbgrel, nel, nute, iel, iad, ima, nncp, iret
    integer(kind=8), pointer :: v_liel(:) => null()
    real(kind=8) :: start, finish
!
    call cpu_time(start)
    call jemarq()
!
!   Champ constant pour construire theta_i dans le te
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
                if (cgTheta%discretization .eq. 'LEGENDRE') then
!
                    call cesexi('C', jcesd, jcesl, ima, 1, 1, 1, iad)
                    iad = abs(iad)
                    zl(jcesl-1+iad) = ASTER_TRUE
                    v_cesv(iad) = i_theta-1
!
                elseif (cgtheta%discretization .eq. 'LINEAIRE') then
!
                    call cesexi('C', jcesd, jcesl, ima, 1, 1, 1, iad)
                    iad = abs(iad)
                    zl(jcesl-1+iad) = ASTER_TRUE
                    if (i_theta .eq. 1) then
                        v_absc(iad) = v_basf(i_theta)
!---------------------- Cas fond fermé
                        if (cgtheta%l_closed) then
                            v_absc(iad) = v_basf(cgTheta%nb_theta_field-1)
                        end if
                    else
                        v_absc(iad) = v_basf(i_theta-1)
                    end if
!
                    call cesexi('C', jcesd, jcesl, ima, 1, 1, 2, iad)
                    iad = abs(iad)
                    zl(jcesl-1+iad) = ASTER_TRUE
                    v_absc(iad) = v_basf(i_theta)
!
                    call cesexi('C', jcesd, jcesl, ima, 1, 1, 3, iad)
                    iad = abs(iad)
                    zl(jcesl-1+iad) = ASTER_TRUE
                    if (i_theta .eq. cgTheta%nb_theta_field) then
                        v_absc(iad) = v_basf(i_theta)
                    else
                        v_absc(iad) = v_basf(i_theta+1)
                    end if
!
                end if
            end do
        end do
!
        if (cgTheta%discretization .eq. 'LEGENDRE') then
!           Conversion du champ par élément simple en champ par éléments
            call cescel(chsdeg, ligrmo, 'CALC_G', 'PDEG', 'NON', &
                        nncp, 'V', chdeg, 'F', iret)
!
            lpain(nchin+1) = 'PDEG'
            lchin(nchin+1) = chdeg
!
        elseif (cgtheta%discretization .eq. 'LINEAIRE') then
!           Conversion du champ par élément simple en champ par éléments
            call cescel(chslag, ligrmo, 'CALC_G', 'PLAG', 'NON', &
                        nncp, 'V', chlag, 'F', iret)
!
            lpain(nchin+1) = 'PLAG'
            lchin(nchin+1) = chlag
        else
            ASSERT(ASTER_FALSE)
        end if
!
        nchin = nchin+1
!
    end if
!
    call jedema()
    call cpu_time(finish)
    cgStat%cgCmpGtheta_disc = cgStat%cgCmpGtheta_disc+finish-start
    cgStat%nb_cgCmpGtheta_disc = cgStat%nb_cgCmpGtheta_disc+1
!
end subroutine
