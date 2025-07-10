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
subroutine cgVerification(cgField, cgTheta, cgStudy, cgStat)
!
    use calcG_type
!
    implicit none
!
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/utmess.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterc/getfac.h"
#include "asterfort/jeveuo.h"
#include "jeveux.h"
!
    type(CalcG_field), intent(in) :: cgField
    type(CalcG_theta), intent(in) :: cgTheta
    type(CalcG_study), intent(in) :: cgStudy
    type(CalcG_stat), intent(inout) :: cgStat
!
! --------------------------------------------------------------------------------------------------
!
!     CALC_G --- Utilities
!
!     Verification of inputs
!
! --------------------------------------------------------------------------------------------------
!
!
    character(len=8) :: model, mesh, typmo, mesh0, nomgd
    aster_logical :: lmodemeca, ldynatrans
    integer(kind=8) :: nexci, nbel, i
    real(kind=8) :: start, finish, dirz, absccur, long
    real(kind=8), pointer :: jvale(:) => null()
!
    call cpu_time(start)
!
    call jemarq()
!
    if (cgField%level_info > 1) then
        call utmess('I', 'RUPTURE3_1')
    end if
!
    call dismoi('MODELE', cgField%result_in, 'RESULTAT', repk=model)
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
    ASSERT(cgTheta%mesh .eq. mesh)
!
    lmodemeca = cgField%isModeMeca()
    ldynatrans = cgField%isDynaTrans()
!
! --- EXCIT is allowed only for MODE_MECA and DYNA_TRANS
!
    call getfac('EXCIT', nexci)
    if (lmodemeca .or. ldynatrans) then
        if (nexci == 0) then
            call utmess('I', 'RUPTURE3_6')
        end if
    else
        if (nexci > 0) then
            call utmess('F', 'RUPTURE3_7')
        end if
    end if

!--- Cas axis : on normalise informe l'utilisateur de la division
!    automatique par 1/R
    call dismoi('MODELISATION', cgStudy%model, 'MODELE', repk=typmo)
    if (typmo(1:4) .eq. 'AXIS') then
        call utmess('I', 'RUPTURE3_10')
    end if
!
!--- Verify the input theta factors field
    if (cgTheta%theta_factors_in) then
        call dismoi('NOM_MAILLA', cgTheta%theta_factors, 'CHAM_NO', repk=mesh0)
        ASSERT(mesh0 .eq. mesh)

        call dismoi('NOM_GD', cgTheta%theta_factors, 'CHAM_NO', repk=nomgd)
        if ((nomgd(1:6) .ne. 'THET_R')) then
            call utmess('F', 'RUPTURE3_5', sk=nomgd)
        end if

        ! Vérifications spécifiques en 3D
        if (cgField%ndim .eq. 3) then
            ! VERIFICATION PRESENCE NB_POINT_FOND
            if (cgTheta%nb_point_fond .ne. 0) then
                ! INTERDICTION D AVOIR NB_POINT_FOND AVEC
                ! DISCTRETISATION =  LEGENDRE
                if (cgTheta%discretization .eq. 'LEGENDRE') then
                    call utmess('F', 'RUPTURE1_73')
                end if
            end if
        end if

        ! Vérifications spécifiques en 3D
        ! Les cmp DIR_Z /ABSC_CUR /LONG doivent 0
        if (cgField%ndim .eq. 2) then
            call dismoi('NB_NO_MAILLA', mesh, 'MAILLAGE', repi=nbel)
            call jeveuo(cgTheta%theta_factors(1:19)//'.VALE', 'L', vr=jvale)
            dirz = 0.0
            absccur = 0.0
            long = 0.0
            do i = 1, nbel
                dirz = dirz+jvale((i-1)*6+4)
                absccur = absccur+jvale((i-1)*6+5)
                long = long+jvale((i-1)*6+6)
            end do
            if (.not. ((dirz .eq. 0.0) .and. (absccur .eq. 0.0) .and. (long .eq. 0.0))) then
                call utmess('F', 'RUPTURE1_76')
            end if
        end if

    end if

    call jedema()
!
    call cpu_time(finish)
    cgStat%cgVerif = finish-start
!
end subroutine
