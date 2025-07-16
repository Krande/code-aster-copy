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
subroutine te0116(nomopt, nomte)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/Behaviour_type.h"
!
    character(len=16), intent(in) :: nomte
    character(len=16), intent(in) :: nomopt
!
! --------------------------------------------------------------------------------------------------
!
! Computing the option REST_ECRO
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbResu = 3
    integer(kind=8) :: codret(nbResu)
    character(len=16), parameter :: resuName(nbResu) = &
                                    (/'TEMP_MINI', 'TEMP_MAXI', 'EPSQ_MINI'/)
    real(kind=8) :: resuVale(nbResu)
    integer(kind=8) :: ipg, npg, nb_vari, ivari, ispg
    integer(kind=8) :: jmate, j_vari_out, j_vari_in, jtimem, jtimep
    character(len=16), pointer :: compor(:) => null()
    character(len=16) :: rela_comp, postIncr
    real(kind=8) :: instp, instm
    real(kind=8) :: tau_inf(1), x0, alpha(1)
    real(kind=8) :: T1, T2, temp, k, dt
    real(kind=8) :: p_in, p_out, ecro_in, ecro_out, valExp
    integer(kind=8) :: iret, i_ecro_init
    integer(kind=8), parameter :: indxEpseq = 1
    real(kind=8) :: epsq_min
    aster_logical :: l_anneal, l_end_anneal
! - Protect against large overflow for exponential
    real(kind=8), parameter :: maxExp = 500.
! - Protect agains division by zero
    real(kind=8), parameter :: toleTemp = 1.
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nomopt .eq. 'REST_ECRO')

! - Initialisation
    x0 = 0.0
    ispg = 1
    p_in = 0.d0
    p_out = 0.d0
    ecro_in = 0.d0
    ecro_out = 0.d0

! - Get informations on current element
    call elrefe_info(fami='RIGI', npg=npg)

! - Get input fields
    call jevech('PMATERC', 'L', jmate)
    call jevech('PCOMPOR', 'L', vk16=compor)
    call jevech('PVARIMR', 'L', j_vari_in)
    call jevech('PINSTPR', 'L', jtimep)
    call jevech('PINSTMR', 'L', jtimem)
!
    rela_comp = compor(RELA_NAME)
    read (compor(NVAR), '(I16)') nb_vari
    postIncr = compor(POSTINCR)

! - Get output field
    call jevech('PVARIPR', 'E', j_vari_out)
!
    if (postIncr .eq. 'REST_ECRO') then
! -     Compute dt
        instm = zr(jtimem)
        instp = zr(jtimep)
        dt = instp-instm

! -     Get T1 et T2 : temperatures of the start and end of annealing!
!       and EPSQ_MIN : minimal strain for annealing to occur
        call rcvalb('RIGI', 1, ispg, &
                    '+', zi(jmate), ' ', 'REST_ECRO', &
                    0, ' ', [0.d0], &
                    nbResu, resuName, resuVale, codret, 2)
        T1 = resuVale(1)
        T2 = resuVale(2)
        epsq_min = resuVale(3)

! -     Get position of cumulative plastic strain of start of annealing (last internal variable)
        i_ecro_init = nb_vari

! -     Modify internal variables
        do ipg = 1, npg
!
            p_in = zr(j_vari_in-1+nb_vari*(ipg-1)+indxEpseq)
            ecro_in = zr(j_vari_in-1+nb_vari*(ipg-1)+i_ecro_init)
!
            call rcvarc(' ', 'TEMP', '+', 'RIGI', ipg, ispg, temp, iret)
            ASSERT(iret .eq. 0)
!
            l_anneal = ((temp .gt. T1) .and. (p_in .gt. epsq_min) .and. (temp .lt. T2))

            ! Je suis très poroche de la borne fin de la restauration
            ! Pour éviter la division par zéro je décide
            l_end_anneal = ASTER_FALSE
            if (abs(temp-T2) .le. toleTemp) then
                l_end_anneal = ASTER_TRUE
            end if
!
            if (.not. l_anneal) then
                p_out = p_in
                ecro_out = p_in
!
            else if (l_anneal) then
                if (l_end_anneal) then
                    p_out = 0.
                    ecro_out = ecro_in
                else
                    x0 = ecro_in
                    call rcvalb('RIGI', ipg, ispg, &
                                '+', zi(jmate), ' ', 'REST_ECRO', &
                                1, 'EPSI', [epsq_min], 1, 'COEF_ECRO', alpha, codret(1), 2)
                    call rcvalb('RIGI', ipg, ispg, &
                                '+', zi(jmate), ' ', 'REST_ECRO', &
                                1, 'EPSI', [epsq_min], 1, 'TAU_INF', tau_inf, codret(1), 2)
                    if (p_in .gt. (x0*tau_inf(1))) then
                        valExp = alpha(1)*(temp-T1)/(T2-temp)
                        if (valExp .ge. maxExp) then
                            p_out = 0.
                        else
                            k = exp(valExp)-1.d0
                            p_out = p_in*(1./(1.+k*dt))+ &
                                    (1-1./(1+k*dt))*(x0*tau_inf(1))
                            p_out = max(p_out, 0.)
                        end if
                    else
                        p_out = p_in
                    end if
                    ecro_out = ecro_in
                end if
            end if
!
            do ivari = 1, nb_vari
                zr(j_vari_out-1+nb_vari*(ipg-1)+ivari) = zr(j_vari_in-1+nb_vari*(ipg-1)+ivari)
            end do
            zr(j_vari_out-1+nb_vari*(ipg-1)+indxEpseq) = p_out
            zr(j_vari_out-1+nb_vari*(ipg-1)+i_ecro_init) = ecro_out

!
        end do
!
    else
!
!     - Internal variables OUT = Internal variables IN (no REST_ECRO)
!
        do ipg = 1, npg
            do ivari = 1, nb_vari
                zr(j_vari_out-1+nb_vari*(ipg-1)+ivari) = zr(j_vari_in-1+nb_vari*(ipg-1)+ivari)
            end do
        end do
    end if
!
end subroutine
