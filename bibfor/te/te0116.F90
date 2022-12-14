! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

subroutine te0116(nomopt, nomte)
!
implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
!
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
    integer :: ipg, npg, nb_vari, ivari, ispg
    integer :: jmate, jcompo, j_vari_out, j_vari_in, jtime, jcarcri
    character(len=16) :: rela_comp
    integer :: nb_res_mx
    parameter (nb_res_mx = 1)
    real(kind=8) :: valres(nb_res_mx), vvalres, check_rest
    integer :: codret(nb_res_mx)
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nomopt.eq.'REST_ECRO')
!
! - Get informations on current element
!
    call elrefe_info(fami = 'RIGI', npg = npg)
!
! - Get input fields
!
    call jevech('PMATERC', 'L', jmate)
    call jevech('PCOMPOR', 'L', jcompo)
    call jevech('PVARIMR', 'L', j_vari_in)
    call jevech('PTEMPSR', 'L', jtime)
    call jevech('PCARCRI', 'L', jcarcri)
!
    rela_comp = zk16(jcompo-1+1)
    read (zk16(jcompo-1+2),'(I16)') nb_vari
!
! - Get output field
!
    call jevech('PVARIPR', 'E', j_vari_out)
!
    check_rest = zr(jcarcri-1+21)
!
    if (check_rest .gt. 0.1) then
!
!     - Modify internal variables
!
        ispg    = 1
        vvalres = zr(jtime)
!
!
        if ((rela_comp.eq.'VMIS_ISOT_LINE').or.(rela_comp.eq.'VMIS_ISOT_TRAC'))  then
            do ipg = 1, npg
!
!             - Evaluate annealing function
!
                call rcvalb('RIGI', ipg, ispg, '+', zi(jmate),&
                            ' '   , 'REST_ECRO', 1, 'INST', [vvalres],&
                            nb_res_mx, 'FONC_MULT', valres, codret, 2)
!
!             - Annealing function bound's checking
!
                if ((valres(1).gt. 1.d0) .or. (valres(1).lt. 0.d0)) then
                    call utmess('F', 'COMPOR1_91', nr=2,valr=valres(1))
                endif
!
                zr(j_vari_out-1+nb_vari*(ipg-1)+1) = zr(j_vari_in-1+nb_vari*(ipg-1)+1)*valres(1)
                zr(j_vari_out-1+nb_vari*(ipg-1)+2) = zr(j_vari_in-1+nb_vari*(ipg-1)+2)
            end do
        elseif (rela_comp.eq.'VMIS_CINE_LINE') then
            do ipg = 1, npg
!
!             - Evaluate annealing function
!
                call rcvalb('RIGI', ipg, ispg, '+', zi(jmate),&
                            ' '   , 'REST_ECRO', 1, 'INST', [vvalres],&
                            nb_res_mx, 'FONC_MULT', valres, codret, 2)
!
!             - Annealing function bound's checking
!
                if ((valres(1).gt. 1.d0) .or. (valres(1).lt. 0.d0)) then
                    call utmess('F', 'COMPOR1_91', nr=2,valr=valres(1))
                endif
!
                do ivari = 1, 6
                    zr(j_vari_out-1+nb_vari*(ipg-1)+ivari) = zr(j_vari_in-1+nb_vari*(ipg-1)+ivari)&
                                                             *valres(1)
                end do
                zr(j_vari_out-1+nb_vari*(ipg-1)+7) = zr(j_vari_in-1+nb_vari*(ipg-1)+7)
            end do
        elseif ((rela_comp.eq.'VMIS_ECMI_LINE').or.(rela_comp.eq.'VMIS_CIN1_CHAB')) then
            do ipg = 1, npg
!
!             - Evaluate annealing function
!
                call rcvalb('RIGI', ipg, ispg, '+', zi(jmate),&
                            ' '   , 'REST_ECRO', 1, 'INST', [vvalres],&
                            nb_res_mx, 'FONC_MULT', valres, codret, 2)
!
!             - Annealing function bound's checking
!
                if ((valres(1).gt. 1.d0) .or. (valres(1).lt. 0.d0)) then
                    call utmess('F', 'COMPOR1_91', nr=2,valr=valres(1))
                endif
!
                zr(j_vari_out-1+nb_vari*(ipg-1)+1) = zr(j_vari_in-1+nb_vari*(ipg-1)+1)*valres(1)
                do ivari = 3, 8
                    zr(j_vari_out-1+nb_vari*(ipg-1)+ivari) = zr(j_vari_in-1+nb_vari*(ipg-1)+ivari)&
                                                             *valres(1)
                end do
                zr(j_vari_out-1+nb_vari*(ipg-1)+2) = zr(j_vari_in-1+nb_vari*(ipg-1)+2)
            end do
        elseif (rela_comp.eq.'VMIS_CIN2_CHAB') then
            do ipg = 1, npg
!
!             - Evaluate annealing function
!
                call rcvalb('RIGI', ipg, ispg, '+', zi(jmate),&
                            ' '   , 'REST_ECRO', 1, 'INST', [vvalres],&
                            nb_res_mx, 'FONC_MULT', valres, codret, 2)
!
!             - Annealing function bound's checking
!
                if ((valres(1).gt. 1.d0) .or. (valres(1).lt. 0.d0)) then
                    call utmess('F', 'COMPOR1_91', nr=2,valr=valres(1))
                endif
!
                zr(j_vari_out-1+nb_vari*(ipg-1)+1) = zr(j_vari_in-1+nb_vari*(ipg-1)+1)*valres(1)
                do ivari = 3, 14
                    zr(j_vari_out-1+nb_vari*(ipg-1)+ivari) = zr(j_vari_in-1+nb_vari*(ipg-1)+ivari)&
                                                             *valres(1)
                end do
                zr(j_vari_out-1+nb_vari*(ipg-1)+2) = zr(j_vari_in-1+nb_vari*(ipg-1)+2)
            end do
        endif
!
    else
!
!     - Internal variables OUT = Internal variables IN (no modifications)
!
        do ipg = 1, npg
            do ivari = 1, nb_vari
                zr(j_vari_out-1+nb_vari*(ipg-1)+ivari) = zr(j_vari_in-1+nb_vari*(ipg-1)+ivari)
            end do
        end do
    endif
!
end subroutine
