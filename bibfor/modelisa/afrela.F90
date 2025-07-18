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

subroutine afrela(coef_real, coef_cplx, dof_name, node_name, repe_type, &
                  repe_defi, nbterm, vale_real, vale_cplx, vale_func, &
                  type_coef, vale_type, epsi, lisrez)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/crelrl.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/juveca.h"
#include "asterfort/utmess.h"
!
!  Person in charge: mickael.abbas at edf.fr
!
    integer(kind=8), intent(in) :: nbterm
    real(kind=8), intent(in) :: coef_real(nbterm)
    complex(kind=8), intent(in) :: coef_cplx(nbterm)
    character(len=8), intent(in) :: dof_name(nbterm)
    character(len=8), intent(in) :: node_name(nbterm)
    integer(kind=8), intent(in) :: repe_type(nbterm)
    real(kind=8), intent(in) :: repe_defi(3, nbterm)
    real(kind=8), intent(in) :: vale_real
    complex(kind=8), intent(in) :: vale_cplx
    character(len=*), intent(in) :: vale_func
    character(len=4), intent(in) :: type_coef
    character(len=4), intent(in) :: vale_type
    real(kind=8), intent(in) :: epsi
    character(len=*), intent(in) :: lisrez
!
! --------------------------------------------------------------------------------------------------
!
! New linear relation
!
!       coef(iterm) . dof_name(iterm) = vale
!
! With:
!       iterm = 1,nbterm
!       coef = coef_real if type_coef = 'REEL'
!       coef = coef_cplx if type_coef = 'COMP'
!       vale = vale_real if vale_type = 'REEL'
!       vale = vale_cplx if vale_type = 'COMP'
!       vale = vale_func if vale_type = 'FONC'
!
! Saving in lisrel object (created if not exist)
!
! --------------------------------------------------------------------------------------------------
!
! In  coef_real : real coefficient
! In  coef_cplx : complex coefficient
! In  dof_name  : name of dof in linear relation
! In  node_name : name of nodes in linear relation
! In  repe_type : type of coordiante system to apply linear relation
!                   if 0 - Global coordinate system
!                   if 2 - Local (2D) coordinate system
!                   if 3 - Local (2D) coordinate system
! In  repe_defi : defintion of local coordinate system
! In  nbterm    : number of terms in linear relation
! In  vale_type : affected value type (real, complex or function)
! In  vale_real : affected value if real
! In  vale_func : affected value if function
! In  vale_cplx : affected value if complex
! In  coef_type : type of coefficient (real or complex)
! In  epsi      : tolerance to detect "zero" coefficient
! In  lisrel    : JEVEUX object list_rela for relations list management
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: vale_real_norm
    complex(kind=8) :: vale_cplx_norm
    integer(kind=8) :: imult
    character(len=8) :: dof_name_tran(3), dof_name_rota(3)
    character(len=19) :: lisrel
    integer(kind=8) :: idbeta, idcoef, idim, idnbre
    integer(kind=8) :: idpoin, idsurc, ifm, ipoint, iret
    integer(kind=8) :: iterm, idirect, lonuti, lveclr, mdim, nbrel0
    integer(kind=8) :: nbrela, nbrmax, nbterr, niv, k
    real(kind=8) :: norm_coef
    aster_logical :: l_rota
    integer(kind=8), pointer :: rlnt(:) => null()
    character(len=8), pointer :: rlno(:) => null()
    character(len=8), pointer :: rldd(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infniv(ifm, niv)
!
! - Initializations
!
    lisrel = lisrez
    dof_name_tran(1) = 'DX'
    dof_name_tran(2) = 'DY'
    dof_name_tran(3) = 'DZ'
    dof_name_rota(1) = 'DRX'
    dof_name_rota(2) = 'DRY'
    dof_name_rota(3) = 'DRZ'
    l_rota = .false.
    vale_real_norm = vale_real
    vale_cplx_norm = vale_cplx
!
! - Information about linear relation before normalization
!
    if (niv .eq. 2) then
        write (ifm, *) ' '
        write (ifm, '(A,I4,A)') ' _RELA IMPRESSION D''UNE RELATION LINEAIRE ENTRE ', nbterm, &
            ' DDLS. (AVANT NORMALISATION DE LA RELATION)'
        do iterm = 1, nbterm
            if (repe_type(iterm) .eq. 0) then
                if (type_coef .eq. 'REEL') then
                    write (ifm, 101) coef_real(iterm), node_name(iterm), dof_name(iterm)
                else if (type_coef .eq. 'COMP') then
                    write (ifm, 103) dble(coef_cplx(iterm)), dimag(coef_cplx(iterm)), &
                        node_name(iterm), dof_name(iterm)
                else
                    ASSERT(.false.)
                end if
            else
                if (type_coef .eq. 'REEL') then
                    write (ifm, 102) coef_real(iterm), node_name(iterm), dof_name(iterm), &
                        (repe_defi(idirect, iterm), idirect=1, repe_type(iterm))
                else if (type_coef .eq. 'COMP') then
                    write (ifm, 104) dble(coef_cplx(iterm)), dimag(coef_cplx(iterm)), &
                        node_name(iterm), dof_name(iterm), &
                        (repe_defi(idirect, iterm), idirect=1, repe_type(iterm))
                else
                    ASSERT(.false.)
                end if
            end if
        end do
        if (vale_type .eq. 'REEL') then
            write (ifm, *) '_RELA = ', vale_real
        else if (vale_type .eq. 'COMP') then
            write (ifm, *) '_RELA = ', vale_cplx
        else if (vale_type .eq. 'FONC') then
            write (ifm, *) '_RELA = ', vale_func
        else
            ASSERT(.false.)
        end if
    end if
!
! - Normalization ratio
!
    if (type_coef .eq. 'REEL') then
        norm_coef = 0.d0
        do iterm = 1, nbterm
            norm_coef = max(norm_coef, abs(coef_real(iterm)))
        end do
        if (norm_coef .eq. 0.d0) then
            call utmess('F', 'CHARGES2_97')
        end if
    else if (type_coef .eq. 'COMP') then
        norm_coef = 0.d0
        do iterm = 1, nbterm
            norm_coef = max(norm_coef, abs(coef_cplx(iterm)))
        end do
        if (norm_coef .eq. 0.d0) then
            call utmess('F', 'CHARGES2_97')
        end if
    else
        ASSERT(.false.)
    end if
!
! - Normalization of values
!
    if (vale_type .eq. 'REEL') then
        vale_real_norm = vale_real_norm/norm_coef
    else if (vale_type .eq. 'COMP') then
        vale_cplx_norm = vale_cplx_norm/norm_coef
    else if (vale_type .eq. 'FONC') then
! ----- Alarm if normalization ratio too much different from 1 ...
        if ((norm_coef .gt. 1.d3) .or. (norm_coef .lt. 1.d-3)) then
            call utmess('A', 'CHARGES2_99')
        end if
! ----- ... but cannot normalize function value !
        norm_coef = 1.d0
    end if
!
! - No <LIST_RELA> object -> creation
!
    call jeexin(lisrel//'.RLCO', iret)
    if (iret .eq. 0) call crelrl(type_coef, vale_type, 'V', lisrel)
!
! - How many linear relations ?
!
    call jeveuo(lisrel//'.RLNR', 'E', idnbre)
    nbrel0 = zi(idnbre)
    nbrela = nbrel0+1
!
! - Initial maximum linear relations number
!
    call jelira(lisrel//'.RLNT', 'LONMAX', nbrmax)
!
! - Length of vectors for all relation terms
!
    call jelira(lisrel//'.RLCO', 'LONMAX', lveclr)
!
! - Real length used
!
    call jeveuo(lisrel//'.RLPO', 'E', idpoin)
    if (nbrel0 .eq. 0) then
        lonuti = 0
    else
        lonuti = zi(idpoin+nbrel0-1)
    end if
!
! - Real number of terms: zero (epsi) terms vanished + active terms (local coordinate system)
!
    nbterr = 0
    do iterm = 1, nbterm
        if (type_coef .eq. 'COMP') then
            if (abs(coef_cplx(iterm)) .gt. epsi) then
                if (repe_type(iterm) .eq. 0) then
                    nbterr = nbterr+1
                else
                    nbterr = nbterr+repe_type(iterm)
                end if
            end if
        else if (type_coef .eq. 'REEL') then
            if (abs(coef_real(iterm)) .gt. epsi) then
                if (repe_type(iterm) .eq. 0) then
                    nbterr = nbterr+1
                else
                    nbterr = nbterr+repe_type(iterm)
                end if
            end if
        else
            ASSERT(.false.)
        end if
    end do
!
! - Increase object size if necessary
!
    if (lonuti+nbterr .ge. lveclr) then
        imult = (lonuti+nbterr)/lveclr+1
        call juveca(lisrel//'.RLCO', imult*lveclr)
        call juveca(lisrel//'.RLDD', imult*lveclr)
        call juveca(lisrel//'.RLNO', imult*lveclr)
    end if
!
! - No enough place for linear relation -> increase objects size
!
    if (nbrela .ge. nbrmax) then
        imult = nbrela/nbrmax+1
        call juveca(lisrel//'.RLBE', imult*nbrmax)
        call juveca(lisrel//'.RLNT', imult*nbrmax)
        call juveca(lisrel//'.RLPO', imult*nbrmax)
        call juveca(lisrel//'.RLSU', imult*nbrmax)
    end if
!
! - Linear relation access
!
    call jeveuo(lisrel//'.RLNR', 'E', idnbre)
    call jeveuo(lisrel//'.RLCO', 'E', idcoef)
    call jeveuo(lisrel//'.RLDD', 'E', vk8=rldd)
    call jeveuo(lisrel//'.RLNO', 'E', vk8=rlno)
    call jeveuo(lisrel//'.RLBE', 'E', idbeta)
    call jeveuo(lisrel//'.RLNT', 'E', vi=rlnt)
    call jeveuo(lisrel//'.RLPO', 'E', idpoin)
    call jeveuo(lisrel//'.RLSU', 'E', idsurc)
!
! - New length
!
    zi(idnbre) = nbrela
    if (nbrel0 .eq. 0) then
        ipoint = 0
    else
        ipoint = zi(idpoin+nbrel0-1)
    end if
    rlnt(nbrela) = nbterr
    if (nbrel0 .eq. 0) then
        zi(idpoin) = nbterr
    else
        zi(idpoin+nbrela-1) = zi(idpoin+nbrel0-1)+nbterr
    end if
!
! - New linear relation affectation
!
    k = 0
!
    if (type_coef .eq. 'COMP') then
        do iterm = 1, nbterm
            if (abs(coef_cplx(iterm)) .gt. epsi) then
                if (repe_type(iterm) .eq. 0) then
!
! ----------------- Global coordinate system
!
                    k = k+1
                    zc(idcoef+ipoint+k-1) = coef_cplx(iterm)/norm_coef
                    rldd(1+ipoint+k-1) = dof_name(iterm)
                    rlno(1+ipoint+k-1) = node_name(iterm)
                else
!
! ----------------- Local coordinate system: rotation or translation ?
!
                    if (dof_name(iterm) .eq. 'DEPL') then
                        l_rota = .false.
                    else if (dof_name(iterm) .eq. 'ROTA') then
                        l_rota = .true.
                    else
                        ASSERT(.false.)
                    end if
!
! ----------------- Change coordinate system
! ----------------- DEPL --> repe_defi(1)*U  + repe_defi(2)*V  + repe_defi(3)*W
! ----------------- ROTA --> repe_defi(1)*RX + repe_defi(2)*RY + repe_defi(3)*RZ
!
                    mdim = repe_type(iterm)
                    do idim = 1, mdim
                        k = k+1
                        zc(idcoef+ipoint+k-1) = coef_cplx(iterm)/norm_coef*repe_defi(idim, iterm)
                        rlno(1+ipoint+k-1) = node_name(iterm)
                        if (.not. l_rota) then
                            rldd(1+ipoint+k-1) = dof_name_tran(idim)
                        else
                            rldd(1+ipoint+k-1) = dof_name_rota(idim)
                        end if
                    end do
                end if
            end if
        end do
!
    else if (type_coef .eq. 'REEL') then
        do iterm = 1, nbterm
            if (abs(coef_real(iterm)) .gt. epsi) then
                if (repe_type(iterm) .eq. 0) then
!
! ----------------- Global coordinate system
!
                    k = k+1
                    zr(idcoef+ipoint+k-1) = coef_real(iterm)/norm_coef
                    rldd(1+ipoint+k-1) = dof_name(iterm)
                    rlno(1+ipoint+k-1) = node_name(iterm)
                else
!
! ----------------- Local coordinate system: rotation or translation ?
!
                    if (dof_name(iterm) .eq. 'DEPL') then
                        l_rota = .false.
                    else if (dof_name(iterm) .eq. 'ROTA') then
                        l_rota = .true.
                    else
                        ASSERT(.false.)
                    end if
!
! ----------------- Change coordinate system
! ----------------- DEPL --> repe_defi(1)*U  + repe_defi(2)*V  + repe_defi(3)*W
! ----------------- ROTA --> repe_defi(1)*RX + repe_defi(2)*RY + repe_defi(3)*RZ
!
                    mdim = repe_type(iterm)
                    do idim = 1, mdim
                        k = k+1
                        zr(idcoef+ipoint+k-1) = coef_real(iterm)/norm_coef*repe_defi(idim, iterm)
                        rlno(1+ipoint+k-1) = node_name(iterm)
                        if (.not. l_rota) then
                            rldd(1+ipoint+k-1) = dof_name_tran(idim)
                        else
                            rldd(1+ipoint+k-1) = dof_name_rota(idim)
                        end if
                    end do
                end if
            end if
        end do
    else
        ASSERT(.false.)
    end if
!
! - Value affectation
!
    if (vale_type .eq. 'REEL') then
        zr(idbeta+nbrela-1) = vale_real_norm
    else if (vale_type .eq. 'COMP') then
        zc(idbeta+nbrela-1) = vale_cplx_norm
    else if (vale_type .eq. 'FONC') then
        zk24(idbeta+nbrela-1) = vale_func
    else
        ASSERT(.false.)
    end if
!
    call jedema()
!
101 format(' _RELA ', e14.7, a10, a10)
103 format(' _RELA ', e14.7, 1x, e14.7, a10, a10)
102 format(' _RELA ', e14.7, a10, a10, 3x, 3(1x, e14.7))
104 format(' _RELA ', e14.7, 1x, e14.7, a10, a10, 3x, 3(1x, e14.7))
!
end subroutine
