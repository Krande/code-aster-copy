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
subroutine dltini(lcrea, nume, result, depini, vitini, &
                  accini, fexini, famini, fliini, neq, &
                  numedd, inchac, ds_energy)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/chpver.h"
#include "asterfort/getvid.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nonlinDSEnergyCreate.h"
#include "asterfort/nonlinDSEnergyInit.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsrusd.h"
#include "asterfort/utmess.h"
#include "asterfort/vtcopy.h"
#include "asterfort/vtcreb.h"
#include "blas/dcopy.h"
!
    character(len=8), intent(in) :: result
    type(NL_DS_Energy), intent(out) :: ds_energy
!
! --------------------------------------------------------------------------------------------------
!
!     CALCUL MECANIQUE TRANSITOIRE PAR INTEGRATION DIRECTE
!     RECUPERATION DES CONDITIONS INITIALES
!
! --------------------------------------------------------------------------------------------------
!
! OUT : LCREA  : CREATION OU NON DU RESULTAT
! OUT : NUME   : NUMERO D'ORDRE DE REPRISE
! OUT : DEPINI : CHAMP DE DEPLACEMENT INITIAL OU DE REPRISE
! OUT : VITINI : CHAMP DE VITESSE INITIALE OU DE REPRISE
! OUT : ACCINI : CHAMP D'ACCELERATION INITIALE OU DE REPRISE
! IN  : NEQ    : NOMBRE D'EQUATIONS
! IN  : NUMEDD : NUMEROTATION DDL
! VAR : INCHAC : CALCUL OU NON DE L'ACCELERATION INITIALE
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: depini(*), vitini(*), accini(*)
    real(kind=8) :: fexini(*), famini(*), fliini(*)
    character(len=24) :: numedd
    aster_logical :: lcrea
    integer(kind=8) :: nume, neq, inchac
    integer(kind=8) :: ire, iret, jvale
    integer(kind=8) :: nai, ndi, ndy, nvi
    integer(kind=8) :: ierr
    character(len=8) :: reuse, dep, vit
    character(len=19) :: champ, cham2
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    lcrea = .true.
!
! - Energy management
    call nonlinDSEnergyCreate(ds_energy)
    call getfac('ENERGIE', iret)
    ds_energy%l_comp = iret .gt. 0
    ds_energy%command = 'DYNA_VIBRA'
    call nonlinDSEnergyInit(result, ds_energy)
!
! - Reuse or not ?
    call getvid('ETAT_INIT', 'RESULTAT', iocc=1, scal=reuse, nbret=ndy)
!
!====
! 3. EN REPRISE
!====
!
    if (ndy .ne. 0) then
!
        call utmess('I', 'DYNAMIQUE_78', sk=reuse)
!
!        --- RECUPERATION DES CHAMPS DEPL VITE ET ACCE ---
        call rsexch(' ', reuse, 'DEPL', nume, champ, &
                    iret)
        if (iret .ne. 0) then
            call utmess('F', 'DYNALINE1_25')
        else
            call jeveuo(champ//'.VALE', 'L', jvale)
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(jvale), b_incx, depini, b_incy)
        end if
        call rsexch(' ', reuse, 'VITE', nume, champ, &
                    iret)
        if (iret .ne. 0) then
            call utmess('F', 'DYNALINE1_26')
        else
            call jeveuo(champ//'.VALE', 'L', jvale)
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(jvale), b_incx, vitini, b_incy)
        end if
        call rsexch(' ', reuse, 'ACCE', nume, champ, &
                    iret)
        if (iret .ne. 0) then
            call utmess('F', 'DYNALINE1_27')
        else
            call jeveuo(champ//'.VALE', 'L', jvale)
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(jvale), b_incx, accini, b_incy)
        end if
        call rsexch(' ', reuse, 'FORC_EXTE', nume, champ, &
                    iret)
        if (iret .eq. 0) then
            call jeveuo(champ//'.VALE', 'L', jvale)
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(jvale), b_incx, fexini, b_incy)
        end if
        call rsexch(' ', reuse, 'FORC_AMOR', nume, champ, &
                    iret)
        if (iret .eq. 0) then
            call jeveuo(champ//'.VALE', 'L', jvale)
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(jvale), b_incx, famini, b_incy)
        end if
        call rsexch(' ', reuse, 'FORC_LIAI', nume, champ, &
                    iret)
        if (iret .eq. 0) then
            call jeveuo(champ//'.VALE', 'L', jvale)
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(jvale), b_incx, fliini, b_incy)
        end if
!
!        --- CREE-T-ON UNE NOUVELLE STRUCTURE ? ---
        if (result .eq. reuse) then
            lcrea = .false.
            call rsrusd(result, nume+1)
        end if
!====
! 4. --- RECUPERATION DES CONDITIONS INITIALES ---
!====
!
    else
        call jeexin(result(1:8)//'           .REFD', ire)
        if (ire .gt. 0) then
            lcrea = .false.
        end if
!
        nume = 0
        call getvid('ETAT_INIT', 'DEPL', iocc=1, scal=champ, nbret=ndi)
        if (ndi .gt. 0) then
            call chpver('F', champ, 'NOEU', 'DEPL_R', ierr)
            inchac = 1
            cham2 = '&&COMDLT.DEPINI'
            call vtcreb(cham2, 'V', 'R', nume_ddlz=numedd, nb_equa_outz=neq)
            call vtcopy(champ, cham2, iret)
            if (iret .ne. 0) then
                call utmess("A", "FIELD0_6", sk='DEPL')
            end if
            call jeveuo(cham2//'.VALE', 'L', jvale)
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(jvale), b_incx, depini, b_incy)
            dep = champ(1:8)
        else
            dep = 'nul'
        end if
!
        call getvid('ETAT_INIT', 'VITE', iocc=1, scal=champ, nbret=nvi)
        if (nvi .gt. 0) then
            call chpver('F', champ, 'NOEU', 'DEPL_R', ierr)
            inchac = 1
            cham2 = '&&COMDLT.VITINI'
            call vtcreb(cham2, 'V', 'R', nume_ddlz=numedd, nb_equa_outz=neq)
            call vtcopy(champ, cham2, iret)
            if (iret .ne. 0) then
                call utmess("F", "FIELD0_6", sk='VITE')
            end if
            call jeveuo(cham2//'.VALE', 'L', jvale)
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(jvale), b_incx, vitini, b_incy)
            vit = champ(1:8)
        else
            vit = 'nul'
        end if
        call utmess('I', 'DYNAMIQUE_79', nk=2, valk=[dep, vit])
!
!
        call getvid('ETAT_INIT', 'ACCE', iocc=1, scal=champ, nbret=nai)
        if (nai .gt. 0) then
            call chpver('F', champ, 'NOEU', 'DEPL_R', ierr)
            inchac = 0
            cham2 = '&&COMDLT.ACCINI'
            call vtcreb(cham2, 'V', 'R', nume_ddlz=numedd, nb_equa_outz=neq)
            call vtcopy(champ, cham2, iret)
            if (iret .ne. 0) then
                call utmess("F", "FIELD0_6", sk='ACCE')
            end if
            call jeveuo(cham2//'.VALE', 'L', jvale)
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(jvale), b_incx, accini, b_incy)
            call utmess('I', 'DYNAMIQUE_84', sk=champ)
        else
            call utmess('I', 'DYNAMIQUE_84', sk='calculee')
        end if
!
    end if
!
    if (ds_energy%l_comp) then
        call utmess('I', 'ETATINIT_5')
    end if
!
    call jedema()
end subroutine
