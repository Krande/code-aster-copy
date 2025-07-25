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

subroutine dtmprep_noli(sd_dtm_)
    implicit none
!
! person_in_charge: hassan.berro at edf.fr
!
! dtmprep_noli : Retreives information regarding 9 types of localized
!            nonlinearities for a transient DYNA_VIBRA calculation
!            on reduced basis (TRAN//GENE)
!   --------------------------------------------------------------------------------------
!       (1)     Stops (chocs)               / CHOC
!       (2)     Anti sismic devices         / ANTI_SISM
!       (3)     Viscous dampers             / DIS_VISC
!       (4)     Nonlinear springs           / DIS_ECRO_TRAC
!       (5)     Buckling                    / FLAMBAGE
!       (6)     Cracked rotor               / ROTOR_FISS
!       (7)     F(V) relationship           / RELA_EFFO_VITE
!       (8)     F(X) relationship           / RELA_EFFO_DEPL
!       (9)     Elastic nonlinear springs   / CHOC_ELAS_TRAC
!   --------------------------------------------------------------------------------------
!
!   Note : Information about these 6 nonlinearity types are read using mdchoc
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/assert.h"
#include "asterfort/dtmcase_coder.h"
#include "asterfort/dtmget.h"
#include "asterfort/dtminivec.h"
#include "asterfort/nlinivec.h"
#include "asterfort/dtmprep_noli_choc.h"
#include "asterfort/dtmprep_noli_flam.h"
#include "asterfort/dtmprep_noli_ants.h"
#include "asterfort/dtmprep_noli_decr.h"
#include "asterfort/dtmprep_noli_dvis.h"
#include "asterfort/dtmprep_noli_rede.h"
#include "asterfort/dtmprep_noli_revi.h"
#include "asterfort/dtmprep_noli_rotf.h"
#include "asterfort/dtmprep_verichoc.h"
#include "asterfort/dtmprep_noli_galet.h"
#include "asterfort/dtmsav.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/nlget.h"
#include "asterfort/nlsav.h"
#include "asterfort/nltype.h"
#include "asterfort/nlvint.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
!   -0.1- Input/output arguments
    character(len=*), intent(in) :: sd_dtm_
!
!   -0.2- Local variables
    integer(kind=8)          :: nbchoc, nbnli, iret, i, nbmode
    integer(kind=8)          :: nbvint, j, nlcase, nbcomp, icomp
    integer(kind=8)          :: nltype_i, ivchoc
    character(len=7) :: casek7
    character(len=8) :: sd_dtm, monmot, sd_nl
    character(len=16):: nltreat_k, nltype_k
    character(len=19):: nomres
!
    real(kind=8), pointer :: basev0(:) => null()
    real(kind=8), pointer :: fext_tmp(:) => null()
!
#define base0(row,col) basev0((col-1)*nbmode+row)
!
!
!   -0.3- Initializations
    call jemarq()
    sd_dtm = sd_dtm_
!
    call getfac('COMPORTEMENT', nbcomp)
    if (nbcomp .eq. 0) then
        call dtmsav(sd_dtm, _NB_NONLI, 1, iscal=0)
        call dtmsav(sd_dtm, _NL_TREAT, 1, iscal=0)
        goto 999
    end if
!
    sd_nl = '&&OP29NL'
    call dtmsav(sd_dtm, _SD_NONL, 1, kscal=sd_nl)
    call nlsav(sd_nl, _MAX_LEVEL, 1, iscal=0)
    call nlsav(sd_nl, _NB_CHOC, 1, iscal=0)
    call nlsav(sd_nl, _NB_FLAMB, 1, iscal=0)
    call nlsav(sd_nl, _NB_ANTSI, 1, iscal=0)
    call nlsav(sd_nl, _NB_DIS_VISC, 1, iscal=0)
    call nlsav(sd_nl, _NB_DIS_ECRO_TRAC, 1, iscal=0)
    call nlsav(sd_nl, _NB_R_FIS, 1, iscal=0)
    call nlsav(sd_nl, _NB_PALIE, 1, iscal=0)
    call nlsav(sd_nl, _NB_REL_FX, 1, iscal=0)
    call nlsav(sd_nl, _NB_REL_FV, 1, iscal=0)
    call nlsav(sd_nl, _NB_DIS_CHOC_ELAS, 1, iscal=0)
!
    do icomp = 1, nbcomp

        call getvtx('COMPORTEMENT', 'RELATION', iocc=icomp, scal=nltype_k)
        do nltype_i = 1, _NL_NB_TYPES
            if (nltype_k .eq. nltype(nltype_i)) goto 5
        end do
        ASSERT(.false.)
5       continue
!
        select case (nltype_i)
!
        case (NL_CHOC)
            call dtmprep_noli_choc(sd_dtm, sd_nl, icomp)
!
        case (NL_BUCKLING)
            call dtmprep_noli_flam(sd_dtm, sd_nl, icomp)
!
        case (NL_ANTI_SISMIC)
            call dtmprep_noli_ants(sd_dtm, sd_nl, icomp)
!
        case (NL_DIS_VISC)
            call dtmprep_noli_dvis(sd_dtm, sd_nl, icomp)
!
        case (NL_DIS_ECRO_TRAC)
            call dtmprep_noli_decr(sd_dtm, sd_nl, icomp)
!
        case (NL_CRACKED_ROTOR)
            call dtmprep_noli_rotf(sd_dtm, sd_nl, icomp)
!
        case (NL_FX_RELATIONSHIP)
            call dtmprep_noli_rede(sd_dtm, sd_nl, icomp)
!
        case (NL_FV_RELATIONSHIP)
            call dtmprep_noli_revi(sd_dtm, sd_nl, icomp)
!
        case (NL_DIS_CHOC_ELAS)
            call dtmprep_noli_galet(sd_dtm, sd_nl, icomp)
!
        case default
            ASSERT(.false.)
        end select
    end do

    call nlvint(sd_nl)
    ! call utimsd(6, 2, .false._1, .true._1, sd_nl, 1, 'V')

    ! ASSERT(.false.)

    call dtmget(sd_dtm, _NB_NONLI, iscal=nbnli)
    call dtmsav(sd_dtm, _NL_TREAT, 1, iscal=0)
    if (nbnli .ne. 0) then

        call dtmget(sd_dtm, _NB_MODES, iscal=nbmode)
        call nlget(sd_nl, _NB_CHOC, iscal=nbchoc)

!       --- Explicit or implicit treatment of choc non-linearities
        if (nbchoc .gt. 0) then
            call getvtx(' ', 'TRAITEMENT_NONL', iocc=1, scal=nltreat_k)
            if (nltreat_k(1:9) .eq. 'IMPLICITE') then

                if (nbchoc .gt. 41) call utmess('F', 'DYNAMIQUE_28', si=41)

                call dtmsav(sd_dtm, _NL_TREAT, 1, iscal=1)
                call dtminivec(sd_dtm, _F_NL_ADD, nbmode)
                call dtminivec(sd_dtm, _IMP_DEPL, nbmode)
                call dtminivec(sd_dtm, _IMP_VITE, nbmode)
                call dtminivec(sd_dtm, _IMP_ACCE, nbmode)
                call dtminivec(sd_dtm, _IMP_FEXT, nbmode)

                nlcase = 0
                call dtmcase_coder(nlcase, casek7)
                call wkvect(sd_dtm//'.PRJ_BAS.'//casek7, 'V V R', nbmode*nbmode, vr=basev0)
                do i = 1, nbmode
                    base0(i, i) = 1.d0
                    do j = i+1, nbmode
                        base0(i, j) = 0.d0
                    end do
                end do

                call nlget(sd_nl, _INTERNAL_VARS, lonvec=nbvint)
                call dtminivec(sd_dtm, _NL_SAVE0, nbvint)
            end if
        end if

!
        call dtmget(sd_dtm, _MULTI_AP, kscal=monmot)
        if (monmot(1:3) .eq. 'OUI') then
            call dtmget(sd_dtm, _CALC_SD, kscal=nomres)
            call jeexin(nomres//'.IPSD', iret)
            if (iret .eq. 0) then
                ASSERT(.false.)
            end if
        end if

        call getfac('VERI_CHOC', ivchoc)
        if (ivchoc .ne. 0) then
            call dtmprep_verichoc(sd_dtm, sd_nl)
        end if

    end if

    call nlinivec(sd_nl, _FEXT_MPI, nbmode, vr=fext_tmp)

999 continue
    call jedema()
end subroutine
