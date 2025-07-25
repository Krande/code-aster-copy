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

subroutine nmiclg(fami, kpg, ksp, option, rela_comp, &
                  imate, epsm, deps, sigm, vim, &
                  sigp, vip, dsde, carcri, codret)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/lcimpl.h"
#include "asterfort/nm1das.h"
#include "asterfort/nm1dci.h"
#include "asterfort/nm1dco.h"
#include "asterfort/nm1dis.h"
#include "asterfort/nm1dpm.h"
#include "asterfort/nmmaba.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
#include "asterfort/verift.h"
#include "asterfort/Behaviour_type.h"
    integer(kind=8) :: imate, kpg, ksp, codret
!
    character(len=16) :: rela_comp, option
    character(len=*) :: fami
    real(kind=8) :: carcri(CARCRI_SIZE)
    real(kind=8) :: vim(*)
    real(kind=8) :: vip(*)
    real(kind=8) :: sigy, sigm, deps, sigp
    integer(kind=8) :: codres(1)
    real(kind=8) :: dsde, epsm
! -------------------------------------------------------------------
!
!    TRAITEMENT DE LA RELATION DE COMPORTEMENT -ELASTOPLASTICITE-
!    ECROUISSAGE ISOTROPE ET CINEMATIQUE- LINEAIRE - VON MISES-
!    POUR UN MODELE CABLE_GAINE (EQUIVALENT A MECA_BARRE)
!
! -------------------------------------------------------------------
! IN  : FAMI  : FAMILLE DE POINT DE GAUSS
!       KPG   : NUMERO DU POINT DE GAUSS DANS LA FAMILLE
!       KSP   : NUMERO DU SOUS POINT
!       OPTION : OPTION DEMANDEE (R_M_T,FULL OU RAPH_MECA)
!       IMATE : POINTEUR MATERIAU CODE
!       EPSM  : DEFORMATION A L'INSTANT MOINS
!       DEPS  : INCREMENT DE DEFORMATION
!       SIGM  : CONTRAINTE A L'INSTANT MOINS
!       VIM   : VARIABLE INTERNE A L'INSTANT MOINS
! OUT : SIGP  : CONTRAINTE A L'INSTANT ACTUEL
!       VIP    : VARIABLE INTERNE A L'INSTANT ACTUEL
!       DSDE   : MATRICE TANGENTE
!       CRILDC : CRITERE LOI DE COMPORTEMENT
!       CODRET : CODE RETOUR
!
!----------VARIABLES LOCALES
!
    integer(kind=8) :: ncstpm
    parameter(ncstpm=13)
    character(len=8) :: nomasl(4)
    real(kind=8) :: cstpm(ncstpm)
    real(kind=8) :: depsth, depsm, tmoins, tplus
    real(kind=8) :: em, ep, dsdem, dsdep
    real(kind=8) :: valres(4), syc, etc, syt, ett, cr, val(1)
    aster_logical :: isot, cine, elas, corr, implex, isotli, pinto, asyml, sans
    data nomasl/'SY_C', 'DC_SIGM_', 'SY_T', 'DT_SIGM_'/
!
!
!----------INITIALISATIONS
!
    elas = .false.
    isot = .false.
    cine = .false.
    corr = .false.
    implex = option(1:16) .eq. 'RIGI_MECA_IMPLEX'
    isotli = .false.
    pinto = .false.
    asyml = .false.
    if (rela_comp .eq. 'ELAS') then
        elas = .true.
    else if ((rela_comp .eq. 'VMIS_ISOT_LINE') .or. (rela_comp .eq. 'VMIS_ISOT_TRAC')) then
        isot = .true.
        if (rela_comp .eq. 'VMIS_ISOT_LINE') then
            isotli = .true.
        end if
    else if (rela_comp .eq. 'VMIS_CINE_LINE') then
        cine = .true.
    else if (rela_comp .eq. 'CORR_ACIER') then
        corr = .true.
    else if (rela_comp .eq. 'PINTO_MENEGOTTO') then
        pinto = .true.
    else if (rela_comp .eq. 'VMIS_ASYM_LINE') then
        asyml = .true.
    else if (rela_comp .eq. 'SANS') then
        sans = .true.
    end if
    if (implex) then
        if ((.not. elas) .and. (.not. isotli)) then
            call utmess('F', 'POUTRE0_49', sk=rela_comp)
        end if
    end if
!
!
! --- CARACTERISTIQUES ELASTIQUES A TMOINS
!
    call rcvalb(fami, kpg, ksp, '-', imate, &
                ' ', 'ELAS', 0, ' ', [0.d0], &
                1, 'E', val, codres, 1)
    em = val(1)
!
! --- CARACTERISTIQUES ELASTIQUES A TPLUS
!
    call rcvalb(fami, kpg, ksp, '+', imate, &
                ' ', 'ELAS', 0, ' ', [0.d0], &
                1, 'E', val, codres, 1)
    ep = val(1)
!
!
    if (isot .and. (.not. implex)) then
        call verift(fami, kpg, ksp, 'T', imate, &
                    epsth_=depsth)
        depsm = deps-depsth
        call nm1dis(fami, kpg, ksp, imate, em, &
                    ep, sigm, depsm, vim, option, &
                    rela_comp, ' ', sigp, vip, dsde)
    else if (cine) then
        call verift(fami, kpg, ksp, 'T', imate, &
                    epsth_=depsth)
        depsm = deps-depsth
        call nm1dci(fami, kpg, ksp, imate, em, &
                    ep, sigm, depsm, vim, option, &
                    ' ', sigp, vip, dsde)
    else if (elas) then
        dsde = ep
        vip(1) = 0.d0
        call verift(fami, kpg, ksp, 'T', imate, &
                    epsth_=depsth)
        sigp = ep*(sigm/em+deps-depsth)
    else if (corr) then
        call nm1dco(fami, kpg, ksp, option, imate, &
                    ' ', ep, sigm, epsm, deps, &
                    vim, sigp, vip, dsde, carcri, &
                    codret)
    else if (implex) then
        call lcimpl(fami, kpg, ksp, imate, em, &
                    ep, sigm, tmoins, tplus, deps, &
                    vim, option, sigp, vip, dsde)
    else if (asyml) then
        call nmmaba(imate, rela_comp, ep, dsde, sigy, &
                    ncstpm, cstpm)
        call rcvalb(fami, 1, 1, '+', imate, &
                    ' ', 'ECRO_ASYM_LINE', 0, ' ', [0.d0], &
                    4, nomasl, valres, codres, 1)
        syc = valres(1)
        etc = valres(2)
        syt = valres(3)
        ett = valres(4)
        cr = 0.d0
        call nm1das(fami, kpg, ksp, ep, syc, &
                    syt, etc, ett, cr, tmoins, &
                    tplus, imate, sigm, deps, vim, &
                    sigp, vip, dsdem, dsdep)
        if (option(1:10) .eq. 'RIGI_MECA_' .or. option(1:9) .eq. 'FULL_MECA') then
            if (option(11:14) .eq. 'ELAS') then
                dsde = ep
            else
                if (option(1:14) .eq. 'RIGI_MECA_TANG') then
                    dsde = dsdem
                else
                    dsde = dsdep
                end if
            end if
        end if
    else if (pinto) then
        call verift(fami, kpg, ksp, 'T', imate, &
                    epsth_=depsth)
        depsm = deps-depsth
        call nmmaba(imate, rela_comp, ep, dsde, sigy, &
                    ncstpm, cstpm)
        call nm1dpm(fami, kpg, ksp, imate, option, &
                    8, ncstpm, cstpm, sigm, vim, &
                    depsm, vip, sigp, dsde)
    else if (sans) then
        sigp = 0.d0
        dsde = 0.d0
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
