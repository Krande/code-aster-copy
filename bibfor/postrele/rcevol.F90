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
subroutine rcevol(typtab, nommat, symax, nbopt, option, lsymm)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8vide.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rcev22.h"
#include "asterfort/rcevfa.h"
#include "asterfort/rcevo0.h"
#include "asterfort/rcevo1.h"
#include "asterfort/rcevo2.h"
#include "asterfort/rcevoa.h"
#include "asterfort/rcevod.h"
#include "asterfort/rcevom.h"
#include "asterfort/rcevse.h"
#include "asterfort/rcevsn.h"
#include "asterfort/rcevsp.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: nbopt
    real(kind=8) :: symax
    character(len=8) :: nommat
    character(len=16) :: typtab, option(*)
    aster_logical, intent(in) :: lsymm
!
!     OPERATEUR POST_RCCM:    TYPE_ANALYSE = 'COMPOSANT'
!                           TYPE_RESU_MECA = 'EVOLUTION'
!
!     ------------------------------------------------------------------
!
    integer(kind=8) :: i, j, n1, nbinti, jinti, nbtran
    real(kind=8) :: para(3), sm
    aster_logical :: lpmpb, lsn, lfatig, flexio, lrocht, lamorc, kemixt
    character(len=8) :: typeke
    character(len=16) :: kinti
    character(len=24) :: cinst, csili, csiex, csno, csne, csneo, csnee, cspo
    character(len=24) :: cspe, cfao, cfae, cnoc, cresu, cresp, intitu, cspto
    character(len=24) :: cspte, cspmo, cspme, cstex, csmex
! DEB ------------------------------------------------------------------
!
! --- VECTEUR DES INSTANTS DEMANDES
!
    cinst = '&&RCEVOL.INSTANTS'
!
! --- VECTEUR DES TRANSITOIRES
!
    cresu = '&&RCEVOL.RESU_MECA'
    cresp = '&&RCEVOL.RESU_PRES'
!
! --- VECTEUR DES CONTRAINTES LINEARISEES AUX EXTREMITES (PMPB, SN)
!
    csili = '&&RCEVOL.SIGM_LINE'
!
! --- VECTEUR DES CONTRAINTES TOTALES AUX EXTREMITES (SP)
!
    csiex = '&&RCEVOL.SIGM_EXTR'
!
! --- VECTEUR DES CONTRAINTES THERMIQUES AUX EXTREMITES (SPTH)
!
    cstex = '&&RCEVOL.STHE_EXTR'
!
! --- VECTEUR DES CONTRAINTES MECANIQUES AUX EXTREMITES (SPMECA)
!
    csmex = '&&RCEVOL.SMEC_EXTR'
!
! --- VECTEUR DES NB_OCCUR
!
    cnoc = '&&RCEVOL.NB_OCCUR'
!
! --- CALCUL DE GRANDEURS A L'ORIGINE ET A L'EXTREMITE
!
    csno = '&&RCEVOL.CALCUL_SN .ORIG'
    csne = '&&RCEVOL.CALCUL_SN .EXTR'
    csneo = '&&RCEVOL.CALCUL_SNE.ORIG'
    csnee = '&&RCEVOL.CALCUL_SNE.EXTR'
    cspo = '&&RCEVOL.CALCUL_SP .ORIG'
    cspe = '&&RCEVOL.CALCUL_SP .EXTR'
    cspto = '&&RCEVOL.CALCUL_SPT.ORIG'
    cspme = '&&RCEVOL.CALCUL_SPT.EXTR'
    cspmo = '&&RCEVOL.CALCUL_SPM.ORIG'
    cspte = '&&RCEVOL.CALCUL_SPM.EXTR'
    cfao = '&&RCEVOL.FATIGUE   .ORIG'
    cfae = '&&RCEVOL.FATIGUE   .EXTR'
!
!     ------------------------------------------------------------------
!                             LES OPTIONS
!     ------------------------------------------------------------------
    lfatig = .false.
    lsn = .false.
    lpmpb = .false.
    flexio = .false.
    lrocht = .false.
    lamorc = .false.
!
    if (lsymm) then
        call utmess('I', 'POSTRCCM_58')
    end if
!
    do i = 1, nbopt
        if (option(i) .eq. 'PM_PB') then
            lpmpb = .true.
        else if (option(i) .eq. 'SN') then
            lsn = .true.
        else if (option(i) .eq. 'FATIGUE_ZH210') then
            lfatig = .true.
            lsn = .true.
        else if (option(i) .eq. 'AMORCAGE') then
            lamorc = .true.
        end if
    end do
!
    if (lamorc .and. (lpmpb .or. lsn .or. lfatig)) then
        call utmess('F', 'POSTRCCM_3')
    end if
!
    kemixt = .false.
    call getvtx(' ', 'TYPE_KE', scal=typeke, nbret=n1)
    if (typeke .eq. 'KE_MIXTE') kemixt = .true.
!
!     ------------------------------------------------------------------
!                      TRAITEMENT DE L'AMORCAGE
!     ------------------------------------------------------------------
!
    if (lamorc) then
        call rcevoa(typtab, nommat)
        goto 999
    end if
!
!     ------------------------------------------------------------------
!                            LE MATERIAU
!     ------------------------------------------------------------------
!
    call rcevo1(nommat, lfatig, sm, para, symax)
!
!     ------------------------------------------------------------------
!                     NOMBRE DE LIGNE A "POST_RCCM"
!     ------------------------------------------------------------------
!
    intitu = '&&RCEVOL.INTITULE'
    call rcevo0(intitu, nbinti, lsn, lfatig, nbtran)
    call jeveuo(intitu, 'L', jinti)
!
!     ------------------------------------------------------------------
!
    do i = 1, nbinti
!
        do j = 1, nbtran
!
            kinti = zk16(jinti-1+nbtran*(i-1)+j)
!
!         --------------------------------------------------------------
!                    TRAITEMENT DU MOT CLE FACTEUR TRANSITOIRE
!         --------------------------------------------------------------
!
            if (lsn .and. .not. lfatig .and. nbtran .gt. 1) then
                call rcev22(nbinti, kinti, j, csili, cinst, &
                            csiex, lfatig, flexio, lrocht, cnoc, &
                            cresu, cresp, lsymm)
            else
                call rcevo2(nbinti, kinti, csili, cinst, csiex, &
                            kemixt, cstex, csmex, lfatig, flexio, &
                            lrocht, cnoc, cresu, cresp, lsymm)
            end if
!
            if (lrocht .and. symax .eq. r8vide()) then
                call utmess('A', 'POSTRCCM_4')
                lrocht = .false.
            end if
!
!         --------------------------------------------------------------
!                          TRAITEMENT DES OPTIONS
!         --------------------------------------------------------------
!
            if (lsn) call rcevsn(csili, cinst, csno, csne, lsymm)
!
            if (flexio) call rcevse(csili, cinst, csneo, csnee, lsymm)
!
            if (lfatig) then
                call rcevsp(csiex, kemixt, cstex, csmex, cinst, &
                            cspo, cspe, cspto, cspte, cspmo, &
                            cspme)
                call rcevfa(nommat, para, sm, cnoc, csno, &
                            csne, cspo, cspe, kemixt, cspto, &
                            cspte, cspmo, cspme, cfao, cfae)
            end if
!
!         --------------------------------------------------------------
!                                 ARCHIVAGE
!         --------------------------------------------------------------
!
            if (typtab .eq. 'VALE_MAX') then
!
                call rcevom(csili, cinst, cnoc, sm, lfatig, &
                            lpmpb, lsn, csno, csne, flexio, &
                            csneo, csnee, cfao, cfae, cspo, &
                            cspe, cresu, kinti, i, j, &
                            lrocht, symax, cresp, kemixt, cspto, &
                            cspte, cspmo, cspme, lsymm, csiex)
!
            else
!
                call rcevod(csili, cinst, cnoc, sm, lfatig, &
                            lpmpb, lsn, csno, csne, flexio, &
                            csneo, csnee, cfao, cfae, cspo, &
                            cspe, cresu, kinti, i, j, &
                            lrocht, symax, cresp, kemixt, cspto, &
                            cspte, cspmo, cspme, lsymm, csiex)
!
            end if
!
            call jedetr(cinst)
            call jedetr(cresu)
            call jedetr(cresp)
            call jedetr(csili)
            call jedetr(csiex)
            call jedetr(cnoc)
            call jedetr(csno)
            call jedetr(csne)
            call jedetr(csneo)
            call jedetr(csnee)
            call jedetr(cspo)
            call jedetr(cspe)
            call jedetr(cfao)
            call jedetr(cfae)
            if (kemixt) then
                call jedetr(cstex)
                call jedetr(csmex)
                call jedetr(cspmo)
                call jedetr(cspme)
                call jedetr(cspto)
                call jedetr(cspte)
            end if
!
        end do
!
    end do
!
    call jedetr(intitu)
!
999 continue
!
end subroutine
