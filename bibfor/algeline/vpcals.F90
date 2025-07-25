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

subroutine vpcals(eigsol, vecrer, vecrei, vecrek, vecvp, &
                  matopa, mxresf, neqact, nblagr, omemax, &
                  omemin, omeshi, solveu, vecblo, veclag, &
                  sigma, npivot, flage, nconv, vpinf, vpmax, mod45b, &
                  vecstb, vecedd, nbddl, vecsdd, nbddl2, csta)
!
! ROUTINE EFFECTUANT LE CALCUL MODAL PARAMETRE DANS EIGSOL PAR LA METHODE DE SORENSEN
! -------------------------------------------------------------------------------------------------
! person_in_charge: olivier.boiteau at edf.fr
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8depi.h"
#include "asterfort/assert.h"
#include "asterfort/infniv.h"
#include "asterfort/freqom.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rectfc.h"
#include "asterfort/rectfr.h"
#include "asterfort/vpbosc.h"
#include "asterfort/vpbost.h"
#include "asterfort/vpecri.h"
#include "asterfort/vplecs.h"
#include "asterfort/vpordi.h"
#include "asterfort/vpordo.h"
#include "asterfort/vpsorc.h"
#include "asterfort/vpsorn.h"
#include "asterfort/vpsor1.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/wpsorc.h"
#include "asterfort/wpsorn.h"
#include "asterfort/wp3vec.h"
#include "asterfort/wp4vec.h"
#include "asterfort/wp5vec.h"
!
! --- INPUT
!
    integer(kind=8), intent(in) :: mxresf, neqact, nblagr, nbddl, nbddl2
    real(kind=8), intent(in)    :: omeshi
    complex(kind=8), intent(in) :: sigma
    character(len=4), intent(in) :: mod45b
    character(len=19), intent(in) :: eigsol, matopa, solveu
    character(len=24), intent(in) :: vecrer, vecrei, vecrek, vecvp, vecblo, veclag
    character(len=24), intent(in) :: vecstb, vecedd, vecsdd
!
! --- OUTPUT
!
    integer(kind=8), intent(out) :: nconv
    real(kind=8), intent(out) :: vpinf, vpmax, csta
    aster_logical, intent(out) :: flage
!
!
! --- INPUT/OUTPUT
!
    integer(kind=8), intent(inout) :: npivot
    real(kind=8), intent(inout) :: omemax, omemin
!
!
! --- VARIABLES LOCALES
!
    integer(kind=8) :: imet, lamor, lmasse, lmatra, lraide, maxitr, nbvect, neq, nfreq
    integer(kind=8) :: lonwl, lselec, lresid, lworkd, lworkl, lworkv, ldsor, laux, lworkr
    integer(kind=8) :: lauc, laur, laul, ldiagr, lsurdr, lprod, lddl, eddl, eddl2
    integer(kind=8) :: nfreq1, izero, mfreq, ifreq, ifm, niv, priram(8)
    integer(kind=8) :: lresui, lresur, lresuk, lvec, redem, jstab
    real(kind=8) :: alpha, quapi2, omecor, precdc, precsh, rbid, rzero, tolsor
    character(len=1) :: appr
    character(len=8) :: method
    character(len=16) :: optiof, stoper, typres
    character(len=19) :: amor, masse, raide
    character(len=24) :: kmetho, k24bid
    aster_logical :: lc, lkr, lns, lpg
!
! -----------------------
! --- CORPS DE LA ROUTINE
! -----------------------
!
!
! ---  INITS.
!
    call jemarq()
    call infniv(ifm, niv)
    if (niv .eq. 2) then
        priram(1) = 2
        priram(2) = 2
        priram(3) = 2
        priram(4) = 2
        priram(5) = 0
        priram(6) = 0
        priram(7) = 0
        priram(8) = 2
    else
        priram(:) = 0
    end if
    quapi2 = r8depi()*r8depi()
    izero = 0
    rzero = 0.d0
    kmetho = 'SORENSEN'
    nconv = 0
    flage = .false.
    call jeveuo(vecrer, 'E', lresur)
    call jeveuo(vecrei, 'E', lresui)
    call jeveuo(vecrek, 'E', lresuk)
    call jeveuo(vecvp, 'E', lvec)
    call jeveuo(vecblo, 'L', lprod)
    call jeveuo(veclag, 'L', lddl)
    if (mod45b(1:4) .EQ. 'SOR1') then
        call jeveuo(vecstb, 'E', jstab)
        call jeveuo(vecedd, 'L', eddl)
        call jeveuo(vecsdd, 'L', eddl2)
    end if
!
! --- LECTURE DES DONNEES DE EIGSOL
!
    call vplecs(eigsol, maxitr_=maxitr, &
                nbvect_=nbvect, nfreq_=nfreq, &
                alpha_=alpha, omecor_=omecor, &
                precdc_=precdc, precsh_=precsh, &
                tolsor_=tolsor, appr_=appr, &
                method_=method, optiof_=optiof, stoper_=stoper, &
                typres_=typres, amor_=amor, masse_=masse, raide_=raide, &
                lc_=lc, lkr_=lkr, lns_=lns, lpg_=lpg)
    ASSERT(method(1:8) .eq. 'SORENSEN')
!
! ---  DESCRIPTEURS MATRICES
    call jeveuo(raide//'.&INT', 'E', lraide)
    neq = zi(lraide+2)
    call jeveuo(masse//'.&INT', 'E', lmasse)
    if (lc) then
        call jeveuo(amor//'.&INT', 'E', lamor)
    else
        lamor = 0
    end if
    call jeveuo(matopa//'.&INT', 'E', lmatra)
!
! --- PRE-ALLOCATIONS MEMOIRE
!
    lonwl = 3*nbvect**2+6*nbvect
    call wkvect('&&VPCALS.SELECT', 'V V L', nbvect, lselec)
!     ------------------------------------------------------------------
!     -------  SORENSEN PB GENERALISE REEL SYMETRIQUE  --------
!     ------------------------------------------------------------------
    if ((.not. lc) .and. lkr .and. (.not. lns)) then
        call wkvect('&&VPCALS.RESID', 'V V R', neq, lresid)
        call wkvect('&&VPCALS.VECT.WORKD', 'V V R', 3*neq, lworkd)
        call wkvect('&&VPCALS.VECT.WORKL', 'V V R', lonwl, lworkl)
        call wkvect('&&VPCALS.VECT.WORKV', 'V V R', 3*nbvect, lworkv)
        call wkvect('&&VPCALS.VAL.PRO', 'V V R', 2*(nfreq+1), ldsor)
        call wkvect('&&VPCALS.VECT.AUX', 'V V R', neq, laux)
!     --------------------------------- ---------------------------------
!     -------  SORENSEN PB GENERALISE COMPLEXE OU REEL NON SYM  --------
!     ------------------------------------------------------------------
    else if ((.not. lc) .and. (lns .or. (.not. lkr))) then
        call wkvect('&&VPCALS.RESID', 'V V C', neq, lresid)
        call wkvect('&&VPCALS.VECT.WORKD', 'V V C', 3*neq, lworkd)
        call wkvect('&&VPCALS.VECT.WORKL', 'V V C', lonwl, lworkl)
        call wkvect('&&VPCALS.VECT.WORKV', 'V V C', 3*nbvect, lworkv)
        call wkvect('&&VPCALS.VAL.PRO', 'V V C', (nfreq+1), ldsor)
        call wkvect('&&VPCALS.VECT.AUX', 'V V C', neq, laux)
        call wkvect('&&VPCALS.VECT.AUR', 'V V R', nbvect, lworkr)
!     ------------------------------------------------------------------
!     -------  SORENSEN PB QUADRATIQUE REEL  SYM  ----------------------
!     -------  APPROCHE REELLE OU IMAGINAIRE      ----------------------
!     ------------------------------------------------------------------
    else if ((lkr .and. lc) .and. (appr .ne. 'C')) then
        call wkvect('&&VPCALS.RESID', 'V V R', 2*neq, lresid)
        call wkvect('&&VPCALS.VECT.WORKD', 'V V R', 6*neq, lworkd)
        call wkvect('&&VPCALS.VECT.AUX', 'V V R', 2*neq, laux)
        call wkvect('&&VPCALS.VECT.AUC', 'V V C', 2*neq*(nbvect+1), lauc)
        call wkvect('&&VPCALS.VECT.AUR', 'V V R', 2*neq*(nbvect+1), laur)
        call wkvect('&&VPCALS.VECT.AUL', 'V V C', neq*(nbvect+1), laul)
        call wkvect('&&VPCALS.VAL.PR', 'V V R', nbvect+1, ldiagr)
        call wkvect('&&VPCALS.VAL.PI', 'V V R', nbvect+1, lsurdr)
        call wkvect('&&VPCALS.VAL.PRO', 'V V R', 2*(nfreq+1), ldsor)
        call wkvect('&&VPCALS.VECT.WORKL', 'V V R', lonwl, lworkl)
        call wkvect('&&VPCALS.VECT.WORKV', 'V V R', 3*nbvect, lworkv)
!     ------------------------------------------------------------------
!     -------  SORENSEN PB QUADRATIQUE REEL, SYM OU NON  ---------------
!     -------  APPROCHE COMPLEXE                         ---------------
!     ------------------------------------------------------------------
    else if ((lkr .and. lc) .and. (appr .eq. 'C')) then
        call wkvect('&&VPCALS.RESID', 'V V C', 2*neq, lresid)
        call wkvect('&&VPCALS.VECT.WORKD', 'V V C', 6*neq, lworkd)
        call wkvect('&&VPCALS.VECT.AUX', 'V V C', 2*neq, laux)
        call wkvect('&&VPCALS.VECT.AUC', 'V V C', 2*neq*(nbvect+1), lauc)
        call wkvect('&&VPCALS.VECT.AUR', 'V V R', 2*neq*(nbvect+1), laur)
        call wkvect('&&VPCALS.VECT.WORKL', 'V V C', lonwl, lworkl)
        call wkvect('&&VPCALS.VECT.WORKV', 'V V C', 3*nbvect, lworkv)
        call wkvect('&&VPCALS.VAL.PRO', 'V V C', 2*(nfreq+1), ldsor)
!     ------------------------------------------------------------------
!     -------  SORENSEN PB QUADRATIQUE COMPLEXE SYM  -------------------
!     -------  APPROCHE COMPLEXE                     -------------------
!     ------------------------------------------------------------------
    else if ((.not. lkr) .and. lc) then
        call wkvect('&&VPCALS.RESID', 'V V C', 2*neq, lresid)
        call wkvect('&&VPCALS.VECT.WORKD', 'V V C', 6*neq, lworkd)
        call wkvect('&&VPCALS.VECT.WORKL', 'V V C', lonwl, lworkl)
        call wkvect('&&VPCALS.VECT.WORKV', 'V V C', 3*nbvect, lworkv)
        call wkvect('&&VPCALS.VAL.PRO', 'V V C', 2*(nfreq+1), ldsor)
        call wkvect('&&VPCALS.VECT.AUX', 'V V C', 2*neq, laux)
        call wkvect('&&VPCALS.VECT.AUC', 'V V C', 2*neq*(nbvect+1), lauc)
        call wkvect('&&VPCALS.VECT.AUR', 'V V R', 2*neq*(nbvect+1), laur)
    else
        ASSERT(.false.)
    end if
!
! --- CAS PARTICULIER DE NMOP45 + FLAMBEMENT SANS MATRICE GEOMETRIQUE
!
    if (mod45b(1:4) .EQ. 'SOR1') then
        redem = 0
        call vpsor1(lmatra, neq, nbvect, nfreq, &
                    tolsor, zr(lvec), zr(lresid), zr(lworkd), zr(lworkl), &
                    lonwl, zl(lselec), zr(ldsor), omeshi, zr(laux), &
                    zr(lworkv), zi(lprod), zi(lddl), zi(eddl), nbddl, neqact, maxitr, &
                    ifm, niv, priram, alpha, omecor, &
                    nconv, flage, solveu, &
                    nbddl2, zi(eddl2), zr(jstab), csta, redem)
        call rectfr(nconv, nconv, omeshi, npivot, nblagr, &
                    zr(ldsor), nfreq+1, zi(lresui), zr(lresur), mxresf)
        nfreq1 = nfreq+1
        call vpbost(typres, nconv, nconv, omeshi, zr(ldsor), &
                    nfreq1, vpinf, vpmax, precdc, method, omecor)
        if (typres(1:9) .eq. 'DYNAMIQUE') &
            call vpordi(1, 0, nconv, zr(lresur+mxresf), zr(lvec), neq, zi(lresui))
!
        do imet = 1, nconv
            zi(lresui-1+mxresf+imet) = izero
            zr(lresur-1+imet) = freqom(zr(lresur-1+mxresf+imet))
!       SI OPTION 'PLUS_GRANDE' : CONVERSION EN VALEUR PHYSIQUE
            if (lpg) zr(lresur-1+imet) = +1.d0/(quapi2*zr(lresur-1+imet))
            zr(lresur-1+2*mxresf+imet) = rzero
            zk24(lresuk-1+mxresf+imet) = kmetho
        end do
        if (typres(1:9) .ne. 'DYNAMIQUE') then
            call vpordo(0, 0, nconv, zr(lresur+mxresf), zr(lvec), neq)
            do imet = 1, nconv
                zr(lresur-1+imet) = freqom(zr(lresur-1+mxresf+imet))
                zi(lresui-1+imet) = imet
            end do
        end if
        goto 99
    end if
!
! --- CALCUL MODAL PROPREMENT DIT
!
    if (.not. lc) then
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     ---------------------  PROBLEME GENERALISE   ---------------------
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!
        if (lkr .and. (.not. lns)) then
!     ------------------------------------------------------------------
!     -------  SORENSEN PB GENERALISE REEL SYMETRIQUE  --------
!     ------------------------------------------------------------------
            call vpsorn(lmasse, lmatra, neq, nbvect, nfreq, &
                        tolsor, zr(lvec), zr(lresid), zr(lworkd), zr(lworkl), &
                        lonwl, zl(lselec), zr(ldsor), omeshi, zr(laux), &
                        zr(lworkv), zi(lprod), zi(lddl), neqact, maxitr, &
                        ifm, niv, priram, alpha, omecor, &
                        nconv, flage, solveu)
            call rectfr(nconv, nconv, omeshi, npivot, nblagr, &
                        zr(ldsor), nfreq+1, zi(lresui), zr(lresur), mxresf)
            nfreq1 = nfreq+1
            call vpbost(typres, nconv, nconv, omeshi, zr(ldsor), &
                        nfreq1, vpinf, vpmax, precdc, method, &
                        omecor)
            if (typres(1:9) .eq. 'DYNAMIQUE') &
                call vpordi(1, 0, nconv, zr(lresur+mxresf), zr(lvec), neq, zi(lresui))
!
            do imet = 1, nconv
                zi(lresui-1+mxresf+imet) = izero
                zr(lresur-1+imet) = freqom(zr(lresur-1+mxresf+imet))
!           SI OPTION 'PLUS_GRANDE' : CONVERSION EN VALEUR PHYSIQUE
                if (lpg) zr(lresur-1+imet) = +1.d0/(quapi2*zr(lresur-1+imet))
                zr(lresur-1+2*mxresf+imet) = rzero
                zk24(lresuk-1+mxresf+imet) = kmetho
            end do
            if (typres(1:9) .ne. 'DYNAMIQUE') then
                call vpordo(0, 0, nconv, zr(lresur+mxresf), zr(lvec), neq)
                do imet = 1, nconv
                    zr(lresur-1+imet) = freqom(zr(lresur-1+mxresf+imet))
                    zi(lresui-1+imet) = imet
                end do
            end if
!
        else if (lns .or. (.not. lkr)) then
!     ------------------------------------------------------------------
!     -------  SORENSEN PB GENERALISE COMPLEXE OU REEL NON SYM  --------
!     ------------------------------------------------------------------
            call vpsorc(lmasse, lmatra, neq, nbvect, nfreq, &
                        tolsor, zc(lvec), zc(lresid), zc(lworkd), zc(lworkl), &
                        lonwl, zl(lselec), zc(ldsor), sigma, zc(laux), &
                        zc(lworkv), zr(lworkr), zi(lprod), zi(lddl), neqact, &
                        maxitr, ifm, niv, priram, abs(alpha), &
                        nconv, flage, solveu)
            npivot = nblagr
            nfreq1 = nfreq+1
            call rectfc(nconv, nconv, sigma, npivot, nblagr, &
                        zc(ldsor), nfreq1, zi(lresui), zr(lresur), nfreq)
            call vpbosc(typres, nconv, nconv, sigma, zc(ldsor), &
                        nfreq1, vpinf, vpmax, precdc, method, &
                        omecor)
            do imet = 1, nconv
                zi(lresui-1+mxresf+imet) = izero
                zr(lresur-1+imet) = freqom(zr(lresur-1+mxresf+imet))
                zk24(lresuk-1+mxresf+imet) = kmetho
            end do
!
        else
            ASSERT(.false.)
        end if
!
    else
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     ---------------------  PROBLEME QUADRATIQUE   --------------------
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!
        if (lkr) then
            select case (appr)
            case ('R', 'I')
!     ------------------------------------------------------------------
!     -------  SORENSEN PB QUADRATIQUE REEL  SYM  ----------------------
!     -------  APPROCHE REELLE OU IMAGINAIRE      ----------------------
!     ------------------------------------------------------------------
                call wpsorn(appr, lmasse, lamor, lmatra, neq, &
                            nbvect, nfreq, tolsor, zc(lvec), zr(lresid), &
                            zr(lworkd), zr(lworkl), lonwl, zl(lselec), zr(ldsor), &
                            zr(lsurdr), zr(ldiagr), sigma, zr(laux), zr(lworkv), &
                            zi(lprod), zi(lddl), neqact, maxitr, ifm, &
                            niv, priram, abs(alpha), nconv, flage, &
                            zr(laur), zc(lauc), zc(laul), solveu)
                nfreq = nconv/2
                call wp3vec(appr, optiof, nfreq, nconv, neq, &
                            sigma, zr(lsurdr), zr(ldiagr), zc(lvec), mxresf, &
                            zi(lresui), zr(lresur), zi(lprod), zc(lauc), omecor)
            case ('C')
!     ------------------------------------------------------------------
!     -------  SORENSEN PB QUADRATIQUE REEL, SYM OU NON  ---------------
!     -------  APPROCHE COMPLEXE                         ---------------
!     ------------------------------------------------------------------
                call wpsorc(lmasse, lamor, lmatra, neq, nbvect, &
                            nfreq, tolsor, zc(lvec), zc(lresid), zc(lworkd), &
                            zc(lworkl), lonwl, zl(lselec), zc(ldsor), sigma, &
                            zc(laux), zc(lworkv), zi(lprod), zi(lddl), neqact, &
                            maxitr, ifm, niv, priram, abs(alpha), &
                            nconv, flage, zc(lauc), zr(laur), solveu)
                nfreq = nconv/2
                call wp4vec(nfreq, nconv, neq, sigma, zc(ldsor), &
                            zc(lvec), mxresf, zi(lresui), zr(lresur), zi(lprod), &
                            zc(lauc), omecor)
            case default
                ASSERT(.false.)
            end select
            do imet = 1, nfreq
                zi(lresui-1+mxresf+imet) = izero
                zr(lresur-1+imet) = freqom(zr(lresur-1+mxresf+imet))
                zk24(lresuk-1+mxresf+imet) = kmetho
            end do
!
        else
!     ------------------------------------------------------------------
!     -------  SORENSEN PB QUADRATIQUE COMPLEXE SYM  -------------------
!     -------  APPROCHE COMPLEXE                     -------------------
!     ------------------------------------------------------------------
            if (lns) then
                ASSERT(.false.)
            end if
            call wpsorc(lmasse, lamor, lmatra, neq, nbvect, &
                        nfreq, tolsor, zc(lvec), zc(lresid), zc(lworkd), &
                        zc(lworkl), lonwl, zl(lselec), zc(ldsor), sigma, &
                        zc(laux), zc(lworkv), zi(lprod), zi(lddl), neqact, &
                        maxitr, ifm, niv, priram, abs(alpha), &
                        nconv, flage, zc(lauc), zr(laur), solveu)
            nfreq = nconv/2
            call wp5vec(nfreq, nconv, neq, zc(ldsor), zc(lvec), &
                        mxresf, zi(lresui), zr(lresur), zc(lauc))
            do imet = 1, nfreq
                zi(lresui-1+mxresf+imet) = izero
                zr(lresur-1+imet) = freqom(zr(lresur-1+mxresf+imet))
                zk24(lresuk-1+mxresf+imet) = kmetho
            end do
        end if
    end if
!
! ---- NOMBRE DE MODES CONVERGES
! ---- SI LE SOLVEUR MODAL A BIEN ACHEVE SON TRAVAIL ON FAIT CETTE AFFEC
! ---- TATION SINON ON NE TIENT COMPTE QUE DES NCONV MODES REELLEMENT CV
99  continue
    if (.not. flage) nconv = nfreq
!
!     ------------------------------------------------------------------
!     -------------------- CORRECTION : OPTION BANDE -------------------
!     ------------------------------------------------------------------
!
! --- SI OPTION BANDE ON NE GARDE QUE LES FREQUENCES DANS LA BANDE
    mfreq = nconv
    if (optiof(1:5) .eq. 'BANDE') then
        if (lc .or. lns .or. .not. lkr) then
            ASSERT(.false.)
        end if
        if ((typres(1:9) .ne. 'DYNAMIQUE') .and. (mod45b(1:4) .ne. 'OP45')) then
            rbid = omemin
            omemin = -omemax
            omemax = -rbid
        end if
        do ifreq = mfreq-1, 0, -1
            if ((zr(lresur+mxresf+ifreq) .gt. omemax) .or. (zr(lresur+mxresf+ifreq) .lt. omemin)) &
                nconv = nconv-1
        end do
        if ((typres(1:9) .ne. 'DYNAMIQUE') .and. (mod45b(1:4) .ne. 'OP45')) then
            rbid = omemin
            omemin = -omemax
            omemax = -rbid
        end if
        if (mfreq .ne. nconv) call utmess('I', 'ALGELINE2_17')
    end if
!
! ---  ON AJUSTE LA VALEUR NFREQ DE LA SD EIGENSOLVER
    if (mod45b(1:4) .eq. 'OP45') then
        call vpecri(eigsol, 'I', 1, k24bid, rbid, nfreq)
    end if

!
!
! --- NETTOYAGE OBJETS TEMPORAIRES
!
    call jedetr('&&VPCALS.SELECT')
    if ((.not. lc) .and. lkr .and. (.not. lns)) then
        call jedetr('&&VPCALS.RESID')
        call jedetr('&&VPCALS.VECT.WORKD')
        call jedetr('&&VPCALS.VECT.WORKL')
        call jedetr('&&VPCALS.VECT.WORKV')
        call jedetr('&&VPCALS.VAL.PRO')
        call jedetr('&&VPCALS.VECT.AUX')
    else if ((.not. lc) .and. (lns .or. (.not. lkr))) then
        call jedetr('&&VPCALS.RESID')
        call jedetr('&&VPCALS.VECT.WORKD')
        call jedetr('&&VPCALS.VECT.WORKL')
        call jedetr('&&VPCALS.VECT.WORKV')
        call jedetr('&&VPCALS.VAL.PRO')
        call jedetr('&&VPCALS.VECT.AUX')
        call jedetr('&&VPCALS.VECT.AUR')
    else if ((lkr .and. lc) .and. (appr .ne. 'C')) then
        call jedetr('&&VPCALS.RESID')
        call jedetr('&&VPCALS.VECT.WORKD')
        call jedetr('&&VPCALS.VECT.AUX')
        call jedetr('&&VPCALS.VECT.AUC')
        call jedetr('&&VPCALS.VECT.AUR')
        call jedetr('&&VPCALS.VECT.AUL')
        call jedetr('&&VPCALS.VAL.PR')
        call jedetr('&&VPCALS.VAL.PI')
        call jedetr('&&VPCALS.VAL.PRO')
        call jedetr('&&VPCALS.VECT.WORKL')
        call jedetr('&&VPCALS.VECT.WORKV')
    else if ((lkr .and. lc) .and. (appr .eq. 'C')) then
        call jedetr('&&VPCALS.RESID')
        call jedetr('&&VPCALS.VECT.WORKD')
        call jedetr('&&VPCALS.VECT.AUX')
        call jedetr('&&VPCALS.VECT.AUC')
        call jedetr('&&VPCALS.VECT.AUR')
        call jedetr('&&VPCALS.VECT.WORKL')
        call jedetr('&&VPCALS.VECT.WORKV')
        call jedetr('&&VPCALS.VAL.PRO')
    else if ((.not. lkr) .and. lc) then
        call jedetr('&&VPCALS.RESID')
        call jedetr('&&VPCALS.VECT.WORKD')
        call jedetr('&&VPCALS.VECT.WORKL')
        call jedetr('&&VPCALS.VECT.WORKV')
        call jedetr('&&VPCALS.VAL.PRO')
        call jedetr('&&VPCALS.VECT.AUX')
        call jedetr('&&VPCALS.VECT.AUC')
        call jedetr('&&VPCALS.VECT.AUR')
    end if
    call jedema()
!
!     FIN DE VPCALS
!
end subroutine
