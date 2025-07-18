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

subroutine vpcalq(eigsol, vecrer, vecrei, vecrek, vecvp, &
                  mxresf, neqact, nblagr, omemax, omemin, &
                  omeshi, vecblo, sigma, npivot, flage, &
                  nconv, vpinf, vpmax)
!
! ROUTINE EFFECTUANT LE CALCUL MODAL PARAMETRE DANS EIGSOL PAR LA METHODE QZ
! -------------------------------------------------------------------------------------------------
! person_in_charge: olivier.boiteau at edf.fr
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8depi.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/freqom.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jerazo.h"
#include "asterfort/rectfc.h"
#include "asterfort/rectfr.h"
#include "asterfort/utmess.h"
#include "asterfort/vpbosc.h"
#include "asterfort/vpbost.h"
#include "asterfort/vpecri.h"
#include "asterfort/vplecs.h"
#include "asterfort/vpordi.h"
#include "asterfort/vpordo.h"
#include "asterfort/vpqzla.h"
#include "asterfort/wkvect.h"
#include "asterfort/wp4vec.h"
!
! --- INPUT
!
    integer(kind=8), intent(in) :: mxresf, neqact, nblagr
    real(kind=8), intent(in) :: omemax, omemin, omeshi
    complex(kind=8), intent(in) :: sigma
    character(len=19), intent(in) :: eigsol
    character(len=24), intent(in) :: vecrer, vecrei, vecrek, vecvp, vecblo
!
! --- OUTPUT
!
    integer(kind=8), intent(out) :: nconv
    real(kind=8), intent(out) :: vpinf, vpmax
    aster_logical, intent(out) :: flage
!
! --- INPUT/OUTPUT
!
    integer(kind=8), intent(inout) :: npivot
!
! --- VARIABLES LOCALES
!
    integer(kind=8) :: imet, lamor, lmasse, lraide, nbvect, neq, nfreq, lprod
    integer(kind=8) :: qrn, qrlwor, qrn2, ilscal, irscal, icscal, ivscal, iiscal
    integer(kind=8) :: lvalpr, iqrn, lqrn, qrar, qrai, qrba, qrvl, kqrn
    integer(kind=8) :: lauc, kqrnr, mfreq, ifreq, izero
    integer(kind=8) :: lresui, lresur, lresuk, lvec
    real(kind=8) :: quapi2, omecor, precdc, rbid, rzero
    character(len=1) :: ktyp
    character(len=8) :: method
    character(len=16) :: optiof, typeqz, typres
    character(len=19) :: amor, masse, raide, numedd
    character(len=24) :: k24bid
    aster_logical :: lc, lkr, lns, lpg
    logical(kind=4), pointer :: bwork(:) => null()
!
! -----------------------
! --- CORPS DE LA ROUTINE
! -----------------------
!
!
! ---  INITS.
!
    call jemarq()
    quapi2 = r8depi()*r8depi()
    izero = 0
    rzero = 0.d0
!
    nconv = 0
    flage = .false.
    call jeveuo(vecrer, 'E', lresur)
    call jeveuo(vecrei, 'E', lresui)
    call jeveuo(vecrek, 'E', lresuk)
    call jeveuo(vecvp, 'E', lvec)
    call jeveuo(vecblo, 'L', lprod)
!
! --- LECTURE DES DONNEES DE EIGSOL
!
    call vplecs(eigsol, nbvect_=nbvect, nfreq_=nfreq, &
                omecor_=omecor, precdc_=precdc, &
                method_=method, optiof_=optiof, typeqz_=typeqz, &
                typres_=typres, amor_=amor, masse_=masse, raide_=raide, &
                lc_=lc, lkr_=lkr, lns_=lns, lpg_=lpg)
    ASSERT(method(1:2) .eq. 'QZ')
    if (lkr) then
        ktyp = 'R'
    else
        ktyp = 'C'
    end if
!
! ---  DESCRIPTEURS MATRICES
    call jeveuo(raide//'.&INT', 'E', lraide)
    neq = zi(lraide+2)
    numedd = ''
    call dismoi('NOM_NUME_DDL', raide, 'MATR_ASSE', repk=numedd)
    call jeveuo(masse//'.&INT', 'E', lmasse)
    if (lc) then
        call jeveuo(amor//'.&INT', 'E', lamor)
    else
        lamor = 0
    end if
!
! --- PRE-ALLOCATIONS MEMOIRE
!
    qrn = nbvect
    qrlwor = 8*qrn
    qrn2 = qrn*qrn
    if (typeqz(1:7) .eq. 'QZ_EQUI') then
        call wkvect('&&VPCALQ.QRLSCALE.WORK', 'V V R', qrn, ilscal)
        call wkvect('&&VPCALQ.QRRSCALE.WORK', 'V V R', qrn, irscal)
        call wkvect('&&VPCALQ.QRRCONDE.WORK', 'V V R', qrn, icscal)
        call wkvect('&&VPCALQ.QRRCONDV.WORK', 'V V R', qrn, ivscal)
        call wkvect('&&VPCALQ.QRI.WORK', 'V V S', qrn+6, iiscal)
        allocate (bwork(qrn))
    end if
    if (lkr .and. (.not. lc) .and. (.not. lns)) then
        call wkvect('&&VPCALQ.QZ.VALPRO', 'V V R', qrn, lvalpr)
        call wkvect('&&VPCALQ.QZ.MATRICEK', 'V V R', qrn2, iqrn)
        call wkvect('&&VPCALQ.QZ.MATRICEM', 'V V R', qrn2, lqrn)
        call wkvect('&&VPCALQ.QZ.ALPHAR', 'V V R', qrn, qrar)
        call wkvect('&&VPCALQ.QZ.ALPHAI', 'V V R', qrn, qrai)
        call wkvect('&&VPCALQ.QZ.BETA', 'V V R', qrn, qrba)
        call wkvect('&&VPCALQ.QZ.VL', 'V V R', qrn, qrvl)
        call wkvect('&&VPCALQ.QZ.WORK', 'V V R', qrlwor, kqrn)
    else
        if (lc) call wkvect('&&VPCALQ.VECT.AUC', 'V V C', qrn2, lauc)
        call wkvect('&&VPCALQ.QZ.VALPRO', 'V V C', qrn, lvalpr)
        call wkvect('&&VPCALQ.QZ.MATRICEK', 'V V C', qrn2, iqrn)
        call wkvect('&&VPCALQ.QZ.MATRICEM', 'V V C', qrn2, lqrn)
        call wkvect('&&VPCALQ.QZ.ALPHA', 'V V C', qrn, qrar)
        call wkvect('&&VPCALQ.QZ.BETA', 'V V C', qrn, qrba)
        call wkvect('&&VPCALQ.QZ.VL', 'V V C', qrn, qrvl)
        call wkvect('&&VPCALQ.QZ.WORK', 'V V C', qrlwor, kqrn)
        call wkvect('&&VPCALQ.QZ.WORKR', 'V V R', qrlwor, kqrnr)
    end if
    call jerazo('&&VPCALQ.QZ.MATRICEK', qrn2, 1)
    call jerazo('&&VPCALQ.QZ.MATRICEM', qrn2, 1)
!
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
        if (lkr .and. (.not. lns)) then
!     ------------------------------------------------------------------
!     -------  QZ PB GENERALISE REEL SYMETRIQUE  -----------------------
!     ------------------------------------------------------------------
            call vpqzla(typeqz, qrn, iqrn, lqrn, qrar, &
                        qrai, qrba, qrvl, lvec, kqrn, &
                        lvalpr, nconv, omecor, ktyp, kqrnr, &
                        neqact, ilscal, irscal, optiof, omemin, &
                        omemax, omeshi, zi(lprod), nfreq, lmasse, &
                        lraide, lamor, numedd, sigma, icscal, &
                        ivscal, iiscal, bwork, flage)
            call rectfr(nconv, nconv, omeshi, npivot, nblagr, &
                        zr(lvalpr), nfreq, zi(lresui), zr(lresur), mxresf)
            call vpbost(typres, nconv, nconv, omeshi, zr(lvalpr), &
                        nfreq, vpinf, vpmax, precdc, method, &
                        omecor)
            if (typres(1:9) .eq. 'DYNAMIQUE') call vpordi(1, 0, nconv, zr(lresur+mxresf), &
                                                          zr(lvec), neq, zi(lresui))
            do imet = 1, nconv
                zi(lresui-1+mxresf+imet) = izero
                zr(lresur-1+imet) = freqom(zr(lresur-1+mxresf+imet))
!           SI OPTION 'PLUS_GRANDE' : CONVERSION EN VALEUR PHYSIQUE
                if (lpg) zr(lresur-1+imet) = +1.d0/(quapi2*zr(lresur-1+imet))
                zr(lresur-1+2*mxresf+imet) = rzero
                zk24(lresuk-1+mxresf+imet) = typeqz
            end do
            if (typres(1:9) .ne. 'DYNAMIQUE') then
                call vpordo(0, 0, nconv, zr(lresur+mxresf), zr(lvec), &
                            neq)
                do imet = 1, nconv
                    zr(lresur-1+imet) = freqom(zr(lresur-1+mxresf+imet))
                    zi(lresui-1+imet) = imet
                end do
            end if
!
        else if ((.not. lkr) .or. lns) then
!     ------------------------------------------------------------------
!     -------  QZ PB GENERALISE COMPLEXE OU REEL NON SYM  --------------
!     ------------------------------------------------------------------
            call vpqzla(typeqz, qrn, iqrn, lqrn, qrar, &
                        qrai, qrba, qrvl, lvec, kqrn, &
                        lvalpr, nconv, omecor, ktyp, kqrnr, &
                        neqact, ilscal, irscal, optiof, omemin, &
                        omemax, omeshi, zi(lprod), nfreq, lmasse, &
                        lraide, lamor, numedd, sigma, icscal, &
                        ivscal, iiscal, bwork, flage)
            npivot = nblagr
            call rectfc(nconv, nconv, sigma, npivot, nblagr, &
                        zc(lvalpr), nfreq, zi(lresui), zr(lresur), nfreq)
            call vpbosc(typres, nconv, nconv, sigma, zc(lvalpr), &
                        nfreq, vpinf, vpmax, precdc, method, &
                        omecor)
            do imet = 1, nconv
                zi(lresui-1+mxresf+imet) = izero
                zr(lresur-1+imet) = freqom(zr(lresur-1+mxresf+imet))
                zk24(lresuk-1+mxresf+imet) = typeqz
            end do
        end if
!
    else
!
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     ---------------------  PROBLEME QUADRATIQUE   --------------------
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!
!
!     ------------------------------------------------------------------
!     -------  QZ PB QUADRATIQUE REEL ET COMPLEXE, SYM OU NON  --------
!     ------------------------------------------------------------------
        call vpqzla(typeqz, qrn, iqrn, lqrn, qrar, &
                    qrai, qrba, qrvl, lvec, kqrn, &
                    lvalpr, nconv, omecor, ktyp, kqrnr, &
                    neqact, ilscal, irscal, optiof, omemin, &
                    omemax, omeshi, zi(lprod), nfreq, lmasse, &
                    lraide, lamor, numedd, sigma, icscal, &
                    ivscal, iiscal, bwork, flage)
        nfreq = nfreq/2
        call wp4vec(nfreq, nconv, neq, sigma, zc(lvalpr), &
                    zc(lvec), mxresf, zi(lresui), zr(lresur), zi(lprod), &
                    zc(lauc), omecor)
        do imet = 1, nfreq
            zi(lresui-1+mxresf+imet) = izero
            zr(lresur-1+imet) = freqom(zr(lresur-1+mxresf+imet))
            zk24(lresuk-1+mxresf+imet) = typeqz
        end do
!
    end if
!
! ---- NOMBRE DE MODES CONVERGES
! ---- SI LE SOLVEUR MODAL A BIEN ACHEVE SON TRAVAIL ON FAIT CETTE AFFEC
! ---- TATION SINON ON NE TIENT COMPTE QUE DES NCONV MODES REELLEMENT CV
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
        do ifreq = mfreq-1, 0, -1
            if ((zr(lresur+mxresf+ifreq) .gt. omemax) .or. (zr(lresur+mxresf+ifreq) .lt. omemin)) &
                nconv = nconv-1
        end do
        if (mfreq .ne. nconv) call utmess('I', 'ALGELINE2_17')
    end if
!
! ---  ON MODIFIE LES VALEURS NFREQ DE LA SD EIGENSOLVER
    call vpecri(eigsol, 'I', 1, k24bid, rbid, &
                nfreq)
!
!
! --- NETTOYAGE OBJETS TEMPORAIRES
!
    if (typeqz(1:7) .eq. 'QZ_EQUI') then
        call jedetr('&&VPCALQ.QRLSCALE.WORK')
        call jedetr('&&VPCALQ.QRRSCALE.WORK')
        call jedetr('&&VPCALQ.QRRCONDE.WORK')
        call jedetr('&&VPCALQ.QRRCONDV.WORK')
        call jedetr('&&VPCALQ.QRI.WORK')
        deallocate (bwork)
    end if
    if (lkr .and. (.not. lc) .and. (.not. lns)) then
        call jedetr('&&VPCALQ.QZ.VALPRO')
        call jedetr('&&VPCALQ.QZ.MATRICEK')
        call jedetr('&&VPCALQ.QZ.MATRICEM')
        call jedetr('&&VPCALQ.QZ.ALPHAR')
        call jedetr('&&VPCALQ.QZ.ALPHAI')
        call jedetr('&&VPCALQ.QZ.BETA')
        call jedetr('&&VPCALQ.QZ.VL')
        call jedetr('&&VPCALQ.QZ.WORK')
    else
        if (lc) call jedetr('&&VPCALQ.VECT.AUC')
        call jedetr('&&VPCALQ.QZ.VALPRO')
        call jedetr('&&VPCALQ.QZ.MATRICEK')
        call jedetr('&&VPCALQ.QZ.MATRICEM')
        call jedetr('&&VPCALQ.QZ.ALPHA')
        call jedetr('&&VPCALQ.QZ.BETA')
        call jedetr('&&VPCALQ.QZ.VL')
        call jedetr('&&VPCALQ.QZ.WORK')
        call jedetr('&&VPCALQ.QZ.WORKR')
    end if
    call jedema()
!
!     FIN DE VPCALQ
!
end subroutine
