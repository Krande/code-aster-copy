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

subroutine vppost(vecrer, vecrei, vecrek, vecvp, nbpark, &
                  nbpari, nbparr, mxresf, nconv, nblagr, &
                  nfreqg, modes, typcon, compex, eigsol, &
                  matopa, matpsc, solveu, vecblo, veclag, &
                  flage, icom1, icom2, mpicou, mpicow, &
                  omemax, omemin, vpinf, vpmax, lcomod, mod45b)
!
! ROUTINE ORGANISANT LES POST-TRAITEMENTS ET NETTOYAGES GENERAUX D'OP0045.
! RQ. ON DETRUITS LES OBJETS GLOBAUX VECBLO, VECLAG ET, SUIVANT LES CAS, VECRIG, SUR BASE VOLATILE.
!     ILS SONT DETRUITS DANS SOLVEU, EIGSOL, VECBLO, VECLAG, VECRER, VECREI, VECREK ET VECVP.
! -------------------------------------------------------------------------------------------------
! person_in_charge: olivier.boiteau at edf.fr
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/detrsd.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/onerrf.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"
#include "asterfort/vpcntl.h"
#include "asterfort/vplecs.h"
#include "asterfort/vpmpi.h"
#include "asterfort/vppara.h"
#include "asterfort/vpwecf.h"
!
!
! --- INPUT
!
    integer(kind=8), intent(in) :: nbpark, nbpari, nbparr, mxresf, nconv, nblagr
    integer(kind=8), intent(in) :: nfreqg
    character(len=4), intent(in) :: mod45b
    character(len=8), intent(in) :: modes
    character(len=16), intent(in) :: typcon, compex
    character(len=19), intent(in) :: eigsol, matopa, matpsc, solveu
    character(len=24), intent(in) :: vecrer, vecrei, vecrek, vecvp, vecblo, veclag
    aster_logical, intent(in) :: flage
!
!
! --- OUTPUT
! None
!
! --- INPUT/OUTPUT
!
    mpi_int, intent(inout) :: mpicou, mpicow
    integer(kind=8), intent(inout) :: icom1, icom2
    real(kind=8), intent(inout) :: omemax, omemin, vpinf, vpmax
    aster_logical, intent(inout) :: lcomod
!
!
! --- VARIABLES LOCALES
!
    integer(kind=8) :: nbpara
    parameter(nbpara=27)
    integer(kind=8) :: nparr, ibid, nbrss, iret, lraide, lmasse, lamor, neq, lddl, lprod
    integer(kind=8) :: lmat(3), lmtpsc, lmatra, ierx, ifm, niv
    integer(kind=8) :: lresui, lresur, lresuk, lvec
    real(kind=8) :: omecor, rbid, precdc, precsh, seuil
    complex(kind=8) :: czero
    character(len=1) :: ktyp, k1blan, ctyp
    character(len=8) :: knega
    character(len=14) :: matra, matrc
    character(len=16) :: k16bid, optiof, optiov, stoper, sturm, typres
    character(len=19) :: amor, masse, raide, k19bid
    character(len=24) :: nopara(nbpara), valk(2)
    aster_logical :: lc, lkr, lns
!
!
! --- DECLARATION DES DATAS
!
!     ------------------------------------------------------------------
    data nopara/&
     &  'NUME_MODE', 'ITER_QR', 'ITER_BATHE',&
     &  'ITER_ARNO', 'ITER_JACOBI', 'ITER_SEPARE',&
     &  'ITER_AJUSTE', 'ITER_INVERSE',&
     &  'NORME', 'METHODE', 'TYPE_MODE',&
     &  'FREQ',&
     &  'OMEGA2', 'AMOR_REDUIT', 'ERREUR',&
     &  'MASS_GENE', 'RIGI_GENE', 'AMOR_GENE',&
     &  'MASS_EFFE_DX', 'MASS_EFFE_DY', 'MASS_EFFE_DZ',&
     &  'FACT_PARTICI_DX', 'FACT_PARTICI_DY', 'FACT_PARTICI_DZ',&
     &  'MASS_EFFE_UN_DX', 'MASS_EFFE_UN_DY', 'MASS_EFFE_UN_DZ'/
!     ------------------------------------------------------------------
!
! -----------------------
! --- CORPS DE LA ROUTINE
! -----------------------
!
!
! --  INITS.
    call jemarq()
    call infniv(ifm, niv)
    call jeveuo(vecrer, 'E', lresur)
    call jeveuo(vecrei, 'E', lresui)
    call jeveuo(vecrek, 'E', lresuk)
    call jeveuo(vecvp, 'E', lvec)
    call jeveuo(vecblo, 'L', lprod)
    call jeveuo(veclag, 'L', lddl)
    czero = dcmplx(0.d0, 0.d0)
    k1blan = ' '
    rbid = 0.d0
!
! --  LECTURE DES DONNEES DE EIGSOL
    call vplecs(eigsol, nbrss_=nbrss, omecor_=omecor, &
                precdc_=precdc, precsh_=precsh, seuil_=seuil, &
                matra_=matra, matrc_=matrc, &
                optiof_=optiof, stoper_=stoper, sturm_=sturm, &
                typres_=typres, amor_=amor, masse_=masse, raide_=raide, &
                lc_=lc, lkr_=lkr, lns_=lns)
    if (lkr) then
        ktyp = 'R'
    else
        ktyp = 'C'
    end if
!
! --  DESCRIPTEURS MATRICES
    call jeveuo(raide//'.&INT', 'E', lraide)
    neq = zi(lraide+2)
    call jeveuo(masse//'.&INT', 'E', lmasse)
    if (lc) then
        call jeveuo(amor//'.&INT', 'E', lamor)
    else
        lamor = 0
    end if
    call jeexin(matpsc//'.&INT', iret)
    if (iret .eq. 0) then
        lmtpsc = 0
    else
        call jeveuo(matpsc//'.&INT', 'E', lmtpsc)
    end if
    call jeexin(matopa//'.&INT', iret)
    if (iret .eq. 0) then
        lmatra = 0
    else
        call jeveuo(matopa//'.&INT', 'E', lmatra)
    end if
!
! -- CALCUL DES PARAMETRES GENERALISES
! --  CALCUL DE LA NORME D'ERREUR SUR LE MODE
! --  STOCKAGE DES VECTEURS PROPRES
! --  PARALLELISME MULTI-NIVEAUX STEP 3
    knega = 'NON'
    nparr = nbparr
    if (typcon(1:9) .eq. 'MODE_ACOU') nparr = 7
!
    if ((.not. lc) .and. lkr .and. (.not. lns)) then
        call vppara(modes, typcon, knega, lraide, lmasse, &
                    lamor, mxresf, neq, nconv, omecor, &
                    zi(lddl), zi(lprod), zr(lvec), [czero], nbpari, &
                    nparr, nbpark, nopara, mod45b, zi(lresui), &
                    zr(lresur), zk24(lresuk), ktyp, lcomod, icom1, &
                    icom2, typres, nfreqg)
    else
        if (lcomod) then
            ASSERT(.false.)
        end if
        call vppara(modes, typcon, knega, lraide, lmasse, &
                    lamor, mxresf, neq, nconv, omecor, &
                    zi(lddl), zi(lprod), [rbid], zc(lvec), nbpari, &
                    nparr, nbpark, nopara, mod45b, zi(lresui), &
                    zr(lresur), zk24(lresuk), ktyp, lcomod, ibid, &
                    ibid, k16bid, ibid)
    end if
!
!
! --  IMPRESSIONS LIEES A LA METHODE
    if ((mod45b(1:4) .eq. 'OP45') .or. (niv .ge. 2)) &
        call vpwecf(k1blan, typres, nconv, mxresf, zi(lresui), &
                    zr(lresur), zk24(lresuk), lamor, ktyp, lns)
    if (mod45b(1:4) .eq. 'OP45') call titre()
!
! --  CONTROLE DE VALIDITE DES MODES CALCULES
    if (sturm(1:3) .eq. 'NON') then
        optiov = ' '
    else
        optiov = optiof
        if (lc .or. (.not. lkr) .or. lns) then
! --  POUR DEBRANCHER LE TEST DE STURM DANS VPCNTL
            optiov = ' '
            valk(1) = matra
            valk(2) = matrc
            call utmess('I', 'ALGELINE2_73', nk=2, valk=valk)
        end if
    end if
!
    lmat(1) = lraide
    lmat(2) = lmasse
    lmat(3) = lmtpsc
! --  SI ON MANIPULE DEUX MATRICES DYNAMIQUES (MATOPA/MATPSC), ON SE DEBARASSE DE CELLE INUTILE
!     (MATRICE + FACTORISEE EVENTUELLE) ET DE SON EVENTUELLE OCCURENCE EXTERNE (MUMPS)
    if ((lmtpsc .ne. lmatra) .and. (lmatra .ne. 0)) then
        call detrsd('MATR_ASSE', matopa)
        call jedetr(matopa(1:19)//'.&INT')
        call jedetr(matopa(1:19)//'.&IN2')
    end if
!
! --  PARALLELISME MULTI-NIVEAUX STEP 4
    call vpmpi(4, lcomod_=lcomod, &
               mpicou_=mpicou, mpicow_=mpicow, &
               omemax_=omemax, omemin_=omemin, vpinf_=vpinf, vpmax_=vpmax)
    if (mod45b(1:4) .eq. 'OP45') then
!
! --  ON PASSE DANS LE MODE "VALIDATION DU CONCEPT EN CAS D'ERREUR"
        call onerrf('EXCEPTION+VALID', k16bid, ibid)
        if (stoper(1:3) .eq. 'OUI') then
            ctyp = 'E'
        else
            ctyp = 'A'
        end if
        call vpcntl(ctyp, modes, optiov, omemin, omemax, &
                    seuil, nconv, zi(lresui), lmat, omecor, &
                    precdc, ierx, vpinf, vpmax, zr(lresur), &
                    zr(lresur+3*mxresf), zr(lresur+mxresf), typres, nblagr, solveu, &
                    nbrss, precsh)
!
        if ((stoper(1:3) .eq. 'OUI') .and. (ierx .ne. 0)) &
            call utmess('Z', 'ALGELINE2_74', num_except=ASTER_SOLVER_ERROR)
!
        if (flage) call utmess('F', 'ALGELINE5_75')
!
!
! --  ON REMET LE MECANISME D'EXCEPTION A SA VALEUR INITIALE
        call onerrf(compex, k16bid, ibid)
    end if
!
! --  DESTRUCTION DE LA MATRICE DYNAMIQUE RESTANTE (VRAI MATPSC DISSOSIEE DE MATOPA OU
! --  MATPSC POINTANT SUR MATOPA D'OU LA RECONSTRUCTION DE NOM CI-DESSOUS
    if (lmtpsc .ne. 0) then
        k19bid = zk24(zi(lmtpsc+1)) (1:19)
        call detrsd('MATR_ASSE', k19bid)
        call jedetr(k19bid//'.&INT')
        call jedetr(k19bid//'.&IN2')
    end if
!
! -- NETTOYAGE DES OBJETS JEVEUX GLOBAUX DE L'OPERATEUR DE LA BASE VOLATILE
    if (mod45b(1:4) .eq. 'OP45') then
        call detrsd('SOLVEUR', solveu)
        call detrsd('EIGENSOLVER', eigsol)
        call jedetr(vecblo)
        call jedetr(veclag)
        call jedetr(vecrer)
        call jedetr(vecrei)
        call jedetr(vecrek)
        call jedetr(vecvp)
    end if

    call jedema()
!
!     FIN DE VPPOST
!
end subroutine
