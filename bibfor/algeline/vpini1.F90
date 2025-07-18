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

subroutine vpini1(eigsol, modes, solveu, typcon, vecblo, &
                  veclag, vecrig, matpsc, matopa, iretr, &
                  nblagr, neqact, npivot, nstoc, omemax, &
                  omemin, omeshi, sigma, mod45, nfreq_calibr_)
! PREPARATION DES ASPECTS NUMERIQUES POUR MODE_ITER_SIMULT: SOLVEUR LINEAIRE, LAGRANGE ET MODES
! RIGIDES, BORNES DE TRAVAIL EFFECTIVES, CALCUL DU NOMBRE DE MODES, FACTO. DE LA MATRICE SHIFTEE,
! DETERMINATION DE LA TAILLE DE L'ESPACE DE PROJECTION.
! RQ1. CODE RETOUR IRETR:
! = -1: PB DS VPDDL
! = -2: PB BANDE VIDE
! = -3: NMOP45, OPTION CALIBRATION, ON NE VEUT QUE LE NBRE DE MODES
! RQ2. ON MODIFIE LES VALEURS NFREQ/NBVECT DE LA SD EIGENSOLVER VIA VPECRI.
! RQ3. ON CREE LES OBJETS GLOBAUX VECBLO, VECLAG ET, SUIVANT LES CAS, VECRIG, SUR BASE VOLATILE.
!      ILS SONT DETRUITS DANS VPPOST POUR LES PREMIERS ET VPCALT POUR LE DERNIER.
! -------------------------------------------------------------------------------------------------
! person_in_charge: olivier.boiteau at edf.fr
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/isnnem.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mtdefs.h"
#include "asterfort/mtdscr.h"
#include "asterfort/omega2.h"
#include "asterfort/rscrsd.h"
#include "asterfort/tldlg2.h"
#include "asterfort/utmess.h"
#include "asterfort/uttcpu.h"
#include "asterfort/vecint.h"
#include "asterfort/vpddl.h"
#include "asterfort/vpecri.h"
#include "asterfort/vpfopc.h"
#include "asterfort/vpfopr.h"
#include "asterfort/vplecs.h"
#include "asterfort/vptabl.h"
#include "asterfort/wkvect.h"
#include "asterfort/wpfopc.h"
#include "asterfort/wpfopr.h"
!
! --- INPUT
!
    character(len=4), intent(in) :: mod45
    character(len=8), intent(in) :: modes
    character(len=16), intent(in) :: typcon
    character(len=19), intent(in) :: eigsol, solveu
    character(len=24), intent(in) :: vecblo, veclag, vecrig
!
! --- OUTPUT
!
    integer(kind=8), intent(out) :: iretr, nblagr, neqact, npivot
    integer(kind=8), intent(out) :: nstoc
    integer(kind=8), optional, intent(out) :: nfreq_calibr_
    real(kind=8), intent(out) :: omemax, omemin, omeshi
    complex(kind=8), intent(out) :: sigma
!
! --- INPUT/OUTPUT
!
    character(len=19), intent(inout) :: matpsc, matopa
!
! --- VARIABLES LOCALES
!
    integer(kind=8) :: ibid, ibid2(2), icoef, ifm, iret, islvi, islvk, krefa, lamor, npiv2(2)
    integer(kind=8) :: lmasse, lmatra, lmtpsc, lraide, nbborn, nbrss, nbvec2, nbvect, nfreq, niv
    integer(kind=8) :: lddl, lprod, nprec, neq, indf
    real(kind=8) :: eps, fmax, fmin, effmax, effmin, freq1, freq2, omecor, precdc, precsh
    real(kind=8) :: rbid, rbid2(2), rtest, rzero
    character(len=1) :: appr
    character(len=8) :: arret, method, chaine
    character(len=9) :: typevp
    character(len=14) :: matra, matrc
    character(len=16) :: modrig, optiof, optiov, typeqz, typres, typcal
    character(len=19) :: amor, masse, raide, tabmod
    character(len=24) :: metres, valk(2), k24bid
    aster_logical :: lc, lkr, lns, lqz, ltabmo
!
! -----------------------
! --- CORPS DE LA ROUTINE
! -----------------------
!
! --  INITS.
    call jemarq()
    call infniv(ifm, niv)
    rzero = 0.d0
    eps = 1.d+4*r8prem()
    iretr = 0
    indf = isnnem()
!
! --  LECTURE DES PARAMETRES MODAUX
    call vplecs(eigsol, nbborn_=nbborn, nbvec2_=nbvec2, &
                nbvect_=nbvect, nbrss_=nbrss, nfreq_=nfreq, &
                omecor_=omecor, freq1_=freq1, freq2_=freq2, &
                precdc_=precdc, precsh_=precsh, &
                appr_=appr, arret_=arret, &
                method_=method, typevp_=typevp, matra_=matra, matrc_=matrc, &
                modrig_=modrig, optiof_=optiof, typcal_=typcal, typeqz_=typeqz, &
                typres_=typres, amor_=amor, masse_=masse, raide_=raide, tabmod_=tabmod, &
                lc_=lc, lkr_=lkr, lns_=lns, lqz_=lqz)
!
! --  DESCRIPTEURS MATRICES
    call mtdscr(raide)
    call jeveuo(raide//'.&INT', 'E', lraide)
    neq = zi(lraide+2)
    call mtdscr(masse)
    call jeveuo(masse//'.&INT', 'E', lmasse)
    if (lc) then
        call mtdscr(amor)
        call jeveuo(amor//'.&INT', 'E', lamor)
    else
        lamor = 0
    end if
    call wkvect(veclag, 'V V I', neq, lddl)
    call wkvect(vecblo, 'V V I', neq, lprod)
!
! --  TRAITEMENTS LAGRANGES
    iretr = 0
    nblagr = -999
    neqact = -999
    call vecint(neq, indf, zi(lddl))
    call vecint(neq, indf, zi(lprod))
    call vpddl(raide, masse, neq, nblagr, ibid, &
               neqact, zi(lddl), zi(lprod), iret)
    if (iret .ne. 0) then
        iretr = -1
        goto 999
    end if
    if (lqz) then
        if (optiof(1:4) .eq. 'TOUT') nfreq = neqact
        if ((typeqz(1:5) .eq. 'QZ_QR') .and. ((nblagr .ne. 0) .or. lns)) then
            valk(1) = matra
            valk(2) = matrc
            call utmess('F', 'ALGELINE5_60', nk=2, valk=valk)
        end if
    end if
!
!
! -- TRAITEMENTS SOLVEUR LINEAIRE
! -- LECTURE DES PARAMETRES SOLVEURS LINEAIRES ET CREATION DE
!    LA SD SOLVEUR ASSOCIEE. CETTE SD SOLVEUR EST LOCALE A L'OPERATEUR. POUR CE CALCUL, C'EST ELLE
!    QUI EST UTILISEE POUR PARAMETREE LE SOLVEUR LINEAIRE, ET NON PAS LA SD SOLVEUR CREE PAR LA
!    CMDE ECLATEE NUME_DDL LORS DE LA CONSTITUTION DES MATRICES.
    call jeveuo(solveu//'.SLVK', 'L', islvk)
    call jeveuo(solveu//'.SLVI', 'E', islvi)
    nprec = zi(islvi)
    metres = zk24(islvk)
    if ((metres(1:4) .ne. 'LDLT') .and. (metres(1:10) .ne. 'MULT_FRONT') .and. &
        (metres(1:5) .ne. 'MUMPS')) call utmess('F', 'ALGELINE5_71')
!
! --  SI ON A BESOIN DE FACTORISER SIMULTANEMENT DEUX MATRICES AVEC LE SOLVEUR MUMPS ON LUI
!     SIGNALE AFIN QU'IL OPTIMISE AU MIEUX LA MEMOIRE POUR CHACUNES D'ELLES.
!     CE N'EST VRAIMENT UTILE QUE SI SOLVEUR/GESTION_MEMOIRE='AUTO'.
! --  CF COMMENTAIRES PLUS LOIN SUR MATOPA/MATPSC.
    if (metres(1:5) .eq. 'MUMPS') then
        if ((lc) .and. (lkr) .and. (appr .eq. 'R')) then
            if (zi(islvi-1+6) .lt. 0) then
! --  PB INITIALISATION DE LA SD_SOLVEUR
                ASSERT(.false.)
            else
                zi(islvi-1+6) = 2
            end if
        end if
    end if
!
! --  DETECTION DES MODES DE CORPS RIGIDE
    nstoc = 0
    if (modrig(1:11) .eq. 'MODE_RIGIDE') then
        call uttcpu('CPU.RESO.1', 'DEBUT', ' ')
        call uttcpu('CPU.RESO.4', 'DEBUT', ' ')
        call tldlg2(lraide, nprec, nstoc, vecrig)
        call uttcpu('CPU.RESO.1', 'FIN', ' ')
        call uttcpu('CPU.RESO.4', 'FIN', ' ')
    end if
!
! --  CONSTRUCTION DES BORNES DE RECHERCHE (FMIN/FMAX, OMEMIN/OMEMAX)
! --  EVENTUELLEMENT VIA UNE TABLE PROVENANT D'INFO_MODE
    if (typres(1:9) .eq. 'DYNAMIQUE') then
        fmin = rzero
        fmax = rzero
        if (nbborn .gt. 0) fmin = freq1
        if (nbborn .gt. 1) fmax = freq2
        if (lc .and. (fmin .lt. 0.d0)) then
            fmin = -fmin
            if (niv .ge. 1) call utmess('I', 'ALGELINE6_10')
        end if
        omemin = omega2(fmin)
        omemax = omega2(fmax)
    else
        omemin = rzero
        omemax = rzero
        if (nbborn .gt. 0) omemin = freq1
        if (nbborn .gt. 1) omemax = freq2
        fmin = omemin
        fmax = omemax
    end if
!
! --  SI AVEC L'OPTION 'BANDE', ON FOURNIT UNE TABLE, CONTROLE ET
! --  LECTURE DES PARAMETRES DE LA TABLE (AU FORMAT INFO_MODE)
! --  GAIN DE TEMPS, ON NE RECALCULE PAS TOUT DS LE VPFOPR SUIVANT
    ltabmo = .false.
    effmin = -999.d0
    effmax = -999.d0
    if (optiof(1:5) .eq. 'BANDE') then
        if (tabmod .ne. '') then
            ltabmo = .true.
            call vptabl(tabmod, typevp, fmin, fmax, precdc, &
                        nfreq, effmin, effmax)
! --  RECUPERATION DES BORNES EFFECTIVES ET CHANGEMENT DES BORNES DE LA BANDE SI NECESSAIRE
            rtest = abs(fmax-effmax)+abs(fmin-effmin)
            if (rtest .gt. eps) then
                valk(1) = tabmod
                call utmess('A', 'ALGELINE2_26', sk=valk(1))
                fmin = effmin
                fmax = effmax
                if (typres(1:9) .eq. 'DYNAMIQUE') then
                    omemin = omega2(fmin)
                    omemax = omega2(fmax)
                else
                    omemin = fmin
                    omemax = fmax
                end if
            end if
        end if
    end if

    if (typres(1:9) .ne. 'DYNAMIQUE') then
        rbid = omemin
        omemin = -omemax
        omemax = -rbid
    end if

!
! --- DETERMINATION D'INFO POUR LE TEST DE STURM ET LES POSITIONS
!     MODALES + CONSTRUCTION DE LA MATRICE DYNAMIQUE/ SA FACTORISEE
!     (STOCKEE DS MATOPA). CAS PARTICULIER DU QEP APPROCHE REEL, POUR
!     LEQUEL ON UTILISE AUSSI UNE MATRICE MATPSC DISTINCTE DE MATOPA.
!     DANS TOUS LES AUTRES CAS, MATPSC=MATOPA.
! --- AVEC QZ, UTILE QUE POUR GENERALISE REEL SYMETRIQUE
!     DANS LES AUTRES CAS (VPFOPC,WPFOPR,WPFOPC) ON SORT DES LE CALCUL
!     DU SHIFT SIGMA
    lmatra = 0
    lmtpsc = 0
    npivot = 0
    sigma = dcmplx(0.d0, 0.d0)
    omeshi = 0.d0
    if (.not. lc) then
! --  PROBLEME GENERALISE REEL SYMETRIQUE
        if (lkr .and. (.not. lns)) then
            call mtdefs(matopa, raide, 'V', 'R')
            call mtdscr(matopa)
            call jeveuo(matopa(1:19)//'.&INT', 'E', lmatra)
! --  POUR EVITER DE REFAIRE LE TEST DE STURM DE PRETTRAITEMENT AVEC L'OPTION 'BANDE'
            if (ltabmo) then
                optiov = 'BANDEA'
            else
                optiov = optiof
            end if
! --  POUR EVITER DE CALCULER LA MATRICE DYNAMIQUE ET SA FACTORISEE
            if ((optiof(1:5) .eq. 'BANDE') .and. (typcal(1:11) .eq. 'CALIBRATION')) then
                optiov = 'BANDEB'
            end if
            call vpfopr(optiov, typres, lmasse, lraide, lmatra, &
                        omemin, omemax, omeshi, nfreq, npiv2, &
                        omecor, precsh, nbrss, nblagr, solveu, &
                        rbid2, ibid2)
            npivot = npiv2(1)
            if (mod45(1:4) .eq. 'OP45') then
                if (nfreq .le. 0) then
                    if (arret(1:3) .eq. 'OUI') then
                        call utmess('Z', 'MODAL_1', num_except=ASTER_SOLVER_ERROR)
                    else
                        nfreq = 1
                        call rscrsd('G', modes, typcon, nfreq)
                        iretr = -2
                        goto 999
                    end if
                end if
            else
                if (nfreq .le. 0) then
                    call utmess('I', 'ALGELINE2_15')
                    iretr = -2
                    goto 999
                else
                    call codent(nfreq, 'G', chaine)
                    call utmess('I', 'ALGELINE2_16', sk=chaine)
                end if
                if (typcal(1:11) .eq. 'CALIBRATION') then
                    iretr = -3
                    nfreq_calibr_ = nfreq
                    goto 999
                end if
            end if
            lmtpsc = lmatra
            matpsc = matopa
        else
! --  PROBLEME GENERALISE COMPLEXE OU REEL NON SYM
            call vpfopc(lmasse, lraide, fmin, sigma, matopa, &
                        raide, lqz, solveu)
            if (.not. lqz) call jeveuo(matopa(1:19)//'.&INT', 'E', lmatra)
        end if
!
    else
!
! --  PROBLEME QUADRATIQUE REEL SYM OU NON SYM
        if (lkr) then
            call wpfopr(lmasse, lamor, lraide, appr, fmin, &
                        sigma, matopa, matpsc, raide, lqz, &
                        solveu)
            if (.not. lqz) then
                call jeveuo(matopa(1:19)//'.&INT', 'E', lmatra)
                call jeexin(matpsc(1:19)//'.&INT', iret)
                if (iret .ne. 0) call jeveuo(matpsc(1:19)//'.&INT', 'E', lmtpsc)
            end if
        else
! --  PROBLEME QUADRATIQUE COMPLEXE SYM
            call wpfopc(lmasse, lamor, lraide, fmin, sigma, &
                        matopa, raide, lqz, solveu)
            if (.not. lqz) call jeveuo(matopa(1:19)//'.&INT', 'E', lmatra)
        end if
    end if

    if ((typres(1:9) .ne. 'DYNAMIQUE') .and. (mod45(1:4) .ne. 'OP45')) then
        rbid = omemin
        omemin = -omemax
        omemax = -rbid
    end if
!
! --- ON BLINDE LES STRUCTURES DE DONNEES DE TYPE MATR_ASSE
! --- ON NE MANIPULE PAR LA SUITE QUE LEUR DESCRIPTEUR
! --- MATOPA --> LMATRA ET MATPSC --> LMTPSC
! --- ON TOUCHE A LEUR .REFA POUR DETRUIRE CORRECTEMENT LES
! --- EVENTUELLES OCCURENCES EXTERNES
    if (lmatra .eq. 0) matopa = ''
    if (lmtpsc .eq. 0) matpsc = ''
    if (lmatra .ne. 0) then
        call jeexin(matopa(1:19)//'.REFA', iret)
        if (iret .ne. 0) then
            call jeveuo(matopa(1:19)//'.REFA', 'E', krefa)
            zk24(krefa-1+7) = solveu
        end if
    end if
    if ((lmtpsc .ne. 0) .and. (lmtpsc .ne. lmatra)) then
        call jeexin(matpsc(1:19)//'.REFA', iret)
        if (iret .ne. 0) then
            call jeveuo(matpsc(1:19)//'.REFA', 'E', krefa)
            zk24(krefa-1+7) = solveu
        end if
    end if
!
!
! --  DETERMINATION DE LA TAILLE DE L'ESPACE DE PROJECTION (NBVECT)
    if (niv .ge. 1) call utmess('I', 'ALGELINE6_11', si=nfreq)
!
! --  CORRECTION DU NOMBRE DE FREQUENCES DEMANDEES
    if (nfreq .gt. neqact) then
        nfreq = neqact
        if (niv .ge. 1) call utmess('I', 'ALGELINE6_12', si=nfreq)
    end if
!
! --  DETERMINATION DE NBVECT (DIMENSION DU SOUS ESPACE)
    icoef = 0
    if (.not. lqz) then
        if (niv .ge. 1) call utmess('I', 'ALGELINE6_13', si=nbvect)
        if (nbvec2 .ne. 0) then
            icoef = nbvec2
        else
            select case (method)
            case ('JACOBI', 'SORENSEN')
                icoef = 2
            case ('TRI_DIAG')
                icoef = 4
            case default
                ASSERT(.false.)
            end select
        end if
        if (nbvect .lt. nfreq) then
            select case (method)
            case ('JACOBI')
                nbvect = min(min(7+nfreq, icoef*nfreq), neqact)
            case ('TRI_DIAG')
                nbvect = min(max(7+nfreq, icoef*nfreq), neqact)
            case ('SORENSEN')
                nbvect = min(max(3+nfreq, icoef*nfreq), neqact)
            case default
                ASSERT(.false.)
            end select
            if (niv .ge. 1) call utmess('I', 'ALGELINE6_14', si=nbvect)
        else
            if (nbvect .gt. neqact) then
                nbvect = neqact
                if (niv .ge. 1) call utmess('I', 'ALGELINE6_15', si=nbvect)
            end if
        end if
    end if
!
! --  TRAITEMENT SPECIFIQUE A SORENSEN
    if ((method .eq. 'SORENSEN') .and. (nbvect-nfreq .le. 2)) then
        if (nfreq .gt. (neqact+2)) nfreq = neqact-2
        nbvect = nfreq+2
    end if
!
! --  TRAITEMENT SPECIFIQUE A QZ
!     AVEC QZ ON A PAS D'ESPACE DE PROJECTION, IL FAUT DONC AFFECTER NBVECT EN DUR
    if (lqz) nbvect = neq
!
! -- CORRECTION DE NBVECT DANS LE CAS QUADRATIQUE
    if (lc) then
        nbvect = 2*nbvect
        nfreq = 2*nfreq
        call utmess('I', 'ALGELINE2_75')
    end if
!
!
! --  ON MODIFIE LES VALEURS NFREQ ET DE NBVECT DE LA SD EIGENSOLVER
    if (nfreq > 200) call utmess('A', 'MODAL_22', si=nfreq)
    call vpecri(eigsol, 'I', 1, k24bid, rbid, nfreq)
    call vpecri(eigsol, 'I', 2, k24bid, rbid, nbvect)
!
999 continue
!
    call jedema()
!
!     FIN DE VPINI1
!
end subroutine
