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

subroutine vpfopr(option, typres, lmasse, lraide, ldynam, &
                  omemin, omemax, omeshi, nbfreq, npivot, &
                  omecor, precsh, nbrssa, nblagr, solveu, &
                  det, idet)
!     DETERMINATION DE SHIFT(S), D'UNE MATRICE SHIFTEE, DE SA FACTORISEE
!     DU NBRE DE PIVOTS NEGATIFS (POUR TEST DE STURM) VOIRE DU NBRE
!     DE FREQ DANS UNE BANDE.
!
!     POUR ETAPE DE PRETRAITEMENTS DE MODE_ITER_SIMULT
!     OPTION='CENTRE' --> OUTPUT: MATRICE SHIFTEE + SA FACTORISEE +
!                         NPIVOT(1) + OMESHI
!                         INPUT: INTENDANCE + OMEMIN
!
!     OPTION='BANDE'  --> OUTPUT: MATRICE SHIFTEE + SA FACTORISEE +
!                         NPIVOT(1) + NBFREQ + OMEMIN/MAX + OMESHI +
!                         AFFICHAGES VPECST
!                         INPUT: INTENDANCE + OMEMIN/MAX
!
!     OPTION='BANDEA'  --> OUTPUT: MATRICE SHIFTEE + SA FACTORISEE +
!                         NPIVOT(1) + OMEMIN/MAX + OMESHI +
!                         AFFICHAGES VPECST ENRICHIS + NBFREQ (JUSTE
!                         POUR AFFICHAGE)
!                         INPUT: INTENDANCE + OMEMIN/MAX + NBFREQ
!
!     OPTION='BANDEB'  --> OUTPUT: NBFREQ + OMEMIN/MAX +
!                         AFFICHAGES VPECST
!                         INPUT: INTENDANCE + OMEMIN/MAX
!
!     OPTION='PLUS_PETITE'/'TOUT' --> OUTPUT: MATRICE SHIFTEE + SA FACTO
!                         RISEE + NPIVOT(1) + OMESHI
!                         INPUT: INTENDANCE
!
!     POUR ETAPE DE POST_TRAITEMENTS DE MODE_ITER_SIMULT
!     OPTION='STURM'  --> OUTPUT: NBFREQ + OMEMIN/MAX + PAS VPECST
!                         AFFICHAGES DEDIES EN AMONT
!                         INPUT: INTENDANCE + OMEMIN/MAX
!
!     POUR INFO_MODE + STURM AVEC LISTE DE FREQ OU CHAR_CRIT
!     PREMIERE BANDE
!     OPTION='STURML1' --> OUTPUT: NBFREQ + OMEMIN/MAX +
!                         AFFICHAGES DEDIE VPECST + NPIVOT(2)
!                         INPUT: INTENDANCE + OMEMIN/MAX
!         ='STURML1P' ...     IDEM + COMM POUR MACRO // (ETAPE INITIALE)
!         ='STURML10/11'...IDEM + COMM POUR MACRO // (ETAPE FINALE)
!
!     BANDES SUIVANTES
!     OPTION='STURMLN' --> OUTPUT: NBFREQ + OMEMIN/MAX +
!                         AFFICHAGES DEDIE VPECST + NPIVOT(2)
!                         INPUT: INTENDANCE + OMEMIN/MAX +
!                         NPIVOT(1)=NPIVOT OMEMIN ET NPIVOT(2)=
!                            NUMERO DE LA BANDE CONSIDEREE
!                         OUTPUT:NPIVOT(2)=NPIVOT OMEMAX.
!           ='STURMLNP' ... IDEM + COMM POUR MACRO //
!     OPTION='STURMLNS' --> OUTPUT: NBFREQ + OMEMIN/MAX +
!                         AFFICHAGES DEDIE VPECST + NPIVOT(2)
!                         INPUT: INTENDANCE + OMEMIN/MAX +
!                         NPIVOT(2)=NPIVOT OMEMAX ET NPIVOT(1)=
!                            NUMERO DE LA BANDE CONSIDEREE
!                         OUTPUT:NPIVOT(1)=NPIVOT OMEMIN.
!           ='STURMLNQ' ... IDEM + COMM POUR MACRO //
!
!     POUR ETAPE PRETRAITEMENT DE MODE_ITER_INV (AJUSTE/SEPARE)
!     OPTION='STURMAD' --> OUTPUT: NBFREQ + OMEMIN/MAX +
!                     AFFICHAGES VPECST + DET(2)/IDET(2) + NPIVOT(2)
!                          INPUT: INTENDANCE + OMEMIN/MAX
!     ------------------------------------------------------------------
! IN  OPTION  : TX : CHOIX DE L'OPTION (PLUS_PETITE,CENTRE,BANDE,STURM)
! IN  TYPRES  : TX : TYPE DU CALCUL (DYNAMIQUE OU FLAMBEMENT)
! IN  LMASSE  : IS : DESCRIPTEUR DE LA MATRICE SECOND MEMBRE
! IN  LRAIDE  : IS : DESCRIPTEUR DE LA MATRICE PREMIER MEMBRE
! IN/OUT LDYNAM :IS : POINTEUR SUR LA FACTORISEE DE LA MATRICE DYNAMIQUE
!                    INDUITE PAR L'OPTION
! IN/OUT OMEMIN : R8 : VALEUR INFERIEURE DE LA BANDE DE RECHERCHE
!                      OU VALEUR DE DEPART POUR LES AUTRES OPTIONS
! IN/OUT OMEMAX : R8 : VALEUR SUPERIEURE DE LA BANDE DE RECHERCHE
!    OUT OMESHI : R8 : VALEUR DU SHIFT  DE LA MATRICE DE TRAVAIL
!    OUT NBFREQ : IS : NOMBRE DE FREQUENCES DANS LA BANDE
! IN/OUT NPIVOT : IS : VECTEUR NOMBRE DE PIVOTS NEGATIFS DE LA MATRICE
!                      DE TRAVAIL FACTORISEE.
!                      ATTENTION PARAMETRE PARFOIS UTILISE
!                      EN INPUT AVEC UN SENS DIFFERENT CF. OPTION.
! IN  OMECOR : R8 : VALEUR DE LA PULSATION AU CARRE DEFINISSANT LES
!                   MODES DE CORPS RIGIDE
! IN  PRECSH : R8 : VALEUR DU DECALAGE DES SHIFTS QUAND LA MATRICE EST
!                   NON INVERSIBLE (CALC_FREQ/PREC_SHIFT)
! IN  NBRSSA : IS : NOMBRE DE DECALAGES DE SHIFTS AUTORISES
! IN  NBLAGR : IS : NOMBRE DE DDLS DE LAGRANGE
! IN  SOLVEU : K19 : SD SOLVEUR POUR PARAMETRER LE SOLVEUR LINEAIRE
! OUT  DET   : R8  : VECTEUR DES DEUX MANTISSES DE DETERMINANT
! OUT  IDET  : IS  : IDEM SUR LES EXPOSANTS
!----------------------------------------------------------------------
!
    implicit none
!
! PARAMETRES D'APPEL
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/asmpi_comm.h"
#include "asterfort/asmpi_barrier.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/freqom.h"
#include "asterfort/infniv.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/vecint.h"
#include "asterfort/vpecst.h"
#include "asterfort/vpstur.h"
#include "asterfort/wkvect.h"
    character(len=*) :: option
    character(len=16) :: typres
    character(len=19) :: solveu
    integer(kind=8) :: lmasse, lraide, ldynam, nbrssa
    real(kind=8) :: omemin, omemax, omeshi, omecor, precsh, det(2)
    integer(kind=8) :: nbfreq, npivot(2), nblagr, idet(2)
!
!
! VARIABLES LOCALES
    mpi_int :: mpicou, mpicow, rang, rangl, nbproc
    character(len=1) :: typep
    character(len=8) :: k8bid
    character(len=16) :: ch16, valk(3)
    character(len=24) :: k24c, k24par
    integer(kind=8) :: niv, ifm, nbessa, ier, nbfmin, nbfmax, ibid, ibande
    integer(kind=8) :: jk24c, jkpar, nbrow, frecou
    real(kind=8) :: valr(2), omgmin, omgmax, omgshi, rbid, prec, omgdec
    aster_logical :: caldet, ldyna
!
    call infniv(ifm, niv)
! MAUVAISE VALEUR DE OPTION
    if ((option .ne. 'CENTRE') .and. (option .ne. 'BANDE') .and. (option .ne. 'BANDEA') .and. &
        (option .ne. 'PLUS_PETITE') .and. (option .ne. 'TOUT') .and. (option .ne. 'STURM') .and. &
        (option .ne. 'STURML1') .and. (option .ne. 'STURML1P') .and. (option .ne. 'STURML10') &
        .and. (option .ne. 'STURML11') .and. (option .ne. 'STURMLN') .and. &
        (option .ne. 'STURMLNP') .and. (option .ne. 'STURMAD') .and. (option .ne. 'STURMLNS') &
        .and. (option .ne. 'STURMLNQ') .and. (option .ne. 'BANDEB')) then
        ASSERT(.false.)
    end if

    det(1) = -9999.d0
    det(2) = -9999.d0
    idet(1) = -9999
    idet(2) = -9999
    ibande = 1
    if (option(1:7) .ne. 'STURMLN') then
        npivot(1) = -9999
        npivot(2) = -9999
    elseif ((option .eq. 'STURMLNS') .or. (option .eq. 'STURMLNQ')) then
        ibande = npivot(1)
        npivot(1) = -9999
    else
        ibande = npivot(2)
        npivot(2) = -9999
    end if
    if (option(1:7) .eq. 'STURML1') ibande = 1
    if (option .eq. 'STURMAD') then
        caldet = .true.
    else
        caldet = .false.
    end if
    valk(1) = 'FREQ'
    valk(2) = 'SEUIL_FREQ'
    valk(3) = 'MATR_RIGI'
    if (typres .eq. 'DYNAMIQUE') then
        ldyna = .true.
    else
        ldyna = .false.
        valk(1) = 'CHAR_CRIT'
        valk(2) = 'SEUIL_CHAR_CRIT'
        if (typres .eq. 'GENERAL') valk(3) = 'MATR_A'
    end if
    if (option(1:5) .ne. 'STURM') then
        if (option(1:5) .eq. 'BANDE') then
            call utmess('I', 'ALGELINE6_41', sk='BANDE')
        else
            call utmess('I', 'ALGELINE6_41', sk=option)
        end if
    end if
!     ------------------------------------------------------------------
!     ------------------------ OPTION CENTRE ---------------------------
!     ------------------------------------------------------------------
!
    if (option .eq. 'CENTRE') then
!
        omgshi = omemin
        nbessa = 0
        prec = precsh
10      continue
        ier = 0
        call vpstur(lraide, omgshi, lmasse, ldynam, rbid, &
                    ibid, npivot(1), ier, solveu, .false._1, &
                    .true._1)
        if (ier .ne. 0) then
            nbessa = nbessa+1
            if (nbessa .le. nbrssa) then
                if (abs(omgshi) .lt. omecor) then
                    omgshi = omecor
                    if (ldyna) then
                        valr(1) = freqom(omgshi)
                    else
                        valr(1) = omgshi
                    end if
                    if (niv .ge. 1) then
                        call utmess('I', 'ALGELINE6_44', sr=valr(1))
                    end if
! --- CE N'EST PLUS LA PEINE DE DECALER, C'EST INUTILE
                    nbessa = nbrssa
                else
                    omgdec = max(omecor, prec*abs(omgshi))
                    omgshi = omgshi+omgdec
                    if (ldyna) then
                        valr(1) = freqom(omgshi)
                    else
                        valr(1) = omgshi
                    end if
                    if (niv .ge. 1) then
                        valr(2) = prec*100.d0
                        call utmess('I', 'ALGELINE6_45', nr=2, valr=valr)
                    end if
                    prec = 2.d0*prec
                end if
                goto 10
            else
                if (ldyna) then
                    valr(1) = freqom(omgshi)
                else
                    valr(1) = omgshi
                end if
                call utmess('F', 'ALGELINE3_81', nk=3, valk=valk, sr=valr(1))
            end if
!
        end if
        omeshi = omgshi
        if (niv .ge. 1) then
            if (ldyna) then
                call utmess('I', 'ALGELINE6_42', sr=freqom(omeshi))
            else
                call utmess('I', 'ALGELINE6_43', sr=-omeshi)
            end if
        end if
!
!     ------------------------------------------------------------------
!     ------------------------ OPTION BANDE* OU STURM** ----------------
!     ------------------------------------------------------------------
!
    else if ((option(1:5) .eq. 'BANDE') .or. (option(1:5) .eq. 'STURM')) &
        then
!
        omgmin = omemin
        if ((option .eq. 'STURMLN') .or. (option .eq. 'STURMLNP')) then
            nbfmin = npivot(1)
        else if ((option .ne. 'BANDEA') .and. (option .ne. 'STURML11')) then
            nbessa = 0
            prec = precsh
21          continue
            ier = 0
            call vpstur(lraide, omgmin, lmasse, ldynam, det(1), &
                        idet(1), npivot(1), ier, solveu, caldet, &
                        .false._1)
            nbfmin = npivot(1)
            if (ier .ne. 0) then
                nbessa = nbessa+1
                if (nbessa .le. nbrssa) then
                    if (abs(omgmin) .lt. omecor) then
                        omgmin = -omecor
                        if (ldyna) then
                            valr(1) = freqom(omgmin)
                        else
                            valr(1) = omgmin
                        end if
                        if (niv .ge. 1) then
                            call utmess('I', 'ALGELINE6_46', sr=valr(1))
                        end if
! --- CE N'EST PLUS LA PEINE DE DECALER, C'EST INUTILE
                        nbessa = nbrssa
                    else
                        omgdec = max(omecor, prec*abs(omgmin))
                        omgmin = omgmin-omgdec
                        if (ldyna) then
                            valr(1) = freqom(omgmin)
                        else
                            valr(1) = omgmin
                        end if
                        if (niv .ge. 1) then
                            valr(2) = prec*100.d0
                            call utmess('I', 'ALGELINE6_47', nr=2, valr=valr)
                        end if
                        prec = 2.d0*prec
                    end if
                    goto 21
                else
                    call utmess('A+', 'ALGELINE3_82')
                    call utmess('A', 'ALGELINE3_84', sk=valk(2))
                end if
            end if
        end if
        omemin = omgmin
        if (omemin .ge. omemax) then
            if (ldyna) then
                call utmess('F', 'ALGELINE3_85')
            else
                call utmess('F', 'ALGELINE3_86')
            end if
        end if
!
        omgmax = omemax
        if ((option .eq. 'STURMLNS') .or. (option .eq. 'STURMLNQ')) then
            nbfmax = npivot(2)
        elseif ((option .ne. 'BANDEA') .and. (option .ne. 'STURML10')) then
            nbessa = 0
            prec = precsh
22          continue
            ier = 0
            call vpstur(lraide, omgmax, lmasse, ldynam, det(2), &
                        idet(2), npivot(2), ier, solveu, caldet, &
                        .false._1)
            nbfmax = npivot(2)
            if (ier .ne. 0) then
                nbessa = nbessa+1
                if (nbessa .le. nbrssa) then
                    if (abs(omgmax) .lt. omecor) then
                        omgmax = omecor
                        if (ldyna) then
                            valr(1) = freqom(omgmax)
                        else
                            valr(1) = omgmax
                        end if
                        if (niv .ge. 1) then
                            call utmess('I', 'ALGELINE6_48', sr=valr(1))
                        end if
! --- CE N'EST PLUS LA PEINE DE DECALER, C'EST INUTILE
                        nbessa = nbrssa
                    else
                        omgdec = max(omecor, prec*abs(omgmax))
                        omgmax = omgmax+omgdec
                        if (ldyna) then
                            valr(1) = freqom(omgmax)
                        else
                            valr(1) = omgmax
                        end if
                        if (niv .ge. 1) then
                            valr(2) = prec*100.d0
                            call utmess('I', 'ALGELINE6_49', nr=2, valr=valr)
                        end if
                        prec = 2.d0*prec
                    end if
                    goto 22
                else
                    call utmess('A+', 'ALGELINE3_83')
                    call utmess('A', 'ALGELINE3_84', sk=valk(2))
                end if
            end if
        end if
        omemax = omgmax
!
!     ------------------------------------------------------------------
!     -------- INFO_MODE OU CALC_MODES SUR PLUSIEURS SOUS-BANDES
!     -------- EN PARALLELE (PART I)
!     ------------------------------------------------------------------
!     --- COMMUNICATION DES PIVOTS POUR LE BON CALCUL DE STURM
        if ((option .eq. 'STURML1P') .or. (option .eq. 'STURMLNP') .or. (option .eq. 'STURML10') &
            .or. (option .eq. 'STURML11') .or. (option .eq. 'STURMLNQ')) then
            call asmpi_comm('GET_WORLD', mpicow)
            call asmpi_comm('GET', mpicou)
            ASSERT(mpicou .ne. mpicow)
!         --- ON REMPLACE LE COMM LOCAL PAR LE COMM WORLD
            call asmpi_info(mpicou, rangl)
            call asmpi_comm('SET', mpicow)
            call asmpi_barrier()
            call asmpi_info(mpicow, rang, nbproc)
!         --- BUFFER DE COM K24C
!         --- K24C(FREQUENCE_COURANTE)=NBFMIN OU MAX
            k24par = '&&OP0032.COULEUR'
            call jeveuo(k24par, 'L', jkpar)
            k24c = '&&VPFOPR.BUFFERMPI'
            nbrow = zi(jkpar+nbproc-1)
            frecou = zi(jkpar+rang)
            call wkvect(k24c, 'V V I', nbrow+1, jk24c)
            call vecint(nbrow+1, 0, zi(jk24c))
            if (option .eq. 'STURML1P') then
                ASSERT(frecou .eq. 1)
                if (rangl .eq. 0) zi(jk24c+1) = nbfmax
!
            else if ((option .eq. 'STURMLNP') .or. (option .eq. 'STURMLNQ')) then
                ASSERT(frecou .gt. 1)
                if (rangl .eq. 0) zi(jk24c+frecou) = nbfmax
!
            else if (option .eq. 'STURML10') then
                ASSERT(frecou .eq. 0)
                if (rangl .eq. 0) zi(jk24c) = nbfmin
!
            else if (option .eq. 'STURML11') then
                ASSERT(frecou .eq. 1)
                if (rangl .eq. 0) zi(jk24c+1) = nbfmax
            end if
!
            call asmpi_comm_vect('MPI_SUM', 'I', nbval=nbrow+1, vi=zi(jk24c))
!
            if ((option .eq. 'STURMLNP') .or. (option .eq. 'STURMLNQ')) then
                nbfmin = zi(jk24c+frecou-1)
            else if ((option .eq. 'STURML10') .or. (option .eq. 'STURML11')) &
                then
                nbfmin = zi(jk24c)
                nbfmax = zi(jk24c+1)
            end if
            call jedetr(k24c)
        end if
!
        k8bid = ' '
        if ((option .eq. 'BANDE') .or. (option .eq. 'BANDEB') .or. (option .eq. 'STURMAD')) then
            typep = 'R'
        else if (option .eq. 'STURM') then
            typep = 'S'
        else if ((option(1:7) .eq. 'STURML1') .or. (option(1:7) &
                                                    .eq. 'STURMLN')) then
            typep = 'D'
            nbfreq = ibande
        else if (option .eq. 'BANDEA') then
            typep = 'F'
        end if
!
        call vpecst(ifm, typres, omgmin, omgmax, nbfmin, &
                    nbfmax, nbfreq, nblagr, typep, k8bid, &
                    0.d0, dcmplx(0.d0, 0.d0))
!
!     ------------------------------------------------------------------
!     -------- INFO_MODE OU CALC_MODES SUR PLUSIEURS SOUS-BANDES
!     -------- EN PARALLELE (PART II)
!     ------------------------------------------------------------------
!     --- SEULS CERTAINS PROCS REMONTENT LES OUTPUTS SINON LA COMM
!     --- EN FIN DE OP0032 VA CUMULER DES INFOS REDONDANTES.
        if ((option .eq. 'STURML1P') .or. (option .eq. 'STURMLNP') .or. (option .eq. 'STURML10') &
            .or. (option .eq. 'STURML11') .or. (option .eq. 'STURMLNQ')) then
            if (rangl .ne. 0) then
                omemin = 0.d0
                omemax = 0.d0
                nbfreq = 0
            end if
            if (option .eq. 'STURML11') then
                omemin = 0.d0
                nbfreq = 0
            else if (option .eq. 'STURML10') then
                omemax = 0.d0
            end if
        end if
!
        omgshi = (omgmax+omgmin)*0.5d0
        if ((option .eq. 'BANDE') .or. (option .eq. 'BANDEA')) then
!          --- CENTRAGE DE L INTERVALLE ---
            nbessa = 0
            prec = precsh
23          continue
            ier = 0
            call vpstur(lraide, omgshi, lmasse, ldynam, rbid, &
                        ibid, npivot(1), ier, solveu, .false._1, &
                        .true._1)
            if (ier .ne. 0) then
                nbessa = nbessa+1
                if (nbessa .le. nbrssa) then
                    if (abs(omgshi) .lt. omecor) then
                        omgshi = omecor
                        if (ldyna) then
                            valr(1) = freqom(omgshi)
                        else
                            valr(1) = omgshi
                        end if
                        if (niv .ge. 1) then
                            call utmess('I', 'ALGELINE6_44', sr=valr(1))
                        end if
! --- CE N'EST PLUS LA PEINE DE DECALER, C'EST INUTILE
                        nbessa = nbrssa
                    else
                        omgdec = max(omecor, prec*abs(omgshi))
                        omgshi = omgshi-omgdec
                        if (ldyna) then
                            valr(1) = freqom(omgshi)
                        else
                            valr(1) = omgshi
                        end if
                        if (niv .ge. 1) then
                            valr(2) = prec*100.d0
                            call utmess('I', 'ALGELINE6_92', nr=2, valr=valr)
                        end if
                        prec = 2.d0*prec
                    end if
                    goto 23
                else
                    if (ldyna) then
                        valr(1) = freqom(omgshi)
                    else
                        valr(1) = omgshi
                    end if
                    call utmess('F', 'ALGELINE3_81', nk=3, valk=valk, sr=valr(1))
                end if
            end if
        end if
        omeshi = omgshi
!
!          --- AFFICHAGE COMMUN ---
        if ((niv .ge. 1) .and. ((option .eq. 'BANDE') .or. (option .eq. 'BANDEA'))) then
            if (ldyna) then
                valr(1) = freqom(omgmin)
                valr(2) = freqom(omgmax)
                call utmess('I', 'ALGELINE6_93', nr=2, valr=valr)
                if (option(1:5) .eq. 'BANDE') then
                    call utmess('I', 'ALGELINE6_42', sr=freqom(omeshi))
                end if
            else
                valr(1) = -omgmax
                valr(2) = -omgmin
                call utmess('I', 'ALGELINE6_94', nr=2, valr=valr)
                if (option(1:5) .eq. 'BANDE') then
                    call utmess('I', 'ALGELINE6_43', sr=-omeshi)
                end if
            end if
        end if
!
!
!     ------------------------------------------------------------------
!     ------------------------ OPTION PLUS_PETITE OU TOUT -------------
!     ------------------------------------------------------------------
!
    else if ((option .eq. 'PLUS_PETITE') .or. (option .eq. 'TOUT')) then
!
        omgshi = 0.d0
        nbessa = 0
        prec = precsh
30      continue
        ier = 0
        call vpstur(lraide, omgshi, lmasse, ldynam, rbid, &
                    ibid, npivot(1), ier, solveu, .false._1, &
                    .true._1)
        if (ier .ne. 0) then
            nbessa = nbessa+1
            if (nbessa .le. nbrssa) then
                if (abs(omgshi) .lt. omecor) then
                    omgshi = -omecor
                    if (ldyna) then
                        valr(1) = freqom(omgshi)
                    else
                        valr(1) = omgshi
                    end if
                    if (niv .ge. 1) then
                        call utmess('I', 'ALGELINE6_48', sr=valr(1))
                    end if
! --- CE N'EST PLUS LA PEINE DE DECALER, C'EST INUTILE
                    nbessa = nbrssa
                else
                    omgdec = max(omecor, prec*abs(omgshi))
                    omgshi = omgshi-omgdec
                    if (ldyna) then
                        valr(1) = freqom(omgshi)
                    else
                        valr(1) = omgshi
                    end if
                    if (niv .ge. 1) then
                        valr(2) = prec*100.d0
                        call utmess('I', 'ALGELINE6_92', nr=2, valr=valr)
                    end if
                    prec = 2.d0*prec
                end if
                goto 30
            else
                if (ldyna) then
                    valr(1) = freqom(omgshi)
                else
                    valr(1) = omgshi
                end if
                call utmess('F', 'ALGELINE3_81', nk=3, valk=valk, sr=valr(1))
            end if
        end if
        omeshi = omgshi
        if (niv .ge. 1) then
            if (ldyna) then
                call utmess('I', 'ALGELINE6_42', sr=freqom(omeshi))
            else
                call utmess('I', 'ALGELINE6_43', sr=-omeshi)
            end if
        end if
!
!     ------------------------------------------------------------------
!     ------------------------ OPTION NON CONNUE -----------------------
!     ------------------------------------------------------------------
!
    else
        ch16 = option
        call utmess('F', 'ALGELINE3_69', sk=ch16)
    end if
!
    if ((niv .ge. 1) .and. (option(1:5) .ne. 'STURM')) write (ifm, 1200)
!
!     -----------------------------FORMAT------------------------------
!
1200 format(72('-'),/)
!     ------------------------------------------------------------------
!
end subroutine
