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
subroutine tensca(tablca, icabl, nbnoca, nbf0, f0, &
                  delta, typrel, trelax, xflu, xret, &
                  ea, rh1000, mu0, fprg, frco, &
                  frli, sa, regl, analy)
    implicit none
!  DESCRIPTION : CALCUL DE LA TENSION LE LONG D'UN CABLE
!  -----------   APPELANT : OP0180 , OPERATEUR DEFI_CABLE_BP
!
!                EN SORTIE ON COMPLETE LA TABLE RESULTAT
!                LES LIGNES COMPLETEES CORRESPONDENT AU DERNIER CABLE
!                LES CASES RENSEIGNEES CORRESPONDENT AU PARAMETRE
!                <TENSION>
!
!  IN     : TABLCA : CHARACTER*19
!                    NOM DE LA TABLE DECRIVANT LES CABLES
!  IN     : ICABL  : INTEGER , SCALAIRE
!                    NUMERO DU CABLE
!  IN     : NBNOCA : INTEGER ,
!                    CONTIENT LE NOMBRE DE NOEUDS DU CABLE ETUDIE
!  IN     : NBF0   : INTEGER , SCALAIRE
!                    NOMBRE D'ANCRAGES ACTIFS DU CABLE (0, 1 OU 2)
!  IN     : F0     : REAL*8 , SCALAIRE
!                    VALEUR DE LA TENSION APPLIQUEE A L'UN OU AUX DEUX
!                    ANCRAGES ACTIFS DU CABLE
!  IN     : DELTA  : REAL*8 , SCALAIRE
!                    VALEUR DU RECUL DE L'ANCRAGE
!  IN     : TYPREL  : CHARACTER*24
!                    TYPE DE RELAXATION UTILISEE
!  IN     : TRELAX : REAL*8 , SCALAIRE
!                    VALEUR DE LA FONCTION CARACTERISANT L'EVOLUTION DE
!                    LA RELAXATION DE L'ACIER DANS LE TEMPS POUR BPEL
!                    OU NOMBRE D'HEURES POUR LA RELAXATION SI ETCC
!                    UTILE SI RELAX = .TRUE.
!  IN     : XFLU   : REAL*8 , SCALAIRE
!                    VALEUR DU TAUX DE PERTE DE TENSION PAR FLUAGE DU
!                    BETON, EN % DE LA TENSION INITIALE
!  IN     : XRET   : REAL*8 , SCALAIRE
!                    VALEUR DU TAUX DE PERTE DE TENSION PAR RETRAIT DU
!                    BETON, EN % DE LA TENSION INITIALE
!  IN     : EA     : REAL*8 , SCALAIRE
!                    VALEUR DU MODULE D'YOUNG DE L'ACIER
!  IN     : RH1000 : REAL*8 , SCALAIRE
!                    VALEUR DE LA RELAXATION A 1000 HEURES EN %
!  IN     : MU0    : REAL*8 , SCALAIRE
!                    VALEUR DU COEFFICIENT DE RELAXATION DE L'ACIER
!                    PRECONTRAINT POUR BPEL
!
!  IN     : FPRG     : REAL*8 , SCALAIRE
!                    VALEUR DE LA CONTRAINTE LIMITE ELASTIQUE DE L'ACIER
!  IN     : FRCO   : REAL*8 , SCALAIRE
!                    VALEUR DU COEFFICIENT DE FROTTEMENT EN COURBE
!                    (CONTACT ENTRE LE CABLE ACIER ET LE MASSIF BETON)
!  IN     : FRLI   : REAL*8 , SCALAIRE
!                    VALEUR DU COEFFICIENT DE FROTTEMENT EN LIGNE
!                    (CONTACT ENTRE LE CABLE ACIER ET LE MASSIF BETON)
!  IN     : SA     : REAL*8 , SCALAIRE
!                    VALEUR DE L'AIRE DE LA SECTION DROITE DU CABLE
!  IN     : REGL   : CHARACTER*4, INDICATION DU REGLEMENT UTILISE
!                    BPEL OU ETCC
!  IN     : ANALY  : CHARACTER*4, TYPE DE CALCUL REALISE
!                    DEFI OU ETCC OU RUPT
!
!
!
!
!-------------------   DECLARATION DES VARIABLES   ---------------------
!
!
! ARGUMENTS
! ---------
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/getvid.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbexip.h"
#include "asterfort/tbexve.h"
#include "asterfort/tensk1.h"
#include "asterfort/tensk2.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    character(len=19) :: tablca
    character(len=4) :: regl, analy
    character(len=24) :: typrel
    integer(kind=8) :: icabl, nbnoca, nbf0
    real(kind=8) :: f0, delta, trelax, xflu, xret, ea, rh1000, mu0, fprg, frco
    real(kind=8) :: frli, sa
!
!
! VARIABLES LOCALES
! -----------------
    real(kind=8), parameter:: prec = 1.d-6

    integer(kind=8) :: ibid, idecno, ino, ipara, jabsc, jalph, jf, nblign
    integer(kind=8) :: nbpara, n1, irt, jtabx, jtaby, nbval
    real(kind=8) :: df, flim, krelax, fi, f2, lg_ref
    complex(kind=8) :: cbid
    aster_logical :: trouv1, trouv2, exi1, exi2
    character(len=3) :: k3b
    character(len=24) :: abscca, alphca
    character(len=8) :: ntable, k8b
    character(len=19) :: newtab
    character(len=24) :: tabx, taby
!
    character(len=24) :: param, parcr(2)
    integer(kind=8), pointer :: tbnp(:) => null()
    character(len=24), pointer :: tblp(:) => null()
    data param/'TENSION                 '/
    data parcr/'ABSC_CURV               ',&
     &                     'ALPHA                   '/
!
!-------------------   DEBUT DU CODE EXECUTABLE    ---------------------
!
    call jemarq()
    cbid = (0.d0, 0.d0)
    ibid = 0
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 1   TRAITEMENT DES CAS PARTICULIERS F0 = 0 OU PAS D'ANCRAGE ACTIF
!     ET VERIFICATION DE LA COHERENCE DES DONNEES
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     NBNO = NBNOCA(ICABL)
!
    call jeveuo(tablca//'.TBNP', 'L', vi=tbnp)
    nblign = tbnp(2)
    idecno = nblign-nbnoca
!
    if ((f0 .eq. 0.0d0) .or. (nbf0 .eq. 0)) then
        do ino = 1, nbnoca
            call tbajli(tablca, 1, param, [ibid], [0.d0], &
                        [cbid], k3b, idecno+ino)
        end do
        goto 999
    end if
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 2   RECUPERATION DE L'ABSCISSE CURVILIGNE ET DE LA DEVIATION ANGULAIRE
!     CUMULEE LE LONG DU CABLE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    nbpara = tbnp(1)
    call jeveuo(tablca//'.TBLP', 'L', vk24=tblp)
    trouv1 = .false.
    trouv2 = .false.
    do ipara = 1, nbpara
        if (tblp(1+4*(ipara-1)) .eq. parcr(1)) then
            trouv1 = .true.
            abscca = tblp(1+4*(ipara-1)+2)
            call jeveuo(abscca, 'L', jabsc)
        end if
        if (tblp(1+4*(ipara-1)) .eq. parcr(2)) then
            trouv2 = .true.
            alphca = tblp(1+4*(ipara-1)+2)
            call jeveuo(alphca, 'L', jalph)
        end if
        if (trouv1 .and. trouv2) goto 30
    end do
!
30  continue
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 3   VERIFICATION DE LA DONNEE DE LA TENSION POUR
!                  MODI_CABLE_ETCC et MODI_CABLE_RUPT
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    call wkvect('&&TENSCA.F', 'V V R', nbnoca, jf)
!
    if ((analy .eq. 'RUPT') .or. (analy .eq. 'ETCC')) then
        call getvid('DEFI_CABLE', 'TENSION', iocc=icabl, scal=ntable, nbret=n1)
        if (n1 .eq. 0) then
            call utmess('F', 'CABLE0_8')
        end if
!
        newtab = ntable
        tabx = '&&TENSCA_TABREF_CURV'
        taby = '&&TENSCA_TABREF_TENS'
!
        call jeexin(newtab//'.TBBA', irt)
        if (irt .eq. 0) then
            call utmess('F', 'UTILITAI4_64')
        end if
!     VERIFICATION DE LA PRESENCE DES BONS PARAMETRES
        call tbexip(newtab, 'ABSC_CURV', exi1, k8b)
        call tbexip(newtab, 'N', exi2, k8b)
!
        if (.not. exi1 .and. .not. exi2) then
            call utmess('F', 'CABLE0_9')
        end if
!
        call tbexve(newtab, 'ABSC_CURV', tabx, 'V', nbval, &
                    k8b)
        call jeveuo(tabx, 'L', jtabx)
        call tbexve(newtab, 'N', taby, 'V', nbval, &
                    k8b)
        call jeveuo(taby, 'L', jtaby)
        if (nbval .ne. nbnoca) then
            call utmess('F', 'CABLE0_10')
        end if
!     ON VERIFIE A MINIMA QUE LES ABSCISSES CURVILIGNES SONT IDENTIQUES
!     (MAIS PAS LES COORDONNEES EXACTES)
        do ino = 1, nbnoca
            ! longueur de reference des mailles attachees aux noeuds courant
            if (ino .eq. 1) then
                lg_ref = abs(zr(jabsc+ino+1-1)-zr(jabsc+ino-1))
            else if (ino .eq. nbnoca) then
                lg_ref = abs(zr(jabsc+ino-1)-zr(jabsc+ino-1-1))
            else
                lg_ref = 0.5d0*(abs(zr(jabsc+ino+1-1)-zr(jabsc+ino-1-1)))
            end if
            if (zr(jtabx+ino-1)-zr(jabsc+ino-1) .ge. prec*lg_ref) then
                call utmess('F', 'CABLE0_11')
            end if
        end do
    end if
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 4   CALCUL DE LA TENSION LE LONG DU CABLE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    if (analy .ne. 'RUPT') then
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    4.1  CAS GENERAL
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! 4.1.1 CALCUL DE LA TENSION LE LONG DU CABLE EN PRENANT EN COMPTE LES
! --- PERTES PAR FROTTEMENT ET PAR RECUL DU(DES) ANCRAGE(S)
!    PAS DE DIFFERENCE ENTRE ETCC ET BPEL
!
        if (nbf0 .eq. 1) then
            call tensk1(icabl, nbnoca, zr(jabsc+idecno), zr(jalph+idecno), f0, &
                        delta, ea, frco, frli, sa, &
                        zr(jf))
        else
            call tensk2(icabl, nbnoca, zr(jabsc+idecno), zr(jalph+idecno), f0, &
                        delta, ea, frco, frli, sa, &
                        zr(jf))
        end if
!
!
! 4.1.2 PRISE EN COMPTE LE CAS ECHEANT DES PERTES DE TENSION PAR
! --- RELAXATION DE L'ACIER
!
!    Verification que le parametre de relaxation n'est pas nul
        if ((typrel .ne. 'SANS') .or. (analy .eq. 'ETCC')) then
            if (rh1000 .le. r8prem()) then
                call utmess('A', 'CABLE0_12')
            end if
        end if
!
        if (typrel .eq. 'BPEL') then
!----------------------------------
!     4.1.2.1 CAS DU BPEL
!-----------------------
            flim = fprg*sa
            krelax = trelax*5.0d-02*rh1000
!
            do ino = 1, nbnoca
                zr(jf+ino-1) = zr(jf+ino-1)*(1.0d0-krelax*(zr(jf+ino-1)/flim-mu0))
            end do
!
        else if (typrel .eq. 'ETCC_DIRECT') then
!
!----------------------------------
!     4.1.2.2   CAS ETCC_DIRECT
!----------------------------------
            flim = fprg*sa
            do ino = 1, nbnoca
                fi = zr(jf+ino-1)
                zr(jf+ino-1) = fi-0.8d0*fi*0.66d-05*rh1000*exp(9.1d0*fi/flim)*(trelax/10&
                               &00.d0)**(0.75d0*(1.d0-(fi/flim)))
!
            end do
        else if (analy .eq. 'ETCC') then
!----------------------------------
!     4.1.2.3   CAS ETCC_INDIRECT (MODI_CABLE_ETCC)
!----------------------------------
            flim = fprg*sa
            do ino = 1, nbnoca
                f2 = zr(jtaby+ino-1)
                zr(jf+ino-1) = zr(jf+ino-1)-0.8d0*0.66d-05*rh1000*exp(9.1d0*f2/flim)*(trel&
                               &ax/1000.d0)**(0.75d0*(1.d0-(f2/flim)))*f2
            end do
!
            call jedetr(tabx)
            call jedetr(taby)
!
!
        end if
!
! 4.1.3 PRISE EN COMPTE LE CAS ECHEANT DES PERTES DE TENSION PAR
! --- FLUAGE ET RETRAIT DU BETON - UNIQUEMENT POUR BPEL
!
        if (regl .eq. 'BPEL') then
!
            if (xflu+xret .ne. 0.0d0) then
                df = (xflu+xret)*f0
                do ino = 1, nbnoca
                    zr(jf+ino-1) = zr(jf+ino-1)-df
                end do
            end if
!
        end if
!
    else
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    4.2  CAS RUPTURE CABLE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! LA TENSION APPLIQUEE EST DIRECTEMENT LA TENSION FOURNIE PAR L'UTILISATEUR
        do ino = 1, nbnoca
            f2 = zr(jtaby+ino-1)
            zr(jf+ino-1) = f2
        end do
!
        call jedetr(tabx)
        call jedetr(taby)
!
    end if
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 5   MISE A JOUR DES OBJETS DE SORTIE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    do ino = 1, nbnoca
        call tbajli(tablca, 1, param, [ibid], zr(jf+ino-1), &
                    [cbid], k3b, idecno+ino)
    end do
!
999 continue
    call jedetr('&&TENSCA.F')
    call jedema()
!
! --- FIN DE TENSCA.
end subroutine
