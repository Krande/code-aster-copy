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
subroutine acgrcr(nbvec, jvectn, jvectu, jvectv, nbordr, &
                  kwork, sompgw, jrwork, tspaq, ipg, &
                  nommet, jvecno, jnorma, forcri, nompar, &
                  vanocr, respc, vnmax)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8pi.h"
#include "asterc/r8prem.h"
#include "asterfort/acplcr.h"
#include "asterfort/anacri.h"
#include "asterfort/fointe.h"
#include "asterfort/fonbpa.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeveuo.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: nbvec, jvectn, jvectu, jvectv, nbordr, kwork
    integer(kind=8) :: sompgw, jrwork, tspaq, ipg, jvecno, jnorma
    character(len=16) :: nommet, forcri
    character(len=8) :: nompar(35)
    real(kind=8) :: respc(24), vnmax(6), vanocr(23)
!
!
! BUT: POUR LA FATIGUE A AMPLITUDE CONSTANTE
!      DETERMINER LE PLAN DES MAX DES TAU_MAX ET CALCULER DES GRANDEURS
!
!
! REMARQUE: CETTE SUBROUTINE EST APPLICABLE POUR UN NOEUD OU IPG EGALE
!           A 1 ET SOMPGW = SOMNOW,JVECPG = JVECNO
! ----------------------------------------------------------------------
! ARGUMENTS :
!     JVECTN  : IN  : ADRESSE DU VECTEUR CONTENANT LES COMPOSANTES DES
!                     VECTEURS NORMAUX.
!     JVECTU  : IN  : ADRESSE DU VECTEUR CONTENANT LES COMPOSANTES DES
!                     VECTEURS u DU PLAN DE CISAILLEMENT.
!     JVECTV  : IN  : ADRESSE DU VECTEUR CONTENANT LES COMPOSANTES DES
!                     VECTEURS v DU PLAN DE CISAILLEMENT.
!     NBORDR  : IN  : NOMBRE DE NUMEROS D'ORDRE.
!     KWORK   : IN  : KWORK = 0 ON TRAITE LA 1ERE MAILLE DU PAQUET DE
!                               MAILLES ;
!                     KWORK = 1 ON TRAITE LA IEME (I>1) MAILLE DU PAQUET
!                               MAILLES.
!     SOMPGW  : IN  : SOMME DES POINTS DE GAUSS DES N MAILLES PRECEDANT
!                     LA MAILLE COURANTE.
!     JRWORK  : IN  : ADRESSE DU VECTEUR DE TRAVAIL CONTENANT
!                     L'HISTORIQUE DES TENSEURS DES CONTRAINTES
!                     ATTACHES A CHAQUE POINT DE GAUSS DES MAILLES
!                     DU <<PAQUET>> DE MAILLES.
!     TSPAQ   : IN  : TAILLE DU SOUS-PAQUET DU <<PAQUET>> DE MAILLES
!                     COURANT.
!     IPG     : IN  : IEME POINT DE GAUSS.
!    NOMMET     IN    NOM DE METHOD D'APPROCHEMENT DE CERCLE ("CERCLE
!                     EXACT" ET "CERCLE APPROCHE")
!    VALA       IN    VALEUR DU PARAMETRE a ASSOCIE AU CRITERE.
!    COEFPA     IN    COEFFICIENT DE PASSAGE CISAILLEMENT - UNIAXIAL.
!   VRSESU      OUT   TABLEAU DES RESULTATS (GRANDEURS ET DOMMAGE).
!                     POUR L'INSTANT, LA DIMENSION DE VRESU EST 24
! ----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
    integer(kind=8) :: i, j, k
    integer(kind=8) :: mnmax(2)
!
    integer(kind=8) :: nparma, praccr(35)
    integer(kind=8) :: ibid, jresun
    integer(kind=8) :: jprof, np, dectau, ipar
    integer(kind=8) :: jdtaum, jtauma, jsgnma, jdsgma
    integer(kind=8) :: jdgama, jgamma, jepnma, jdenma
    integer(kind=8) :: jdgpma, jgapma, jeppma, jdepma
!
    real(kind=8) :: epsilo, pi
    real(kind=8) :: valpar(35), valpu(35)
    real(kind=8) :: grcrma(2), gcmax
    character(len=16) :: typch, nomcr
    character(len=24) :: chnom, cbid
    character(len=8) :: nompf(35)
    aster_logical :: lbid, crsigm, crepst, crepse, crepsp, rayon
!
    epsilo = 10*r8prem()
    pi = r8pi()
!
    typch = 'PERIODIQUE'
    nomcr = 'FORMULE_CRITERE'
!
    call anacri(nomcr, forcri, typch, 'NON', praccr, &
                lbid, crsigm, crepst, crepse, crepsp)
!
    rayon = .false.
    if ((praccr(24) .eq. 1) .or. (praccr(25) .eq. 1) .or. (praccr(32) .eq. 1)) then
        rayon = .true.
    end if
! initialisation
    do i = 1, 35
        valpar(i) = 0.0d0
    end do
!
! Récuperer les paramètres qui ne dépendent pas de plan
    do i = 7, 23
        valpar(i) = vanocr(i)
    end do
!
!
    call wkvect('&&ACGRCR.RESU_N', 'V V I', nbvec, jresun)
!
    if (crsigm) then
!
        call wkvect('&&ACGRCR.DTAU_MAX', 'V V R', nbvec, jdtaum)
        call wkvect('&&ACGRCR.TAUMAX', 'V V R', nbvec, jtauma)
        call wkvect('&&ACGRCR.SGNMAX', 'V V R', nbvec, jsgnma)
        call wkvect('&&ACGRCR.DSNMAX', 'V V R', nbvec, jdsgma)
!
        dectau = 0
        call acplcr(nbvec, jvectn, jvectu, jvectv, nbordr, &
                    kwork, sompgw, jrwork, tspaq, ipg, &
                    dectau, nommet, jvecno, jnorma, rayon, &
                    jresun, jdtaum, jtauma, jsgnma, jdsgma)
!
    end if
!
    if (crepst) then
        call wkvect('&&ACGRCR.DGAMAX', 'V V R', nbvec, jdgama)
        call wkvect('&&ACGRCR.GAMMAX', 'V V R', nbvec, jgamma)
        call wkvect('&&ACGRCR.EPNMAX', 'V V R', nbvec, jepnma)
        call wkvect('&&ACGRCR.DENMAX', 'V V R', nbvec, jdenma)
!
        dectau = 6
        call acplcr(nbvec, jvectn, jvectu, jvectv, nbordr, &
                    kwork, sompgw, jrwork, tspaq, ipg, &
                    dectau, nommet, jvecno, jnorma, rayon, &
                    jresun, jdgama, jgamma, jepnma, jdenma)
!
    end if
!
    if (crepsp) then
        call wkvect('&&ACGRCR.DGPMAX', 'V V R', nbvec, jdgpma)
        call wkvect('&&ACGRCR.GAPMAX', 'V V R', nbvec, jgapma)
        call wkvect('&&ACGRCR.EPPMAX', 'V V R', nbvec, jeppma)
        call wkvect('&&ACGRCR.DEPMAX', 'V V R', nbvec, jdepma)
!
        dectau = 12
        call acplcr(nbvec, jvectn, jvectu, jvectv, nbordr, &
                    kwork, sompgw, jrwork, tspaq, ipg, &
                    dectau, nommet, jvecno, jnorma, rayon, &
                    jresun, jdgpma, jgapma, jeppma, jdepma)
!
    end if
!
!
! 3/ CDU 1ER MAX DES DELTA_TAU ET DU VECTEUR NORMAL ASSOCIE
!
    grcrma(1) = r8prem()
    grcrma(2) = r8prem()
    mnmax(1) = 1
    mnmax(2) = 1
!
    gcmax = 0.0d0
!
!  RECUPERER LES NOMS DE PARAMETRES FOURNIS PAR L'UTILISATEUR
!  NOMBRE DE PARAMETRES DISPONIBLES
!
    nparma = 35
!
!
    chnom(20:24) = '.PROL'
    chnom(1:19) = forcri
!
    call jeveuo(chnom, 'L', jprof)
    call fonbpa(forcri, zk24(jprof), cbid, nparma, np, &
                nompf)
!
!
    do i = 1, nbvec
        if (crsigm) then
            valpar(24) = zr(jdtaum+i-1)
            valpar(26) = zr(jdsgma+i-1)
            valpar(28) = zr(jtauma+i-1)
            valpar(30) = zr(jsgnma+i-1)
        end if
!
        if (crepst) then
!! POUR ENGEERING STRAIN
            valpar(25) = zr(jdgama+i-1)*2
            valpar(27) = zr(jdenma+i-1)
            valpar(29) = zr(jgamma+i-1)*2
            valpar(31) = zr(jepnma+i-1)
        end if
        if (crepsp) then
            valpar(32) = zr(jdgpma+i-1)*2
            valpar(33) = zr(jdepma+i-1)
            valpar(34) = zr(jgapma+i-1)*2
            valpar(35) = zr(jeppma+i-1)
        end if
!
        do j = 1, np
            do ipar = 1, nparma
                if (nompf(j) .eq. nompar(ipar)) then
                    valpu(j) = valpar(ipar)
                    goto 30
                end if
            end do
30          continue
        end do
        !
        call fointe('F', forcri, np, nompf, valpu, &
                    gcmax, ibid)
!
        if (gcmax .gt. epsilo) then
            if ((gcmax-grcrma(1)) .gt. epsilo) then
                grcrma(2) = grcrma(1)
                mnmax(2) = mnmax(1)
                grcrma(1) = gcmax
                mnmax(1) = i
!
            end if
            if ((abs(grcrma(2)-grcrma(1)) .lt. epsilo) .and. (i .ne. mnmax(1))) then
                grcrma(2) = gcmax
                mnmax(2) = i
            end if
        end if
    end do
!!!RECUPER DES VALUERS DE LA GRANDEUR CRITIQUE
    do k = 1, 2
        i = mnmax(k)
!
        if (crsigm) then
            valpar(24) = zr(jdtaum+i-1)
            valpar(26) = zr(jdsgma+i-1)
            valpar(28) = zr(jtauma+i-1)
            valpar(30) = zr(jsgnma+i-1)
        end if
!
!! *2 POUR ENGEERING STRAIN
        if (crepst) then
            valpar(25) = zr(jdgama+i-1)*2
            valpar(27) = zr(jdenma+i-1)
            valpar(29) = zr(jgamma+i-1)*2
            valpar(31) = zr(jepnma+i-1)
        end if
!
        if (crepsp) then
            valpar(32) = zr(jdgpma+i-1)*2
            valpar(33) = zr(jdepma+i-1)
            valpar(34) = zr(jgapma+i-1)*2
            valpar(35) = zr(jeppma+i-1)
        end if
!
        do j = 1, 12
            respc(j+(k-1)*12) = valpar(23+j)
        end do
!
!!!ORIENTATION DU PLAN CRITIQUE
        vnmax(3*(k-1)+1) = zr(jvectn+(i-1)*3)
        vnmax(3*(k-1)+2) = zr(jvectn+(i-1)*3+1)
        vnmax(3*(k-1)+3) = zr(jvectn+(i-1)*3+2)
!
    end do
!
    if (crsigm) then
        call jedetr('&&ACGRCR.TAUMAX')
        call jedetr('&&ACGRCR.SGNMAX')
        call jedetr('&&ACGRCR.DSNMAX')
    end if
!
    if (crepst) then
        call jedetr('&&ACGRCR.DGAMAX')
        call jedetr('&&ACGRCR.GAMMAX')
        call jedetr('&&ACGRCR.EPNMAX')
        call jedetr('&&ACGRCR.DENMAX')
    end if
!
    if (crepsp) then
        call jedetr('&&ACGRCR.DGPMAX')
        call jedetr('&&ACGRCR.GAPMAX')
        call jedetr('&&ACGRCR.EPPMAX')
        call jedetr('&&ACGRCR.DEPMAX')
    end if
!
    call jedetr('&&ACGRCR.DTAU_MAX')
    call jedetr('&&ACGRCR.RESU_N')
!
end subroutine
