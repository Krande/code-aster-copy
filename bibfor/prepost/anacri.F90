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
subroutine anacri(nomcri, nomfor, typcha, impgrd, paract, &
                  fordef, crsigm, crepst, crepse, crepsp)
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/fonbpa.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
    character(len=3) :: impgrd
    character(len=16) :: nomcri, nomfor, typcha
    integer(kind=8) :: paract(35)
    aster_logical :: fordef, crsigm, crepst, crepse, crepsp
!
! ---------------------------------------------------------------------
! BUT: ANALYSER LE CRITERE POUR DETERMINER LES GRANDEURS NECESSAIARES
!                                A CALCUCLER
! ---------------------------------------------------------------------
! ARGUMENTS:
! NOMCRI   IN    K16: NOM DU CRITERE D'ENDOMMAGEMENT PAR FATIGUE.
! NOMFOR   IN    K16: LE NOM DE FORMULE DE GRNADUER EQUIVALENTE
! IMPGRD   IN    K3 : 'OUI' : IMPRIMER LES GRANDEURS A CALCULER
!                     'NON':  PAS IMPRIMER
! TYPCHA   IN    K16: TYPE DE CHARGEMENT (PERIODIQUE OU NON).
! PARACT   OUT   REAL: INDICATEUR DU GRANDEUR ACTIVE
!                      PARACT(K) = 1: K-IEME GRANDEUR EST ACTIVE
!
! FORDEF  LOGICAL : 'OUI' POUR LA PROJECTION DE L'HISTOIRE
!       DE DEFORMATION CISSAILLEMEMENT ET 'NON' POUR LA PROJECTION DE
!                   L'HISTOIRE DE CONTRAINET CISSAILLEMEMENT
! CRSIGM  LOGICAL : HISTOIRE DE CONTRAINTE NECESSAIRE
! CREPST  LOGICAL : HISTOIRE DE DEFORMATION TOTALE NECESSAIRE
! CREPSE  LOGICAL : HISTOIRE DE DEFORMATION TOTALE NECESSAIRE
! CREPSP  LOGICAL : HISTOIRE DE DEFORMATION PLASTIQUE NECESSAIRE
!
!-----------------------------------------------------------------------
    integer(kind=8) :: ip, id, nparma, nparm2, jprof, np, l
    character(len=2) :: nomty1(35), nomty2(30)
    character(len=8) :: nompa1(35), nompa2(30), nompf(35)
    character(len=24) :: chnom, cbid
    aster_logical :: grdexi
!
!     ---------------------------------------------------------------
    data nompa1/'DTAUMA', 'PHYDRM', 'NORMAX', 'NORMOY',&
     &                  'EPNMAX', 'EPNMOY', 'DEPSPE', 'EPSPR1',&
     &                  'SIGNM1', 'DENDIS', 'DENDIE', 'APHYDR',&
     &                  'MPHYDR', 'DSIGEQ', 'SIGPR1', 'EPSNM1',&
     &                  'INVA2S', 'DSITRE', 'DEPTRE', 'EPSPAC',&
     &                  'RAYSPH', 'AMPCIS', 'DEPSEE',&
     &                  'DTAUCR', 'DGAMCR', 'DSINCR', 'DEPNCR', &
     &                  'MTAUCR', 'MGAMCR', 'MSINCR', 'MEPNCR', &
     &                  'DGAMPC', 'DEPNPC', 'MGAMPC', 'MEPNPC'/
!     ---------------------------------------------------------------
!     ---------------------------------------------------------------
!      C = CONTRAINTE, T = DEF TOTALE, E = DEF ELAS, P = DEF PLAS
!
!
    data nomty1/'CC', 'CC', 'CC', 'CC',&
     &                  'TT', 'TT', 'PP', 'TT',&
     &                  'CT', 'CP', 'CE', 'CC',&
     &                  'CC', 'CC', 'CC', 'TC',&
     &                  'TT', 'CC', 'TT', 'PP',&
     &                  'CC', 'CC', 'EE',&
     &                  'CC', 'TT', 'CC', 'TT',&
     &                  'CC', 'TT', 'CC', 'TT',&
     &                  'PP', 'PP', 'PP', 'PP'/
!     ---------------------------------------------------------------
!     ---------------------------------------------------------------
    data nompa2/'TAUPR_1', 'TAUPR_2', 'SIGN_1', 'SIGN_2',&
     &                 'PHYDR_1', 'PHYDR_2', 'EPSPR_1', 'EPSPR_2',&
     &                 'SIPR1_1', 'SIPR1_2', 'EPSN1_1', 'EPSN1_2',&
     &                 'ETPR1_1', 'ETPR1_2', 'SITN1_1', 'SITN1_2',&
     &                 'EPPR1_1', 'EPPR1_2', 'SIPN1_1', 'SIPN1_2',&
     &                 'SIGEQ_1', 'SIGEQ_2', 'ETEQ_1', 'ETEQ_2',&
     &                 'EPEQ_1', 'EPEQ_2', 'INVJ2_1', 'INVJ2_2',&
     &                 'SITRE_1', 'SITRE_2'/
!       -------------------------------------------------------------
    data nomty2/'CC', 'CC', 'CC', 'CC',&
     &                 'CC', 'CC', 'TT', 'TT',&
     &                 'CC', 'CC', 'TC', 'TC',&
     &                 'TT', 'TT', 'CT', 'CT',&
     &                 'PP', 'PP', 'CP', 'CP',&
     &                 'CC', 'CC', 'TT', 'TT',&
     &                 'PP', 'PP', 'TT', 'TT',&
     &                 'CC', 'CC'/
!
!
!-----------------------------------------------------------------------
!234567                                                              012
!
    call jemarq()
!
!
    nparm2 = 30
!
! NOMBRE MAX DE PARAMETRES DISPONIBLES
    nparma = 35
!
!     INITIALISATION
    do ip = 1, nparma
        paract(ip) = 0
    end do
!
    fordef = .false.
!
    if (nomcri(1:7) .eq. 'FORMULE') then
!
! RECUPERER LES NOMS DE PARAMETRES FOURNIS PAR L'UTILISATEUR
        chnom(20:24) = '.PROL'
        chnom(1:19) = nomfor
!
        call jeveuo(chnom, 'L', jprof)
        call fonbpa(nomfor, zk24(jprof), cbid, nparma, np, &
                    nompf)
!
! VERIFIER QUE LE NOM DE GRANDEUR A CALCULER EST BON
        if (typcha(1:14) .eq. 'NON_PERIODIQUE') then
            do id = 1, np
                grdexi = .false.
                do ip = 1, nparm2
                    if (nompf(id) .eq. nompa2(ip)) then
                        grdexi = .true.
                        paract(ip) = 1
                    end if
                end do
                if (.not. grdexi) then
                    call utmess('F', 'FATIGUE1_91', sk=nompf(id))
                end if
!
                if (nompf(id) (1:3) .eq. 'EPS') then
                    fordef = .true.
                    do ip = 1, np
                        if (nompf(ip) (1:3) .eq. 'TAU') then
                            call utmess('F', 'FATIGUE1_92')
                        end if
                    end do
                end if
                if (nompf(id) (1:3) .eq. 'TAU') then
                    do ip = 1, np
                        if (nompf(ip) (1:3) .eq. 'EPS') then
                            call utmess('F', 'FATIGUE1_92')
                        end if
                    end do
                end if
            end do
!
        end if
!
        if (typcha(1:10) .eq. 'PERIODIQUE') then
            do id = 1, np
                grdexi = .false.
                do ip = 1, nparma
                    if (nompf(id) .eq. nompa1(ip)) then
                        grdexi = .true.
                        paract(ip) = 1
                    end if
                end do
!
                if (.not. grdexi) then
                    call utmess('F', 'FATIGUE1_91', sk=nompf(id))
                end if
!
            end do
        end if
!
    end if
!
    if (nomcri(1:14) .eq. 'MATAKE_MODI_AC') then
        paract(1) = 1
        paract(3) = 1
        paract(4) = 1
        paract(5) = 1
        paract(6) = 1
    end if
!
    if (nomcri(1:16) .eq. 'DANG_VAN_MODI_AC') then
        paract(1) = 1
        paract(2) = 1
        paract(4) = 1
        paract(5) = 1
        paract(6) = 1
    end if
!
    if (nomcri(1:14) .eq. 'MATAKE_MODI_AV') then
        paract(1) = 1
        paract(2) = 1
        paract(3) = 1
        paract(4) = 1
    end if
!
    if (nomcri(1:16) .eq. 'DANG_VAN_MODI_AV') then
        paract(1) = 1
        paract(2) = 1
        paract(5) = 1
        paract(6) = 1
    end if
!
    if (nomcri(1:16) .eq. 'FATESOCI_MODI_AV') then
        paract(3) = 1
        paract(4) = 1
        paract(7) = 1
        paract(8) = 1
    end if
!
    if (nomcri(1:11) .eq. 'VMIS_TRESCA') then
        paract(14) = 1
        paract(18) = 1
    end if
!
! DANS POST_FATIGUE
    if (nomcri(1:9) .eq. 'CROSSLAND') then
        paract(2) = 1
        paract(22) = 1
    end if
!
    if (nomcri(1:12) .eq. 'PAPADOPOULOS') then
        paract(2) = 1
        paract(21) = 1
    end if
! POUR OPERATEUR POST_FATIGUE
!
! ANALYSER LES HISTORES NECESSAIRE
    crsigm = .false.
    crepst = .false.
    crepse = .false.
    crepsp = .false.
!
    do ip = 1, nparma
        if (paract(ip) .eq. 1) then
            do l = 1, 2
                if (typcha .eq. 'PERIODIQUE') then
                    if (nomty1(ip) (l:l) .eq. 'C') then
                        crsigm = .true.
                    end if
                    if (nomty1(ip) (l:l) .eq. 'T') then
                        crepst = .true.
                    end if
                    if (nomty1(ip) (l:l) .eq. 'E') then
                        crepse = .true.
                    end if
                    if (nomty1(ip) (l:l) .eq. 'P') then
                        crepsp = .true.
                    end if
                else
                    if (nomty2(ip) (l:l) .eq. 'C') then
                        crsigm = .true.
                    end if
                    if (nomty2(ip) (l:l) .eq. 'T') then
                        crepst = .true.
                    end if
                    if (nomty2(ip) (l:l) .eq. 'E') then
                        crepse = .true.
                    end if
                    if (nomty2(ip) (l:l) .eq. 'P') then
                        crepsp = .true.
                    end if
                end if
            end do
        end if
    end do
!
! IMPRIMER DES INFO
    if (impgrd .eq. 'OUI') then
        write (6, *) 'CRITERE AMORCAGE A UTILISER ==>', nomcri
        write (6, *) ' '
        write (6, *) 'LES GRANDEURS A CALCULER : '
        do ip = 1, nparma
            if (paract(ip) .eq. 1) then
!
                if (typcha .eq. 'PERIODIQUE') then
                    write (6, *) '    ', nompa1(ip)
                    write (6, *) ' '
                else
                    write (6, *) '    ', nompa2(ip)
                    write (6, *) ' '
                end if
            end if
!
        end do
!
        write (6, *) 'HISTOIRES DE CHARGEMENT DOIVENT CONSISTER :'
!
        if (crsigm) then
            write (6, *) '    CONTRAINTE'
        end if
!
        if (crepst) then
            write (6, *) '    DEFORMATION TOTALE'
        end if
!
        if (crepse) then
            write (6, *) '    DEFORMATION ELASTIQUE'
        end if
!
        if (crepsp) then
            write (6, *) '    DEFORMATION PLASTIQUE'
        end if
        write (6, *) ' '
!
        if (crepse) then
            write (6, *) 'ON NOTE: DEFORMATION ELASTIQUE = DEFORMATION'//&
                &'TOTALE - DEFORMATION PLASTIQUE'
            if (.not. crepst) then
                write (6, *) 'LE CHARGEMENT DOIT CONSISTER EN PLUS:'//&
                    &'DEFORMATION TOTALE (OBLIGATOIRE)'
            end if
!
            if (.not. crepsp) then
                write (6, *) 'LE CHARGEMENT DOIT CONSISTER EN PLUS:'//&
                    &'DEFORMATION PLASTIQUE (OPTIONEL)'
            end if
!
        end if
    end if
!
    call jedema()
end subroutine
