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

subroutine amumps(action, kxmps, rsolu, vcine, nbsol, &
                  iret, impr, ifmump, prepos, pcentp)
!
!
    implicit none
!--------------------------------------------------------------
! OBJET: DRIVER EN MODE REEL DE LA RESOLUTION DE SYSTEMES LINEAIRES
!        VIA MUMPS (EN SIMPLE PRECISION POUR MUMPS UNIQUEMENT)
!
! IN : ACTION :
!     /'PRERES'  : POUR DEMANDER LA FACTORISATION
!     /'RESOUD'  : POUR DEMANDER LA DESCENTE/REMONTEE
!     /'DETR_MAT/OCC': POUR DEMANDER LA DESTRUCTION DE L'INSTANCE MUMPS
!                  ASSOCIEE A UNE MATRICE
!
! IN : KXMPS (I)   : INDICE DE L'INSTANCE MUMPS DANS SMPS
! VAR: RSOLU (R)   : EN ENTREE : VECTEUR SECOND MEMBRE (REEL)
!                    EN SORTIE : VECTEUR SOLUTION (REEL)
!            (SI ACTION=RESOUD)
! IN : VCINE (K19) : NOM DU CHAM_NO DE CHARGEMENT CINEMATIQUE
!            (SI ACTION=RESOUD)
! OUT : IRET (I) : CODE_RETOUR :
!            0 : OK
!            1 : ERREUR (DANS LE CAS OU MUMPS EST UTILISE EN PRE_COND)
!            2 : MATRICE NUMERIQUEMENT SINGULIERE
! IN  : NBSOL  : NBRE DE SYSTEMES A RESOUDRE
! IN  : IMPR,IFMUMP : PARAMETRES POUR SORTIE FICHIER MATRICE CF AMUMPH
! IN : PREPOS (LOG) : SI .TRUE. ON FAIT LES PRE ET POSTTRAITEMENTS DE
!           MISE A L'ECHELLE DU RHS ET DE LA SOLUTION (MRCONL) ET DE LA
!           PRISE EN COMPTE DES AFFE_CHAR_CINE (CSMBGG).
!           SI .FALSE. ON NE LES FAIT PAS (PAR EXEMPLE EN MODAL).
! IN  : PCENTP VECTEUR D'ENTIER GERE PAR AMUMPH POUR PARAMETRER LES
!                STRATEGIES D'ADAPTATION EN CAS DE PB PCENT_PIVOT
!---------------------------------------------------------------
! person_in_charge: olivier.boiteau at edf.fr
!
#include "asterc/matfpe.h"
#include "asterf_types.h"
#include "asterf.h"
#include "asterfort/amumpi.h"
#include "asterfort/amumpm.h"
#include "asterfort/amumpp.h"
#include "asterfort/amumpt.h"
#include "asterfort/amumpu.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
#include "mumps/smumps.h"
#include "asterfort/isParallelMatrix.h"
!
    character(len=*) :: action
    character(len=14) :: impr
    character(len=19) :: vcine, nosolv
    integer(kind=8) :: iret, nbsol, kxmps, ifmump, pcentp(2)
    real(kind=8) :: rsolu(*)
    aster_logical :: prepos
!
#ifdef ASTER_HAVE_MUMPS
#include "asterf_mumps.h"
    type(smumps_struc), pointer :: smpsk => null()
    integer(kind=8) :: rang, nbproc, niv, ifm, ibid, ietdeb, ifactm, nbfact
    integer(kind=8) :: ietrat, nprec, ifact, iaux, vali(4), pcpi
    character(len=1) :: rouc, type, prec
    character(len=3) :: matd
    character(len=5) :: etam, klag2
    character(len=8) :: ktypr
    character(len=12) :: usersm, k12bid
    character(len=14) :: nonu
    character(len=19) :: nomat
    character(len=24) :: kmonit(12), k24aux, posttrait
    real(kind=8) :: epsmax, valr(2), rctdeb, temps(6), epsmat
    complex(kind=8) :: cbid(1)
    aster_logical :: lquali, ldist, lresol, lmd, lbid, lpreco, lbis, lpb13, ldet
    aster_logical :: lopfac, l_parallel_matrix, lbloc
    character(len=24), pointer :: refa(:) => null()
    character(len=24), pointer :: slvk(:) => null()
    real(kind=8), pointer :: slvr(:) => null()
    integer(kind=8), pointer :: slvi(:) => null()
    call jemarq()
!
!       ------------------------------------------------
!        INITS
!       ------------------------------------------------
! --- ON DESACTIVE LA LEVEE D'EXCEPTION FPE DANS LA BIBLIOTHEQUE MKL
! --  CAR CES EXCEPTIONS NE SONT PAS JUSTIFIEES
    call matfpe(-1_8)
    call infdbg('SOLVEUR', ifm, niv)
!
! --- PARAMETRE POUR IMPRESSION FICHIER
    lresol = ((impr(1:3) .eq. 'NON') .or. (impr(1:9) .eq. 'OUI_SOLVE'))
!
! --- TYPE DE SYSTEME: REEL OU COMPLEXE
    type = 'S'
    ASSERT(kxmps .gt. 0)
    ASSERT(kxmps .le. nmxins)
    nomat = nomats(kxmps)
    nosolv = nosols(kxmps)
    nonu = nonus(kxmps)
    etam = etams(kxmps)
    rouc = roucs(kxmps)
    prec = precs(kxmps)
    ASSERT((rouc .eq. 'R') .and. (prec .eq. 'S'))
    smpsk => smps(kxmps)
    iret = 0
!
! --- MUMPS PARALLELE DISTRIBUE ?
    call jeveuo(nomat//'.REFA', 'L', vk24=refa)
    ldist = (refa(11) .ne. 'MPI_COMPLET')
    rang = smpsk%myid
    nbproc = smpsk%nprocs
!
! --- MATRICE ASTER DISTRIBUEE ?
    call dismoi('MATR_DISTRIBUEE', nomat, 'MATR_ASSE', repk=matd)
    lmd = matd .eq. 'OUI'
!
! --- MATRICE ASTER HPC ?
    l_parallel_matrix = isParallelMatrix(nomat)
!
    lquali = .false.
!
    lbloc = .false.
    call jeexin(nosolv//'.SLVK', ibid)
    if (ibid .ne. 0) then
        call jeveuo(nosolv//'.SLVK', 'L', vk24=slvk)
! ---   ANALYSE PAR BLOCS
!       PAS ENCORE ETENDU AU MODE DISTRIBUE
        lbloc = (((slvk(5) (1:3) .eq. 'FR+') .or. (slvk(5) (1:3) .eq. 'LR+') .or. &
                  (slvk(5) (1:4) .eq. 'AUTO')) .and. (.not. l_parallel_matrix) .and. (.not. lmd))
    end if
!
    if (action(1:5) .ne. 'DETR_') then
        call jeveuo(nosolv//'.SLVK', 'E', vk24=slvk)
        call jeveuo(nosolv//'.SLVR', 'L', vr=slvr)
        call jeveuo(nosolv//'.SLVI', 'E', vi=slvi)
!
! --- L'UTILISATEUR VEUT-IL UNE ESTIMATION DE LA QUALITE DE LA SOL ?
! --- => LQUALI
        epsmax = slvr(2)
        posttrait = slvk(11)
        lquali = (epsmax .gt. 0.d0)
!
! --- POUR "ELIMINER" LE 2EME LAGRANGE :
! --- OPTION DEBRANCHEE SI CALCUL DE DETERMINANT
        klag2 = slvk(6) (1:5)
        lbis = klag2(1:5) .eq. 'LAGR2'
!
! --- TRES PROBABLEMENT COMMANDE FACTORISER (POSTTRAITEMENTS
! --- INITIALISE A 'XXXX'). ON NE DETRUIRA RIEN A L'ISSU DE LA
! --- FACTO, AU CAS OU UN OP. RESOUDRE + RESI_RELA>0 SUIVRAIT
        if (slvk(11) (1:4) .eq. 'XXXX') then
            lopfac = .true.
        else
            lopfac = .false.
        end if
!
! --- TYPE DE RESOLUTION
        ktypr = slvk(3) (1:8)
!
! --- PARAMETRE NPREC
        nprec = slvi(1)
!
! --- MUMPS EST-IL UTILISE COMME PRECONDITIONNEUR ?
! --- SI OUI, ON DEBRANCHE LES ALARMES ET INFO (PAS LES UTMESS_F)
        lpreco = slvk(8) (1:3) .eq. 'OUI'
!
! --- FILTRAGE DE LA MATRICE DONNEE A MUMPS (UNIQUEMENT NON LINEAIRE)
        epsmat = slvr(1)
!
! --- STRATEGIE MEMOIRE POUR MUMPS
        usersm = slvk(9) (1:12)
        nbfact = slvi(6)
    end if
!
! --- POUR MONITORING
    call amumpt(0, kmonit, temps, rang, nbproc, &
                kxmps, lquali, type, ietdeb, ietrat, &
                rctdeb, ldist)
!
!     ------------------------------------------------
!     ------------------------------------------------
    if (action(1:6) .eq. 'PRERES') then
!     ------------------------------------------------
!     ------------------------------------------------
!
!       ------------------------------------------------
!        INITIALISATION DE L'OCCURENCE MUMPS KXMPS:
!       ------------------------------------------------
        call amumpi(0, lquali, ldist, kxmps, type, l_parallel_matrix, lbid)
        call smumps(smpsk)

        rang = smpsk%myid
        nbproc = smpsk%nprocs
!
!       --------------------------------------------------------------
!        CHOIX ICNTL VECTEUR DE PARAMETRES POUR MUMPS (ANALYSE+FACTO):
!       --------------------------------------------------------------
        call amumpi(2, lquali, ldist, kxmps, type, l_parallel_matrix, lbloc)
!
!       -----------------------------------------------------
!       CALCUL DU DETERMINANT PART I ?
!       -----------------------------------------------------
        ldet = .false.
        if (slvi(5) .eq. 1) then
            if ((niv .ge. 2) .and. (lbis) .and. (.not. lpreco)) then
                call utmess('I', 'FACTOR_88')
            end if
            slvk(6) = 'NON'
            klag2 = 'NON'
            lbis = .false.
            ldet = .true.
        end if
!
!       ------------------------------------------------
!        REMPLISSAGE DE LA MATRICE MUMPS :
!       ------------------------------------------------
        call amumpt(1, kmonit, temps, rang, nbproc, &
                    kxmps, lquali, type, ietdeb, ietrat, &
                    rctdeb, ldist)
        call amumpm(ldist, kxmps, kmonit, impr, ifmump, &
                    klag2, type, lmd, epsmat, ktypr, &
                    lpreco, l_parallel_matrix, lbloc)
!
!       -----------------------------------------------------
!       CONSERVE-T-ON LES FACTEURS OU NON ?
!       -----------------------------------------------------
        if (slvi(4) .eq. 1) then
            smpsk%icntl(31) = 1
        end if
!
!       ------------------------------------------------
!        ANALYSE MUMPS:
!       ------------------------------------------------
!       INITIALISATIONS POUR ANALYSE+FACTO+CORRECTION EVENTUELLE
        ifact = 0
        lpb13 = .false.
        if (usersm(1:4) .eq. 'AUTO') then
            ifactm = pcentp(1)
        else
            ifactm = 1
        end if
!
10      continue
        call amumpt(2, kmonit, temps, rang, nbproc, &
                    kxmps, lquali, type, ietdeb, ietrat, &
                    rctdeb, ldist)
        smpsk%job = 1
        call smumps(smpsk)
        call amumpt(4, kmonit, temps, rang, nbproc, &
                    kxmps, lquali, type, ietdeb, ietrat, &
                    rctdeb, ldist)
!
!       ------------------------------------------------
!        GESTION ERREURS ET MENAGE ASTER:
!       ------------------------------------------------
        if (smpsk%infog(1) .eq. 0) then
!              -- C'EST OK
        else if ((smpsk%infog(1) .eq. -5) .or. (smpsk%infog(1) .eq. -7)) then
            call utmess('F', 'FACTOR_64')
        else if (smpsk%infog(1) .eq. -6) then
            iret = 2
            goto 99
        else if (smpsk%infog(1) .eq. -38) then
            call utmess('F', 'FACTOR_91')
        else if (smpsk%infog(1) .eq. -51) then
            call utmess('F', 'FACTOR_92')
        else
            iaux = smpsk%infog(1)
            if (iaux .lt. 0) then
                call utmess('F', 'FACTOR_55', si=iaux)
            else
                if ((.not. lpreco) .and. (niv .ge. 2)) then
                    call utmess('A', 'FACTOR_55', si=iaux)
                end if
            end if
        end if
        if ((slvk(4) .ne. 'AUTO') .and. (smpsk%icntl(7) .ne. smpsk%infog(7)) .and. &
            (.not. lpreco) .and. (smpsk%infog(32) .eq. 1) .and. (niv .ge. 2)) then
            call utmess('I', 'FACTOR_50', sk=slvk(4))
        end if
        if ((slvk(4) .ne. 'AUTO') .and. (smpsk%icntl(29) .ne. smpsk%infog(7)) .and. &
            (.not. lpreco) .and. (smpsk%infog(32) .eq. 2) .and. (niv .ge. 2)) then
            call utmess('I', 'FACTOR_50', sk=slvk(4))
        end if
!
!       -----------------------------------------------------
!        CHOIX DE LA STRATEGIE MUMPS POUR LA GESTION MEMOIRE
!       -----------------------------------------------------
        if (.not. lpb13) call amumpu(1, 'S', kxmps, usersm, ibid, &
                                     lbid, nbfact)
!
! ---   ON SORT POUR REVENIR A AMUMPH ET DETRUIRE L'OCCURENCE MUMPS
! ---   ASSOCIEE
        if (usersm(1:4) .eq. 'EVAL') goto 99
!
!       -----------------------------------------------------
!       CALCUL DU DETERMINANT PART II ?
!       -----------------------------------------------------
        if (ldet) smpsk%icntl(33) = 1
!
!       ------------------------------------------------
!        FACTORISATION NUMERIQUE MUMPS:
!       ------------------------------------------------
!
! --- SI GESTION_MEMOIRE='AUTO'
! --- ON TENTE PLUSIEURS (PCENTP(1)) FACTORISATIONS NUMERIQUES EN
! --- MULTIPLIANT, A CHAQUE ECHEC, L'ANCIEN PCENT_PIVOT PAR PCENTP(2)
! --- VOIRE EN PASSANT EN OOC (EN DERNIER RESSORT).
! --- AUTO-ADAPTATION DU PARAMETRAGE SOLVEUR/PCENT_PIVOT:
! --- ON MODIFIE LE PARAMETRE DANS LA SD_SOLVEUR A LA VOLEE POUR NE
! --- PAS PERDRE DE TEMPS LA PROCHAINE FOIS. CETTE VALEUR N'EST VALABLE
! --- QUE DANS L'OPERATEUR CONSIDERE.
! --- ON FAIT LA MEME CHOSE EN CAS DE PB D'ALLOCATION MEMOIRE (INFOG=-13
! --- CELA PEUT ETRE DU A UN ICNTL(23) MAL ESTIME
!
        smpsk%job = 2
        if (lresol) then
            pcpi = smpsk%icntl(14)
            do ifact = 1, ifactm
                call smumps(smpsk)
                iaux = smpsk%infog(1)
!
! --- TRAITEMENT CORRECTIF ICNTL(14)
                if ((iaux .eq. -8) .or. (iaux .eq. -9) .or. (iaux .eq. -14) &
                    .or. (iaux .eq. -15) .or. (iaux .eq. -17) .or. (iaux .eq. -20)) then
                    if (ifact .eq. ifactm) then
! ---  ICNTL(14): PLUS DE NOUVELLE TENTATIVE POSSIBLE
                        if (lpreco) then
!                 -- MUMPS EST APPELE COMME PRECONDITIONNEUR
!                 -- ON SORT AVEC UN CODE RETOUR NON NUL
                            iret = 1
                            goto 99
                        else
                            vali(1) = ifactm
                            vali(2) = pcpi
                            vali(3) = smpsk%icntl(14)
                            call utmess('F', 'FACTOR_53', ni=3_8, vali=vali)
                        end if
                    else
! ---  ICNTL(14): ON MODIFIE DES PARAMETRES POUR LA NOUVELLE TENTATIVE ET ON REVIENT A L'ANALYSE
                        smpsk%icntl(14) = smpsk%icntl(14)*to_mumps_int(pcentp(2))
                        slvi(2) = smpsk%icntl(14)
                        if ((niv .ge. 2) .and. (.not. lpreco)) then
                            vali(1) = smpsk%icntl(14)/pcentp(2)
                            vali(2) = smpsk%icntl(14)
                            vali(3) = ifact
                            vali(4) = ifactm
                            call utmess('I', 'FACTOR_58', ni=4_8, vali=vali)
                        end if
! --- DERNIERE CHANCE: ON RAJOUTE L'OOC
                        if (ifact .eq. (ifactm-1)) then
                            lpb13 = .true.
                            smpsk%icntl(23) = 0
                            smpsk%icntl(22) = 1
                        end if
                        ifactm = max(ifactm-ifact, 1_8)
                        goto 10
                    end if
!
! --- TRAITEMENT CORRECTIF ICNTL(23)
! --- CE N'EST UTILE QU' UNE FOIS D'OU LE CONTROLE DE LPB13
                else if (((iaux .eq. -13) .or. (iaux .eq. -19)) .and. (.not. lpb13)) then
! ---  ICNTL(23): ON MODIFIE DES PARAMETRES POUR LA NOUVELLE TENTATIVE ET ON REVIENT A L'ANALYSE
                    if ((niv .ge. 2) .and. (.not. lpreco)) then
                        vali(1) = smpsk%icntl(23)
                        call utmess('I', 'FACTOR_85', si=vali(1))
                    end if
                    lpb13 = .true.
                    smpsk%icntl(23) = 0
                    smpsk%icntl(22) = 1
                    ifactm = max(ifactm-ifact, 1_8)
                    goto 10
                else
! ---  SORTIE STANDARD SANS ERREUR
                    exit
                end if
            end do
        end if
!
! ---  AFFICHAGE DE CONTROLE
        if (niv .ge. 2) then
            write (ifm, *)
            write (ifm, *)&
     &      '<AMUMPS> FACTO. NUM. - NBRE TENTATIVES/MAX: ', ifact, ifactm
        end if
        call amumpt(6, kmonit, temps, rang, nbproc, &
                    kxmps, lquali, type, ietdeb, ietrat, &
                    rctdeb, ldist)
!
!       ------------------------------------------------
!        GESTION ERREURS ET MENAGE ASTER (SAUF ERREUR ICNTL(14/23)
!           TRAITEE EN AMONT):
!       ------------------------------------------------
        valr(1) = (smpsk%infog(13)*100.d0)/smpsk%n
        if (valr(1) .gt. 10.0d0) then
            if ((niv .ge. 2) .and. (.not. lpreco)) then
                call utmess('I', 'FACTOR_73')
            end if
        end if
        if (smpsk%infog(1) .eq. 0) then
!              -- C'EST OK
        else if (smpsk%infog(1) .eq. -10) then
            iret = 2
            goto 99
        else if (smpsk%infog(1) .eq. -13) then
            call utmess('F', 'FACTOR_54')
        else if (smpsk%infog(1) .eq. -37) then
            call utmess('F', 'FACTOR_65')
        else if (smpsk%infog(1) .eq. -90) then
            call utmess('F', 'FACTOR_66')
        else
            iaux = smpsk%infog(1)
            if (iaux .lt. 0) then
                call utmess('F', 'FACTOR_55', si=iaux)
            else
                if ((.not. lpreco) .and. (niv .ge. 2)) then
                    call utmess('A', 'FACTOR_55', si=iaux)
                end if
            end if
        end if
!
!       ------------------------------------------------
!        DETECTION DE SINGULARITE SI NECESSAIRE:
!       ------------------------------------------------
        call amumpu(2, 'S', kxmps, k12bid, nprec, &
                    lresol, ibid)
!
!       ------------------------------------------------
!        RECUPERATION DU DETERMINANT SI NECESSAIRE:
!       ------------------------------------------------
        call amumpu(4, 'S', kxmps, k12bid, ibid, &
                    lbid, ibid)
!
!       ON SOULAGE LA MEMOIRE JEVEUX DES QUE POSSIBLE D'OBJETS MUMPS
!       INUTILES
        if (((rang .eq. 0) .and. (.not. ldist)) .or. (ldist) .or. (l_parallel_matrix)) then
            if (.not. (lquali) .and. (posttrait(1:4) .ne. 'MINI') .and. .not. lopfac) then
                if (ldist .or. l_parallel_matrix) then
                    deallocate (smpsk%a_loc, stat=ibid)
                    deallocate (smpsk%irn_loc, stat=ibid)
                    deallocate (smpsk%jcn_loc, stat=ibid)
                else
                    deallocate (smpsk%a, stat=ibid)
                    deallocate (smpsk%irn, stat=ibid)
                    deallocate (smpsk%jcn, stat=ibid)
                end if
            end if
        end if
!
!     ------------------------------------------------
!     ------------------------------------------------
    else if (action(1:6) .eq. 'RESOUD') then
!     ------------------------------------------------
!     ------------------------------------------------
!
!       ------------------------------------------------
!        PRETRAITEMENTS ASTER DU/DES SECOND(S) MEMBRE(S) :
!       ------------------------------------------------
        call amumpt(7, kmonit, temps, rang, nbproc, &
                    kxmps, lquali, type, ietdeb, ietrat, &
                    rctdeb, ldist)
        call amumpp(0, nbsol, kxmps, ldist, type, &
                    impr, ifmump, lbis, rsolu, cbid, &
                    vcine, prepos, lpreco, l_parallel_matrix)
!
!       --------------------------------------------------------------
!        CHOIX ICNTL VECTEUR DE PARAMETRES POUR MUMPS (SOLVE):
!       --------------------------------------------------------------
        call amumpi(3, lquali, ldist, kxmps, type, l_parallel_matrix, lbid)
!
!       ------------------------------------------------
!        RESOLUTION MUMPS :
!       ------------------------------------------------
        call amumpt(8, kmonit, temps, rang, nbproc, &
                    kxmps, lquali, type, ietdeb, ietrat, &
                    rctdeb, ldist)
        smpsk%job = 3
        if (lresol) call smumps(smpsk)

        call amumpt(10, kmonit, temps, rang, nbproc, &
                    kxmps, lquali, type, ietdeb, ietrat, &
                    rctdeb, ldist)
!
!       ------------------------------------------------
!        GESTION ERREURS ET MENAGE ASTER:
!       ------------------------------------------------
        if (smpsk%infog(1) .eq. 0) then
!              -- C'EST OK
        else if ((smpsk%infog(1) .eq. 8) .and. (lquali)) then
            iaux = smpsk%infog(10)
            if ((.not. lpreco) .and. (niv .ge. 2)) then
                call utmess('A', 'FACTOR_62', si=iaux)
            end if
        else if (smpsk%infog(1) .lt. 0) then
            iaux = smpsk%infog(1)
            call utmess('F', 'FACTOR_55', si=iaux)
        else if (smpsk%infog(1) .eq. 4) then
!          -- PERMUTATION DE COLONNES, SMPSK%JCN MODIFIE VOLONTAIREMENT
!          -- PAR MUMPS. IL NE FAUT DONC PAS LE MANIPULER TEL QUE
!          -- PAS GRAVE POUR ASTER.
        else
            iaux = smpsk%infog(1)
            if ((.not. lpreco) .and. (iaux .ne. 2) .and. (niv .ge. 2)) then
                call utmess('A', 'FACTOR_55', si=iaux)
            end if
        end if
! --- CONTROLE DE L'ERREUR SUR LA SOLUTION :
        if ((lquali) .and. (posttrait(1:4) .ne. 'MINI')) then
            if (smpsk%rinfog(9) .gt. epsmax) then
                valr(1) = smpsk%rinfog(9)
                valr(2) = epsmax
                call utmess('F', 'FACTOR_57', nr=2_8, valr=valr)
            end if
        end if
!
!       ------------------------------------------------
!        POST-TRAITEMENTS ASTER DE/DES (LA) SOLUTION(S) :
!       ------------------------------------------------
        call amumpp(2, nbsol, kxmps, ldist, type, &
                    impr, ifmump, lbis, rsolu, cbid, &
                    vcine, prepos, lpreco, l_parallel_matrix)
!
!       ------------------------------------------------
!        AFFICHAGE DU MONITORING :
!       ------------------------------------------------
        call amumpt(12, kmonit, temps, rang, nbproc, &
                    kxmps, lquali, type, ietdeb, &
                    ietrat, rctdeb, ldist)
!
!     ------------------------------------------------
!     ------------------------------------------------
    else if (action(1:5) .eq. 'DETR_') then
!     ------------------------------------------------
!     ------------------------------------------------
!
!       ------------------------------------------------
!        MENAGE ASTER ET MUMPS:
!       ------------------------------------------------
        if (nomats(kxmps) .ne. ' ') then
            if (((rang .eq. 0) .and. (.not. ldist)) .or. (ldist) .or. (l_parallel_matrix)) then
                if (ldist .or. l_parallel_matrix) then
                    deallocate (smpsk%a_loc, stat=ibid)
                    deallocate (smpsk%irn_loc, stat=ibid)
                    deallocate (smpsk%jcn_loc, stat=ibid)
                else
                    deallocate (smpsk%a, stat=ibid)
                    deallocate (smpsk%irn, stat=ibid)
                    deallocate (smpsk%jcn, stat=ibid)
                end if
            end if
            if ((rang .eq. 0) .and. lbloc) then
                deallocate (smpsk%blkptr, stat=ibid)
            end if
            etams(kxmps) = ' '
            nonus(kxmps) = ' '
            nomats(kxmps) = ' '
            nosols(kxmps) = ' '
            roucs(kxmps) = ' '
            precs(kxmps) = ' '
            smpsk%job = -2
            call smumps(smpsk)
! NETTOYAGE OBJETS AUXILIAIRES AU CAS OU
            k24aux = '&&TAILLE_OBJ_MUMPS'
            call jeexin(k24aux, ibid)
            if (ibid .ne. 0) call jedetr(k24aux)
!
            k24aux = '&&AMUMP.PIVNUL'
            call jeexin(k24aux, ibid)
            if (ibid .ne. 0) call jedetr(k24aux)
!
            k24aux = '&&AMUMP.DETERMINANT'
            call jeexin(k24aux, ibid)
            if (ibid .ne. 0) call jedetr(k24aux)
        end if
    end if
!
!     -- ON REACTIVE LA LEVEE D'EXCEPTION
99  continue
    call matfpe(1_8)
    call jedema()
!
#endif
end subroutine
