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
! aslint: disable=W1501
!
! person_in_charge: vinicius.alves-fernandes at edf.fr
! contributor: cyril.borely at setec.com
!
subroutine difoncalc(tirela, raidTang, vloc, vpara, nbVloc, nbPara, iret, nbdecp, errmax)
!
! --------------------------------------------------------------------------------------------------
!  IN
!     tirela   : Tir élastique en entrée
!     nbVloc   : nombre de paramètre de la loi locale
!     vloc     : Paramètre de la loi locale
!     raidTang : Raideur tangeante
!     nbPara   : nombre variables locales
!     valvarloc: Variables locales en entrée
!     nbdecp   : Nombres d'itérations dans la boucle de convergence
!     errmax   : Erreur relative autorisée (indexée sur Vélas)
!
!  OUT
!     tirela   : Tir élastique corrigé en sortie
!     valvarloc: Variables locales modifiées en sortie
!     iret     : code de retour du matériau, on le met à 1 quand les convergences sont dépassées
! --------------------------------------------------------------------------------------------------
!
!    Contenu de vpara
!    data vpara /'LONG_X','LONG_Y',&
!                 'PHI','COHESION','RAID_GLIS',&
!                 'GAMMA_REFE','CP_SERVICE','CP_ULTIME','DEPL_REFE',&
!                 'RAID_CP_X','GAMMA_CP_X','RAID_CP_Y','GAMMA_CP_Y',&
!                 'RAID_CP_RX','GAMMA_CP_RX','RAID_CP_RY','GAMMA_CP_RY','DECOL'/
!
!    Contenu de vloc
!    data vloc ('DASLX','DASLY','DASLZ','DASLRX','DASLRY','DACPX','DACPY','DACPZ',
!              'DACPRX','DACPRY','DAQSLX','DAQSLY','DAQCPX','DAQCPY','DARCP','DAQCPRX','DAQCPRY',
!              'DAUTRE','DABSSLX','DABSSLY','DABSCPZ'),
! --------------------------------------------------------------------------------------------------
!
    implicit none
!
#include "asterc/r8prem.h"
#include "asterfort/mgauss.h"
!
    integer(kind=8) :: nbVloc, nbPara, iret, nbdecp
    real(kind=8) :: vpara(nbPara), tirela(6), raidTang(6), vloc(nbVloc), errmax
!
!   compteur du nombre d'itération
    integer(kind=8) :: Ite, IteCrit, ii, jj
!   nombre de critères à vérifier
    integer(kind=8) :: nbDDL, nbDDLii
    integer(kind=8) :: nbcomb, comppt, comptcombi
!
!   fcp     : valeur réelle du coefficient de sécuirté à la capacité portante
!   fslid   : idem pour le glissement
!   fsini   : vecteur qui sauvegarde les critère de rupture normé
!       fs (1) =fslid(H,Mx,My)
!       fs (2) =fslid(-H,Mx,My)
!       fs (3) =fslid(H,-Mx,My)
!       fs (4) =fslid(H,Mx,-My)
!       fs (5) =fCP(H,Mx,My)
!       fs (6) =fCP(-H,Mx,My)
!       fs (7) =fCP(H,-Mx,My)
!       fs (8) =fCP(H,Mx,-My)
!   fsret   : les écoulements plastiques correspodants
    real(kind=8) :: fcp, fslid, fsini(8), fsret(8), fsite(8), fslbis, fcpbis, fsitenow(8)
!
!   sHslIni     : force horizontale initiale en glissement
!   sMxslIni    : moment x initial en glissement
!   sMyslIni    : moment y initial en glissement
!   sHCPini     : force horizontale initiale en CP
!   sMxCPIni    : moment x initial en CP
!   sMyCPIni    : moment y initial en CP
    real(kind=8) :: sHslIni(2), sMxslIni, sMyslIni, sHCPini(2), sMxCPIni, sMyCPIni
!   la même chose mais à l'itération de calcul
    real(kind=8) :: sHslIte, sMxslIte, sMyslIte, sHCPIte, sMxCPIte, sMyCPIte
!   valeur avec le signe des moments et de l'effort horizontal pour le calcul des
!       différents critères
    real(kind=8) :: valH, valMx, valMy, valMxPabs, valMyPabs
!   paramètres important du macroéléemnt renommé
    real(kind=8) :: Lx, Ly, phi, Velas, error, Cinter
!   Force horizontales et excentrement pour les deux critères
!       avec les écrouissages cinématiques inclus
    real(kind=8) :: Hsl, exsl, eysl, HCP, exCP, eyCP
!
!   Fite    : à chaque itération de calcul la valeur de la force après ou avant correction (à voir)
!   dFcor   : la correction du retour élastique de la force
!   QCPite  : à chaque itération de calcul la valeur de l'écrouissage CP
!   dQCP    : l'incrément de l'écrouissage CP'
!   Qslite  : à chaque itération de calcul la valeur de l'écrouissage sl
!   dQsl    : l'incrément de l'écrouissage sl
!   dUpCP   : incrément de dépalcement plastiques en CP
!   dUpSL   : incrément de dépalcement plastiques en sliding
    real(kind=8) :: Fite(5), dFcor(5), QCPite(5), dQCP(5), Qslite(2), dQsl(2), dUpCP(5), dUpsl(5)
!   valeur de calcul intermédiaire
    real(kind=8) :: valint, det
!
!   HHCP        : matrice diagonale qui relie écrouissage CP et déplacement plastique
!   Hslid       : matrice diagonale qui relie écrouissage slid et déplacement plastique
!   KtangPetit  : matrice tangeante élastique réduite à 5 termes
    real(kind=8) :: HHCP(5), Hslid(2), KtangPetit(5)
!   matrices globales des dérivées du critères
    real(kind=8) :: dfdF(8, 5), dfdQCP(4, 5), dfdQsl(4, 2)
    real(kind=8) :: MatAinvTot(8, 8)
!   numérotation pour savoir quels critères sont concernés H et H- pour le transfert à
!       la matrice de raideur
    real(kind=8) :: etatCP, etatG, signeHHCP3(5), fssomme
!
    integer(kind=8), allocatable :: combss(:, :), VectVrai(:), vectPass(:)
!   MatAinverser : la matrice à inverser pour obtenir les avancements
!   fsInverse    : les deltas(finaux)
    real(kind=8), allocatable :: MatAinverser(:, :), fsInverse(:)
!
!   critère de convergence et listes des critères à regarder
    logical :: NoConv, criteres(8)
!   critere de convergence pour savoir si on a regardé les bons critères
    logical :: NoConvCrit
!   valeur de calcul intermédiaire
    logical :: LogiInt
!
! --------------------------------------------------------------------------------------------------
! initialisation de paramètre de calcul pour simplicité des écritures
!
    ! longueur de la fondation selon x
    Lx = vpara(1)
    ! longueur de la fondation selon y
    Ly = vpara(2)
    ! angle phi en rad (donné en degrés)
    phi = vpara(3)/45.0*Datan(1.D0)
    ! Capacité portant où on estime que la fondation est dans le domaine linéaire
    Velas = vpara(7)
    Cinter = vpara(4)
    ! erreur sur le critère autorisée ramenée à Vélastique
    error = Velas*errmax
    !
    ! initilisation des paramètres pour débuggage
    etatG = 0.0
    etatCP = 0.0
    allocate (fsInverse(1))
    fsInverse(:) = 0.0
    allocate (combss(1, 1))
    combss(:, :) = 0
    fsite(:) = 0.d0
    !
    ! calcul de fCP et fslid
    ! force horizontale pour le glissement
    Hsl = ((tirela(1)-vloc(11))**2+(tirela(2)-vloc(12))**2)**0.5
    ! la gestion de Fz (tirela(3)) en traction -Fz < Velas/100 est gérée dans une autre routine
    ! excentrement x pour le glissement
    exsl = abs(tirela(5))/(-tirela(3))
    ! excentrement y pour le glissement
    eysl = abs(tirela(4))/(-tirela(3))
!   calcul du critère de glissement
    fslid = Hsl+tirela(3)*tan(phi)-Cinter*Lx*Ly*(1.0-2.0*exsl/Lx)*(1.0-2.0*eysl/Ly)
    ! force horizontale pour le CP
    HCP = ((tirela(1)-vloc(13))**2+(tirela(2)-vloc(14))**2)**0.5
!   la gestion de Fz (tirela(3)) en traction -Fz < Velas/100 est gérée dans une autre routine
    ! excentrement x pour le CP
    exCP = abs(tirela(5)-vloc(17))/(-tirela(3))
    ! excentrement y pour le CP
    eyCP = abs(tirela(4)-vloc(16))/(-tirela(3))
!   calcul du critère de CP
    fCP = -tirela(3)-(Velas+vloc(15))*((1.0+HCP/tirela(3))**3)*(1.0-2.0*exCP/Lx)*(1.0-2.0*eyCP/Ly)
    if ((fCP .LT. error) .AND. (fslid .LT. error)) then
        ! on regarde si on est dans le domaine élastique et on gère la suite facilement en fait
        ! il n'y a rien à faire excepté la gestion de dautre
        vloc(18) = 1.0
        ! la variable interne montre que l'on est en élasticité'
    else
        ! on doit traiter la plasticité
        ! enregistrement des paramètres de signes initiaux normés à un
        valint = ((tirela(1)-vloc(11))**2+(tirela(2)-vloc(12))**2)**0.5
        if (valint .LE. r8prem()) then
            sHslIni(1) = 0.0
            sHslIni(2) = 0.0
        else
            sHslIni(1) = (tirela(1)-vloc(11))/valint
            sHslIni(2) = (tirela(2)-vloc(12))/valint
        end if
        sMxslIni = signe(tirela(4), r8prem()*Ly)
        sMyslIni = signe(tirela(5), r8prem()*Lx)
        valint = ((tirela(1)-vloc(13))**2+(tirela(2)-vloc(14))**2)**0.5
        if (valint .LE. r8prem()) then
            sHCPini(1) = 0.0
            sHCPini(2) = 0.0
        else
            sHCPini(1) = (tirela(1)-vloc(13))/valint
            sHCPini(2) = (tirela(2)-vloc(14))/valint
        end if
        sMxCPIni = signe(tirela(4)-vloc(16), r8prem()*Ly)
        sMyCPIni = signe(tirela(5)-vloc(17), r8prem()*Lx)
        ! Calcul préalable des paramètres hors boucle
        ! Calcul de HHCP telle que dqCP=HHCP*dupCP, HHCP(3) manquant
        ! car dépend si le déplacement plastique CP est vers le bas ou le haut
        HHCP(1) = vpara(10)*exp(-1.0*vpara(11)*vloc(21))
        HHCP(2) = vpara(12)*exp(-1.0*vpara(13)*vloc(21))
        HHCP(4) = vpara(14)*exp(-1.0*vpara(15)*vloc(21))
        HHCP(5) = vpara(16)*exp(-1.0*vpara(17)*vloc(21))
        ! Calcul de Hslid telle que dqs=Hslid*dupslid
        Hslid(1) = vpara(5)*exp(-1.0*vpara(6)*vloc(19))
        Hslid(2) = vpara(5)*exp(-1.0*vpara(6)*vloc(20))
        if (vpara(9) .GT. r8prem()) then
            ! traitement du cas DEPL_REFE=0 pas d'écrouissage cinématique' pour éviter
            ! une division par zéro
            HHCP(3) = vpara(8)*vpara(9)/(vpara(9)+vloc(21))**2
        else
            HHCP(3) = 0.0
        end if
        signeHHCP3 = 1.0
        ! sera utilisé pour calculer le signe de HHCP(3)
        ! lancement de la boucle de retour sur le critère
        ! initialisation des paramètres initiaux de la boucle
        ! initialisation du nombre d'itération
        Ite = 0
        ! logical de sortie de boucle initialisé
        NoConv = .TRUE.
        do ii = 1, 5
            ! Initialisation de la force qui va être corrigée au tir élastique
            Fite(ii) = tirela(ii)
            ! Initialisation de la correction à 0
            dFcor(ii) = 0.0
            ! Initialisation de l'écrouissage en CP qui va être modifié à l'écrousissage initial
            QCPite(ii) = vloc(12+ii)
            ! initilisation de l'incrément d'écrouissage en CP à 0
            dQCP(ii) = 0.0
            ! stockage des rigidité tangeante élastique (avec décollement) dans un vecteur
            !   à 5 composantes pour des calculs simplifiés
            KtangPetit(ii) = raidTang(ii)
        end do
        do ii = 1, 2
            ! Initialisation de l'écrouissage en glissement qui va être modifié
            !   à l'écrousissage initial
            Qslite(ii) = vloc(10+ii)
            ! initilisation de l'incrément d'écrouissage en glissement à 0
            dQsl(ii) = 0.0
        end do
        do while ((Noconv) .and. (Ite .LE. nbdecp))
            ! limitation de la boucle à nombre d'itérations
            ! vérification du changement de signe ou pas 0 est considéré même signe à error près
            valint = (Fite(1)-Qslite(1))*sHslIni(1)+(Fite(2)-Qslite(2))*sHslIni(2)
            sHslIte = signeExcl(valint, r8prem())
            sMxslIte = signeExcl(sMxslIni*Fite(4), r8prem()*Ly)
            sMyslIte = signeExcl(sMyslIni*Fite(5), r8prem()*Lx)
            valint = (Fite(1)-QCPite(1))*sHCPIni(1)+(Fite(2)-QCPite(2))*sHCPIni(2)
            sHCPIte = signeExcl(valint, r8prem())
            sMxCPIte = signeExcl(sMxCPIni*(Fite(4)-QCPIte(4)), r8prem()*Ly)
            sMyCPIte = signeExcl(sMyCPIni*(Fite(5)-QCPIte(5)), r8prem()*Lx)
            ! calcul des différents df/dF et du terme résultant en
            ! ici en glissement
            do ii = 1, 4
                ! création des Variables intermédiaires de signe des efforts et moment
                if (ii .EQ. 2) then
                    ! le critere et sa dérivée en ii=2 : glissement avec H de signe opposé
                    valH = -sHslIte
                    ! valH signe à mettre devant Fh pour le calcul
                else
                    valH = sHslIte
                end if
                if (ii .EQ. 3) then
                    ! le critere et sa dérivée en ii=3 : glissement avec Mx de signe opposé
                    valMx = -sMxslIte
                    ! valMx signe à mettre devant Mx pour le calcul
                    valMxPabs = -sMxslIni
                    ! valMxPabs signe de Mx pour le calcul
                else
                    valMx = sMxslIte
                    valMxPabs = sMxslIni
                end if
                if (ii .EQ. 4) then
                    ! le critere et sa dérivée en ii=4 : glissement avec My de signe opposé
                    valMy = -sMyslIte
                    ! valMy signe à mettre devant My pour le calcul
                    valMyPabs = -sMyslIni
                    ! valMyPabs signe de My pour le calcul
                else
                    valMy = sMyslIte
                    valMyPabs = sMyslIni
                end if
                ! calcul de la force horizontale
                Hsl = ((Fite(1)-Qslite(1))**2+(Fite(2)-Qslite(2))**2)**0.5
                ! traitement de dfdF(ii,1) et dfdF(ii,2) si la force horizontale est nulle
                if (Hsl .LE. r8prem()) then
                    dfdF(ii, 1) = sHslIni(1)*valH
                    dfdF(ii, 2) = sHslIni(2)*valH
                else
                    dfdF(ii, 1) = (Fite(1)-Qslite(1))/Hsl*valH
                    dfdF(ii, 2) = (Fite(2)-Qslite(2))/Hsl*valH
                end if
                ! calcul de dfdF(ii,3) dfdF(ii,4) dfdF(ii,5), la gestion d'un Fz
                !   en traction est déjà gérée'
                if ((abs(Fite(5)) .LT. (-Lx*Fite(3)/2.0)) &
                    .AND. (abs(Fite(4)) .LT. (-Ly*Fite(3))/2.0)) then
                    dfdF(ii, 3) = tan(phi)+2.0*Cinter*Lx*Ly*valMy* &
                                  abs(Fite(5))/(Lx*Fite(3)**2)* &
                                  (1.0+2.0*valMx*abs(Fite(4))/(Ly*Fite(3)))
                    dfdF(ii, 3) = dfdF(ii, 3)+2.0*Cinter*Lx*Ly*valMx* &
                                  abs(Fite(4))/(Ly*Fite(3)**2)* &
                                  (1.0+2.0*valMy*abs(Fite(5))/(Lx*Fite(3)))
                    dfdF(ii, 4) = -2.0*valMxPabs*Cinter*Lx*Ly* &
                                  (1.0+2.0*valMy*abs(Fite(5))/(Lx*Fite(3)))/(Ly*Fite(3))
                    dfdF(ii, 5) = -2.0*valMyPabs*Cinter*Lx*Ly* &
                                  (1.0+2.0*valMx*abs(Fite(4))/(Ly*Fite(3)))/(Lx*Fite(3))
                    fsite(ii) = valH*Hsl+Fite(3)*tan(phi)-Cinter*Lx*Ly* &
                                (1.0+2.0*valMy*abs(Fite(5))/(Lx*Fite(3)))* &
                                (1.0+2.0*valMx*abs(Fite(4))/(Ly*Fite(3)))
                else
                    ! si les moments sont trop forts, on n'a plus de cohésion
                    dfdF(ii, 3) = tan(phi)
                    dfdF(ii, 4) = 0.0
                    dfdF(ii, 5) = 0.0
                end if
                ! calcul des dérivées du critère en glissement par rapport aux écrouissages
                !   qui ne sont que l'opposé des dérivées par rapport aux force
                !   (car écrouissage cinématique)
                dfdQsl(ii, 1) = -dfdF(ii, 1)
                dfdQsl(ii, 2) = -dfdF(ii, 2)
                ! on introduit dans le calcul du critère la correction en fonction de
                !   l'itération pour calcul avant l'inversion'
                fsini(ii) = fsite(ii)+(vloc(11)-Qslite(1))*dfdQsl(ii, 1)+ &
                            (vloc(12)-Qslite(2))*dfdQsl(ii, 2)
                do jj = 1, 5
                    fsini(ii) = fsini(ii)+(tirela(jj)-Fite(jj))*dfdF(ii, jj)
                end do
            end do
            ! même pratique mais pour le critère CP
            do ii = 1, 4
                if (ii .EQ. 2) then
                    valH = -sHCPIte
                else
                    valH = sHCPIte
                end if
                if (ii .EQ. 3) then
                    valMx = -sMxCPIte
                    valMxPabs = -sMxCPIni
                else
                    valMx = sMxCPIte
                    valMxPabs = sMxCPIni
                end if
                if (ii .EQ. 4) then
                    valMy = -sMyCPIte
                    valMyPabs = -sMyCPIni
                else
                    valMy = sMyCPIte
                    valMyPabs = sMyCPIni
                end if
                HCP = ((Fite(1)-QCPite(1))**2+(Fite(2)-QCPite(2))**2)**0.5
                valint = -3.0*(Velas+QCPite(3))* &
                         (1.0+2.0*valMy*abs(Fite(5)-QCPite(5))/(Lx*Fite(3)))* &
                         (1.0+2.0*valMx*abs(Fite(4)-QCPite(4))/(Ly*Fite(3)))* &
                         (1.0+valH*HCP/Fite(3))**2
                if (HCP .LT. r8prem()) then
                    dfdF(ii+4, 1) = sHCPIni(1)*valH*valint/Fite(3)
                    dfdF(ii+4, 2) = sHCPIni(2)*valH*valint/Fite(3)
                else
                    dfdF(ii+4, 1) = (Fite(1)-QCPIte(1))*valH*valint/HCP/Fite(3)
                    dfdF(ii+4, 2) = (Fite(2)-QCPIte(2))*valH*valint/HCP/Fite(3)
                end if
                dfdF(ii+4, 3) = -1.0+(Velas+QCPIte(3))*2.0*valMx* &
                                abs(Fite(4)-QCPIte(4))/(Ly*Fite(3)**2)* &
                                (1.0+valH*HCP/Fite(3))**3* &
                                (1.0+2.0*valMy*abs(Fite(5)-QCPIte(5))/(Lx*Fite(3)))
                dfdF(ii+4, 3) = dfdF(ii+4, 3)+(Velas+QCPIte(3))*2.0*valMy* &
                                abs(Fite(5)-QCPIte(5))/(Lx*Fite(3)**2)*(1.0+valH*HCP/Fite(3))**3* &
                                (1.0+2.0*valMx*abs(Fite(4)-QCPIte(4))/(Ly*Fite(3)))
                dfdF(ii+4, 3) = dfdF(ii+4, 3)+(Velas+QCPIte(3))*3.0*valH*HCP/(Fite(3)**2)* &
                                (1.0+valH*HCP/Fite(3))**2* &
                                (1.0+2.0*valMx*abs(Fite(4)-QCPIte(4))/(Ly*Fite(3)))* &
                                (1.0+2.0*valMy*abs(Fite(5)-QCPIte(5))/(Lx*Fite(3)))
                dfdF(ii+4, 4) = -2.0*(Velas+QCPIte(3))/(Ly*Fite(3))*(1.0+valH*HCP/Fite(3))**3* &
                                (1.0+2.0*valMy*abs(Fite(5)-QCPIte(5))/(Lx*Fite(3)))*valMxPabs
                dfdF(ii+4, 5) = -2.0*(Velas+QCPIte(3))/(Lx*Fite(3))*(1.0+valH*HCP/Fite(3))**3* &
                                (1.0+2.0*valMx*abs(Fite(4)-QCPIte(4))/(Ly*Fite(3)))*valMyPabs
                ! la gestion d'un Fz en traction est déjà gérée'
                dfdQCP(ii, :) = -dfdF(ii+4, :)
                dfdQCP(ii, 3) = -(1.0+valH*HCP/Fite(3))**3* &
                                (1.0+2.0*valMx*abs(Fite(4)-QCPIte(4))/(Ly*Fite(3)))* &
                                (1.0+2.0*valMy*abs(Fite(5)-QCPIte(5))/(Lx*Fite(3)))
                fsite(ii+4) = -Fite(3)-(Velas+QCPIte(3))*(1.0+valH*HCP/Fite(3))**3* &
                              (1.0+2.0*valMx*abs(Fite(4)-QCPIte(4))/(Ly*Fite(3)))* &
                              (1.0+2.0*valMy*abs(Fite(5)-QCPIte(5))/(Lx*Fite(3)))
                fsini(ii+4) = fsite(ii+4)
                do jj = 1, 5
                    fsini(ii+4) = fsini(ii+4)+(tirela(jj)-Fite(jj))* &
                                  dfdF(ii+4, jj)+(vloc(12+jj)-QCPite(jj))*dfdQCP(ii, jj)
                end do
            end do
            ! les matrices globales sont déjà calculées, maintenant on crée le vecteur
            !       de booléen pour savoir quel critère calculer
            ! si criteres en glissement dépassé par rapport à r_prem ici
            if (fslid .GE. error) then
                ! glissement de base toujours vérifié
                criteres(1) = .TRUE.
                ! si la cohésion est non nulle alors pas de double critère possible
                if (abs(Cinter) .LT. errmax) then
                    do ii = 2, 4
                        LogiInt = .True.
                        do jj = 1, 5
                            LogiInt = LogiInt .and. (abs(dfdF(1, jj)-dfdF(ii, jj)) .LT. errmax)
                        end do
                        ! on regarde si il y a une différence entre les dérivées du critère
                        !   en glissement, si non, on ne considère pas le critère (glissement
                        !   avec H, Mx ou My de signe opposé)
                        criteres(ii) = .NOT. (LogiInt)
                    end do
                else
                    criteres(2:4) = .FALSE.
                end if
            else
                criteres(:4) = .FALSE.
            end if
            ! on réalise la même chose sur le critère en CP sauf qu'il n'y a pas
            !   de contraintes avec la cohésion'
            if (fCP .GE. error) then
                criteres(5) = .TRUE.
                do ii = 6, 8
                    LogiInt = .True.
                    do jj = 1, 5
                        LogiInt = LogiInt .and. (abs(dfdF(5, jj)-dfdF(ii, jj)) .LT. errmax)
                    end do
                    criteres(ii) = .NOT. (LogiInt)
                end do
            else
                criteres(5:) = .FALSE.
            end if
            ! On lance la boucle de convergence sur l'ensemble des critères
            IteCrit = 0
            NoConvCrit = .TRUE.
            ! calcul de la matrice avec tous les éléments de calcul de la matrice à inverser
            ! (dfdF)*K*(dfDF)-dfdQ*I*dfDF
            do ii = 1, 8
                do jj = 1, 8
                    MatAinvTot(ii, jj) = sum(dfdF(ii, :)*KtangPetit(:)*dfDF(jj, :))
                end do
            end do
            do ii = 1, 4
                do jj = 1, 4
                    ! on regarde si il faut rajouter le terme d'écrouissage sliding
                    MatAinvTot(ii, jj) = MatAinvTot(ii, jj)-dfdQsl(ii, 1)*Hslid(1)*dfdF(jj, 1)- &
                                         dfdQsl(ii, 2)*Hslid(2)*dfdF(jj, 2)
                end do
            end do
            do ii = 5, 8
                do jj = 5, 8
                    ! on regarde si il faut rajouter le terme d'écrousissage CP
                    if (dfdF(ii, 3) .GT. r8prem()) then
                        signeHHCP3(3) = +1.0
                    else
                        signeHHCP3(3) = -1.0
                    end if
                    MatAinvTot(ii, jj) = MatAinvTot(ii, jj)- &
                                         sum(dfdQCP(ii-4, :)*HHCP(:)*signeHHCP3(:)*dfdF(jj, :))
                end do
            end do
            ! on teste toutes les combinaisons et on enlève celle la plus haute
            NoConvCrit = .TRUE.
            if (criteres(1) .and. criteres(5)) then
                nbDDL = 0
                do jj = 1, 8
                    if (criteres(jj)) then
                        nbDDL = nbDDL+1
                    end if
                end do
                comppt = nbDDL-2
                if (allocated(VectVrai)) then
                    deallocate (VectVrai)
                end if
                allocate (vectVrai(comppt))
                comppt = 0
                do jj = 2, 4
                    if (criteres(jj)) then
                        comppt = comppt+1
                        vectVrai(comppt) = jj
                    end if
                end do
                do jj = 6, 8
                    if (criteres(jj)) then
                        comppt = comppt+1
                        vectVrai(comppt) = jj
                    end if
                end do
                NoConvCrit = .true.
                do while (NoConvCrit .and. (comppt >= 0))
                    if (allocated(MatAinverser)) then
                        deallocate (MatAinverser)
                    end if
                    if (allocated(fsInverse)) then
                        deallocate (fsInverse)
                    end if
                    if (allocated(vectPass)) then
                        deallocate (vectPass)
                    end if
                    if (allocated(combss)) then
                        deallocate (combss)
                    end if
                    allocate (MatAinverser(comppt+2, comppt+2), fsInverse(comppt+2), &
                              vectPass(comppt+2))
                    ! on parcours toutes les combinaisons
                    combss = combinaisons(nbDDL-2, comppt)
                    nbcomb = ubound(combss, 1)
                    comptcombi = 0
                    ! on affecte le bon côté ! cas final où on a fait le double critère avec
                    !   aucun autres critères inversé'
                    if ((sHslIte .lt. -r8prem()) .and. (comppt .EQ. 0)) then
                        ! on vérifie que le critère reste bien orienté : exemple on peut de H- à H+
                        !   en plein milieu de l'itération il convient alors de changer le
                        !   critère vérifié sur le double critère'
                        vectPass(1) = 2
                    else
                        vectPass(1) = 1
                    end if
                    ! cas final où on a fait le double critère avec aucun autres critères inversé'
                    if ((sHCPIte .lt. -r8prem()) .and. (comppt .EQ. 0)) then
                        ! on vérifie que le critère reste bien orienté : exemple on peut de H- à H+
                        !   en plein milieu de l'itération il convient alors de changer le
                        !   critère vérifié sur le double critère'
                        vectPass(2) = 6
                    else
                        vectPass(2) = 5
                    end if
                    do while (NoConvCrit .and. (comptcombi < nbcomb))
                        comptcombi = comptcombi+1
                        do ii = 1, comppt
                            vectPass(2+ii) = vectVrai(combss(comptcombi, ii))
                        end do
                        fssomme = 0
                        do ii = 1, 2+comppt
                            fsInverse(ii) = fsini(vectPass(ii))
                            fssomme = fssomme+abs(fsInverse(ii))
                            do jj = 1, 2+comppt
                                MatAinverser(ii, jj) = MatAinvTot(vectPass(ii), vectPass(jj))
                            end do
                        end do
                        fssomme = fssomme/(2+comppt)
                        ! on résout le système les écoulement locaux sont stocké dans fsInverse
                        call mgauss('NCVP', MatAinverser, fsInverse, 2+comppt, 2+comppt, &
                                    1, det, iret)
                        ! si pas de problème d'inversion
                        if (iret .EQ. 0) then
                            ! On regarde si tous les écoulements sont positifs, si oui, il
                            !   faut garder l'ensemble de ces DDLs
                            !   + regarder le cas des inversions iréalistes fs trop gros
                            NoConvCrit = .false.
                            do ii = 1, 2+comppt
                                NoConvCrit = NoConvCrit .or. &
                                             (fsInverse(ii) .LE. (-r8prem())) .or. &
                                             (abs(fsInverse(ii)) .GT. fssomme)
                            end do
                        end if
                    end do
                    comppt = comppt-1
                end do
                if (NoConvCrit) then
                    ! on est en critere de CP seul
                    if (fsInverse(1) .LE. (-r8prem())) then
                        ! on désactive le critère de glissement
                        criteres(1:4) = .False.
                    else if (fsInverse(2) .LE. (-r8prem())) then
                        ! on désactive le critère de CP
                        criteres(5:8) = .False.
                    end if
                else
                    ! sinon on remet a jour les critères convergés
                    criteres = .False.
                    do ii = 1, 2+comppt+1
                        criteres(vectPass(ii)) = .True.
                    end do
                end if
            end if
            ! fin de la boucle du double criteres
            ! boucle dans le cas critère de glissement seul
            if (NoConvCrit .and. criteres(1) .and. (.not. (criteres(5)))) then
                nbDDL = 0
                do jj = 1, 4
                    if (criteres(jj)) then
                        nbDDL = nbDDL+1
                    end if
                end do
                comppt = nbDDL-1
                if (allocated(VectVrai)) then
                    deallocate (VectVrai)
                end if
                allocate (vectVrai(comppt))
                comppt = 0
                do jj = 2, 4
                    if (criteres(jj)) then
                        comppt = comppt+1
                        vectVrai(comppt) = jj
                    end if
                end do
                !
                NoConvCrit = .true.
                do while (NoConvCrit .and. (comppt >= 0))
                    if (allocated(MatAinverser)) then
                        deallocate (MatAinverser)
                    end if
                    if (allocated(fsInverse)) then
                        deallocate (fsInverse)
                    end if
                    if (allocated(vectPass)) then
                        deallocate (vectPass)
                    end if
                    if (allocated(combss)) then
                        deallocate (combss)
                    end if
                    allocate (MatAinverser(comppt+1, comppt+1), fsInverse(comppt+1), &
                              vectPass(comppt+1))
                    ! on parcours toutes les combinaisons
                    combss = combinaisons(nbDDL-1, comppt)
                    nbcomb = ubound(combss, 1)
                    comptcombi = 0
                    vectPass(1) = 1
                    do while (NoConvCrit .and. (comptcombi < nbcomb))
                        comptcombi = comptcombi+1
                        do ii = 1, comppt
                            vectPass(1+ii) = vectVrai(combss(comptcombi, ii))
                        end do
                        do ii = 1, 1+comppt
                            fsInverse(ii) = fsini(vectPass(ii))
                            do jj = 1, 1+comppt
                                MatAinverser(ii, jj) = MatAinvTot(vectPass(ii), vectPass(jj))
                            end do
                        end do
                        ! on résout le système les écoulement locaux sont stocké dans fsInverse
                        call mgauss('NCVP', MatAinverser, fsInverse, 1+comppt, 1+comppt, &
                                    1, det, iret)
                        if (iret .EQ. 0) then
                            ! si pas de problème d'inversion
                            ! On regarde si tous les écoulements sont positifs, si oui,
                            !   il faut garder l'ensemble de ces DDLs'
                            NoConvCrit = .false.
                            do ii = 1, 1+comppt
                                NoConvCrit = NoConvCrit .or. (fsInverse(ii) .LE. (-r8prem()))
                            end do
                        end if
                    end do
                    comppt = comppt-1
                end do
                ! cas final sans convergence implique problème
                if (NoConvCrit) then
                    iret = 1
                    goto 999
                else
                    ! sinon on remet a jour les critères convergés
                    criteres = .False.
                    do ii = 1, 1+comppt+1
                        criteres(vectPass(ii)) = .True.
                    end do
                end if
            end if
            ! fin de la boucle du critère glissement
            ! boucle dans le cas critère de CP seul
            if (NoConvCrit .and. criteres(5) .and. (.not. (criteres(1)))) then
                nbDDL = 0
                do jj = 1, 4
                    if (criteres(jj+4)) then
                        nbDDL = nbDDL+1
                    end if
                end do
                comppt = nbDDL-1
                if (allocated(VectVrai)) then
                    deallocate (VectVrai)
                end if
                allocate (vectVrai(comppt))
                comppt = 0
                do jj = 2, 4
                    if (criteres(jj+4)) then
                        comppt = comppt+1
                        vectVrai(comppt) = jj+4
                    end if
                end do
                NoConvCrit = .true.
                do while (NoConvCrit .and. (comppt >= 0))
                    if (allocated(MatAinverser)) then
                        deallocate (MatAinverser)
                    end if
                    if (allocated(fsInverse)) then
                        deallocate (fsInverse)
                    end if
                    if (allocated(vectPass)) then
                        deallocate (vectPass)
                    end if
                    if (allocated(combss)) then
                        deallocate (combss)
                    end if
                    allocate (MatAinverser(comppt+1, comppt+1), &
                              fsInverse(comppt+1), vectPass(comppt+1))
                    ! on parcours toutes les combinaisons
                    combss = combinaisons(nbDDL-1, comppt)
                    nbcomb = ubound(combss, 1)
                    comptcombi = 0
                    vectPass(1) = 5
                    do while (NoConvCrit .and. (comptcombi < nbcomb))
                        comptcombi = comptcombi+1
                        do ii = 1, comppt
                            vectPass(1+ii) = vectVrai(combss(comptcombi, ii))
                        end do
                        do ii = 1, 1+comppt
                            fsInverse(ii) = fsini(vectPass(ii))
                            do jj = 1, 1+comppt
                                MatAinverser(ii, jj) = MatAinvTot(vectPass(ii), vectPass(jj))
                            end do
                        end do
                        ! on résout le système les écoulement locaux sont stocké dans fsInverse
                        call mgauss('NCVP', MatAinverser, fsInverse, 1+comppt, 1+comppt, &
                                    1, det, iret)
                        ! si pas de problème d'inversion
                        if (iret .EQ. 0) then
                            ! On regarde si tous les écoulements sont positifs, si oui,
                            ! il faut garder l'ensemble de ces DDLs'
                            NoConvCrit = .false.
                            do ii = 1, 1+comppt
                                NoConvCrit = NoConvCrit .or. (fsInverse(ii) .LE. (-r8prem()))
                            end do
                        end if
                    end do
                    comppt = comppt-1
                end do
                ! cas final sans convergence implique problème
                if (NoConvCrit) then
                    iret = 1
                    goto 999
                else
                    ! sinon on remet a jour les critères convergés
                    criteres = .False.
                    do ii = 1, 1+comppt+1
                        criteres(vectPass(ii)) = .True.
                    end do
                end if
            end if
            ! fin de la boucle du critère CP
            !
            fsret = 0
            do ii = 1, ubound(fsInverse, 1)
                fsret(vectPass(ii)) = fsInverse(ii)
            end do
            ! récupération des écoulements et cacluls des termes d'itérations
            ! les écoulements sont stockés dans fs, on calcule Fcor, dQcp ite et dQs
            dFcor = KtangPetit*matmul(transpose(dfdF), fsret)
            !
            dUpCP = matmul(transpose(dfdF(5:, :)), fsret(5:))
            dUpsl = matmul(transpose(dfdF(:4, :)), fsret(:4))
            dQCP = HHCP*dUpCP
            dQsl(1) = HSLID(1)*dUpsl(1)
            dQsl(2) = HSLID(2)*dUpsl(2)
            do ii = 1, 5
                Fite(ii) = tirela(ii)-dFcor(ii)
                QCPIte(ii) = vloc(12+ii)+dQCP(ii)
            end do
            ! on regarde l'évolution du dépalcement plastique et recalcul HHCP et R
            if (vpara(9) .GT. r8prem()) then
                QCPite(3) = vpara(8)*(vloc(21)+abs(dUpCP(3)))/(vpara(9)+vloc(21)+abs(dUpCP(3)))
            else
                QCPite(3) = 0.0
            end if
            if (abs(dUpCP(3)) .GT. r8prem()) then
                ! HHCP(3) est calculé comme un opérateur sécant et non tangeant à la
                !   formulation de l'écrouissage isotrope'
                HHCP(3) = abs((QCPite(3)-vloc(15))/dUpCP(3))
            end if
            Qslite(1) = vloc(11)+dQsl(1)
            Qslite(2) = vloc(12)+dQsl(2)
            ! on s'assure que l'on n'a pas explosé les critères de convergences dans la boucle
            if ((abs(Fite(4)-QCPite(4)) .GT. -Fite(3)*Ly/2.0) .OR. &
                (abs(Fite(5)-QCPite(5)) .GT. -Fite(3)*Lx/2.0) .OR. &
                (((Fite(1)-QCPIte(1))**2+(Fite(2)-QCPite(2))**2)**0.5 .GT. -Fite(3))) then
                iret = 1
                goto 999
            end if
            ! vérification de convergence
            nbDDLii = 0
            NoConv = .FALSE.
            ! recalcul de fCP et fsl
            fslbis = ((Fite(1)-Qslite(1))**2+(Fite(2)-Qslite(2))**2)**0.5+ &
                     Fite(3)*tan(phi)-Cinter*Lx*Ly*(1.0+2.0*abs(Fite(5))/(Lx*Fite(3)))* &
                     (1.0+2.0*abs(Fite(4))/(Ly*Fite(3)))
            fCPbis = -Fite(3)-(Velas+QCPIte(3))*(1.0+((Fite(1)-QCPIte(1))**2+ &
                                                      (Fite(2)-QCPIte(2))**2)**0.5/Fite(3))**3* &
                     (1.0+2.0*abs(Fite(4)-QCPIte(4))/(Ly*Fite(3)))* &
                     (1.0+2.0*abs(Fite(5)-QCPIte(5))/(Lx*Fite(3)))
            if (fslbis .GT. error) then
                NoConv = .TRUE.
                criteres(1) = .TRUE.
            end if
            if (fCPbis .GT. error) then
                NoConv = .TRUE.
                criteres(5) = .TRUE.
            end if
            ! recalcul de fsite obligatoire avec les mêmes signes qu'avant'
            valint = (Fite(1)-Qslite(1))*sHslIni(1)+(Fite(2)-Qslite(2))*sHslIni(2)
            sHslIte = signeExcl(valint, r8prem())
            sMxslIte = signeExcl(sMxslIni*Fite(4), r8prem()*Ly)
            sMyslIte = signeExcl(sMyslIni*Fite(5), r8prem()*Lx)
            valint = (Fite(1)-QCPite(1))*sHCPIni(1)+(Fite(2)-QCPite(2))*sHCPIni(2)
            sHCPIte = signeExcl(valint, r8prem())
            sMxCPIte = signeExcl(sMxCPIni*(Fite(4)-QCPIte(4)), r8prem()*Ly)
            sMyCPIte = signeExcl(sMyCPIni*(Fite(5)-QCPIte(5)), r8prem()*Lx)
            !
            fsitenow(1) = sHslIte*((Fite(1)-Qslite(1))**2+ &
                                   (Fite(2)-Qslite(2))**2)**0.5+Fite(3)*tan(phi)- &
                          Cinter*Lx*Ly*(1.0+2.0*sMyslIte*abs(Fite(5))/(Lx*Fite(3)))* &
                          (1.0+2.0*sMxslIte*abs(Fite(4))/(Ly*Fite(3)))
            fsitenow(2) = -sHslIte*((Fite(1)-Qslite(1))**2+(Fite(2)-Qslite(2))**2)**0.5+ &
                          Fite(3)*tan(phi)-Cinter*Lx*Ly* &
                          (1.0+2.0*sMyslIte*abs(Fite(5))/(Lx*Fite(3)))* &
                          (1.0+2.0*sMxslIte*abs(Fite(4))/(Ly*Fite(3)))
            fsitenow(3) = sHslIte*((Fite(1)-Qslite(1))**2+(Fite(2)-Qslite(2))**2)**0.5+ &
                          Fite(3)*tan(phi)-Cinter*Lx*Ly*(1.0+2.0*abs(Fite(5))/(Lx*Fite(3)))* &
                          (1.0-2.0*abs(Fite(4))/(Ly*Fite(3)))
            fsitenow(4) = sHslIte*((Fite(1)-Qslite(1))**2+(Fite(2)-Qslite(2))**2)**0.5+ &
                          Fite(3)*tan(phi)-Cinter*Lx*Ly*(1.0-2.0*abs(Fite(5))/(Lx*Fite(3)))* &
                          (1.0+2.0*sMxslIte*abs(Fite(4))/(Ly*Fite(3)))
            fsitenow(5) = -Fite(3)- &
                          (Velas+QCPIte(3))* &
                          (1.0+sHCPIte*((Fite(1)-QCPIte(1))**2+ &
                                        (Fite(2)-QCPIte(2))**2)**0.5/Fite(3))**3* &
                          (1.0+2.0*sMxCPIte*abs(Fite(4)-QCPIte(4))/(Ly*Fite(3)))* &
                          (1.0+2.0*sMyCPIte*abs(Fite(5)-QCPIte(5))/(Lx*Fite(3)))
            fsitenow(6) = -Fite(3)- &
                          (Velas+QCPIte(3))* &
                          (1.0-sHCPIte*((Fite(1)-QCPIte(1))**2+ &
                                        (Fite(2)-QCPIte(2))**2)**0.5/Fite(3))**3* &
                          (1.0+2.0*sMxCPIte*abs(Fite(4)-QCPIte(4))/(Ly*Fite(3)))* &
                          (1.0+2.0*sMyCPIte*abs(Fite(5)-QCPIte(5))/(Lx*Fite(3)))
            fsitenow(7) = -Fite(3)- &
                          (Velas+QCPIte(3))* &
                          (1.0+sHCPIte*((Fite(1)-QCPIte(1))**2+ &
                                        (Fite(2)-QCPIte(2))**2)**0.5/Fite(3))**3* &
                          (1.0-2.0*sMxCPIte*abs(Fite(4)-QCPIte(4))/(Ly*Fite(3)))* &
                          (1.0+2.0*sMyCPIte*abs(Fite(5)-QCPIte(5))/(Lx*Fite(3)))
            fsitenow(8) = -Fite(3)- &
                          (Velas+QCPIte(3))* &
                          (1.0+sHCPIte*((Fite(1)-QCPIte(1))**2+ &
                                        (Fite(2)-QCPIte(2))**2)**0.5/Fite(3))**3* &
                          (1.0+2.0*sMxCPIte*abs(Fite(4)-QCPIte(4))/(Ly*Fite(3)))* &
                          (1.0-2.0*sMyCPIte*abs(Fite(5)-QCPIte(5))/(Lx*Fite(3)))
            do ii = 1, 8
                if (criteres(ii)) then
                    ! on vérifie pour les critères dépassés par le tir élastique la bonne
                    !   convergence sur le critère
                    nbDDLii = nbDDLii+1
                    if (abs(fsitenow(ii)) .GT. error) then
                        NoConv = .TRUE.
                    end if
                end if
            end do
            ! après convergence, on regarde les autres critères
            if (.not. (NoConv)) then
                do ii = 1, 8
                    ! pour les autres critères, on vérifie juste qu'ils ne sont pas dépassés
                    !   suites aux convergences précédentes
                    if (.not. (criteres(ii))) then
                        if ((fsite(ii)) .GT. error) then
                            NoConv = .TRUE.
                            if ((ii .EQ. 1)) then
                                criteres(ii) = .TRUE.
                                fslid = fsite(ii)
                            end if
                            if ((ii .EQ. 5)) then
                                criteres(ii) = .TRUE.
                                fCP = fsite(ii)
                            end if
                        end if
                    end if
                end do
            end if
            !on sort de la boucle de convergence
            Ite = Ite+1
        end do
        ! vérification de convergence
        if (NoConv) then
            iret = 1
            goto 999
        end if
        ! on mets alors à jour les force et paramètres internes en sortie
        do ii = 1, 5
            tirela(ii) = Fite(ii)
            vloc(12+ii) = QCPIte(ii)
            vloc(ii) = vloc(ii)+dUpsl(ii)
            vloc(ii+5) = vloc(ii+5)+dUpCP(ii)
        end do
        vloc(11) = Qslite(1)
        vloc(12) = Qslite(2)
        vloc(19) = vloc(19)+abs(dUpsl(1))
        vloc(20) = vloc(20)+abs(dUpsl(2))
        vloc(21) = vloc(21)+abs(dUpCP(3))
        ! définition pour le calcul de matrice des critères dépassé dans vloc(18) selon une
        !       numérotation précise
        ! alors glissement, on initialise etatG à 0 au cas où
        etatG = 0.0
        if (criteres(1)) then
            if (criteres(2)) then
                if (criteres(3)) then
                    if (criteres(4)) then
                        ! H, Mx et My concomittant
                        etatG = 8.0
                    else
                        ! H et Mx concomittant
                        etatG = 5.0
                    end if
                else
                    if (criteres(4)) then
                        ! H et My concomittant
                        etatG = 6.0
                    else
                        ! H concomittant seul
                        etatG = 2.0
                    end if
                end if
            else
                if (criteres(3)) then
                    if (criteres(4)) then
                        ! Mx et My concomittant
                        etatG = 7.0
                    else
                        ! Mx concomittant seul
                        etatG = 3.0
                    end if
                else
                    if (criteres(4)) then
                        !My concomittant seul
                        etatG = 4.0
                    else
                        etatG = 1.0
                    end if
                end if
            end if
        end if
        ! alors CP, on initialise etatCP à 0 au cas où
        etatCP = 0.0
        if (criteres(5)) then
            if (criteres(6)) then
                if (criteres(7)) then
                    if (criteres(8)) then
                        ! H, Mx et My concomittant
                        etatCP = 8.0
                    else
                        ! H et Mx concomittant
                        etatCP = 5.0
                    end if
                else
                    if (criteres(8)) then
                        ! H et My concomittant
                        etatCP = 6.0
                    else
                        ! H concomittant seul
                        etatCP = 2.0
                    end if
                end if
            else
                if (criteres(7)) then
                    if (criteres(8)) then
                        ! Mx et My concomittant
                        etatCP = 7.0
                    else
                        ! Mx concomittant seul
                        etatCP = 3.0
                    end if
                else
                    if (criteres(8)) then
                        !My concomittant seul
                        etatCP = 4.0
                    else
                        etatCP = 1.0
                    end if
                end if
            end if
        end if
        if (criteres(1)) then
            !double critère
            if (criteres(5)) then
                vloc(18) = 100.0+10.0*etatCP+etatG
            else
                !critère glissement seul
                vloc(18) = 10.0+etatG
            end if
        else
            ! CP seul
            if (criteres(5)) then
                vloc(18) = 20.0+etatCP
            else
                ! élastique mais dans ce cas on y arrive pas là
                vloc(18) = 1.0
            end if
        end if
        if (criteres(1) .and. criteres(6) .and. (.not. criteres(5))) then
            vloc(18) = 111.0
        end if
        if (criteres(5) .and. criteres(2) .and. (.not. criteres(1))) then
            vloc(18) = 111.0
        end if
        if (criteres(6) .and. (.not. criteres(5)) .and. criteres(2) .and. (.not. criteres(1))) then
            vloc(18) = 111.0
        end if
    end if

999 continue
    !
    if (allocated(combss)) then
        deallocate (combss)
    end if
    if (allocated(VectVrai)) then
        deallocate (VectVrai)
    end if
    if (allocated(vectPass)) then
        deallocate (vectPass)
    end if
    if (allocated(MatAinverser)) then
        deallocate (MatAinverser)
    end if
    if (allocated(fsInverse)) then
        deallocate (fsInverse)
    end if
    !
contains

! fonction rapide de signe
!   attention ici est égal à 0 pour a=0 contrairement à sign de fortran avec erreur
    function signe(a, error)
        real(kind=8) :: a, error
        integer(kind=8) :: signe
        if (a .LT. -error) then
            signe = -1
        else if (a .GT. error) then
            signe = +1
        else
            signe = 0
        end if
    end function

! fonction rapide de signe
!   attention ici est égal à 1 pour a=0 contrairement à sign de fortran avec erreur
    function signeExcl(a, error)
        real(kind=8) :: a, error
        integer(kind=8) :: signeExcl
        if (a .LT. -error) then
            signeExcl = -1
        else if (a .GT. error) then
            signeExcl = +1
        else
            signeExcl = 1
        end if
    end function

! Fonction qui fait la liste des combinaisons
    function combinaisons(nbDDL, compt)
        integer(kind=8) :: nbDDL, compt, jj, nbComb, comptloc, kk, ll, mm, nn
        integer(kind=8), allocatable :: combinaisons(:, :)
        nbComb = 1
        do jj = 1, compt
            nbComb = nbComb*(nbDDL+1-jj)/jj
        end do
        allocate (combinaisons(nbComb, compt))
        !
        if (compt == 1) then
            do jj = 1, nbDDL
                combinaisons(jj, 1) = jj
            end do
        else if (compt == 2) then
            comptloc = 0
            do jj = 1, nbDDL
                do kk = jj+1, nbDDL
                    comptloc = comptloc+1
                    combinaisons(comptloc, 1) = jj
                    combinaisons(comptloc, 2) = kk
                end do
            end do
        else if (compt == 3) then
            comptloc = 0
            do jj = 1, nbDDL
                do kk = jj+1, nbDDL
                    do ll = kk+1, nbDDL
                        comptloc = comptloc+1
                        combinaisons(comptloc, 1) = jj
                        combinaisons(comptloc, 2) = kk
                        combinaisons(comptloc, 3) = ll
                    end do
                end do
            end do
        else if (compt == 4) then
            comptloc = 0
            do jj = 1, nbDDL
                do kk = jj+1, nbDDL
                    do ll = kk+1, nbDDL
                        do mm = ll+1, nbDDL
                            comptloc = comptloc+1
                            combinaisons(comptloc, 1) = jj
                            combinaisons(comptloc, 2) = kk
                            combinaisons(comptloc, 3) = ll
                            combinaisons(comptloc, 4) = mm
                        end do
                    end do
                end do
            end do
        else if (compt == 5) then
            comptloc = 0
            do jj = 1, nbDDL
                do kk = jj+1, nbDDL
                    do ll = kk+1, nbDDL
                        do mm = ll+1, nbDDL
                            do nn = mm+1, nbDDL
                                comptloc = comptloc+1
                                combinaisons(comptloc, 1) = jj
                                combinaisons(comptloc, 2) = kk
                                combinaisons(comptloc, 3) = ll
                                combinaisons(comptloc, 4) = mm
                                combinaisons(comptloc, 5) = nn
                            end do
                        end do
                    end do
                end do
            end do
        else if (compt == 6) then
            do jj = 1, 6
                combinaisons(1, jj) = jj
            end do
        else
            combinaisons(:, :) = 0
        end if
    end function

end subroutine
