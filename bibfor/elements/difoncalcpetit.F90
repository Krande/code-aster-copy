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
subroutine difoncalcpetit(tirela, raidTang, vloc, vpara, nbVloc, &
                          nbPara, iret, nbdecp, errmax, calculNormal, &
                          calculPetitH)
!
! --------------------------------------------------------------------------------------------------
!
!  IN
!     tirela   : Tir élastique en entrée
!     nbVloc   : nombre de paramètre de la loi locale
!     vloc     : Paramètre de la loi locale
!     raidTang : Raideur tangeante
!     nbPara   : nombre variables locales
!     valvarloc: Variables locales en entrée
!     nbdecp   : Nombres d'itérations dans la boucle de convergence
!     errmax   : Erreur relative autorisée (indexée sur Vélas)
!     calculPetitH : critère de dépassement du calcul
!
!  OUT
!     tirela   : Tir élastique corrigé en sortie
!     valvarloc: Variables locales modifiées en sortie
!     iret     : code de retour du matériau, on le met à 1 quand les convergences sont dépassées
!     calculNormal : vérification si le calcul à tourner normalement
!
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
!
! --------------------------------------------------------------------------------------------------
!
    implicit none
!
#include "asterc/r8prem.h"
#include "asterfort/mgauss.h"
!
    integer(kind=8) :: nbVloc, nbPara, iret, nbdecp
    real(kind=8) :: vpara(nbPara), tirela(6), raidTang(6), vloc(nbVloc), errmax
    logical :: calculNormal, calculPetitH
!   compteur du nombre d'itération
    integer(kind=8) :: Ite, IteCrit, ii, jj, compt(3)
!   nombre de critères à vérifier et compteurs
    integer(kind=8) :: nbDDL, nbDDLii, nbDDLjj
!
!   fcpH : valeur réelle du coefficient de sécurité à la capacité portante sur H seulement
!   fCPMx : idem pour Mx
!   fCPMy : idem pour My
!   fslid : idem pour le glissement
!   fsret : les écoulements plastiques correspodants
!   fsini : vecteur qui sauvegarde les critère de rupture normé
!       fs (1)  =fslid(H,Mx,My)
!       fs (2)  =fslid(-H,Mx,My)
!       fs (3)  =fslid(H,-Mx,My)
!       fs (4)  =fslid(H,Mx,-My)
!       fs (5)  =fCPH(H,Mx,My)
!       fs (6)  =fCPH(-H,Mx,My)
!       fs (7)  =fCPH(H,-Mx,My)
!       fs (8)  =fCPH(H,Mx,-My)
!       fs (9)  =fCPMx(H,Mx,My)
!       fs (10) =fCPMx(-H,Mx,My)
!       fs (11) =fCPMx(H,-Mx,My)
!       fs (12) =fCPMx(H,Mx,-My)
!       fs (13) =fCPMx(H,Mx,My)
!       fs (14) =fCPMx(-H,Mx,My)
!       fs (15) =fCPMx(H,-Mx,My)
!       fs (16) =fCPMx(H,Mx,-My)
    real(kind=8) :: fCPh, fCPMx, fCPMy, fslid, fsini(16), fsret(16), fsite(16)
!
!   sHslIni : force horizontale initiale en glissement
!   sMxslIni : moment x initial en glissement
!   sMyslIni : moment y initial en glissement
!   sHCPini : force horizontale initiale en CP
!   sMxCPIni : moment x initial en CP
!   sMyCPIni : moment y initial en CP
    real(kind=8) :: sHslIni(2), sMxslIni, sMyslIni, sHCPini(2), sMxCPIni, sMyCPIni
!   la même chose mais à l'itération de calcul
    real(kind=8) :: sHslIte, sMxslIte, sMyslIte, sHCPIte, sMxCPIte, sMyCPIte
!   valeur avec le signe des moments et de l'effort horizontal pour les différents critères
    real(kind=8) :: valH, valMx, valMy, valMxPabs, valMyPabs
!   paramètres important du macroélément renommé
    real(kind=8) :: Lx, Ly, phi, Velas, error, Cinter
!   Force horizontales et excentrement pour les deux critères
!       avec les écrouissages cinématiques inclus
    real(kind=8) :: Hsl, exsl, eysl, HCP
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
!   HHCP matrice diagonale qui relie écrouissage CP et déplacement plastique
!   Hslid matrice diagonale qui relie écrouissage slid et déplacement plastique
!   KtangPetit matrice tangeante élastique réduite à 5 termes
    real(kind=8) :: HHCP(5), Hslid(2), KtangPetit(5)
!   matrices globales des dérivées du critères
    real(kind=8) :: dfdF(16, 5), dfdQCP(12, 5), dfdQsl(4, 2), MatAinvTot(16, 16)
!   MatAinverser : la matrice à inverser pour obtenir les écoulements plastiques
!   fsInverse : après inversion du système contient les écoulements plastiques finaux
    real(kind=8), allocatable :: MatAinverser(:, :), fsInverse(:)
!   numérotation pour savoir quels critères sont concernés H et H- (idem avec les moments)
!   pour le transfert à la matrice de raideur
    real(kind=8) :: etatCPh, etatCPMx, etatCPMy, etatG, etatfactor
!   factor : 0.1 ou 0.01 , facteur de linéarisation déterminé dans difonabb
    real(kind=8) :: factor
!
!   calculS     : pour savoir si on calcule le critère en glissement
!   calculH     : pour savoir si on calcule le critère en CP linéarisé selon H
!   calculMx    : pour savoir si on calcule le critère en CP linéarisé selon Mx
!   calculMy    : pour savoir si on calcule le critère en CP linéarisé selon My
    logical :: calculS, calculH, calculMx, calculMy
!   valeur de calcul intermédiaire, paramètre regardant si on est en traction pure
    logical :: LogiInt, TracPure
!   critère de convergence et listes des critères à regarder
    logical :: NoConv, criteres(16)
!   critere de convergence pour savoir si on a regardé les bons critères
    logical :: NoConvCrit, NoConvCritInt
!
! -------------------------------------------------------------------------------------------------
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

    ! initialisation de factor
    if (calculPetitH) then
        factor = 0.01
    else
        factor = 0.10
    end if
    ! vérification que les moments et efforts sont bien suffisamment faibles
    !   par rapport au facteur utilisé
    if ((((tirela(1)-vloc(13))**2+(tirela(2)-vloc(14))**2)**0.5) .GE. &
        (factor*(Velas+vloc(15))-error)) then
        iret = 1
        goto 999
    end if
    if (abs(tirela(4)-vloc(16))/Ly .GE. (factor*(Velas+vloc(15))/2.0-error)) then
        iret = 1
        goto 999
    end if
    if (abs(tirela(5)-vloc(17))/Lx .GE. (factor*(Velas+vloc(15))/2.0-error)) then
        iret = 1
        goto 999
    end if
    ! calcul de fCPs et fslid
    ! force horizontale pour le glissement
    Hsl = ((tirela(1)-vloc(11))**2+(tirela(2)-vloc(12))**2)**0.5
    ! calcul de l'excentrement si Fz négatif'
    if (tirela(3) .LT. -error) then
        ! excentrement x pour le glissement
        exsl = min(Ly/2.0, abs(tirela(5))/(-tirela(3)))
        ! excentrement y pour le glissement
        eysl = min(abs(tirela(4))/(-tirela(3)), Lx/2.0)
    else
        ! gestion de l'excentrement si Fz positif alors l'excentrement est soit nul ou max'
        ! à voir si on doit comparer à error (dimensionné) ou à r8prem()
        if (abs(tirela(4)) .LT. error*Ly) then
            eysl = 0.0
        else
            eysl = Ly/2.0
        end if
        if (abs(tirela(5)) .LT. error*Lx) then
            exsl = 0.0
        else
            exsl = Lx/2.0
        end if
    end if
    ! calcul du critère de glissement
    fslid = Hsl+tirela(3)*tan(phi)-Cinter*Lx*Ly*(1-2*exsl/Ly)*(1-2*eysl/Lx)
    ! force horizontale pour le CP
    HCP = ((tirela(1)-vloc(13))**2+(tirela(2)-vloc(14))**2)**0.5
    ! calcul des critères linéarisés de CP
    fCPH = (1.0-(factor/((1.0-2.0*abs(tirela(4)-vloc(16))/ &
                          (factor*Ly*(Velas+vloc(15))))* &
                         (1.0-2.0*abs(tirela(5)-vloc(17))/ &
                          (factor*Lx*(Velas+vloc(15))))))**(1.0/3.0))*tirela(3)+HCP
    fCPMx = (1.0-factor/((1.0-HCP/(factor*(Velas+vloc(15))))**3* &
                         (1.0-2.0*abs(tirela(5)-vloc(17))/ &
                          (factor*Lx*(Velas+vloc(15))))))*Ly*tirela(3)+ &
            2.0*abs(tirela(4)-vloc(16))
    fCPMy = (1.0-factor/((1.0-HCP/(factor*(Velas+vloc(15))))**3* &
                         (1.0-2.0*abs(tirela(4)-vloc(16))/ &
                          (factor*Ly*(Velas+vloc(15))))))*Lx*tirela(3)+ &
            2.0*abs(tirela(5)-vloc(17))
    ! on regarde en plus si on est en traction pure
    TracPure = ((abs(((tirela(1)-vloc(13))**2+(tirela(2)-vloc(14))**2)**0.5) .LE. error) &
                .AND. (abs(tirela(4)-vloc(16)) .LE. error*Ly) &
                .AND. (abs(tirela(5)-vloc(17)) .LE. error*Lx))
    !
    if ((fCPH .LT. r8prem()) .AND. (fCPMx .LT. r8prem()) .AND. &
        (fCPMy .LT. r8prem()) .AND. (fslid .LT. r8prem())) then
        ! on regarde si on est dans le domaine élastique et on gère la suite facilement
        !   en fait il n'y a rien à faire excepté la gestion de dautre
        ! la variable interne montre que l'on est en élasticité
        vloc(18) = 1.0
        ! sinon on calcule en f=0
        calculNormal = .FALSE.
    else
        ! on doit traiter la plasticité
        ! on stocke dans des variables boolean quels criteres est dépassé
        calculS = (fslid .GT. r8prem())
        calculH = (fCPH .GE. r8prem()) .AND. ((abs(HCP) .GT. r8prem()) .OR. TracPure)
        calculMx = ((fCPMx .GE. r8prem()*Ly) .AND. (abs(tirela(4)-vloc(16)) .GT. r8prem()*Ly))
        calculMy = ((fCPMy .GE. r8prem()*Ly) .AND. (abs(tirela(5)-vloc(17)) .GT. r8prem()*Lx))
        ! enregistrement des paramètres de signes initiaux normés à un
        valint = ((tirela(1)-vloc(11))**2+(tirela(2)-vloc(12))**2)**0.5
        if (valint .LT. r8prem()) then
            sHslIni(1) = 0.0
            sHslIni(2) = 0.0
        else
            sHslIni(1) = (tirela(1)-vloc(11))/valint
            sHslIni(2) = (tirela(2)-vloc(12))/valint
        end if
        sMxslIni = signe(tirela(4), error*Ly)
        sMyslIni = signe(tirela(5), error*Lx)
        valint = ((tirela(1)-vloc(13))**2+(tirela(2)-vloc(14))**2)**0.5
        if (valint .LT. error) then
            sHCPini(1) = 0.0
            sHCPini(2) = 0.0
        else
            sHCPini(1) = (tirela(1)-vloc(13))/valint
            sHCPini(2) = (tirela(2)-vloc(14))/valint
        end if
        sMxCPIni = signe(tirela(4)-vloc(16), error*Ly)
        sMyCPIni = signe(tirela(5)-vloc(17), error*Lx)
        !   Calcul préalable des paramètres hors boucle
        !   Calcul de HHCP telle que dqCP=HHCP*dupCP, HHCP(3) manquant
        !   car dépend si le déplacement plastique CP est vers le bas ou le haut
        HHCP(1) = vpara(10)*exp(-1.0*vpara(11)*vloc(21))
        HHCP(2) = vpara(12)*exp(-1.0*vpara(13)*vloc(21))
        HHCP(4) = vpara(14)*exp(-1.0*vpara(15)*vloc(21))
        HHCP(5) = vpara(16)*exp(-1.0*vpara(17)*vloc(21))
        ! Calcul de Hslid telle que dqs=Hslid*dupslid
        Hslid(1) = vpara(5)*exp(-1.0*vpara(6)*vloc(19))
        Hslid(2) = vpara(5)*exp(-1.0*vpara(6)*vloc(20))
        ! Traitement particulier de HHCP(3) qui va dépendre du signe dfCP/dFz
        !   (en soit dans la calcul en petit critère le signe devrait toujorus être positif)
        ! traitement du cas DEPL_REFE=0 pas d'écrouissage cinématique' pour éviter
        !   une division par zéro
        if (vpara(9) .GT. r8prem()) then
            HHCP(3) = vpara(8)*vpara(9)/(vpara(9)+vloc(21))**2
        else
            HHCP(3) = 0.0
        end if
        ! lancement de la boucle de convergence sur le/les critère(s)
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
            ! initilisation de l'incrément d'écrouissage en CP à 0'
            dQCP(ii) = 0.0
            ! stockage des rigidité tangeante élastique (avec décollement) dans un vecteur
            !   à 5 composantes pour des calculs simplifiés
            KtangPetit(ii) = raidTang(ii)
        end do
        do ii = 1, 2
            ! Initialisation de l'écrouissage en glissement qui va être modifié
            !   à l'écrousissage initial
            Qslite(ii) = vloc(10+ii)
            ! initilisation de l'incrément d'écrouissage en glissement à 0'
            dQsl(ii) = 0.0
        end do
        ! limitation de la boucle à nombre d'itérations
        do while ((Noconv) .and. (Ite .LE. nbdecp))
            ! vérification du changement de signe ou pas au cours de l'itération' 0 est
            !   considéré même signe à error près
            valint = (Fite(1)-Qslite(1))*sHslIni(1)+(Fite(2)-Qslite(2))*sHslIni(2)
            sHslIte = signeExcl(valint, error)
            sMxslIte = signeExcl(sMxslIni*Fite(4), error*Ly)
            sMyslIte = signeExcl(sMyslIni*Fite(5), error*Lx)
            valint = (Fite(1)-QCPite(1))*sHCPIni(1)+(Fite(2)-QCPite(2))*sHCPIni(2)
            sHCPIte = signeExcl(valint, error)
            sMxCPIte = signeExcl(sMxCPIni*(Fite(4)-QCPIte(4)), error*Ly)
            sMyCPIte = signeExcl(sMyCPIni*(Fite(5)-QCPIte(5)), error*Lx)
            ! calcul des différents df/dF et du terme résultant en
            ! ici en glissement
            do ii = 1, 4
                if (calculS) then
                    ! création des Variables intermédiaires de signe des efforts et moment
                    ! le critere et sa dérivée en ii=2 correspond au glissement
                    !   avec H de signe opposé
                    if (ii .EQ. 2) then
                        ! valH signe à mettre devant Fh pour le calcul
                        valH = -sHslIte
                    else
                        valH = sHslIte
                    end if
                    ! le critere et sa dérivée en ii=3 correspond au glissement avec
                    !   Mx de signe opposé
                    if (ii .EQ. 3) then
                        ! valMx signe à mettre devant Mx pour le calcul
                        valMx = -sMxslIte
                        ! valMxPabs signe de Mx pour le calcul
                        valMxPabs = -sMxslIni
                    else
                        valMx = sMxslIte
                        valMxPabs = sMxslIni
                    end if
                    ! le critere et sa dérivée en ii=4 correspond au glissement avec
                    !   My de signe opposé
                    if (ii .EQ. 4) then
                        ! valMy signe à mettre devant My pour le calcul
                        valMy = -sMyslIte
                        ! valMyPabs signe de My pour le calcul
                        valMyPabs = -sMyslIni
                    else
                        valMy = sMyslIte
                        valMyPabs = sMyslIni
                    end if
                    ! calcul de la force horizontale
                    Hsl = ((Fite(1)-Qslite(1))**2+(Fite(2)-Qslite(2))**2)**0.5
                    ! traitement particulier de dfdF(ii,1) et dfdF(ii,2) si la force
                    !   horizontale est nulle
                    if (Hsl .LT. r8prem()) then
                        dfdF(ii, 1) = sHslIni(1)*valH
                        dfdF(ii, 2) = sHslIni(2)*valH
                    else
                        dfdF(ii, 1) = (Fite(1)-Qslite(1))/Hsl*valH
                        dfdF(ii, 2) = (Fite(2)-Qslite(2))/Hsl*valH
                    end if
                    ! calcul de dfdF(ii,3) dfdF(ii,4) dfdF(ii,5)
                    ! dans le cas d'unr force verticale en compression
                    if (Fite(3) .LT. -error) then
                        ! attention si les moments ne sont pas trops forts
                        if ((abs(Fite(5)) .LT. (-Lx*Fite(3)/2.0)) .AND. &
                            (abs(Fite(4)) .LT. (-Ly*Fite(3)/2.0))) then
                            dfdF(ii, 3) = tan(phi)+2.0*Cinter*Lx*Ly*valMy*abs(Fite(5))/ &
                                          (Lx*Fite(3)**2)*(1.0+2.0*valMx*abs(Fite(4))/(Ly*Fite(3)))
                            dfdF(ii, 3) = dfdF(ii, 3)+2.0*Cinter*Lx*Ly*valMx* &
                                          abs(Fite(4))/(Ly*Fite(3)**2)* &
                                          (1.0+2.0*valMy*abs(Fite(5))/(Lx*Fite(3)))
                            dfdF(ii, 4) = -2.0*valMxPabs*Cinter*Lx*Ly* &
                                          (1.0+2.0*valMy*abs(Fite(5))/(Lx*Fite(3)))/(Ly*Fite(3))
                            dfdF(ii, 5) = -2.0*valMyPabs*Cinter*Lx*Ly* &
                                          (1.0+2.0*valMx*abs(Fite(4))/(Ly*Fite(3)))/(Lx*Fite(3))
                            fsite(ii) = valH*Hsl+Fite(3)*tan(phi)- &
                                        Cinter*Lx*Ly*(1.0+2.0*valMy*abs(Fite(5))/(Lx*Fite(3)))* &
                                        (1.0+2.0*valMx*abs(Fite(4))/(Ly*Fite(3)))
                        else
                            ! si les moments sont trop forts, on annule la composante liée à
                            !   la cohésion ! sinon subdvision
                            iret = 1
                            goto 999
                        end if
                    else if ((abs(Fite(5)) .LT. error) .AND. (abs(Fite(4)) .LT. error)) then
                        ! si force verticale en traction
                        ! on regarde si les moments sont nuls ou pas pour prise en compte de la
                        !   cohséion dans le calcul du citère
                        dfdF(ii, 3) = tan(phi)
                        dfdF(ii, 4) = 0
                        dfdF(ii, 5) = 0
                        fsite(ii) = valH*Hsl+Fite(3)*tan(phi)-Cinter*Lx*Ly
                    else
                        ! sinon subdvision
                        dfdF(ii, 3) = tan(phi)
                        dfdF(ii, 4) = 0
                        dfdF(ii, 5) = 0
                        fsite(ii) = valH*Hsl+Fite(3)*tan(phi)
                    end if
                    ! calcul des dérivées du critère en glissement par rapport aux écrouissages
                    ! qui ne sont que l'opposé des dérivées par rapport aux forces
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
                else
                    ! si pas de critère de glissement à calculer, on mets tout à 0
                    dfdF(ii, :) = 0.0
                    fsite(ii) = 0.0
                    fsini(ii) = 0.0
                    dfdQsl(ii, :) = 0.0
                end if
            end do
            !
            ! calculs des dérivées et terme résultant pour le critère CP
            do ii = 1, 4
                ! création des Variables intermédiaires de signe des efforts et moment
                ! le critere et sa dérivée en ii=2 c'est le CP linéarisé avec H de signe opposé
                if (ii .EQ. 2) then
                    ! valH signe à mettre devant Fh pour le calcul
                    valH = -sHCPIte
                else
                    valH = sHCPIte
                end if
                ! le critere et sa dérivée en ii=3 c'est le CP linéarisé avec Mx de signe opposé
                if (ii .EQ. 3) then
                    ! valMx signe à mettre devant Mx pour le calcul
                    valMx = -sMxslIte
                    ! valMxPabs signe de Mx pour le calcul
                    valMxPabs = -sMxslIni
                else
                    valMx = sMxCPIte
                    valMxPabs = sMxCPIni
                end if
                ! le critere et sa dérivée en ii=3 c'est le CP linéarisé avec My de signe opposé
                if (ii .EQ. 4) then
                    ! valMy signe à mettre devant My pour le calcul
                    valMy = -sMyCPIte
                    ! valMyPabs signe de My pour le calcul
                    valMyPabs = -sMyCPIni
                else
                    valMy = sMyCPIte
                    valMyPabs = sMyCPIni
                end if
                ! calcul de la force horizontale prenant en compte les écrouissage cinématiques
                HCP = ((Fite(1)-QCPite(1))**2+(Fite(2)-QCPite(2))**2)**0.5
                ! Calcul pour le critère linéarisé selon H de CP (si nécessaire)
                if (calculH) then
                    ! traitement particulier de dfdF(5:8,1) et dfdF(5:8,2) si la force
                    !   horizontale est nulle
                    if (HCP .LT. r8prem()) then
                        dfdF(ii+4, 1) = sHCPIni(1)*valH
                        dfdF(ii+4, 2) = sHCPIni(2)*valH
                    else
                        dfdF(ii+4, 1) = (Fite(1)-QCPIte(1))/HCP*valH
                        dfdF(ii+4, 2) = (Fite(2)-QCPIte(2))/HCP*valH
                    end if
                    ! ici comme le critère est linéarisé, il n'est pas nécessaire de
                    !   distinguer le cas Fz nul
                    ! calcul de dfdF(5:8,3:5)
                    dfdF(ii+4, 3) = 1.0-(factor/((1.0-2.0*valMx*abs(Fite(4)-QCPIte(4))/ &
                                                  (factor*Ly*(Velas+QCPIte(3))))* &
                                                 (1.0-2.0*valMy*abs(Fite(5)-QCPIte(5))/ &
                                                  (factor*Lx*(Velas+QCPIte(3))))))**(1.0/3.0)
                    dfdF(ii+4, 4) = -2.0*Fite(3)*valMxPabs/ &
                                    (3.0*Lx*(Velas+QCPIte(3))* &
                                     (1.0-2.0*valMy*abs(Fite(5)-QCPIte(5))/ &
                                      (factor*Lx*(Velas+QCPIte(3))))**(1.0/3.0)* &
                                     (1.0-2.0*valMx*abs(Fite(4)-QCPIte(4))/ &
                                      (factor*Ly*(Velas+QCPIte(3))))**(4.0/3.0)*factor**(2.0/3.0))
                    dfdF(ii+4, 5) = -2.0*Fite(3)*valMyPabs/ &
                                    (3.0*Ly*(Velas+QCPIte(3))* &
                                     (1.0-2.0*valMy*abs(Fite(5)-QCPIte(5))/ &
                                      (factor*Lx*(Velas+QCPIte(3))))**(4.0/3.0)* &
                                     (1.0-2.0*valMx*abs(Fite(4)-QCPIte(4))/ &
                                      (factor*Ly*(Velas+QCPIte(3))))**(1.0/3.0)*factor**(2.0/3.0))
                    ! calcul des dérivées du critère en CP linéarisé en H par rapport aux
                    !   écrouissages qui ne sont que l'opposé des dérivées par rapport aux
                    !   forces (car écrouissage cinématique)
                    dfdQCP(ii, :) = -dfdF(ii+4, :)
                    ! sauf dfdQCP(5:8,3) qui est un écrouissage isotrope
                    dfdQCP(ii, 3) = -2.0*Fite(3)*(valMx*abs(Fite(4)-QCPIte(4))*Ly* &
                                                  (1.0-2.0*valMy*abs(Fite(5)-QCPIte(5))/ &
                                                   (factor*Lx*(Velas+QCPite(3))))+ &
                                                  valMy*abs(Fite(5)-QCPIte(5))*Lx* &
                                                  (1.0-2.0*valMx*abs(Fite(4)-QCPIte(4))/ &
                                                   (factor*Ly*(Velas+QCPite(3)))))
                    dfdQCP(ii, 3) = dfdQCP(ii, 3)/(3.0*Ly*Lx*(Velas+QCPite(3))**2* &
                                                   (1.0-2.0*valMy*abs(Fite(5)-QCPIte(5))/ &
                                                    (factor*Lx*(Velas+QCPite(3))))**(4.0/3.0)* &
                                                   (1.0-2.0*valMx*abs(Fite(4)-QCPIte(4))/ &
                                                    (factor*Ly*(Velas+QCPite(3))))**(4.0/3.0)* &
                                                   factor**(2.0/3.0))
                    ! calcul du critère linéarisé en H du CP
                    fsite(ii+4) = (1.0-(factor/((1.0-2.0*valMy*abs(Fite(5)-QCPIte(5))/ &
                                                 (factor*Lx*(Velas+QCPite(3))))* &
                                                (1.0-2.0*valMx*abs(Fite(4)-QCPIte(4))/ &
                                                 (factor*Ly*(Velas+QCPite(3))))))**(1.0/3.0))* &
                                  Fite(3)+HCP*valH
                    fsini(ii+4) = fsite(ii+4)
                    ! on introduit dans le calcul du critère la correction en fonction de
                    !   l'itération pour calcul avant l'inversion
                    do jj = 1, 5
                        fsini(ii+4) = fsini(ii+4)+(tirela(jj)-Fite(jj))*dfdF(ii+4, jj)+ &
                                      (vloc(12+jj)-QCPite(jj))*dfdQCP(ii, jj)
                    end do
                else
                    ! si pas de critère de CP linérisé selon H à calculer, on mets tout à 0
                    dfdF(ii+4, :) = 0.0
                    fsite(ii+4) = 0.0
                    fsini(ii+4) = 0.0
                    dfdQCP(ii, :) = 0.0
                    ! on calcule quand même le critère linéarisé en H de CP pour des
                    !   vérifications ultérieures
                    fsite(5) = (1.0-(factor/((1.0-2.0*abs(Fite(5)-QCPIte(5))/ &
                                              (factor*Lx*(Velas+QCPite(3))))* &
                                             (1.0-2.0*abs(Fite(4)-QCPIte(4))/ &
                                              (factor*Ly*(Velas+QCPite(3))))))**(1.0/3.0)) &
                               *Fite(3)+HCP
                end if
                ! même chose mais critère CP linéarisé en MX
                if (calculMx) then
                    valint = -3.0*Ly*Fite(3)/((Velas+QCPite(3))* &
                                              (1.0-valH*HCP/(factor*(Velas+QCPite(3))))**4* &
                                              (1.0-2.0*valMy*abs(Fite(5)-QCPIte(5))/ &
                                               (factor*Lx*(Velas+QCPite(3)))))
                    if (HCP .LT. error) then
                        dfdF(ii+8, 1) = sHCPIni(1)*valint*valH
                        dfdF(ii+8, 2) = sHCPIni(2)*valint*valH
                    else
                        dfdF(ii+8, 1) = (Fite(1)-QCPIte(1))/HCP*valH*valint
                        dfdF(ii+8, 2) = (Fite(2)-QCPIte(2))/HCP*valH*valint
                    end if
                    dfdF(ii+8, 3) = 1.0-factor/((1.0-valH*HCP/(factor*(Velas+QCPIte(3))))**3* &
                                                (1.0-2.0*valMy*abs(Fite(5)-QCPIte(5))/ &
                                                 (factor*Lx*(Velas+QCPIte(3)))))
                    dfdF(ii+8, 4) = 2.0*valMxPabs
                    dfdF(ii+8, 5) = -2.0*Ly*Fite(3)*valMyPabs/ &
                                    (Lx*(Velas+QCPIte(3))* &
                                     (1.0-valH*HCP/(factor*(Velas+QCPite(3))))**3* &
                                     (1.0-2.0*valMy*abs(Fite(5)-QCPIte(5))/ &
                                      (factor*Lx*(Velas+QCPite(3)))))
                    dfdQCP(ii+4, :) = -dfdF(ii+8, :)
                    dfdQCP(ii+4, 3) = 2.0*valMy*abs(Fite(5)-QCPIte(5))*Ly*Fite(3)/ &
                                      (Lx*(Velas+QCPITe(3))**2* &
                                       (1.0-2.0*valMy*abs(Fite(5)-QCPIte(5))/ &
                                        (factor*Lx*(Velas+QCPIte(3))))**2* &
                                       (1.0-valH*HCP/(factor*(Velas+QCPIte(3))))**3)
                    dfdQCP(ii+4, 3) = dfdQCP(ii+4, 3)+ &
                                      3.0*valH*HCP*Fite(3)/( &
                                      (Velas+QCPITe(3))**2* &
                                      (1.0-2.0*valMy*(Fite(5)-QCPIte(5))/ &
                                       (factor*Lx*(Velas+QCPIte(3))))* &
                                      (1.0-valH*HCP/(factor*(Velas+QCPIte(3))))**4)
                    fsite(ii+8) = (1.0-factor/((1.0-2.0*valMy*abs(Fite(5)-QCPIte(5))/ &
                                                (factor*Lx*(Velas+QCPite(3))))* &
                                               (1.0-valH*HCP/(factor*(Velas+QCPite(3))))**3))* &
                                  Ly*Fite(3)+2.0*valMx*abs(Fite(4)-QCPIte(4))
                    fsini(ii+8) = fsite(ii+8)
                    do jj = 1, 5
                        fsini(ii+8) = fsini(ii+8)+(tirela(jj)-Fite(jj))*dfdF(ii+8, jj)+ &
                                      (vloc(12+jj)-QCPite(jj))*dfdQCP(ii+4, jj)
                    end do
                else
                    dfdF(ii+8, :) = 0.0
                    fsite(ii+8) = 0.0
                    fsini(ii+8) = 0.0
                    dfdQCP(ii+4, :) = 0.0
                    fsite(9) = (1.0-factor/((1.0-2.0*abs(Fite(5)-QCPIte(5))/ &
                                             (factor*Lx*(Velas+QCPite(3))))* &
                                            (1.0-HCP/(factor*(Velas+QCPite(3))))**3))*Ly*Fite(3)+ &
                               2.0*abs(Fite(4)-QCPIte(4))
                end if
                !  même chose mais critère CP linéarisé en MY
                if (calculMy) then
                    valint = -3.0*Ly*Fite(3)/( &
                             (Velas+QCPite(3))* &
                             (1.0-valH*HCP/(factor*(Velas+QCPite(3))))**4* &
                             (1.0-2.0*valMx*abs(Fite(4)-QCPIte(4))/(factor*Ly*(Velas+QCPite(3)))))
                    if (HCP .LT. error) then
                        dfdF(ii+12, 1) = sHCPIni(1)*valint*valH
                        dfdF(ii+12, 2) = sHCPIni(2)*valint*valH
                    else
                        dfdF(ii+12, 1) = (Fite(1)-QCPIte(1))/HCP*valH*valint
                        dfdF(ii+12, 2) = (Fite(2)-QCPIte(2))/HCP*valH*valint
                    end if
                    dfdF(ii+12, 3) = 1.0-factor/( &
                                     (1.0-valH*HCP/(factor*(Velas+QCPIte(3))))**3* &
                                     (1.0-2.0*valMx*abs(Fite(4)-QCPIte(4))/ &
                                      (factor*Ly*(Velas+QCPIte(3)))))
                    dfdF(ii+12, 5) = 2.0*valMyPabs
                    dfdF(ii+12, 4) = -2.0*Lx*Fite(3)*valMxPabs/( &
                                     Ly*(Velas+QCPIte(3))* &
                                     (1.0-valH*HCP/(factor*(Velas+QCPite(3))))**3* &
                                     (1.0-2.0*valMx*abs(Fite(4)-QCPIte(4))/ &
                                      (factor*Ly*(Velas+QCPite(3)))))
                    dfdQCP(ii+8, :) = -dfdF(ii+12, :)
                    dfdQCP(ii+8, 3) = 2.0*valMx*abs(Fite(4)-QCPIte(4))*Lx*Fite(3)/ &
                                      (Ly*(Velas+QCPITe(3))**2* &
                                       (1.0-2.0*valMx*abs(Fite(4)-QCPIte(4))/ &
                                        (factor*Ly*(Velas+QCPIte(3))))**2* &
                                       (1.0-valH*HCP/(factor*(Velas+QCPIte(3))))**3)
                    dfdQCP(ii+8, 3) = dfdQCP(ii+8, 3)+3.0*valH*HCP*Fite(3)/( &
                                      (Velas+QCPITe(3))**2* &
                                      (1.0-2.0*valMx*abs(Fite(4)-QCPIte(4))/ &
                                       (factor*Ly*(Velas+QCPIte(3))))* &
                                      (1.0-valH*HCP/(factor*(Velas+QCPIte(3))))**4)
                    fsite(ii+12) = (1.0-factor/((1.0-2.0*valMx*abs(Fite(4)-QCPIte(4))/ &
                                                 (factor*Ly*(Velas+QCPite(3))))* &
                                                (1.0-valH*HCP/(factor*(Velas+QCPite(3))))**3)) &
                                   *Lx*Fite(3)+ &
                                   2.0*valMy*abs(Fite(5)-QCPIte(5))
                    fsini(ii+12) = fsite(ii+12)
                    do jj = 1, 5
                        fsini(ii+12) = fsini(ii+12)+(tirela(jj)-Fite(jj))* &
                                       dfdF(ii+12, jj)+(vloc(12+jj)-QCPite(jj))*dfdQCP(ii+8, jj)
                    end do
                else
                    dfdF(ii+12, :) = 0.0
                    fsite(ii+12) = 0.0
                    fsini(ii+12) = 0.0
                    dfdQCP(ii+8, :) = 0.0
                    fsite(13) = (1.0-factor/( &
                                 (1.0-2.0*abs(Fite(4)-QCPIte(4))/ &
                                  (factor*Ly*(Velas+QCPite(3))))* &
                                 (1.0-HCP/(factor*(Velas+QCPite(3))))**3))* &
                                Lx*Fite(3)+ &
                                2.0*abs(Fite(5)-QCPIte(5))
                end if
            end do
            ! les matrices globales sont déjà calculées, maintenant on crée le vecteur
            !   de booléen pour savoir quel critère calculer
            ! si criteres en glissement dépassé
            if (calculS) then
                ! glissement de base toujours vérifié
                criteres(1) = .TRUE.
                ! si la cohésion est non nulle alors pas de double critère possible
                if (abs(Cinter*Lx*Ly) .LT. error) then
                    do ii = 2, 4
                        LogiInt = (abs( &
                                   abs(sum(dfdF(ii, :)*dfdF(1, :))**2/ &
                                       (sum(dfdF(ii, :)**2)*sum(dfdF(1, :)**2))) &
                                   -1.0) .LT. errmax)
                        ! on regarde si il y a une différence entre les dérivées du critère en
                        !   glissement, si non, on ne considère pas le critère
                        !   (glissement avec H, Mx ou My de signe opposé)
                        criteres(ii) = .NOT. (LogiInt)
                    end do
                else
                    criteres(2:4) = .FALSE.
                end if
            else
                criteres(:4) = .FALSE.
            end if
            ! si criteres en CP linéarisé selon H dépassé
            if (calculH) then
                ! CP linéarisé selon H  de base toujours vérifié
                criteres(5) = .TRUE.
                ! on ne considère que le critère CP linéarisé selon H avec H de signe
                !   opposé (Mx et My de signe opposé sont compris dans les autres linéarisation)
                criteres(7:8) = .FALSE.
                ! on regarde si il y a une différence entre les dérivées du critère en CP
                !   linéarisé selon H avec ou sans H de signe opposé,
                !   si non, on ne considère pas le critère
                LogiInt = (abs(abs(sum(dfdF(6, :)*dfdF(5, :))**2/ &
                                   (sum(dfdF(6, :)**2)*sum(dfdF(5, :)**2)))-1.0) .LT. errmax)
                criteres(6) = .NOT. (LogiInt)
            else
                criteres(5:8) = .FALSE.
            end if
            ! si criteres en CP linéarisé selon Mx dépassé
            if (calculMx) then
                ! CP linéarisé selon Mx de base toujours vérifié
                criteres(9) = .TRUE.
                ! on ne considère que le critère CP linéarisé selon Mx avec Mx de signe
                !   opposé (H et My de signe opposé sont compris dans les autres linéarisation)
                criteres(10:12) = .FALSE.
                ! on regarde si il y a une différence entre les dérivées du critère en CP
                !   linéarisé selon Mx avec ou sans Mx de signe opposé,
                !   si non, on ne considère pas le critère
                LogiInt = (abs(abs(sum(dfdF(11, :)*dfdF(9, :))**2/ &
                                   (sum(dfdF(11, :)**2)*sum(dfdF(9, :)**2)))-1.0) .LT. errmax)
                criteres(11) = .NOT. (LogiInt)
            else
                criteres(9:12) = .FALSE.
            end if
            if (calculMy) then
                ! si criteres en CP linéarisé selon My dépassé
                ! CP linéarisé selon My de base toujours vérifié
                criteres(13) = .TRUE.
                ! on ne considère que le critère CP linéarisé selon My avec My de signe
                !   opposé (H et Mx de signe opposé sont compris dans les autres linéarisation)
                criteres(14:15) = .FALSE.
                ! on regarde si il y a une différence entre les dérivées du critère en CP
                !   linéarisé selon My avec ou sans My de signe opposé,
                !   si non, on ne considère pas le critère
                LogiInt = (abs(abs(sum(dfdF(16, :)*dfdF(13, :))**2/ &
                                   (sum(dfdF(16, :)**2)*sum(dfdF(13, :)**2)))-1.0) .LT. errmax)
                criteres(16) = .NOT. (LogiInt)
            else
                criteres(13:16) = .FALSE.
            end if
            ! On lance la boucle de convergence sur l'ensemble des critères
            IteCrit = 0
            NoConvCrit = .TRUE.
            ! calcul de la matrice avec tous les éléments de calcul de la matrice à inverser
            ! (dfdF)*K*(dfDF)-dfdQ*I*dfDF
            do ii = 1, 16
                do jj = 1, 16
                    MatAinvTot(ii, jj) = sum(dfdF(ii, :)*KtangPetit(:)*dfDF(jj, :))
                end do
            end do

            do ii = 1, 4
                do jj = 1, 4
                    ! on regarde si il faut rajouter le terme d'écrousissage sliding
                    MatAinvTot(ii, jj) = MatAinvTot(ii, jj)-dfdQsl(ii, 1)*Hslid(1)*dfdF(jj, 1)- &
                                         dfdQsl(ii, 2)*Hslid(2)*dfdF(jj, 2)
                end do
            end do
            do ii = 5, 16
                do jj = 5, 16
                    ! on regarde si il faut rajouter le terme d'écrousissage CP
                    MatAinvTot(ii, jj) = MatAinvTot(ii, jj)-sum(dfdQCP(ii-4, :)*HHCP(:)*dfdF(jj, :))
                end do
            end do

            do ii = 1, 4
                ! calcul du nombre de critères dépassés avec signe opposé (H Mx et My)
                !   pour obtention des ddimension du système à inverser
                nbDDL = 0
                do jj = 1, 4
                    if (criteres((ii-1)*4+jj)) then
                        nbDDL = nbDDL+1
                    end if
                    ! on regarde combien de criteres on a
                end do
                if (nbDDL .GE. 2) then
                    ! on teste toutes les combinaisons du plus grand nombres de DDL au plus petit
                    !   pour garder le plus grand degré de liberté (si plus d'un degré de liberté)
                    ! il s'agit de la variable qui définit quand on a le plus grand nombre de DDL
                    NoConvCritInt = .True.
                    if (nbDDL .EQ. 4) then
                        !on crée le système local à inverser
                        allocate (MatAinverser(nbDDL, nbDDL), fsInverse(nbDDL))
                        MatAinverser = MatAinvTot((4*ii-3):(4*ii), (4*ii-3):(4*ii))
                        fsInverse = fsini((4*ii-3):(4*ii))
                        ! on résout le système les écoulement locaux sont stocké dans fsInverse
                        call mgauss('NCVP', MatAinverser, fsInverse, nbDDL, nbDDL, 1, det, iret)
                        if (iret .EQ. 0) then
                            ! si pas de problème d'inversion
                            ! si pas de problème d'inversion On regarde si tous les
                            !   écoulements sont positifs,
                            !   si oui,il faut garder l'ensemble de ces DDLs'
                            NoConvCritInt = .NOT. ((fsInverse(1) .GT. -r8prem()) .AND. &
                                                   (fsInverse(2) .GT. -r8prem()) .AND. &
                                                   (fsInverse(3) .GT. -r8prem()) .AND. &
                                                   (fsInverse(4) .GT. -r8prem()))
                        end if
                        ! desallocation pour utilisation des mêmes variables
                        !   avec une dimension différente
                        deallocate (MatAinverser, fsInverse)
                    end if
                    ! si une combinaison de 4 criteres s'avère insatisfaisante ou inexistante,
                    !   on regarde pour des dimensions plus petites (3)'
                    if (NoConvCritInt .AND. (nbDDL .GE. 3)) then
                        compt(1) = 1
                        ! on tourne sur les 3 combinaisons de criteres possible
                        !   (à savoir qu'on inclus forcement 4*(ii-1)+1 qui le critère
                        !   sans signe opposé
                        do jj = 1, 3
                            if (NoConvCritInt) then
                                ! on créee le sous-système  sur les combinaisons où criteres
                                !   est vrai pour les 3 DDLs
                                allocate (MatAinverser(3, 3), fsInverse(3))
                                compt(2) = jj+1
                                compt(3) = mod(jj, 3)+2
                                do nbDDLii = 1, 3
                                    if (.NOT. (criteres(4*(ii-1)+compt(nbDDLii)))) then
                                        goto 998
                                    end if
                                    fsInverse(nbDDLii) = fsini(4*(ii-1)+compt(nbDDLii))
                                    do nbDDLjj = 1, 3
                                        MatAinverser(nbDDLii, nbDDLjj) = &
                                            MatAinvTot(4*(ii-1)+compt(nbDDLii), &
                                                       4*(ii-1)+compt(nbDDLjj))
                                    end do
                                end do
                                ! inversion du sous-sytèmes
                                call mgauss('NCVP', MatAinverser, fsInverse, 3, 3, 1, det, iret)
                                ! si pas de problème d'inversion'
                                if (iret .EQ. 0) then
                                    ! On regarde si tous les écoulements sont positifs, si oui,
                                    !   il faut garder l'ensemble de ces DDLs'
                                    NoConvCritInt = .NOT. ((fsInverse(1) .GT. -r8prem()) .AND. &
                                                           (fsInverse(2) .GT. -r8prem()) .AND. &
                                                           (fsInverse(3) .GT. -r8prem()))
                                    if (.NOT. NoConvCritInt) then
                                        criteres(4*(ii-1)+mod(jj+1, 3)+2) = .FALSE.
                                        ! si on a trouvé la combinaison qui marche avec
                                        !   3 criteres on annule le dernier
                                    end if
                                end if
                                ! desallocation pour utilisation des mêmes variables avec
                                !   une dimension différente
998                             continue
                                deallocate (MatAinverser, fsInverse)
                            end if
                        end do
                    end if
                    ! si une combinaison de 3 criteres s'avère insatisfaisante ou inexistante,
                    !   on regarde pour des dimensions plus petites (2) selon les mêmes
                    !   schémas que précédemment
                    if (NoConvCritInt .AND. (nbDDL .GE. 2)) then
                        compt(1) = 1
                        do jj = 1, 3
                            if (NoConvCritInt) then
                                allocate (MatAinverser(2, 2), fsInverse(2))
                                compt(2) = jj+1
                                do nbDDLii = 1, 2
                                    if (.NOT. (criteres(4*(ii-1)+compt(nbDDLii)))) then
                                        goto 997
                                    end if
                                    fsInverse(nbDDLii) = fsini(4*(ii-1)+compt(nbDDLii))
                                    do nbDDLjj = 1, 2
                                        MatAinverser(nbDDLii, nbDDLjj) = &
                                            MatAinvTot(4*(ii-1)+compt(nbDDLii), &
                                                       4*(ii-1)+compt(nbDDLjj))
                                    end do
                                end do
                                call mgauss('NCVP', MatAinverser, fsInverse, 2, 2, 1, det, iret)
                                if (iret .EQ. 0) then
                                    NoConvCritInt = .NOT. ((fsInverse(1) .GT. -r8prem()) .AND. &
                                                           (fsInverse(2) .GT. -r8prem()))
                                    if (.NOT. NoConvCritInt) then
                                        criteres(4*(ii-1)+mod(jj, 3)+2) = .FALSE.
                                        criteres(4*(ii-1)+mod(jj+1, 3)+2) = .FALSE.
                                    else
                                        criteres(4*(ii-1)+jj+1) = .FALSE.
                                    end if
                                end if
997                             continue
                                deallocate (MatAinverser, fsInverse)
                            end if
                        end do
                    end if
                end if
            end do
            ! gestion particulière des multicombinaisons en CP
            if (criteres(5) .AND. criteres(6)) then
                ! si à la fois les criteres CP linéarisé sur H et avec H en signe opposé
                if (criteres(9) .AND. criteres(11)) then
                    ! on regarde si le critère linéarisé en Mx avec signe négatif est présent
                    criteres(11) = .FALSE.
                    ! Si oui on annule ce dernier car redondance avec les 3 autres critères
                end if
                if (criteres(13) .AND. criteres(16)) then
                    ! on regarde si le critère linéarisé en My avec signe négatif est présent
                    criteres(16) = .FALSE.
                    ! Si oui on annule ce dernier car redondance avec les 3 autres critères
                end if
            end if
            if (criteres(9) .AND. criteres(11)) then
                ! si à la fois les criteres CP linéarisé sur Mx et avec M en signe opposé
                if (criteres(13) .AND. criteres(16)) then
                    ! on regarde si le critère linéarisé en My avec signe négatif est présent
                    criteres(16) = .FALSE.
                    ! Si oui on annule ce dernier car redondance avec les 3 autres critères
                end if
            end if
            ! gestion particulière de double combinaisons
            if (criteres(2)) then
                ! si on a le glissement avec H négatif
                if (criteres(6)) then
                    ! on regarde si on a le critère CP linéarisé sur H avec H négatif
                    if (dfdF(6, 3) .GT. dfdF(1, 3)) then
                        ! on garde le critère avec la plus petite pente car c'est celui
                        !   qui est dépassé en premier
                        criteres(6) = .FALSE.
                    else
                        criteres(2) = .FALSE.
                    end if
                end if
                ! on annule les composante avec signe opposé des autres linéarisation
                !   du critère CP car redondante
                criteres(11) = .FALSE.
                criteres(16) = .FALSE.
            end if
            ! dans les cas particuliers sans écrouissage en CP, on ne peut pas avoir de double
            !   critère (glissement et CPh)
            if (criteres(1) .AND. criteres(5) .AND. (abs(vpara(5)) .LE. r8prem()) .AND. &
                ((abs(vpara(10)) .LE. r8prem()) .OR. (abs(tirela(1)) .LE. r8prem())) .AND. &
                ((abs(vpara(12)) .LE. r8prem()) .OR. (abs(tirela(2)) .LE. r8prem()))) then
                ! on garde alors critère à la pente la plus petite
                if (dfdF(5, 3) .GT. dfdF(1, 3)) then
                    criteres(5) = .FALSE.
                else
                    criteres(1) = .FALSE.
                end if
            end if
            ! dans les cas particuliers en H =0 on ne peut pas avoir de double
            !   critère (glissement et CPh)
            if (criteres(1) .AND. criteres(5) .AND. &
                (((Fite(1)-Qslite(1))**2+(Fite(2)-Qslite(2))**2)**0.5 .LE. r8prem()) .AND. &
                (((Fite(1)-QCPite(1))**2+(Fite(2)-QCPite(2))**2)**0.5 .LE. r8prem())) then
                ! on garde alors critère à la pente la plus petite
                if (dfdF(5, 3) .GT. dfdF(1, 3)) then
                    criteres(5) = .FALSE.
                else
                    criteres(1) = .FALSE.
                end if
            end if
            ! les degrés de libertés sont calculés, on lance le calcul d'inversion de
            !   l'ensemble des critères restants
            do while ((NoConvCrit) .and. (IteCrit .LE. nbdecp))
                ! jusqu'à un certain nombre d'itéraation ou qu'on ait convergé avec des
                !   écoulements positifs sur l'ensemble des critère'
                nbDDL = 0
                do ii = 1, 16
                    ! comptage du noimbre de degré de libertés
                    if (criteres(ii)) then
                        nbDDL = nbDDL+1
                    end if
                end do
                if (nbDDL .EQ. 0) then
                    ! si pas de degré de liberté, il y a un problème car on est alors en élasticité
                    iret = 1
                    goto 999
                end if
                allocate (MatAinverser(nbDDL, nbDDL), fsInverse(nbDDL))
                nbDDLii = 0
                ! On extrait la sous matrice que l'on inversera
                do ii = 1, 16
                    if (criteres(ii)) then
                        nbDDLii = nbDDLii+1
                        fsInverse(nbDDLii) = fsini(ii)
                        nbDDLjj = 0
                        do jj = 1, 16
                            if (criteres(jj)) then
                                nbDDLjj = nbDDLjj+1
                                MatAinverser(nbDDLii, nbDDLjj) = MatAinvTot(ii, jj)
                            end if
                        end do
                    end if
                end do
                ! On inverse la sous matrice et le vecteur fsinverse qui recueillera les lambdas
                call mgauss('NCVP', MatAinverser, fsInverse, nbDDL, nbDDL, 1, det, iret)
                if (iret .NE. 0) then
                    goto 999
                end if
                nbDDLii = 0
                NoConvCrit = .FALSE.
                ! on regarde si on converge sur le critère
                do ii = 1, 16
                    if (criteres(ii)) then
                        nbDDLii = nbDDLii+1
                        ! on regarde si on a un écoulement négatif
                        if (fsInverse(nbDDLii) .LE. -r8prem()) then
                            ! dans ce cas on considère que ce critère n'est en réalité pas
                            !   franchi et on reprend le calcul sans ce critère'
                            NoConvCrit = .TRUE.
                            criteres(ii) = .FALSE.
                        end if
                    end if
                end do
                ! si convergence les écoulements sont stockés dans fsret en mettant 0
                !   pour les écoulements non concernés
                if (.NOT. (NoConvCrit)) then
                    nbDDLii = 0
                    do ii = 1, 16
                        if (criteres(ii)) then
                            nbDDLii = nbDDLii+1
                            fsret(ii) = fsInverse(nbDDLii)
                        else
                            fsret(ii) = 0.0
                        end if
                    end do
                end if
                ! maintenant on peut désallouer les matrices de travail intermédiaires (bouclées)
                deallocate (fsInverse, MatAinverser)
                IteCrit = IteCrit+1
            end do
            ! récupération des écoulements et cacluls des termes d'itérations
            ! les écoulements sont stockés dans fs, on calcule Fcor, dQcp ite et dQs
            dFcor = KtangPetit*matmul(transpose(dfdF), fsret)
            dUpCP = matmul(transpose(dfdF(5:, :)), fsret(5:))
            dUpsl = matmul(transpose(dfdF(:4, :)), fsret(:4))
            dQCP = HHCP*dUpCP
            dQsl(1) = HSLID(1)*dUpsl(1)
            dQsl(2) = HSLID(2)*dUpsl(2)
            do ii = 1, 5
                Fite(ii) = tirela(ii)-dFcor(ii)
                QCPIte(ii) = vloc(12+ii)+dQCP(ii)
            end do
            ! on regarde l'évolution du déplacement plastique et recalcule R
            if (vpara(9) .GT. r8prem()) then
                QCPite(3) = vpara(8)*(vloc(21)+abs(dUpCP(3)))/(vpara(9)+vloc(21)+abs(dUpCP(3)))
            else
                QCPite(3) = 0.0
            end if
            if (abs(dUpCP(3)) .GT. error/raidTang(3)) then
                ! HHCP(3) est calculé comme un opérateur sécant et non tangeant à la
                !   formulation de l'écrouissage isotrope
                HHCP(3) = (QCPite(3)-vloc(15))/dUpCP(3)
            end if
            Qslite(1) = vloc(11)+dQsl(1)
            Qslite(2) = vloc(12)+dQsl(2)

            ! vérification de convergence
            nbDDLii = 0
            NoConv = .FALSE.
            do ii = 1, 16
                if (criteres(ii)) then
                    ! on vérifie pour les critères dépassés par le tir élastique la bonne
                    !   convergence sur le critère
                    nbDDLii = nbDDLii+1
                    if (abs(fsite(ii)) .GT. error) then
                        NoConv = .TRUE.
                    end if
                else
                    ! pour les autres critères, on vérifie juste qu'ils ne sont pas dépassés
                    !   suites aux convergences précédentes
                    if ((fsite(ii)) .GT. error) then
                        ! si ces critères sont dépassés, on relance la boucle de clacul en
                        !   considérant ces critères supplémentaires dans la convergence
                        if (ii .EQ. 1) then
                            calculS = .TRUE.
                            criteres(1) = .TRUE.
                            NoConv = .TRUE.
                        end if
                        if ((ii .EQ. 5) .AND. (.NOT. criteres(9)) .AND. (.NOT. criteres(13))) then
                            calculH = .TRUE.
                            criteres(5) = .TRUE.
                            NoConv = .TRUE.
                        end if
                        if ((ii .EQ. 9) .AND. (.NOT. criteres(5)) .AND. (.NOT. criteres(13))) then
                            calculMx = .TRUE.
                            criteres(9) = .TRUE.
                            NoConv = .TRUE.
                        end if
                        if ((ii .EQ. 13) .AND. (.NOT. criteres(5)) .AND. (.NOT. criteres(9))) then
                            calculMy = .TRUE.
                            criteres(13) = .TRUE.
                            NoConv = .TRUE.
                        end if
                    end if
                end if
            end do
            ! on sort de la boucle de convergence en incrémentant
            Ite = Ite+1
        end do
        if (NoConv) then
            ! vérification de convergence
            iret = 1
            goto 999
        end if
        ! on vérifie qu'on reste bien sur la plage de valider de la linéarisation'
        if ((-Fite(3)/(Velas+QCPite(3))) .GT. (factor)) then
            calculNormal = .TRUE.
            goto 999
        else
            ! sinon on effectue le calcul sans linéarisation
            calculNormal = .FALSE.
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
        !
        ! définition pour le calcul de matrice des critères dépassé dans vloc(18)
        !   selon une numérotation précise
        etatG = 0.0
        ! alors glissement
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
        etatCPH = 0.0
        ! alors CP
        if (criteres(5)) then
            if (criteres(6)) then
                if (criteres(7)) then
                    if (criteres(8)) then
                        ! H, Mx et My concomittant
                        etatCPH = 8.0
                    else
                        ! H et Mx concomittant
                        etatCPH = 5.0
                    end if
                else
                    if (criteres(8)) then
                        ! H et My concomittant
                        etatCPH = 6.0
                    else
                        ! H concomittant seul
                        etatCPH = 2.0
                    end if
                end if
            else
                if (criteres(7)) then
                    if (criteres(8)) then
                        ! Mx et My concomittant
                        etatCPH = 7.0
                    else
                        ! Mx concomittant seul
                        etatCPH = 3.0
                    end if
                else
                    if (criteres(8)) then
                        !My concomittant seul
                        etatCPH = 4.0
                    else
                        etatCPH = 1.0
                    end if
                end if
            end if
        end if
        etatCPMx = 0.0
        ! alors CP
        if (criteres(9)) then
            if (criteres(10)) then
                if (criteres(11)) then
                    if (criteres(12)) then
                        ! H, Mx et My concomittant
                        etatCPMx = 8.0
                    else
                        ! H et Mx concomittant
                        etatCPMx = 5.0
                    end if
                else
                    if (criteres(12)) then
                        ! H et My concomittant
                        etatCPMx = 6.0
                    else
                        ! H concomittant seul
                        etatCPMx = 2.0
                    end if
                end if
            else
                if (criteres(11)) then
                    if (criteres(12)) then
                        ! Mx et My concomittant
                        etatCPMx = 7.0
                    else
                        ! Mx concomittant seul
                        etatCPMx = 3.0
                    end if
                else
                    if (criteres(12)) then
                        !My concomittant seul
                        etatCPMx = 4.0
                    else
                        etatCPMx = 1.0
                    end if
                end if
            end if
        end if
        etatCPMy = 0.0
        ! alors CP
        if (criteres(13)) then
            if (criteres(14)) then
                if (criteres(15)) then
                    if (criteres(16)) then
                        ! H, Mx et My concomittant
                        etatCPMy = 8.0
                    else
                        ! H et Mx concomittant
                        etatCPMy = 5.0
                    end if
                else
                    if (criteres(16)) then
                        ! H et My concomittant
                        etatCPMy = 6.0
                    else
                        ! H concomittant seul
                        etatCPMy = 2.0
                    end if
                end if
            else
                if (criteres(15)) then
                    if (criteres(16)) then
                        ! Mx et My concomittant
                        etatCPMy = 7.0
                    else
                        ! Mx concomittant seul
                        etatCPMy = 3.0
                    end if
                else
                    if (criteres(16)) then
                        ! My concomittant seul
                        etatCPMy = 4.0
                    else
                        etatCPMy = 1.0
                    end if
                end if
            end if
        end if
        if (calculPetitH) then
            etatfactor = 1.0
        else
            etatfactor = 2.0
        end if
        vloc(18) = etatfactor*10000.0+etatCPMy*1000.0+etatCPMx*100.0+etatCPH*10.0+etatG
    end if
    !
999 continue
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

end subroutine
