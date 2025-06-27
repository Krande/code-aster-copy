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
subroutine difondmatpetit(tirela, raidTang, vloc, vpara, nbVloc, nbPara, klr, errmax, dulmat, iret)
!
! --------------------------------------------------------------------------------------------------
!  IN
!     tirela    : Tir élastique en entrée
!     nbVloc    : nombre de paramètre de la loi locale
!     vloc      : Valeur des Paramètres de la loi locale
!     raidTang  : Raideur tangeante
!     nbPara    : nombre variables locales
!     vpara     : Variables locales en entrée
!     errmax    : erreur
!     dulmat    : déplacement à l'incrément pour repérer l'état de H
!
!  OUT
!     klr      : matrice symétrique de retour
!     iret     : code de retour du matériau, on le met à 1 quand les convergences sont dépassées
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
    use linalg_ops_module, only: as_matmul
    implicit none
!
#include "asterc/r8prem.h"
#include "asterfort/mgauss.h"
!
    integer(kind=8)      :: nbVloc, nbPara, iret
    real(kind=8) :: vpara(nbPara), tirela(6), raidTang(6), vloc(nbVloc), klr(78), errmax, dulmat(4)
!
! --------------------------------------------------------------------------------------------------
!
!   nombre et compteur de critères à vérifier
    integer(kind=8) :: nbDDL, nbDDLii, nbDDLjj
!   Compteurs
    integer(kind=8) :: compt, ii, jj
!   numérotation pour savoir quels critères sont concernés H et H- pour le transfert
    integer(kind=8) :: etatG, etatHCP, etatMxCP, etatMyCP
!   valeur avec le signe des moments et de l'effort horizontal pour le calcul
!     des différents critères
    real(kind=8) :: valH, valMx, valMy, valMxPabs, valMyPabs
!   paramètres important du macroélément renommé
    real(kind=8) :: Lx, Ly, phi, Velas, error, Cinter
!   Force horizontales et excentrement pour les deux critères
!     avec les écrouissages cinématiques inclus
    real(kind=8) :: Hsl, HCP, Fxint, Fyint
!   valeur de calcul intermédiaire
    real(kind=8) :: valint, det
!
!   HHCP matrice diagonale qui relie écrouissage CP et déplacement plastique
!   Hslid matrice diagonale qui relie écrouissage slid et déplacement plastique
    real(kind=8) :: HHCP(5), Hslid(2)
!   matrices globales des dérivées du critères
    real(kind=8) :: dfdF(16, 5), dfdQCP(12, 5), dfdQsl(4, 2)
!
!   MatAinverser    : la matrice à inverser pour obtenir les avancements
!   fsInverse       : les deltas(finaux)
!   dfOrdo          : récupération des dérivés des critères uniquement utiles au calcul
    real(kind=8), allocatable :: MatAinverser(:, :), fsInverse(:, :), dfOrdo(:, :)
    integer(kind=8), allocatable :: assoddl(:)
!   La partie plastique de la raideur
    real(kind=8) :: Kplastique(5, 5)
!   factor  de linéarisation,reste dans des division euclidienne
    real(kind=8) :: factor, reste
!
!   calculS     : pour savoir si on calcule le critère en glissement
!   calculH     : pour savoir si on calcule le critère en CP linéarisé selon H
!   calculMx    : pour savoir si on calcule le critère en CP linéarisé selon Mx
!   calculMy    : pour savoir si on calcule le critère en CP linéarisé selon My
    logical :: calculS, calculH, calculMx, calculMy
!   listes des critères à regarder
    logical :: criteres(16)
!
! --------------------------------------------------------------------------------------------------
!   initialisation de paramètre de calcul pour simplicité des écritures
!
    iret = 0
    ! longueur de la fondation selon x
    Lx = vpara(1)
    ! longueur de la fondation selon y
    Ly = vpara(2)
    ! tangeante de l'angle phi
    phi = vpara(3)/45.0*Datan(1.D0)
    ! Capacité portant où on estime que la fondation est dans le domaine linéaire
    Velas = vpara(7)
    Cinter = vpara(4)
    ! erreur sur le critère autorisée ramenée à Vélastique
    error = Velas*errmax
    !
    ! recupération des critères dépassés sur la base de vloc(18) qui a une numérotation codée
    ! extraction de la valeur du facteur de linéarisation
    if (vloc(18) .LT. 9999.0) then
        ! on ne devrait pas se trouver dans cette subroutine
        iret = 1
        goto 999
    else if (vloc(18) .GT. 19999.0) then
        factor = 0.10
    else
        factor = 0.01
    end if
    ! récupération des états
    ! en glissement
    etatG = nint(mod(vloc(18), 10.0))
    criteres(1:4) = .FALSE.
    if (.NOT. (etatG .EQ. 0)) then
        criteres(1) = .True.
        if (etatG .EQ. 2) then
            criteres(2) = .True.
        else if (etatG .EQ. 3) then
            criteres(3) = .True.
        else if (etatG .EQ. 4) then
            criteres(4) = .True.
        else if (etatG .EQ. 5) then
            criteres(2:3) = .True.
        else if (etatG .EQ. 6) then
            criteres(2) = .True.
            criteres(4) = .True.
        else if (etatG .EQ. 7) then
            criteres(3:4) = .True.
        else if (etatG .EQ. 8) then
            criteres(2:4) = .True.
        end if
    end if
    !
    calculS = criteres(1)
    ! en critère CP selon H
    reste = (vloc(18)-etatG)/10.0
    etatHCP = nint(mod(reste, 10.0))
    criteres(5:8) = .FALSE.
    if (.NOT. (etatHCP .EQ. 0)) then
        criteres(5) = .True.
        if (etatHCP .EQ. 2) then
            criteres(6) = .True.
        else if (etatHCP .EQ. 3) then
            criteres(7) = .True.
        else if (etatHCP .EQ. 4) then
            criteres(8) = .True.
        else if (etatHCP .EQ. 5) then
            criteres(6:7) = .True.
        else if (etatHCP .EQ. 6) then
            criteres(6) = .True.
            criteres(8) = .True.
        else if (etatHCP .EQ. 7) then
            criteres(7:8) = .True.
        else if (etatHCP .EQ. 8) then
            criteres(6:8) = .True.
        end if
    end if
    calculH = criteres(5)
    ! en critère CP selon Mx
    reste = (reste-etatHCP)/10
    etatMxCP = nint(mod(reste, 10.0))
    criteres(9:12) = .FALSE.
    if (.NOT. (etatMxCP .EQ. 0)) then
        criteres(9) = .True.
        if (etatMxCP .EQ. 2) then
            criteres(10) = .True.
        else if (etatMxCP .EQ. 3) then
            criteres(11) = .True.
        else if (etatMxCP .EQ. 4) then
            criteres(12) = .True.
        else if (etatMxCP .EQ. 5) then
            criteres(10:11) = .True.
        else if (etatMxCP .EQ. 6) then
            criteres(10) = .True.
            criteres(12) = .True.
        else if (etatMxCP .EQ. 7) then
            criteres(11:12) = .True.
        else if (etatMxCP .EQ. 8) then
            criteres(10:12) = .True.
        end if
    end if
    calculMx = criteres(9)
    ! en critère CP selon My
    reste = (reste-etatMxCP)/10.0
    etatMyCP = nint(mod(reste, 10.0))
    criteres(13:16) = .FALSE.
    if (.NOT. (etatMyCP .EQ. 0)) then
        criteres(13) = .True.
        if (etatMyCP .EQ. 2) then
            criteres(14) = .True.
        else if (etatMyCP .EQ. 3) then
            criteres(15) = .True.
        else if (etatMyCP .EQ. 4) then
            criteres(16) = .True.
        else if (etatMyCP .EQ. 5) then
            criteres(14:15) = .True.
        else if (etatMyCP .EQ. 6) then
            criteres(14) = .True.
            criteres(16) = .True.
        else if (etatMyCP .EQ. 7) then
            criteres(15:16) = .True.
        else if (etatMyCP .EQ. 8) then
            criteres(14:16) = .True.
        end if
    end if
    calculMy = criteres(13)
    ! calcul des éléments df/dF
    do ii = 1, 4
        ! ici pour le glissement
        if (calculS) then
            ! création des Variables intermédiaires de signe des efforts et moment
            ! le critere et sa dérivée en ii=2 correspond au glissement avec H de signe opposé
            if (ii .EQ. 2) then
                valH = -1.0
            else
                valH = 1.0
            end if
            ! le critere et sa dérivée en ii=3 correspond au glissement avec Mx de signe opposé
            if (ii .EQ. 3) then
                ! valMx signe à mettre devant Mx pour le calcul
                valMx = -1.0
                ! valMxPabs signe de Mx pour le calcul
                valMxPabs = -signe(tirela(4), error*Ly)
            else
                valMx = 1.0
                valMxPabs = signe(tirela(4), error*Ly)
            end if
            ! le critere et sa dérivée en ii=4 correspond au glissement avec My de signe opposé
            if (ii .EQ. 4) then
                ! valMy signe à mettre devant My pour le calcul
                valMy = -1.0
                ! valMyPabs signe de My pour le calcul
                valMyPabs = -signe(tirela(5), error*Lx)
            else
                valMy = 1.0
                valMyPabs = signe(tirela(5), error*Lx)
            end if
            ! calcul de la force horizontale
            Hsl = ((tirela(1)-vloc(11))**2+(tirela(2)-vloc(12))**2)**0.5
            if (Hsl .LT. r8prem()) then
                ! traitement particulier de dfdF(ii,1) et dfdF(ii,2) si la force
                ! horizontale est nulle
                if (criteres(1) .AND. criteres(2)) then
                    ! mais dans le cas où on est au point 0000
                    ! il faut bien créé le double critère mais de la bonne façon
                    ! on recalcule le déplacement
                    Fxint = dulmat(1)/raidtang(1)
                    Fyint = dulmat(2)/raidtang(2)
                    if (((Fxint)**2+(Fyint)**2)**0.5 .LT. r8prem()) then
                        dfdF(ii, 1) = 0
                        dfdF(ii, 2) = 0
                        criteres(2) = .False.
                    else
                        dfdF(ii, 1) = (Fxint)/(((Fxint)**2+(Fyint)**2)**0.5)*valH
                        dfdF(ii, 2) = (Fyint)/(((Fxint)**2+(Fyint)**2)**0.5)*valH
                    end if
                else
                    dfdF(ii, 1) = 0
                    dfdF(ii, 2) = 0
                end if
                ! pas de possibilité de double critère en H négatif
                criteres(2) = .FALSE.
            else
                dfdF(ii, 1) = (tirela(1)-vloc(11))/Hsl*valH
                dfdF(ii, 2) = (tirela(2)-vloc(12))/Hsl*valH
            end if
            ! calcul de dfdF(ii,3) dfdF(ii,4) dfdF(ii,5)
            ! dans le cas d'une force verticale en compression
            if (tirela(3) .LT. -error) then
                ! attention si les moments ne sont pas trops forts
                if ((abs(tirela(5)) .LT. (-Lx*tirela(3)/2.0)) .AND. &
                    (abs(tirela(4)) .LT. (-Ly*tirela(3)/2.0))) then
                    dfdF(ii, 3) = tan(phi)+2.0*Cinter*Lx*Ly*valMy*abs(tirela(5))/ &
                                  (Lx*tirela(3)**2)*(1.0+2.0*valMx*abs(tirela(4))/(Ly*tirela(3)))
                    dfdF(ii, 3) = dfdF(ii, 3)+2.0*Cinter*Lx*Ly*valMx* &
                                  abs(tirela(4))/(Ly*tirela(3)**2)* &
                                  (1.0+2.0*valMy*abs(tirela(5))/(Lx*tirela(3)))
                    dfdF(ii, 4) = &
                        -2.0*valMxPabs*Cinter*Lx*Ly* &
                        (1+2*valMy*abs(tirela(5))/(Lx*tirela(3)))/(Ly*tirela(3))
                    dfdF(ii, 5) = &
                        -2.0*valMyPabs*Cinter*Lx*Ly* &
                        (1+2*valMx*abs(tirela(4))/(Ly*tirela(3)))/(Lx*tirela(3))
                else
                    ! si les moments sont trop forts, on annule la composante liée à la cohésion
                    dfdF(ii, 3) = tan(phi)
                    dfdF(ii, 4:5) = 0.0
                end if
            else
                ! si force verticale en traction
                dfdF(ii, 3) = tan(phi)
                dfdF(ii, 4:5) = 0.0
            end if
            ! calcul des dérivées du critère en glissement par rapport aux écrouissages
            !   qui ne sont que l'opposé des dérivées par rapport aux force
            !   car écrouissage cinématique
            dfdQsl(ii, 1) = -dfdF(ii, 1)
            dfdQsl(ii, 2) = -dfdF(ii, 2)
        else
            ! si pas de critère de glissement à calculer, on mets tout à 0
            dfdF(ii, :) = 0.0
            dfdQsl(ii, :) = 0.0
        end if
    end do
    ! calcul pour le critère de CP
    do ii = 1, 4
        ! création des Variables intermédiaires de signe des efforts et moment (idem glissement)
        if (ii .EQ. 2) then
            valH = -1.0
        else
            valH = 1.0
        end if
        if (ii .EQ. 3) then
            valMx = -1.0
            ! gestion des doubles criteres en 0
            if (abs(tirela(4)-vloc(16)) .LT. error*Ly) then
                if (criteres(5) .AND. criteres(7)) then
                    valMxPabs = -1
                else
                    valMxPabs = 0
                end if
            else
                valMxPabs = -signe(tirela(4)-vloc(16), error*Ly)
            end if
            valMxPabs = -signe(tirela(4)-vloc(16), error*Ly)
        else
            valMx = 1.0
            ! gestion des doubles criteres en 0
            if (abs(tirela(4)-vloc(16)) .LT. error*Ly) then
                if (criteres(5) .AND. criteres(7)) then
                    valMxPabs = 1
                else
                    valMxPabs = 0
                end if
            else
                valMxPabs = signe(tirela(4)-vloc(16), error*Ly)
            end if
        end if
        if (ii .EQ. 4) then
            valMy = -1.0
            ! gestion des doubles criteres en 0
            if (abs(tirela(5)-vloc(17)) .LT. error*Lx) then
                if (criteres(5) .AND. criteres(8)) then
                    valMyPabs = -1
                else
                    valMyPabs = 0
                end if
            else
                valMyPabs = -signe(tirela(5)-vloc(17), error*Lx)
            end if
        else
            valMy = 1.0
            ! gestion des doubles criteres en 0
            if (abs(tirela(5)-vloc(17)) .LT. error*Lx) then
                if (criteres(5) .AND. criteres(8)) then
                    valMyPabs = 1
                else
                    valMyPabs = 0
                end if
            else
                valMyPabs = signe(tirela(5)-vloc(17), error*Lx)
            end if
        end if
        ! calcul de la force horizontale prenant en compte les écrouissage cinématiques
        HCP = ((tirela(1)-vloc(13))**2+(tirela(2)-vloc(14))**2)**0.5
        if (calculH) then
            ! Calcul pour le critère linéarisé selon H de CP (si nécessaire)
            if (HCP .LT. r8prem()) then
                if (criteres(5) .AND. criteres(6)) then
                    ! mais dans le cas où on est au point 0000
                    ! il faut bien créé le double critère mais de la bonne façon
                    ! on recalcule le déplacement
                    Fxint = dulmat(1)/raidtang(1)
                    Fyint = dulmat(2)/raidtang(2)
                    if (((Fxint)**2+(Fyint)**2)**0.5 .LT. r8prem()) then
                        dfdF(ii+4, 1) = 0
                        dfdF(ii+4, 2) = 0
                        criteres(6) = .False.
                    else
                        dfdF(ii+4, 1) = (Fxint)/(((Fxint)**2+(Fyint)**2)**0.5)*valH
                        dfdF(ii+4, 2) = (Fyint)/(((Fxint)**2+(Fyint)**2)**0.5)*valH
                    end if
                else
                    ! traitement particulier de dfdF(5:8,1) et dfdF(5:8,2) si la force
                    ! horizontale est nulle
                    dfdF(ii+4, 1) = 0.0
                    dfdF(ii+4, 2) = 0.0
                    ! suppression dede l'inconnu en HCP
                    criteres(6) = .False.
                end if
            else
                dfdF(ii+4, 1) = (tirela(1)-vloc(13))/HCP*valH
                dfdF(ii+4, 2) = (tirela(2)-vloc(14))/HCP*valH
            end if
            ! Le critère est linéarisé, il n'est pas nécessaire de distinguer le cas Fz nul
            ! calcul de dfdF(5:8,3:5)
            dfdF(ii+4, 3) = 1.0-(factor/((1.0-2.0*valMx*abs(tirela(4)-vloc(16))/ &
                                          (factor*Ly*(Velas+vloc(15))))* &
                                         (1.0-2.0*valMy*abs(tirela(5)-vloc(17))/ &
                                          (factor*Lx*(Velas+vloc(15))))))**(1.0/3.0)
            dfdF(ii+4, 4) = &
                -2.0*tirela(3)*valMxPabs/ &
                (3*Lx*(Velas+vloc(15))* &
                 (1.0-2.0*valMy*abs(tirela(5)-vloc(17))/ &
                  (factor*Lx*(Velas+vloc(15))))**(1.0/3.0)* &
                 (1.0-2.0*valMx*abs(tirela(4)-vloc(16))/ &
                  (factor*Ly*(Velas+vloc(15))))**(4.0/3.0)*factor**(2.0/3.0))
            dfdF(ii+4, 5) = &
                -2.0*tirela(3)*valMyPabs/ &
                (3.0*Ly*(Velas+vloc(15))* &
                 (1.0-2.0*valMy*abs(tirela(5)-vloc(17))/ &
                  (factor*Lx*(Velas+vloc(15))))**(4.0/3.0)* &
                 (1.0-2.0*valMx*abs(tirela(4)-vloc(16))/ &
                  (factor*Ly*(Velas+vloc(15))))**(1.0/3.0)*factor**(2.0/3.0))
            ! calcul des dérivées du critère en CP linéarisé en H par rapport aux écrouissages
            !   qui ne sont que l'opposé des dérivées par rapport aux forces
            !   car écrouissage cinématique
            dfdQCP(ii, :) = -dfdF(ii+4, :)
            ! sauf dfdQCP(5:8,3) qui est un écrouissage isotrope
            dfdQCP(ii, 3) = &
                -2.0*tirela(3)* &
                (valMx*abs(tirela(4)-vloc(16))*Ly* &
                 (1.0-2.0*valMy*abs(tirela(5)-vloc(17))/(factor*Lx*(Velas+vloc(15))))+ &
                 valMy*abs(tirela(5)-vloc(17))*Lx* &
                 (1.0-2.0*valMx*abs(tirela(4)-vloc(16))/(factor*Ly*(Velas+vloc(15)))))
            dfdQCP(ii, 3) = &
                dfdQCP(ii, 3)/(3.0*Ly*Lx*(Velas+vloc(15))**2* &
                               (1.0-2.0*valMy*abs(tirela(5)-vloc(17))/ &
                                (factor*Lx*(Velas+vloc(15))))**(4.0/3.0)* &
                               (1.0-2.0*valMx*abs(tirela(4)-vloc(16))/ &
                                (factor*Ly*(Velas+vloc(15))))**(4.0/3.0)*factor**(2.0/3.0))
        else
            ! si pas de critère de CP linérisé selon H à calculer, on mets tout à 0
            dfdF(ii+4, :) = 0.0
            dfdQCP(ii, :) = 0.0
        end if
        ! même chose mais critère CP linéarisé selon Mx
        if (calculMx) then
            if (abs(tirela(4)-vloc(16)) .LT. error) then
                criteres(11) = .FALSE.
            end if
            valint = -3.0*Ly*tirela(3)/ &
                     ((Velas+vloc(15))*(1.0-valH*HCP/(factor*(Velas+vloc(15))))**4* &
                      (1.0-2.0*valMy*abs(tirela(5)-vloc(17))/(factor*Lx*(Velas+vloc(15)))))
            if (HCP .LT. r8prem()) then
                if (criteres(9) .AND. criteres(10)) then
                    ! mais dans le cas où on est au point 0000
                    ! il faut bien créé le double critère mais de la bonne façon
                    ! on recalcule le déplacement
                    Fxint = dulmat(1)/raidtang(1)
                    Fyint = dulmat(2)/raidtang(2)
                    if (((Fxint)**2+(Fyint)**2)**0.5 .LT. r8prem()) then
                        dfdF(ii+8, 1) = 0
                        dfdF(ii+8, 2) = 0
                        criteres(10) = .False.
                    else
                        dfdF(ii+8, 1) = (Fxint)/(((Fxint)**2+(Fyint)**2)**0.5)*valH
                        dfdF(ii+8, 2) = (Fyint)/(((Fxint)**2+(Fyint)**2)**0.5)*valH
                    end if
                else
                    dfdF(ii+8, 1) = 0.0
                    dfdF(ii+8, 2) = 0.0
                    criteres(10) = .False.
                end if
            else
                dfdF(ii+8, 1) = (tirela(1)-vloc(13))/HCP*valH*valint
                dfdF(ii+8, 2) = (tirela(2)-vloc(14))/HCP*valH*valint
            end if
            dfdF(ii+8, 3) = 1.0-factor/((1-valH*HCP/(factor*(Velas+vloc(15))))**3* &
                                        (1.0-2.0*valMy*abs(tirela(5)-vloc(17))/ &
                                         (factor*Lx*(Velas+vloc(15)))))
            dfdF(ii+8, 4) = 2.0*valMxPabs
            dfdF(ii+8, 5) = &
                -2.0*Ly*tirela(3)* &
                valMyPabs/(Lx*(Velas+vloc(15))*(1-valH*HCP/(factor*(Velas+vloc(15))))**3* &
                           (1.0-2.0*valMy*abs(tirela(5)-vloc(17))/(factor*Lx*(Velas+vloc(15)))))
            dfdQCP(ii+4, :) = -dfdF(ii+8, :)
            dfdQCP(ii+4, 3) = &
                2.0*valMy*abs(tirela(5)-vloc(17))*Ly*tirela(3)/ &
                (Lx*(Velas+vloc(15))**2* &
                 (1.0-2.0*valMy*abs(tirela(5)-vloc(17))/(factor*Lx*(Velas+vloc(15))))**2* &
                 (1.0-valH*HCP/(factor*(Velas+vloc(15))))**3)
            dfdQCP(ii+4, 3) = &
                dfdQCP(ii+4, 3)+3.0*valH*HCP* &
                tirela(3)/((Velas+vloc(15))**2* &
                           (1.0-2.0*valMy*(tirela(5)-vloc(17))/(factor*Lx*(Velas+vloc(15))))* &
                           (1.0-valH*HCP/(factor*(Velas+vloc(15))))**4)
        else
            dfdF(ii+8, :) = 0.0
            dfdQCP(ii+4, :) = 0.0
        end if
        ! même chose mais critère CP linéarisé selon Mx
        if (calculMy) then
            if (abs(tirela(5)-vloc(17)) .LT. error) then
                criteres(16) = .FALSE.
            end if
            valint = -3.0*Ly*tirela(3)/ &
                     ((Velas+vloc(15))*(1.0-valH*HCP/(factor*(Velas+vloc(15))))**4* &
                      (1.0-2.0*valMx*abs(tirela(4)-vloc(16))/(factor*Ly*(Velas+vloc(15)))))
            if (HCP .LT. r8prem()) then
                if (criteres(13) .AND. criteres(14)) then
                    ! mais dans le cas où on est au point 0000
                    ! il faut bien créé le double critère mais de la bonne façon
                    ! on recalcule le déplacement
                    Fxint = dulmat(1)/raidtang(1)
                    Fyint = dulmat(2)/raidtang(2)
                    if (((Fxint)**2+(Fyint)**2)**0.5 .LT. r8prem()) then
                        dfdF(ii+12, 1) = 0
                        dfdF(ii+12, 2) = 0
                        criteres(14) = .False.
                    else
                        dfdF(ii+12, 1) = (Fxint)/(((Fxint)**2+(Fyint)**2)**0.5)*valH
                        dfdF(ii+12, 2) = (Fyint)/(((Fxint)**2+(Fyint)**2)**0.5)*valH
                    end if
                else
                    dfdF(ii+12, 1) = 0.0
                    dfdF(ii+12, 2) = 0.0
                    criteres(14) = .False.
                end if
            else
                dfdF(ii+12, 1) = (tirela(1)-vloc(13))/HCP*valH*valint
                dfdF(ii+12, 2) = (tirela(2)-vloc(15))/HCP*valH*valint
            end if
            dfdF(ii+12, 3) = &
                1-factor/((1-valH*HCP/(factor*(Velas+vloc(15))))**3* &
                          (1.0-2.0*valMx*abs(tirela(4)-vloc(16))/(factor*Ly*(Velas+vloc(15)))))
            dfdF(ii+12, 5) = 2.0*valMyPabs
            dfdF(ii+12, 4) = &
                -2.0*Lx*tirela(3)* &
                valMxPabs/(Ly*(Velas+vloc(15))*(1-valH*HCP/(factor*(Velas+vloc(15))))**3* &
                           (1.0-2.0*valMx*abs(tirela(4)-vloc(16))/(factor*Ly*(Velas+vloc(15)))))
            dfdQCP(ii+8, :) = -dfdF(ii+12, :)
            dfdQCP(ii+8, 3) = &
                2.0*valMx*abs(tirela(4)-vloc(16))*Lx* &
                tirela(3)/ &
                (Ly*(Velas+vloc(15))**2* &
                 (1.0-2.0*valMx*abs(tirela(4)-vloc(16))/(factor*Ly*(Velas+vloc(15))))**2* &
                 (1.0-valH*HCP/(factor*(Velas+vloc(15))))**3)
            dfdQCP(ii+8, 3) = &
                dfdQCP(ii+8, 3)+3.0*valH*HCP*tirela(3)/ &
                ((Velas+vloc(15))**2* &
                 (1.0-2.0*valMx*abs(tirela(4)-vloc(16))/(factor*Ly*(Velas+vloc(15))))* &
                 (1.0-valH*HCP/(factor*(Velas+vloc(15))))**4)
        else
            dfdF(ii+12, :) = 0.0
            dfdQCP(ii+8, :) = 0.0
        end if
    end do
    ! calcul des degré de liberté
    nbDDL = 0
    ! comptage du nombre de degré de libertés
    do ii = 1, 16
        if (criteres(ii)) then
            nbDDL = nbDDL+1
        end if
    end do
    ! si pas de degré de liberté, il y a un problème car on est alors en élasticité
    if (nbDDL .EQ. 0) then
        goto 999
    end if
    !
    ! Calcul préalable des paramètres hors boucle
    ! Calcul de HHCP telle que dqCP=HHCP*dupCP, HHCP(3) traité à aprt
    !   car dépend si le déplacement plastique CP est vers le bas ou le haut
    HHCP(1) = vpara(10)*exp(-1.0*vpara(11)*vloc(21))
    HHCP(2) = vpara(12)*exp(-1.0*vpara(13)*vloc(21))
    HHCP(4) = vpara(14)*exp(-1.0*vpara(15)*vloc(21))
    HHCP(5) = vpara(16)*exp(-1.0*vpara(17)*vloc(21))
    ! traitement du cas DEPL_REFE=0 pas d'écrouissage cinématique pour éviter une division
    !   par zéro pas besoin de déterminer le signe car on dans les Fz petits
    if (vpara(9) .GT. r8prem()) then
        HHCP(3) = vpara(8)*vpara(9)/(vpara(9)+vloc(21))**2
    else
        HHCP(3) = 0.0
    end if
    ! Calcul de Hslid telle que dqs=Hslid*dupslid
    Hslid(1) = vpara(5)*exp(-1.0*vpara(6)*vloc(19))
    Hslid(2) = vpara(5)*exp(-1.0*vpara(6)*vloc(20))
    !
    ! création de la matrice A inverser par la suite
    allocate (MatAinverser(nbDDL, nbDDL), fsInverse(nbDDL, nbDDL), &
              dfOrdo(nbDDL, 5), assoddl(nbDDL))
    ! matrice identite
    forall (ii=1:nbDDL, jj=1:nbDDL) fsInverse(ii, jj) = (ii/jj)*(jj/ii)

    nbDDLii = 0
    do ii = 1, 16
        if (criteres(ii)) then
            nbDDLii = nbDDLii+1
            assoddl(nbDDLii) = ii
            dfOrdo(nbDDLii, :) = dfdF(ii, :)
        end if
    end do

    do nbDDLii = 1, nbDDL
        ii = assoddl(nbDDLii)
        do nbDDLjj = 1, nbDDL
            jj = assoddl(nbDDLjj)
            if ((ii .le. 4) .and. (jj .le. 4)) then
                MatAinverser(nbDDLii, nbDDLjj) = &
                    sum(dfdF(ii, :)*raidTang(:5)*dfdF(jj, :))- &
                    dfdQsl(ii, 1)*Hslid(1)*dfdF(jj, 1)-dfdQsl(ii, 2)*Hslid(2)*dfdF(jj, 2)
            else if ((ii .ge. 5) .and. (jj .ge. 5)) then
                MatAinverser(nbDDLii, nbDDLjj) = &
                    sum(dfdF(ii, :)*raidTang(:5)*dfdF(jj, :))- &
                    sum(dfdQCP(ii-4, :)*HHCP(:)*dfdF(jj, :))
            else
                MatAinverser(nbDDLii, nbDDLjj) = sum(dfdF(ii, :)*raidTang(:5)*dfdF(jj, :))
            end if
        end do
    end do

    ! on inverse la matrice entièrement
    call mgauss('NCVP', MatAinverser, fsInverse, nbDDL, nbDDL, nbDDL, det, iret)
    if (iret .NE. 0) then
        goto 999
    end if
    ! on reconstruit alors la matrice tangente plastique
    Kplastique = as_matmul(as_matmul(transpose(dfOrdo), fsInverse), dfOrdo)
    compt = 0
    do ii = 1, 5
        do jj = 1, ii
            Kplastique(ii, jj) = Kplastique(ii, jj)*raidTang(ii)*raidTang(jj)
            compt = compt+1
            klr(compt) = klr(compt)-Kplastique(ii, jj)
        end do
    end do

999 continue
    !
    if (allocated(MatAinverser)) then
        deallocate (MatAinverser)
    end if
    if (allocated(fsInverse)) then
        deallocate (fsInverse)
    end if
    if (allocated(dfOrdo)) then
        deallocate (dfOrdo)
    end if
    if (allocated(assoddl)) then
        deallocate (assoddl)
    end if
    !
contains

! fonction rapide de signe
!   attention ici est égal à 0 pour a=0 contrairement à sign de fortran
    function signe(a, error)
        real(kind=8) :: a, signe, error
        if (a .LT. -error) then
            signe = -1.0
        else if (a .GT. error) then
            signe = +1.0
        else
            signe = 0.0
        end if
    end function

end subroutine
