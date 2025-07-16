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
! person_in_charge: vinicius.alves-fernandes at edf.fr
! contributor: cyril.borely at setec.com
!
subroutine difondmat(tirela, raidTang, vloc, vpara, nbVloc, nbPara, klr, errmax, dulMat, iret)
!
! --------------------------------------------------------------------------------------------------
!  IN
!     tirela   : Tir élastique en entrée
!     nbVloc   : nombre de paramètre de la loi locale
!     vloc     : Paramètre de la loi locale
!     raidTang : Raideur tangeante
!     nbPara   : nombre variables locales
!     valvarloc: Variables locales en entrée
!     errmax   : erreur
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
    real(kind=8) :: vpara(nbPara), tirela(6), raidTang(6), vloc(nbVloc), klr(78), errmax, dulMat(4)
!
! --------------------------------------------------------------------------------------------------
!
!   compteur du nombre d'itération
    integer(kind=8) :: ii, jj, compt
!   nombre de critères à vérifier
    integer(kind=8) :: nbDDL, nbDDLii, nbDDLjj
!
!   valeur avec le signe des moments et de l'effort horizontal pour les différents critères
    real(kind=8) :: valH, valMx, valMy, valMxPabs, valMyPabs
!   paramètres important du macroélément renommé
    real(kind=8) :: Lx, Ly, phi, Velas, error, Cinter
!   Force horizontales et excentrement pour les deux critères
!       (avec les écrouissages cinématiques inclus)
    real(kind=8) :: Hsl, HCP
!   listes des critères à regarder
    logical :: criteres(8)
!   valeur de calcul intermédiaire
    real(kind=8) :: valint, det
!
!   HHCP        : matrice diagonale qui relie écrouissage CP et déplacement plastique
!   Hslid       : matrice diagonale qui relie écrouissage slid et déplacement plastique
!   signeHHCP3  : contient le signe de HHCP(3)
    real(kind=8) :: HHCP(5), Hslid(2), signeHHCP3(5)
!   matrices globales des dérivées du critères
    real(kind=8) :: dfdF(8, 5), dfdQCP(4, 5), dfdQsl(4, 2)
!
!   MatAinverser    : la matrice à inverser pour obtenir les avancements
!   fsInverse       : les deltas(finaux)
!   dfOrdo          : récupération des dérivés des critères uniquement utiles au calcul
    real(kind=8), allocatable :: MatAinverser(:, :), fsInverse(:, :), dfOrdo(:, :)
    integer(kind=8), allocatable :: assoddl(:)
    real(kind=8) :: etatCP, etatG
!   numérotation pour savoir quels critères sont concernés H et H- pour le
!   transfert à la matrice de raideur
!   la partie plastique de la raideur
    real(kind=8) :: Kplastique(5, 5)
!
! --------------------------------------------------------------------------------------------------
!
! initialisation de paramètre de calcul pour simplicité des écritures
!
    iret = 0
    ! longueur de la fondation selon x
    Lx = vpara(1)
    ! longueur de la fondation selon y
    Ly = vpara(2)
    ! tangente de l'angle phi
    phi = vpara(3)/45.0*Datan(1.D0)
    ! Capacité portant où on estime que la fondation est dans le domaine linéaire
    Velas = vpara(7)
    Cinter = vpara(4)
    ! erreur sur le critère autorisée ramenée à Vélastique
    error = Velas*errmax
    !
    ! récupération des critères dépassés sur la base de vloc(18) qui a une numérotation codée
    if (vloc(18) .LT. 1.0) then
        ! élastique, rien besoin de faire
        goto 999
    end if
    ! double critère
    if (vloc(18) .GT. 99.0) then
        criteres(1) = .TRUE.
        criteres(5) = .TRUE.
        etatG = mod(vloc(18), 10.0)
        etatCP = (vloc(18)-100.0-etatG)/10.0
        ! etat CP seul
    else if (vloc(18) .GT. 19.0) then
        criteres(1) = .FALSE.
        criteres(5) = .TRUE.
        etatG = 1.0
        etatCP = (vloc(18)-20.0)
        ! etat SL seul
    else if (vloc(18) .GT. 9.0) then
        criteres(1) = .True.
        criteres(5) = .False.
        etatG = vloc(18)-10.0
        etatCP = 1.0
    else
        criteres(1) = .False.
        criteres(5) = .False.
        etatG = 1.0
        etatCP = 1.0
    end if
    ! récupération des états
    ! en glissement
    if (etatG .LT. 1.5) then
        criteres(2) = .False.
        criteres(3) = .False.
        criteres(4) = .False.
    else if (etatG .LT. 2.5) then
        criteres(2) = .True.
        criteres(3) = .False.
        criteres(4) = .False.
    else if (etatG .LT. 3.5) then
        criteres(2) = .False.
        criteres(3) = .True.
        criteres(4) = .False.
    else if (etatG .LT. 4.5) then
        criteres(2) = .False.
        criteres(3) = .False.
        criteres(4) = .True.
    else if (etatG .LT. 5.5) then
        criteres(2) = .True.
        criteres(3) = .True.
        criteres(4) = .False.
    else if (etatG .LT. 6.5) then
        criteres(2) = .True.
        criteres(3) = .False.
        criteres(4) = .True.
    else if (etatG .LT. 7.5) then
        criteres(2) = .False.
        criteres(3) = .True.
        criteres(4) = .True.
    else
        criteres(2) = .True.
        criteres(3) = .True.
        criteres(4) = .True.
    end if
    ! en critère CP
    if (etatCP .LT. 1.5) then
        criteres(6) = .False.
        criteres(7) = .False.
        criteres(8) = .False.
    else if (etatCP .LT. 2.5) then
        criteres(6) = .True.
        criteres(7) = .False.
        criteres(8) = .False.
    else if (etatCP .LT. 3.5) then
        criteres(6) = .False.
        criteres(7) = .True.
        criteres(8) = .False.
    else if (etatCP .LT. 4.5) then
        criteres(6) = .False.
        criteres(7) = .False.
        criteres(8) = .True.
    else if (etatCP .LT. 5.5) then
        criteres(6) = .True.
        criteres(7) = .True.
        criteres(8) = .False.
    else if (etatCP .LT. 6.5) then
        criteres(6) = .True.
        criteres(7) = .False.
        criteres(8) = .True.
    else if (etatCP .LT. 7.5) then
        criteres(6) = .False.
        criteres(7) = .True.
        criteres(8) = .True.
    else
        criteres(6) = .True.
        criteres(7) = .True.
        criteres(8) = .True.
    end if
    !
    ! comptage du nombre de critères à considérer
    nbDDL = 0
    do ii = 1, 8
        if (criteres(ii)) then
            nbDDL = nbDDL+1
        end if
    end do
    if (nbDDL .EQ. 0) then
        goto 999
    end if
    ! Calcul préalable des paramètres hors boucle
    ! Calcul de HHCP telle que dqCP=HHCP*dupCP, HHCP(3) manquant
    !   car dépend si le déplacement plastique CP est vers le bas ou le haut
    HHCP(1) = vpara(10)*exp(-1.0*vpara(11)*vloc(21))
    HHCP(2) = vpara(12)*exp(-1.0*vpara(13)*vloc(21))
    HHCP(4) = vpara(14)*exp(-1.0*vpara(15)*vloc(21))
    HHCP(5) = vpara(16)*exp(-1.0*vpara(17)*vloc(21))
    ! Calcul de HHCP telle que dqs=Hslid*dupslid
    Hslid(1) = vpara(5)*exp(-1.0*vpara(6)*vloc(19))
    Hslid(2) = vpara(5)*exp(-1.0*vpara(6)*vloc(20))
    !
    ! On calcule la matrice des dfdF dfdQCP et dfdQsl
    do ii = 1, 4
        ! Création des Variables intermédiaires de signe des efforts et moment
        if (ii .EQ. 2) then
            valH = -1.0
        else
            valH = 1.0
        end if
        if (ii .EQ. 3) then
            valMx = -1.0
            valMxPabs = -signe(tirela(4), error*Ly)
        else
            valMx = 1.0
            valMxPabs = signe(tirela(4), error*Ly)
        end if
        if (ii .EQ. 4) then
            valMy = -1.0
            valMyPabs = -signe(tirela(5), error*Lx)
        else
            valMy = 1.0
            valMyPabs = signe(tirela(5), error*Lx)
        end if
        ! Calcul des dérivées relative au glissement
        Hsl = ((tirela(1)-vloc(11))**2.0+(tirela(2)-vloc(12))**2.0)**0.5
        if (Hsl .LT. r8prem()) then
            ! cas spécial de H =0
            dfdF(ii, 1) = 0
            dfdF(ii, 2) = 0
            ! pas de possibilité de double critère en H négatif
            criteres(2) = .FALSE.
        else
            dfdF(ii, 1) = (tirela(1)-vloc(11))/Hsl*valH
            dfdF(ii, 2) = (tirela(2)-vloc(12))/Hsl*valH
        end if
        ! la gestion d'un Fz en traction est déjà gérée' on calcule dfdF(1:4,3:5)
        if ((abs(tirela(5)) .LT. (-Lx*tirela(3))) .AND. (abs(tirela(4)) .LT. (-Ly*tirela(3)))) then
            dfdF(ii, 3) = tan(phi)+2.0*Cinter*Lx*Ly*valMy*abs(tirela(5))/(Lx*tirela(3)**2)* &
                          (1.0+2.0*valMx*abs(tirela(4))/(Ly*tirela(3)))
            dfdF(ii, 3) = dfdF(ii, 3)+2.0*Cinter*Lx*Ly*valMx*abs(tirela(4))/(Ly*tirela(3)**2)* &
                          (1.0+2.0*valMy*abs(tirela(5))/(Lx*tirela(3)))
            dfdF(ii, 4) = -2.0*valMxPabs*Cinter*Lx*Ly* &
                          (1.0+2.0*valMy*abs(tirela(5))/(Lx*tirela(3)))/(Ly*tirela(3))
            dfdF(ii, 5) = -2.0*valMyPabs*Cinter*Lx*Ly* &
                          (1.0+2.0*valMx*abs(tirela(4))/(Ly*tirela(3)))/(Lx*tirela(3))
        else
            dfdF(ii, 3) = tan(phi)
            dfdF(ii, 4:5) = 0.0
        end if
        dfdQsl(ii, 1) = -dfdF(ii, 1)
        dfdQsl(ii, 2) = -dfdF(ii, 2)
    end do
    ! calcul des dérivées relative au CP
    do ii = 1, 4
        ! création des Variables intermédiaires de signe des efforts et moment
        if (ii .EQ. 2) then
            valH = -1.0
        else
            valH = 1.0
        end if
        if (ii .EQ. 3) then
            valMx = -1.0
            valMxPabs = -signe(tirela(4)-vloc(16), error*Ly)
            if (nint(valMxPabs) == 0) then
                if (dulMat(3) >= error*Ly/raidTang(4)) then
                    valMxPabs = -1.0
                else if (dulMat(3) <= -error*Ly/raidTang(4)) then
                    valMxPabs = 1.0
                end if
            end if
        else
            valMx = 1.0
            valMxPabs = signe(tirela(4)-vloc(16), error*Ly)
            if (nint(valMxPabs) == 0) then
                if (dulMat(3) >= error*Ly/raidTang(4)) then
                    valMxPabs = 1.0
                else if (dulMat(3) <= -error*Ly/raidTang(4)) then
                    valMxPabs = -1.0
                end if
            end if

        end if
        if (ii .EQ. 4) then
            valMy = -1.0
            valMyPabs = -signe(tirela(5)-vloc(17), error*Lx)
            if (nint(valMyPabs) == 0) then
                if (dulMat(4) >= error*Lx/raidTang(5)) then
                    valMyPabs = -1.0
                else if (dulMat(4) <= -error*Lx/raidTang(5)) then
                    valMyPabs = 1.0
                end if
            end if
        else
            valMy = 1.0
            valMyPabs = signe(tirela(5)-vloc(17), error*Lx)
            if (nint(valMyPabs) == 0) then
                if (dulMat(4) >= error*Lx/raidTang(5)) then
                    valMyPabs = 1.0
                else if (dulMat(4) <= -error*Lx/raidTang(5)) then
                    valMyPabs = -1.0
                end if
            end if
        end if
        HCP = ((tirela(1)-vloc(13))**2.0+(tirela(2)-vloc(14))**2.0)**0.5
        valint = -3.0*(Velas+vloc(15))*(1.0+2.0*valMy*abs(tirela(5)-vloc(17))/(Lx*tirela(3)))* &
                 (1.0+2.0*valMx*abs(tirela(4)-vloc(16))/(Ly*tirela(3)))*(1+valH*HCP/tirela(3))**2
        ! cas spécial de H =0
        if (HCP .LT. r8prem()) then
            dfdF(ii+4, 1) = 0.0
            dfdF(ii+4, 2) = 0.0
            ! suppression dede l'inconnu en HCP
            criteres(6) = .False.
        else
            dfdF(ii+4, 1) = (tirela(1)-vloc(13))*valH*valint/HCP/tirela(3)
            dfdF(ii+4, 2) = (tirela(2)-vloc(14))*valH*valint/HCP/tirela(3)
        end if
        dfdF(ii+4, 3) = -1.0+ &
                        (Velas+vloc(15))*2.0*valMx* &
                        abs(tirela(4)-vloc(16))/(Ly*tirela(3)**2)* &
                        (1.0+valH*HCP/tirela(3))**3* &
                        (1.0+2.0*valMy*abs(tirela(5)-vloc(17))/(Lx*tirela(3)))
        dfdF(ii+4, 3) = dfdF(ii+4, 3)+ &
                        (Velas+vloc(15))*2.0*valMy*abs(tirela(5)-vloc(17))/(Lx*tirela(3)**2)* &
                        (1.0+valH*HCP/tirela(3))**3* &
                        (1.0+2.0*valMx*abs(tirela(4)-vloc(16))/(Ly*tirela(3)))
        dfdF(ii+4, 3) = dfdF(ii+4, 3)+(Velas+vloc(15))*3.0*valH*HCP/(tirela(3)**2)* &
                        (1.0+valH*HCP/tirela(3))**2* &
                        (1.0+2.0*valMx*abs(tirela(4)-vloc(16))/(Ly*tirela(3)))* &
                        (1.0+2.0*valMy*abs(tirela(5)-vloc(17))/(Lx*tirela(3)))
        dfdF(ii+4, 4) = -2.0*(Velas+vloc(15))/(Ly*tirela(3))*(1+valH*HCP/tirela(3))**3* &
                        (1.0+2.0*valMy*abs(tirela(5)-vloc(17))/(Lx*tirela(3)))*valMxPabs
        dfdF(ii+4, 5) = -2.0*(Velas+vloc(15))/(Lx*tirela(3))*(1+valH*HCP/tirela(3))**3* &
                        (1.0+2.0*valMx*abs(tirela(4)-vloc(16))/(Ly*tirela(3)))*valMyPabs
        ! la gestion d'un Fz en traction est déjà gérée'
        dfdQCP(ii, :) = -dfdF(ii+4, :)
        dfdQCP(ii, 3) = -(1+valH*HCP/tirela(3))**3.0* &
                        (1.0+2.0*valMx*abs(tirela(4)-vloc(16))/(Ly*tirela(3)))* &
                        (1.0+2.0*valMy*abs(tirela(5)-vloc(17))/(Lx*tirela(3)))
    end do
    ! gestion de HHCP(3)
    ! traitement du cas DEPL_REFE=0 pas d'écrouissage cinématique' pour éviter une
    !   division par zéro
    ! pas besoin de déterminer le signe car on dans les Fz petits
    if (vpara(9) .LE. r8prem()) then
        HHCP(3) = 0.0
    else
        HHCP(3) = vpara(8)*vpara(9)/(vpara(9)+vloc(21))**2.0
    end if
    ! création de la matrice de signe de HHCP(3) qui dépendra de dfdF(ii,3)
    signeHHCP3 = 1.0
    !
    ! création de la matrice A inverser par la suite
    ! comptage du nombre de critères à considérer
    nbDDL = 0
    do ii = 1, 8
        if (criteres(ii)) then
            nbDDL = nbDDL+1
        end if
    end do
    if (nbDDL .EQ. 0) then
        goto 999
    end if
    allocate (MatAinverser(nbDDL, nbDDL), fsInverse(nbDDL, nbDDL), dfOrdo(nbDDL, 5), assoddl(nbDDL))
    ! matrice identite
    forall (ii=1:nbDDL, jj=1:nbDDL) fsInverse(ii, jj) = (ii/jj)*(jj/ii)

    nbDDLii = 0
    do ii = 1, 8
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
                MatAinverser(nbDDLii, nbDDLjj) = sum(dfdF(ii, :)*raidTang(:5)*dfdF(jj, :))- &
                                                 dfdQsl(ii, 1)*Hslid(1)*dfdF(jj, 1)- &
                                                 dfdQsl(ii, 2)*Hslid(2)*dfdF(jj, 2)
            else if ((ii .ge. 5) .and. (jj .ge. 5)) then
                if (dfdF(ii, 3) .gt. r8prem()) then
                    signeHHCP3(3) = +1.0
                else
                    signeHHCP3(3) = -1.0
                end if
                MatAinverser(nbDDLii, nbDDLjj) = sum(dfdF(ii, :)*raidTang(:5)*dfdF(jj, :))- &
                                                 sum(dfdQCP(ii-4, :)*HHCP(:)* &
                                                     signeHHCP3(:)*dfdF(jj, :))
            else
                MatAinverser(nbDDLii, nbDDLjj) = sum(dfdF(ii, :)*raidTang(:5)*dfdF(jj, :))
            end if
        end do
    end do

    !
    ! on inverse la matrice entièrement
    call mgauss('NCVP', MatAinverser, fsInverse, nbDDL, nbDDL, nbDDL, det, iret)
    ! vérification d'inversion
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
    !
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

! fonction rapide de signe ( attention ici est égal à 0 pour a=0 contrairement à sign de fortran)
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
