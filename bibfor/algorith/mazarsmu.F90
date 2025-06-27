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

subroutine mazarsmu(option, epsela, deps, dimloc, mazars, &
                    varm, varp, sigp, dsidep)
!
! --------------------------------------------------------------------------------------------------
!
!           Loi de Mazars unilatérale :
!               en 1D
!               en contraintes planes
!               en 3D
!
! --------------------------------------------------------------------------------------------------
!
    use tenseur_meca_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/sgmxve.h"
!
    integer(kind=8)             :: dimloc
    character(len=16)   :: option
    real(kind=8)        :: epsela(6), deps(6)
    real(kind=8)        :: varm(*), varp(*), sigp(*), dsidep(6, 6), mazars(*)
!
! --------------------------------------------------------------------------------------------------
!
!  IN :
!     option    : FULL_MECA RAPH_MECA RIGI_MECA_TANG
!     epsela    : déformation totale
!     deps      : incrément de déformation
!     dimloc    : dimension du problème : 6=3D, 4=CP, 1=1D
!     mazars    : les coefficients de la loi
!     varm      : variables internes à l'instant moins
!
!  OUT :
!     varp      : variables internes a l'instant plus
!     sigp      : contrainte a l'instant plus
!     dsidep    : matrice tangente
!
! --------------------------------------------------------------------------------------------------
!
!     Variables internes
!        1  -> icels  : critère sigma
!        2  -> icelu  : critère epsi
!        3  -> idomm  : endommagement
!        4  -> iepsqt : valeur de epseqt de traction
!        5  -> iepsqc : valeur de epseqc de compression
!        6  -> irsigm : facteur de triaxialité en contrainte
!        7  -> itemp  : température maximale atteinte par le matériau
!        8  -> idissd : dissipation d'endommagement
!
! --------------------------------------------------------------------------------------------------
!   Index des variables internes
!     icels  = 1, icelu  = 2, idomm = 3, iepsqt = 4
!     iepsqc = 5, irsigm = 6, itemp = 7, idissd = 8
    integer(kind=8), parameter :: icels = 1, icelu = 2, idomm = 3, iepsqt = 4
    integer(kind=8), parameter :: iepsqc = 5, irsigm = 6, idissd = 8
! --------------------------------------------------------------------------------------------------
!   Index des coefficients de la loi
!     iepsd0   = 1, ik       = 2, iac    = 3, ibc    = 4, iat = 5, ibt = 6
!     isigmlim = 7, iepsilim = 8, iepsc0 = 9, iepst0 = 10
!     iyoung   = 11, inu     = 12
    integer(kind=8), parameter :: ik = 2, iac = 3, ibc = 4, iat = 5, ibt = 6
    integer(kind=8), parameter :: isigmlim = 7, iepsilim = 8, iepsc0 = 9, iepst0 = 10
    integer(kind=8), parameter :: iyoung = 11, inu = 12
!
    aster_logical :: rigi, resi
!
    real(kind=8), parameter :: grdexp = 200.0
    real(kind=8), parameter :: domlim = 0.99999
!
    real(kind=8) :: epst0, kk, ac, bc, at, bt, ee, nu, sgels, epelu
    real(kind=8) :: aa, bb, rr, epsc0, ddomdy1d
!
    real(kind=8) :: sigeq, Epsic, Epsit, Iepsi, Jepsi, sigseuil
    real(kind=8) :: epspri(3), sigpri(3), epseq, epsloc(6), sigploc(6), gr, Y0, Yt, Yc
    real(kind=8) :: trsiga, trsigt, lambda, deuxmu, YYplus, YYmoins, DOMplus, DOMmoins
    real(kind=8) :: rtemp
!
    real(kind=8), parameter :: NormEpsi = 1.0E-06
    type(tenseur2)   :: TEpsi, TSigma
    type(basepropre) :: TEpsiPrin
!
! --------------------------------------------------------------------------------------------------
!
!   RIGI_MECA_TANG ->        DSIDEP        -->  RIGI
!   FULL_MECA      ->  SIGP  DSIDEP  VARP  -->  RIGI  RESI
!   RAPH_MECA      ->  SIGP          VARP  -->        RESI
!
!   Gestion à faire après l'appel
!       rigi : est toujours vrai. Le module tangent est toujours retourné
!       sigp : est toujours calculée.
!   rigi = (option(1:4).eq.'RIGI' .or. option(1:4).eq.'FULL')
    resi = (option(1:4) .eq. 'RAPH' .or. option(1:4) .eq. 'FULL')
    rigi = .true.
!
! --------------------------------------------------------------------------------------------------
!   Caractéristiques matériaux
    epst0 = mazars(iepst0)
    at = mazars(iat)
    bt = mazars(ibt)
    epsc0 = mazars(iepsc0)
    ac = mazars(iac)
    bc = mazars(ibc)
    sgels = mazars(isigmlim)
    epelu = mazars(iepsilim)
    kk = mazars(ik)
    ee = mazars(iyoung)
    nu = mazars(inu)
!
    lambda = ee*nu/(1.0+nu)/(1.0-2.0*nu)
    deuxmu = ee/(1.0+nu)
!
! --------------------------------------------------------------------------------------------------
!   Déformation actuelle
    epsloc(1:6) = 0.0
    epsloc(1:dimloc) = epsela(1:dimloc)+deps(1:dimloc)
    if (dimloc == 4) then
        ! On est en contraintes planes.
        epsloc(3) = -nu*(epsloc(1)+epsloc(2))/(1.0-nu)
        ! Le cisaillement est mis comme il faut, pour utiliser le 'vrai' calcul tensoriel
        epsloc(4) = epsloc(4)/sqrt(2.d0)
    else if (dimloc == 6) then
        !  On est en 3D
        ! Le cisaillement est mis comme il faut, pour utiliser le 'vrai' calcul tensoriel
        epsloc(4) = epsloc(4)/sqrt(2.d0)
        epsloc(5) = epsloc(5)/sqrt(2.d0)
        epsloc(6) = epsloc(6)/sqrt(2.d0)
    else if (dimloc == 1) then
        ! On est en 1D
        epsloc(2) = -nu*epsloc(1)
        epsloc(3) = -nu*epsloc(1)
        ! Le cisaillement n'existe pas
    end if
! --------------------------------------------------------------------------------------------------
!   Repère principal de EPS
    if (dimloc == 1) then
        epspri = [epsloc(1), epsloc(2), epsloc(3)]
    else
        TEpsi = epsloc
        ! TEpsiPrin : tenseur propre, base propre, ...
        TEpsiPrin = basevecteurpropre(TEpsi, NormEpsi)
        ! Valeurs propres
        epspri = [TEpsiPrin%valep%vect(1), TEpsiPrin%valep%vect(2), TEpsiPrin%valep%vect(3)]
    end if
    Iepsi = epspri(1)+epspri(2)+epspri(3)
! --------------------------------------------------------------------------------------------------
!   Calcul des contraintes élastiques, dans le repère principal de déformation
    if (dimloc == 1) then
        sigpri(:) = 0.0
        sigpri(1) = ee*epspri(1)
    else
        sigpri(:) = deuxmu*epspri(:)+lambda*Iepsi
    end if
! --------------------------------------------------------------------------------------------------
!   Calcul des critères en contraintes
    trsigt = max(0.0, sigpri(1))+max(0.0, sigpri(2))+max(0.0, sigpri(3))
    trsiga = abs(sigpri(1))+abs(sigpri(2))+abs(sigpri(3))
! --------------------------------------------------------------------------------------------------
!   Domaine où l'endommagement ne progresse pas (autour de 0)
!           En 1D     : dépend uniquement de Exx (1 c'est la limite)
!           En 2D, 3D : dépend des déformations principales
!           ==> Choix de mettre 1.0E-06 sur la contrainte PIC (la plus petite entre Trac et Comp)
    sigseuil = abs(ee*min(epst0, epsc0)*1.0E-06)
!   si trsiga < sigseuil : La somme des contraintes principales est donc petite
!       Lorsque l'on fait :
!           Sigp = sigpri*(1.0-DOMplus) Toutes les composantes du Sigp sont petites
!           Quelle que soit la valeur de "rr dans [0, 1]" cela n'a pas d'influence
!           Remplacer (trsiga > r8prem()) par (trsiga > sigseuil)
!
!   Calcul de R : 0 en compression pure, 1 en traction pure
    rr = 0.0
    if (trsiga > sigseuil) then
        rr = trsigt/trsiga
    else if (trsigt > 0.0) then
        rr = 1.0
    end if
    rr = min(max(0.0, rr), 1.0)
! --------------------------------------------------------------------------------------------------
!   Endommagement précédent dans la direction du chargement
!       Récupération des : epsqt, epsqc
    DOMmoins = 0.0
    Yt = max(epst0, varm(iepsqt))
    Yc = max(epsc0, varm(iepsqc))
!
    YYmoins = rr*Yt+(1.0-rr)*Yc
    Y0 = rr*epst0+(1.0-rr)*epsc0
    if (YYmoins > Y0) then
        ! calcul de l'endommagement
        aa = 2.0*(rr*rr)*(at-2.0*kk*at+ac)-rr*(at-4.0*kk*at+3.0*ac)+ac
        gr = rr**(rr**2-2.0*rr+2.0)
        bb = gr*bt+(1.0-gr)*bc
        ! Éviter que le calcul plante dans l'évaluation de exp(rtemp) si rtemp est trop grand
        rtemp = bb*(YYmoins-Y0)
        DOMmoins = 1.0-Y0*(1.0-aa)/YYmoins
        if (rtemp .le. grdexp) DOMmoins = DOMmoins-aa*exp(-rtemp)
        DOMmoins = min(max(DOMmoins, 0.0), domlim)
    end if
    DOMplus = DOMmoins
! --------------------------------------------------------------------------------------------------
    ddomdy1d = 0.0
!   Calcul des critères en déformation
    Jepsi = (epspri(1)-epspri(2))**2+(epspri(1)-epspri(3))**2+(epspri(2)-epspri(3))**2
    Jepsi = sqrt(Jepsi*0.5)
!   Fonctions seuils
    Epsit = (Iepsi/(1.0-2.0*nu)+Jepsi/(1.0+nu))/2.0
    Epsic = (Iepsi/(1.0-2.0*nu)+6.0*Jepsi/(1.0+nu))/5.0
!   Max sur tout le chargement
    Epsit = max(varm(iepsqt), Epsit)
    Epsic = max(varm(iepsqc), Epsic)
!   Sup par rapport au seuil d'endommagement initial
    Yt = max(epst0, Epsit)
    Yc = max(epsc0, Epsic)
!   Le YY actuel dans la direction de chargement
    YYplus = rr*Yt+(1.0-rr)*Yc
    Y0 = rr*epst0+(1.0-rr)*epsc0
! --------------------------------------------------------------------------------------------------
    if (YYplus > Y0) then
        ! calcul de l'endommagement
        aa = 2.0*(rr*rr)*(at-2.0*kk*at+ac)-rr*(at-4.0*kk*at+3.0*ac)+ac
        gr = rr**(rr**2-2.0*rr+2.0)
        bb = gr*bt+(1.0-gr)*bc
        ! Éviter que le calcul plante dans l'évaluation de exp(rtemp) si rtemp est trop grand
        rtemp = bb*(YYplus-Y0)
        DOMplus = 1.0-Y0*(1.0-aa)/YYplus
        ddomdy1d = Y0*(1.0-aa)/YYplus**2
        if (rtemp .le. grdexp) then
            DOMplus = DOMplus-aa*exp(-rtemp)
            ddomdy1d = ddomdy1d+aa*bb*exp(-rtemp)
        end if
        if (DOMplus .le. DOMmoins) then
            DOMplus = DOMmoins
            ddomdy1d = 0.0
        else if (DOMplus .ge. domlim) then
            DOMplus = domlim
            ddomdy1d = 0.0
        end if
    end if
!
!   Calcul des contraintes dans le repère initial
!       Passage base propre vers base initiale
    if (dimloc == 1) then
        sigp(1) = sigpri(1)*(1.0-DOMplus)
    else
        TSigma = [sigpri(1), sigpri(2), sigpri(3)]
        TSigma = VersBaseInitiale(TEpsiPrin, TSigma)
        sigploc = TSigma
        sigp(1:dimloc) = sigploc(1:dimloc)*(1.0-DOMplus)
    end if
!   Mise à jour des variables internes
    if (resi) then
        ! Correspond aux critères ELS, ELU dans le cas non-linéaire 1D
        epseq = sqrt(epspri(1)**2+epspri(2)**2+epspri(3)**2)
        sigeq = sqrt(sigpri(1)**2+sigpri(2)**2+sigpri(3)**2)
        varp(icels) = sigeq*sgmxve(3, sigpri)*(1.0-DOMplus)/sgels
        varp(icelu) = epseq*sgmxve(3, epspri)/epelu
        !
        varp(idomm) = DOMplus
        ! rr : [ 0 , ... , 1 ] <==> [ comp , cisa , trac ]
        if (rr .gt. 0.0) varp(iepsqt) = Epsit
        if (rr .lt. 1.0) varp(iepsqc) = Epsic
        varp(irsigm) = rr
        varp(idissd) = 0.0
    end if
!
! --------------------------------------------------------------------------------------------------
!   Si rigi : On retourne la matrice élastique endommagée
    if (rigi) then
        dsidep(:, :) = 0.0
        ! L'influence de l'endommagement est limité sur la tangente (peut-être plus d'itération)
        lambda = lambda*(1.0-DOMplus)
        deuxmu = deuxmu*(1.0-DOMplus)
        ! En contraintes planes ou pas, les termes 14 24 41 42 sont toujours nuls.
        if (dimloc == 4) then
            ! Les termes non nuls qui servent en CP : 11 12 21 22 44
            !     Dans dktnli on utilise            : 11 12 21 22 44 14 24 41 42
            dsidep(1, 1) = (1.0-DOMplus)*ee/(1.0-nu*nu)
            dsidep(1, 2) = nu*dsidep(1, 1)
            dsidep(2, 1) = dsidep(1, 2)
            dsidep(2, 2) = dsidep(1, 1)
            dsidep(4, 4) = deuxmu
        else if (dimloc == 6) then
            dsidep(1, 1) = deuxmu+lambda
            dsidep(1, 2) = lambda
            dsidep(1, 3) = lambda
            dsidep(2, 1) = lambda
            dsidep(2, 2) = deuxmu+lambda
            dsidep(2, 3) = lambda
            dsidep(3, 1) = lambda
            dsidep(3, 2) = lambda
            dsidep(3, 3) = deuxmu+lambda
            dsidep(4, 4) = deuxmu
            dsidep(5, 5) = deuxmu
            dsidep(6, 6) = deuxmu
        else if (dimloc == 1) then
            dsidep(1, 1) = ee*((1.0-DOMplus)-abs(epspri(1))*ddomdy1d)
        end if
    end if
end subroutine
