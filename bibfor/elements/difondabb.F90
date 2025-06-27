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
subroutine difondabb(for_discret, iret)
!
! --------------------------------------------------------------------------------------------------
!
! IN    for_discret : voir l'appel
! OUT   iret        : code retour
!
! --------------------------------------------------------------------------------------------------
!
    use te0047_type
    implicit none
!
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/diraidklv.h"
#include "asterfort/diklvraid.h"
#include "asterfort/infdis.h"
#include "asterfort/jevech.h"
#include "asterfort/difondmat.h"
#include "asterfort/difondmatpetit.h"
#include "asterfort/difoncalc.h"
#include "asterfort/difoncalcpetit.h"
#include "asterfort/tecael.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
#include "asterfort/utpslg.h"
#include "asterfort/utpvlg.h"
#include "blas/dcopy.h"
#include "asterfort/Behaviour_type.h"
!
    type(te0047_dscr), intent(in) :: for_discret
    integer(kind=8), intent(out) :: iret
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: imat, jdc, irep, neq, ii, ifono, icontp, icontm, icarcr
    integer(kind=8) :: iadzi, iazk24
!   paramètres des fonctions bizarres
    integer(kind=8) :: nbpar, kpg, spt, icodma, nbrePara, nbVarloc, nbdecp
!
    parameter(nbrePara=18, nbVarloc=21)
    integer(kind=8) :: codre1(nbrePara)
!
    real(kind=8) :: r8bid, klv(78), fl(12), valre1(nbrePara), dulMat(4)
    real(kind=8) :: valpar, valvarloc(nbVarloc), fzzz, mxxx, myyy, rNLx, rNLy
    real(kind=8) :: exxx, eyyy, tirela(6), raidTang(6), lxxx, lyyy, errmax, error, HCP
!
    character(len=8) :: k8bid, nompar, fami, poum
    character(len=24) :: messak(6)
    character(len=16) :: nomre1(nbrePara)
    character(len=16), pointer :: compor(:) => null()
!
!   calculPetit : pour savoir si on lance le calcul de linéarisation
!   PetitH      : basé sur 0.01Fz ou 0.2 Fz
    logical :: calculNormal, calculPetit, calculPetitH
!   increment : pour savoir si il y a un incrément de chargement dans le modèle
    logical :: increment
    blas_int :: b_incx, b_incy, b_n
!                 1234567890123456   1234567890123456   1234567890123456   1234567890123456
    data nomre1/'LONG_X          ', 'LONG_Y          ', 'PHI             ', 'COHESION        ', &
        'RAID_GLIS       ', 'GAMMA_REFE      ', 'CP_SERVICE      ', 'CP_ULTIME       ', &
        'DEPL_REFE       ', 'RAID_CP_X       ', 'GAMMA_CP_X      ', 'RAID_CP_Y       ', &
        'GAMMA_CP_Y      ', 'RAID_CP_RX      ', 'GAMMA_CP_RX     ', 'RAID_CP_RY      ', &
        'GAMMA_CP_RY     ', 'DECOLLEMENT     '/
!
! --------------------------------------------------------------------------------------------------
!
    iret = 0
! récupération des paramètres de comportement du calcul : pointeur icompo
    call jevech('PCOMPOR', 'L', vk16=compor)
! vérification du bon matériau, on regarde dans irep l'indication de la matrice
    call infdis('REPK', irep, r8bid, k8bid)
! Seulement en 3D TR et en poi1
    if ((for_discret%nomte(1:11) .ne. 'MECA_DIS_TR') .or. (for_discret%nno .ne. 1)) then
        messak(1) = for_discret%nomte
        messak(2) = for_discret%option
        messak(3) = compor(INCRELAS)
        messak(4) = compor(RELA_NAME)
        call tecael(iadzi, iazk24)
        messak(5) = zk24(iazk24-1+3)
        call utmess('F', 'DISCRETS_23', nk=5, valk=messak)
    end if
! Seulement en repère local : irep = 2
    if (irep .ne. 2) then
        messak(1) = for_discret%nomte
        messak(2) = for_discret%option
        messak(3) = compor(INCRELAS)
        messak(4) = compor(RELA_NAME)
        call tecael(iadzi, iazk24)
        messak(5) = zk24(iazk24-1+3)
        call utmess('F', 'DISCRETS_5', nk=5, valk=messak)
    end if
! récupération de la matrice de rigidité avec le pointeur jdc
    call jevech('PCADISK', 'L', jdc)
! récupération de l'emplacement des variables internes matériau avec le pointeur icontm
    call jevech('PVARIMR', 'L', icontm)
! stockage en local de ces variables internes
    do ii = 1, nbVarloc
        valvarloc(ii) = zr(icontm-1+ii)
    end do
! récupération de l'emplacement des paramètres matériaux avec le même pointeur icontm
    call jevech('PMATERC', 'L', icontm)
! définition des paramètres de rcvalb
    valre1(:) = 0.d0; valpar = 0.d0; nbpar = 0; nompar = ' '
    fami = 'FPG1'; kpg = 1; spt = 1; poum = '+'
    icodma = zi(icontm)
! récupération des paramètres matériaux dans valre1, comme tous les paramètres ont une valeur
! obligatoire ou, facultatif par défaut, pas de vérification supplémentaire
    call rcvalb(fami, kpg, spt, poum, icodma, &
                ' ', 'FONDA_SUPERFI', nbpar, nompar, [valpar], &
                nbrePara, nomre1, valre1, codre1, 1)
! vérification que Lx et Ly soit strictement positif (ou supérieur à r8prem()),
! Vélas peut être nul : plastification directe
    lxxx = valre1(1)
    lyyy = valre1(2)
    if ((lxxx .le. r8prem()) .or. (lyyy .le. r8prem())) then
        messak(1) = for_discret%nomte
        messak(2) = for_discret%option
        messak(3) = compor(INCRELAS)
        messak(4) = compor(RELA_NAME)
        call tecael(iadzi, iazk24)
        messak(5) = zk24(iazk24-1+3)
        call utmess('F', 'DISCRETS_44', nk=5, valk=messak)
    end if
! modification de Vult en DVult=Vult-Vélas positif avec vérification
! On ne peut pas avoir Vult=0.0 et Vélas=0
    if ((valre1(7) .le. r8prem()) .and. (valre1(8) .le. r8prem())) then
        messak(1) = for_discret%nomte
        messak(2) = for_discret%option
        messak(3) = compor(INCRELAS)
        messak(4) = compor(RELA_NAME)
        call tecael(iadzi, iazk24)
        messak(5) = zk24(iazk24-1+3)
        call utmess('F', 'DISCRETS_45', nk=5, valk=messak)
    end if
    valre1(8) = max(valre1(8)-valre1(7), 0.0)
! le repère est forcement local (voir plus haut) on copie donc directement
! la matrice (symmétrique) dans klv
    b_n = to_blas_int(for_discret%nbt)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, zr(jdc), b_incx, klv, b_incy)
! Récupère les termes diagonaux de la matrice de raideur dans raidTang
    call diraidklv(for_discret%nomte, raidTang, klv)
    !
! récupération de la force à l'instant précédent avec pointeur icontm dans le repère local
    call jevech('PCONTMR', 'L', icontm)
! récupération des paramètres d'intégration
    call jevech('PCARCRI', 'L', icarcr)
! nombre d'itérations maxi (ITER_INTE_MAXI) qui doit être renseigné avec la commande de
!      subdivision
    nbdecp = int(zr(icarcr))
! tolérance de convergence (RESI_INTE)
    errmax = zr(icarcr+2)
! erreur dans la même dimension des critère de palsticité (multiplié par Vélas))
    error = errmax*valre1(7)
    !
! on regarde si DECOL est renseigné, si oui alors on calcul le décollement si non
! les facteur multi sont égaux à 1
! si = 1 alors décollement  =2 sinon
    if (valre1(18) .lt. 1.50) then
        fzzz = zr(icontm-1+3)
        mxxx = zr(icontm-1+4)
        myyy = zr(icontm-1+5)
! on calcule la modification de raideur élastique due au décollement/excentrement
! si fzzz <=0.0 alors on a un problème, on suppose le cas Fz=0.0 M=0.0
! pas de modification de la raideur,
! si un moment est présent alors la raideur en rotation est nulle pour éviter
        if (fzzz .GE. (-r8prem())) then
            if (abs(mxxx) .le. error*Lyyy) then
                rNLx = 1.0
            else
                rNLx = 0.0
            end if
            if (abs(myyy) .le. error*Lxxx) then
                rNLy = 1.00
            else
                rNLy = 0.0
            end if
        else
            eyyy = abs(mxxx/fzzz)
            exxx = abs(myyy/fzzz)
            if (eyyy .lt. lyyy/6.0) then
                rNLx = 1.0
            else
                rNLx = 27.0/8.0*((1.0-2.0*eyyy/lyyy)**3.0)
            end if
            if (exxx .lt. lxxx/6.0) then
                rNLy = 1.0
            else
                rNLy = 27.0/8.0*((1.0-2.0*exxx/lxxx)**3.0)
            end if
        end if
    else
! cas sans décollement
        rNLx = 1.0
        rNLy = 1.0
    end if
! modification des raideurs diagonales en rotations
    raidTang(4) = raidTang(4)*rNLx
    raidTang(5) = raidTang(5)*rNLy
! recréation de klv diagonal
    call diklvraid(for_discret%nomte, klv, raidTang)
! stockage de la matrice tangente dans klv
    neq = for_discret%nno*for_discret%nc
! au lieu de reformer la matrice totale non symétrique (comme la dimension est fixée neq=6)
! on calcul l'incrément de force 'à la main' et on regarde si il y a un increment
    increment = .FALSE.
    do ii = 1, neq
! calcul de l'incrément élastique
        fl(ii) = raidTang(ii)*for_discret%dul(ii)
! calcul de la force élastique
        tirela(ii) = fl(ii)+zr(icontm-1+ii)
! increment  = increment .OR. (tirela(ii) .NE. zr(icontm-1+ii))
        increment = increment .OR. (abs(fl(ii)) .ge. r8prem())
    end do
! force verticale trop petite on linéarise le critère ce qui force un traitement particulier
! si calculNormal reste vrai alors on doit calculer le critère non linéarisé
    calculNormal = .TRUE.
! deux cas de linéarisation :
! si on dépasse 1 % du Vélas -> linéarisation automatique
    if (tirela(3) .GE. (-0.01*(valre1(7)+valvarloc(15)))) then
! indique (si vrai) qu'on linéarise à 0.01 Vélas
        calculPetit = .TRUE.
        if (tirela(3) .lt. -r8prem()) then
            calculPetitH = .NOT. ( &
                           ( &
                           abs( &
                           tirela(4)-valvarloc(16))/(-lyyy/2.0*tirela(3)) .GT. 0.8) .OR. (abs(ti&
                           &rela(5)-valvarloc(17))/(-lxxx/2.0*tirela(3) &
                           ) .GT. 0.8 &
                           ) &
                           )
! si les moments sont trop important on utilise la même linéarisation que le cas 10%
!   du Vélas, n'a plus aucun sens si la force verticale est en traction'
        else
            calculPetitH = .TRUE.
        end if
    else if (tirela(3) .GE. (-0.1*(valre1(7)+valvarloc(15)))) then
! sinon on regarde si on dépasse 10 % du Vélas -> linéarisation si moment fort
        calculPetit = ( &
                      abs( &
                      tirela(4)-valvarloc(16))/(-lyyy/2.0*tirela(3)) .GT. 0.8) .OR. (abs(tirela(5&
                      &)-valvarloc(17))/(-lxxx/2.0*tirela(3) &
                      ) .GT. 0.8 &
                      )
        calculPetitH = .FALSE.
    else
! pas de calcul de linéarisation
        calculPetit = .FALSE.
    end if
! on effectue le calcul de l'incrément de chargement avec critère linéarisé si les conditions
! sont remplis (Vélas suffisamment petit)'
    if (calculPetit) then
        if (increment) then
            call difoncalcpetit(tirela, raidTang, valvarloc, valre1, nbVarloc, &
                                nbrePara, iret, nbdecp, errmax, calculNormal, &
                                calculPetitH)
! si il n'y a pas d'incrément on passe direct à la suite en mettant calculNormal .FALSE.
! on évite ainsi le cas 0, 0, 0 calculé dans le calculateur normal
        else
            calculNormal = .FALSE.
        end if
    end if
! on vérifie la convergence
    if (iret .NE. 0) then
        goto 999
    end if
! si besoin de faire le calcul d'incrément de chargement sans linéarisation du critère
    if (calculNormal) then
! on lance le calcul de la force corrigée
        HCP = ((tirela(1)-valvarloc(13))**2.0+(tirela(2)-valvarloc(14))**2.0)**0.5
! on regarde si les moments ou les forces horizontales ne sont pas trop importants
        if ((HCP .lt. (-tirela(3))) .and. &
            (abs(tirela(4)-valvarloc(16)) .lt. (-lyyy/2.0*tirela(3))) .and. &
            (abs(tirela(5)-valvarloc(17)) .lt. (-lxxx/2.0*tirela(3)))) then
            call difoncalc(tirela, raidTang, valvarloc, valre1, nbVarloc, &
                           nbrePara, iret, nbdecp, errmax)
        else
            iret = 1
        end if
    end if
! on vérifie la convergence
    if (iret .NE. 0) then
        goto 999
    end if
! on lance le calcul de la matrice tangente corrigée
! calcul de dulmat qui stocke le déplacement du matériau
    dulMat(1) = for_discret%dul(1)
    dulMat(2) = for_discret%dul(2)
    dulMat(3) = for_discret%dul(4)
    dulMat(4) = for_discret%dul(5)
    if (valvarloc(18) .GT. 9999.0) then
! sur la variable locale 18 (dautre) on stocke l'état du macroélément, si >=10000
!   alors on est dans le cas de la linéarisation'
        call difondmatpetit(tirela, raidTang, valvarloc, valre1, nbVarloc, &
                            nbrePara, klv, errmax, dulMat, iret)
    else
        call difondmat(tirela, raidTang, valvarloc, valre1, nbVarloc, &
                       nbrePara, klv, errmax, dulMat, iret)
    end if
    !
! Sortie : Matrice tangente
    if (for_discret%lMatr) then
! on récupère le pointeur d'écriture de la matrice de raideur tangente
! (dans le repère globale) dans imat
        call jevech('PMATUUR', 'E', imat)
! on est forcement en 3D, utpslg retourne sur le pointeur imat,
! la matrice tangente préalablement changée de repère par pgl
        call utpslg(for_discret%nno, for_discret%nc, for_discret%pgl, klv, zr(imat))
    end if
! Sortie : Contrainte
    if (for_discret%lSigm) then
! récupération du pointeur d'écriture des forces internes à cet incrément (icontp)
        call jevech('PCONTPR', 'E', icontp)
! écriture des forces internes en sortie à cet incrément
        do ii = 1, neq
            zr(icontp-1+ii) = tirela(ii)
        end do
    end if
! Sortie : Efforts
    if (for_discret%lVect) then
!  récupération du pointeur d'écriture de la force globale à cet incrément (ifono)
        call jevech('PVECTUR', 'E', ifono)
! on est forcement en 3D, utpvlg retourne sur le pointeur ifono,
! la force globale préalablement changée de repère par pgl
        call utpvlg(for_discret%nno, for_discret%nc, for_discret%pgl, tirela, zr(ifono))
    end if
! Sortie : Paramètres internes
    if (for_discret%lVari) then
! récupération du pointeur d'écriture des variables internes à cet incrément (icontm)
        call jevech('PVARIPR', 'E', icontm)
! mise à jour des variables internes
        do ii = 1, nbVarloc
            zr(icontm+ii-1) = valvarloc(ii)
        end do
    end if
999 continue
end subroutine
