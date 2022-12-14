! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
!
type(te0047_dscr), intent(in) :: for_discret
integer, intent(out)          :: iret
!
! --------------------------------------------------------------------------------------------------
!
    integer :: imat, jdc, irep, neq, ii, ifono, icontp, icontm,icompo,icarcr
    integer :: iadzi, iazk24
!   param??tres des fonctions bizarres
    integer :: nbpar,kpg,spt,icodma,nbrePara,nbVarloc,nbdecp
!
    parameter (nbrePara = 18,nbVarloc = 21)
    integer      :: codre1(nbrePara)
!
    real(kind=8) :: r8bid, klv(78), fl(12),valre1(nbrePara),dulMat(4)
    real(kind=8) :: valpar,valvarloc(nbVarloc),fzzz,mxxx,myyy,rNLx,rNLy
    real(kind=8) :: exxx,eyyy,tirela(6),raidTang(6),lxxx,lyyy,errmax,error,HCP
!
    character(len=8)    :: k8bid, nompar, fami, poum
    character(len=24)   :: messak(6)
    character(len=16)   :: nomre1(nbrePara)
!
!   calculPetit : pour savoir si on lance le calcul de lin??arisation
!   PetitH      : bas?? sur 0.01Fz ou 0.2 Fz
    logical :: calculNormal,calculPetit,calculPetitH
!   increment : pour savoir si il y a un incr??ment de chargement dans le mod??le
    logical :: increment
!                 1234567890123456   1234567890123456   1234567890123456   1234567890123456
    data nomre1 /'LONG_X          ','LONG_Y          ','PHI             ','COHESION        ', &
                 'RAID_GLIS       ','GAMMA_REFE      ','CP_SERVICE      ','CP_ULTIME       ', &
                 'DEPL_REFE       ','RAID_CP_X       ','GAMMA_CP_X      ','RAID_CP_Y       ', &
                 'GAMMA_CP_Y      ','RAID_CP_RX      ','GAMMA_CP_RX     ','RAID_CP_RY      ', &
                 'GAMMA_CP_RY     ','DECOLLEMENT     '/
!
! --------------------------------------------------------------------------------------------------
!
    iret = 0
    ! r??cup??ration des param??tres de comportement du calcul : pointeur icompo
    call jevech('PCOMPOR', 'L', icompo)
    ! v??rification du bon mat??riau, on regarde dans irep l'indication de la matrice
    call infdis('REPK', irep, r8bid, k8bid)
    ! Seulement en 3D TR et en poi1
    if ((for_discret%nomte(1:11).ne.'MECA_DIS_TR').or.(for_discret%nno.ne.1) ) then
        messak(1) = for_discret%nomte
        messak(2) = for_discret%option
        messak(3) = zk16(icompo+3)
        messak(4) = zk16(icompo)
        call tecael(iadzi, iazk24)
        messak(5) = zk24(iazk24-1+3)
        call utmess('F', 'DISCRETS_23', nk=5, valk=messak)
    endif
    ! Seulement en rep??re local : irep = 2
    if (irep .ne. 2) then
        messak(1) = for_discret%nomte
        messak(2) = for_discret%option
        messak(3) = zk16(icompo+3)
        messak(4) = zk16(icompo)
        call tecael(iadzi, iazk24)
        messak(5) = zk24(iazk24-1+3)
        call utmess('F', 'DISCRETS_5', nk=5, valk=messak)
    endif
    ! r??cup??ration de la matrice de rigidit?? avec le pointeur jdc
    call jevech('PCADISK', 'L', jdc)
    ! r??cup??ration de l'emplacement des variables internes mat??riau avec le pointeur icontm
    call jevech('PVARIMR', 'L', icontm)
    ! stockage en local de ces variables internes
    do ii=1, nbVarloc
        valvarloc(ii)=zr(icontm-1+ii)
    enddo
    ! r??cup??ration de l'emplacement des param??tres mat??riaux avec le m??me pointeur icontm
    call jevech('PMATERC', 'L', icontm)
    ! d??finition des param??tres de rcvalb
    valre1(:) = 0.d0 ; valpar = 0.d0; nbpar = 0; nompar = ' '
    fami = 'FPG1'; kpg = 1; spt = 1; poum = '+'
    icodma=zi(icontm)
    ! r??cup??ration des param??tres mat??riaux dans valre1, comme tous les param??tres ont une valeur
    ! obligatoire ou, facultatif par d??faut, pas de v??rification suppl??mentaire
    call rcvalb(fami, kpg, spt, poum, icodma, ' ', 'FONDA_SUPERFI', nbpar, nompar, [valpar],&
                nbrePara, nomre1, valre1, codre1, 1)
    ! v??rification que Lx et Ly soit strictement positif (ou sup??rieur ?? r8prem()),
    ! V??las peut ??tre nul : plastification directe
    lxxx = valre1(1)
    lyyy = valre1(2)
    if ((lxxx.le.r8prem()) .or. (lyyy.le.r8prem()) ) then
        messak(1) = for_discret%nomte
        messak(2) = for_discret%option
        messak(3) = zk16(icompo+3)
        messak(4) = zk16(icompo)
        call tecael(iadzi, iazk24)
        messak(5) = zk24(iazk24-1+3)
        call utmess('F', 'DISCRETS_44', nk=5, valk=messak)
    endif
    ! modification de Vult en DVult=Vult-V??las positif avec v??rification
    ! On ne peut pas avoir Vult=0.0 et V??las=0
    if ((valre1(7).le.r8prem()) .and. (valre1(8).le.r8prem()) ) then
        messak(1) = for_discret%nomte
        messak(2) = for_discret%option
        messak(3) = zk16(icompo+3)
        messak(4) = zk16(icompo)
        call tecael(iadzi, iazk24)
        messak(5) = zk24(iazk24-1+3)
        call utmess('F', 'DISCRETS_45', nk=5, valk=messak)
    endif
    valre1(8)=max(valre1(8)-valre1(7),0.0)
    ! le rep??re est forcement local (voir plus haut) on copie donc directement
    ! la matrice (symm??trique) dans klv
    call dcopy(for_discret%nbt, zr(jdc), 1, klv, 1)
    ! R??cup??re les termes diagonaux de la matrice de raideur dans raidTang
    call diraidklv(for_discret%nomte,raidTang,klv)
    !
    ! r??cup??ration de la force ?? l'instant pr??c??dent avec pointeur icontm dans le rep??re local
    call jevech('PCONTMR', 'L', icontm)
    ! r??cup??ration des param??tres d'int??gration
    call jevech('PCARCRI', 'L', icarcr)
    ! nombre d'it??rations maxi (ITER_INTE_MAXI) qui doit ??tre renseign?? avec la commande de
    !      subdivision
    nbdecp = int(zr(icarcr))
    ! tol??rance de convergence (RESI_INTE_RELA)
    errmax = zr(icarcr+2)
    ! erreur dans la m??me dimension des crit??re de palsticit?? (multipli?? par V??las))
    error=errmax*valre1(7)
    !
    ! on regarde si DECOL est renseign??, si oui alors on calcul le d??collement si non
    ! les facteur multi sont ??gaux ?? 1
    ! si = 1 alors d??collement  =2 sinon
    if (valre1(18) .lt. 1.50) then
        fzzz = zr(icontm-1+3)
        mxxx = zr(icontm-1+4)
        myyy = zr(icontm-1+5)
        ! on calcule la modification de raideur ??lastique due au d??collement/excentrement
        ! si fzzz <=0.0 alors on a un probl??me, on suppose le cas Fz=0.0 M=0.0
        ! pas de modification de la raideur,
        ! si un moment est pr??sent alors la raideur en rotation est nulle pour ??viter
        if (fzzz .GE. (-r8prem())) then
            if (abs(mxxx).le.error*Lyyy) then
                rNLx=1.0
            else
                rNLx=0.0
            endif
            if (abs(myyy).le.error*Lxxx) then
                rNLy=1.00
            else
                rNLy=0.0
            endif
        else
            eyyy=abs(mxxx/fzzz)
            exxx=abs(myyy/fzzz)
            if (eyyy .lt. lyyy/6.0) then
                rNLx=1.0
            else
                rNLx=27.0/8.0*((1.0-2.0*eyyy/lyyy)**3.0)
            endif
            if (exxx .lt. lxxx/6.0) then
                rNLy=1.0
            else
                rNLy=27.0/8.0*((1.0-2.0*exxx/lxxx)**3.0)
            endif
        endif
    else
        ! cas sans d??collement
        rNLx = 1.0
        rNLy = 1.0
    endif
    ! modification des raideurs diagonales en rotations
    raidTang(4)=raidTang(4)*rNLx
    raidTang(5)=raidTang(5)*rNLy
    ! recr??ation de klv diagonal
    call diklvraid(for_discret%nomte, klv, raidTang)
    ! stockage de la matrice tangente dans klv
    neq = for_discret%nno*for_discret%nc
    ! au lieu de reformer la matrice totale non sym??trique (comme la dimension est fix??e neq=6)
    ! on calcul l'incr??ment de force '?? la main' et on regarde si il y a un increment
    increment=.FALSE.
    do ii=1,neq
        ! calcul de l'incr??ment ??lastique
        fl(ii)     = raidTang(ii)*for_discret%dul(ii)
        ! calcul de la force ??lastique
        tirela(ii) = fl(ii) + zr(icontm-1+ii)
        ! increment  = increment .OR. (tirela(ii) .NE. zr(icontm-1+ii))
        increment  = increment .OR. (abs(fl(ii)) .ge. r8prem() )
    enddo
    ! force verticale trop petite on lin??arise le crit??re ce qui force un traitement particulier
    ! si calculNormal reste vrai alors on doit calculer le crit??re non lin??aris??
    calculNormal=.TRUE.
    ! deux cas de lin??arisation :
    ! si on d??passe 1 % du V??las -> lin??arisation automatique
    if (tirela(3) .GE. (-0.01*(valre1(7)+valvarloc(15)))) then
        ! indique (si vrai) qu'on lin??arise ?? 0.01 V??las
        calculPetit=.TRUE.
        if (tirela(3).lt.-r8prem()) then
            calculPetitH = .NOT.((abs(tirela(4)-valvarloc(16))/(-lyyy/2.0*tirela(3)).GT.0.8) &
                                  .OR.(abs(tirela(5)-valvarloc(17))/(-lxxx/2.0*tirela(3)).GT.0.8))
            ! si les moments sont trop important on utilise la m??me lin??arisation que le cas 10%
            !   du V??las, n'a plus aucun sens si la force verticale est en traction'
        else
            calculPetitH=.TRUE.
        endif
    else if (tirela(3) .GE. (-0.1*(valre1(7)+valvarloc(15)))) then
        ! sinon on regarde si on d??passe 10 % du V??las -> lin??arisation si moment fort
        calculPetit = (abs(tirela(4)-valvarloc(16))/(-lyyy/2.0*tirela(3)).GT.0.8) &
                       .OR.(abs(tirela(5)-valvarloc(17))/(-lxxx/2.0*tirela(3)).GT.0.8)
        calculPetitH=.FALSE.
    else
        ! pas de calcul de lin??arisation
        calculPetit=.FALSE.
    endif
    ! on effectue le calcul de l'incr??ment de chargement avec crit??re lin??aris?? si les conditions
    ! sont remplis (V??las suffisamment petit)'
    if (calculPetit) then
        if(increment)then
            call difoncalcpetit(tirela,raidTang,valvarloc,valre1,nbVarloc, &
                                nbrePara,iret,nbdecp,errmax,calculNormal,calculPetitH)
        ! si il n'y a pas d'incr??ment on passe direct ?? la suite en mettant calculNormal .FALSE.
        ! on ??vite ainsi le cas 0, 0, 0 calcul?? dans le calculateur normal
        else
            calculNormal=.FALSE.
        endif
    endif
    ! on v??rifie la convergence
    if (iret .NE. 0) then
        goto 999
    endif
    ! si besoin de faire le calcul d'incr??ment de chargement sans lin??arisation du crit??re
    if (calculNormal) then
        ! on lance le calcul de la force corrig??e
        HCP=((tirela(1)-valvarloc(13))**2.0+(tirela(2)-valvarloc(14))**2.0)**0.5
        ! on regarde si les moments ou les forces horizontales ne sont pas trop importants
        if ( (HCP.lt.(-tirela(3))).and.(abs(tirela(4)-valvarloc(16)).lt.(-lyyy/2.0*tirela(3))) &
              .and. (abs(tirela(5)-valvarloc(17)) .lt.(-lxxx/2.0*tirela(3))) ) then
            call difoncalc(tirela,raidTang,valvarloc,valre1,nbVarloc,nbrePara,iret,nbdecp,errmax)
        else
            iret=1
        endif
    endif
    ! on v??rifie la convergence
    if (iret .NE. 0) then
        goto 999
    endif
    ! on lance le calcul de la matrice tangente corrig??e
    ! calcul de dulmat qui stocke le d??placement du mat??riau
    dulMat(1)=for_discret%dul(1)
    dulMat(2)=for_discret%dul(2)
    dulMat(3)=for_discret%dul(4)
    dulMat(4)=for_discret%dul(5)
    if (valvarloc(18).GT. 9999.0) then
        ! sur la variable locale 18 (dautre) on stocke l'??tat du macro??l??ment, si >=10000
        !   alors on est dans le cas de la lin??arisation'
        call difondmatpetit(tirela,raidTang,valvarloc,valre1,nbVarloc,nbrePara, &
                            klv,errmax,dulMat,iret)
    else
        call difondmat(tirela,raidTang,valvarloc,valre1,nbVarloc,nbrePara,klv,errmax,dulMat,iret)
    endif
    !
    ! Sortie : Matrice tangente
    if ( for_discret%lMatr ) then
        ! on r??cup??re le pointeur d'??criture de la matrice de raideur tangente
        ! (dans le rep??re globale) dans imat
        call jevech('PMATUUR', 'E', imat)
        ! on est forcement en 3D, utpslg retourne sur le pointeur imat,
        ! la matrice tangente pr??alablement chang??e de rep??re par pgl
        call utpslg(for_discret%nno, for_discret%nc, for_discret%pgl, klv, zr(imat))
    endif
    ! Sortie : Contrainte
    if ( for_discret%lSigm ) then
        ! r??cup??ration du pointeur d'??criture des forces internes ?? cet incr??ment (icontp)
        call jevech('PCONTPR', 'E', icontp)
        ! ??criture des forces internes en sortie ?? cet incr??ment
        do ii=1,neq
            zr(icontp-1+ii) = tirela(ii)
        enddo
    endif
    ! Sortie : Efforts
    if ( for_discret%lVect ) then
        !  r??cup??ration du pointeur d'??criture de la force globale ?? cet incr??ment (ifono)
        call jevech('PVECTUR', 'E', ifono)
        ! on est forcement en 3D, utpvlg retourne sur le pointeur ifono,
        ! la force globale pr??alablement chang??e de rep??re par pgl
        call utpvlg(for_discret%nno, for_discret%nc, for_discret%pgl, tirela, zr(ifono))
    endif
    ! Sortie : Param??tres internes
    if ( for_discret%lVari ) then
        ! r??cup??ration du pointeur d'??criture des variables internes ?? cet incr??ment (icontm)
        call jevech('PVARIPR', 'E', icontm)
        ! mise ?? jour des variables internes
        do ii = 1, nbVarloc
            zr(icontm+ii-1) = valvarloc(ii)
        enddo
    endif
999 continue
end subroutine
