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
! person_in_charge: jean-luc.flejou at edf.fr
!
subroutine dis_contact_frot(for_discret, iret)
!
! --------------------------------------------------------------------------------------------------
!
!     RELATION DE COMPORTEMENT "DIS_CHOC" : DYNAMIQUE AVEC OU SANS FROTTEMENT
!                                         : STATIQUE SANS FROTTEMENT
!
! --------------------------------------------------------------------------------------------------
!
! IN    for_discret : voir l'appel
! OUT   iret        : code retour
!
! --------------------------------------------------------------------------------------------------
!
    use te0047_type, only: te0047_dscr
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/diraidklv.h"
#include "asterfort/diklvraid.h"
#include "asterfort/ldc_dis_contact_frot.h"
#include "asterfort/infdis.h"
#include "asterfort/jevech.h"
#include "asterfort/pmavec.h"
#include "asterfort/rcadlv.h"
#include "asterfort/rcvala.h"
#include "asterfort/rk5adp.h"
#include "asterfort/tecach.h"
#include "asterfort/ut2mgl.h"
#include "asterfort/ut2mlg.h"
#include "asterfort/ut2vgl.h"
#include "asterfort/ut2vlg.h"
#include "asterfort/utmess.h"
#include "asterfort/utpsgl.h"
#include "asterfort/utpslg.h"
#include "asterfort/utpvgl.h"
#include "asterfort/utpvlg.h"
#include "asterfort/vecma.h"
#include "blas/dcopy.h"
!
    type(te0047_dscr), intent(in) :: for_discret
    integer(kind=8), intent(out) :: iret
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: jdc, irep, imat, ivarim, ii, ivitp, idepen, iviten, igeom, ivarip
    integer(kind=8) :: iretlc, ifono, imatsym, icarcr, iiter, iterat
    integer(kind=8) :: icontm, icontp
!
    real(kind=8) :: klc(for_discret%neq*for_discret%neq), klv(for_discret%nbt)
    real(kind=8) :: dvl(for_discret%neq), dpe(for_discret%neq), dve(for_discret%neq)
    real(kind=8) :: fl(for_discret%neq)
    real(kind=8) :: force(3), raide(6)
    real(kind=8) :: r8bid
    character(len=8) :: k8bid
    aster_logical :: Prediction, Dynamique
! --------------------------------------------------------------------------------------------------
    integer(kind=8), parameter :: nbre1 = 9
    real(kind=8) :: valre1(nbre1)
    integer(kind=8) :: codre1(nbre1)
    character(len=8) :: nomre1(nbre1)
    integer(kind=8) :: nbpar
    real(kind=8) :: valpar
    character(len=8) :: nompar
    integer(kind=8) :: jadre1, jcodre1
!
    data nomre1/'RIGI_NOR', 'RIGI_TAN', 'AMOR_NOR', 'AMOR_TAN', 'COULOMB', &
        'DIST_1', 'DIST_2', 'JEU', 'CONTACT'/
! --------------------------------------------------------------------------------------------------
!   Pour l'intégration de la loi de comportement
    real(kind=8) :: temps0, temps1, dtemps
!   Paramètres de la loi :        Kn       Kt       mu       cn
    integer(kind=8), parameter :: ikn = 1, ikt = 2, imu = 3, icn = 4
!   Paramètres de la loi :        ct       jeu,      ky,    kz
    integer(kind=8), parameter :: ict = 5, ijeu = 6, iky = 7, ikz = 8
    integer(kind=8), parameter :: nbpara = 8
    real(kind=8) :: ldcpar(nbpara)
    integer(kind=8) :: ldcpai(2)
    character(len=8) :: ldcpac(1)
!   Équations du système
    integer(kind=8), parameter :: nbequa = 14
    real(kind=8) :: y0(nbequa), dy0(nbequa), resu(nbequa*2), errmax
    integer(kind=8) :: nbdecp
!   Variables internes
    integer(kind=8), parameter :: nbvari = 9, nbcorr = 8, idebut = 9
    integer(kind=8) :: Correspond(nbcorr)
    real(kind=8) :: varmo(nbvari), varpl(nbvari)
!
    character(len=32) :: messak(3)
! --------------------------------------------------------------------------------------------------
    integer(kind=8) :: nbout
    real(kind=8) :: xl(6), xd(3), rignor, rigtan, coulom, deplac, evoljeu0, evoljeu1, xjeu
    real(kind=8) :: LgDiscret, Dist12, inst0, inst1
    blas_int :: b_incx, b_incy, b_n
! --------------------------------------------------------------------------------------------------
!
    iret = 0
!   Paramètres en entrée
    call jevech('PCADISK', 'L', jdc)
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PCONTMR', 'L', icontm)
    call jevech('PMATERC', 'L', imat)
!   on recupere le no de l'iteration de newton
    call jevech('PITERAT', 'L', iiter)
    iterat = zi(iiter)
!
    call infdis('REPK', irep, r8bid, k8bid)
!   irep = 1 = matrice en repère global ==> passer en local
    if (irep .eq. 1) then
        if (for_discret%ndim .eq. 3) then
            call utpsgl(for_discret%nno, for_discret%nc, for_discret%pgl, zr(jdc), klv)
        else if (for_discret%ndim .eq. 2) then
            call ut2mgl(for_discret%nno, for_discret%nc, for_discret%pgl, zr(jdc), klv)
        end if
    else
        b_n = to_blas_int(for_discret%nbt)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(jdc), b_incx, klv, b_incy)
    end if
!   Récupération des termes diagonaux : raide = klv(i,i)
    call diraidklv(for_discret%nomte, raide, klv)
!
!   Champ de vitesse
    Dynamique = ASTER_FALSE
    call tecach('ONO', 'PVITPLU', 'L', iretlc, iad=ivitp)
    if (iretlc .eq. 0) then
        if (for_discret%ndim .eq. 3) then
            call utpvgl(for_discret%nno, for_discret%nc, for_discret%pgl, zr(ivitp), dvl)
        else if (for_discret%ndim .eq. 2) then
            call ut2vgl(for_discret%nno, for_discret%nc, for_discret%pgl, zr(ivitp), dvl)
        end if
        Dynamique = ASTER_TRUE
    else
        dvl(:) = 0.0d0
    end if
!   Champ de déplacement d'entrainement
    call tecach('ONO', 'PDEPENT', 'L', iretlc, iad=idepen)
    if (iretlc .eq. 0) then
        if (for_discret%ndim .eq. 3) then
            call utpvgl(for_discret%nno, for_discret%nc, for_discret%pgl, zr(idepen), dpe)
        else if (for_discret%ndim .eq. 2) then
            call ut2vgl(for_discret%nno, for_discret%nc, for_discret%pgl, zr(idepen), dpe)
        end if
    else
        dpe(:) = 0.0d0
    end if
!   Champ de vitesse d'entrainement
    call tecach('ONO', 'PVITENT', 'L', iretlc, iad=iviten)
    if (iretlc .eq. 0) then
        if (for_discret%ndim .eq. 3) then
            call utpvgl(for_discret%nno, for_discret%nc, for_discret%pgl, zr(iviten), dve)
        else if (for_discret%ndim .eq. 2) then
            call ut2vgl(for_discret%nno, for_discret%nc, for_discret%pgl, zr(iviten), dve)
        end if
    else
        dve(:) = 0.d0
    end if
!
!   Variables internes
    call jevech('PVARIMR', 'L', ivarim)
    do ii = 1, nbvari
        varmo(ii) = zr(ivarim+ii-1)
        varpl(ii) = varmo(ii)
    end do
!
!   Temps + et - , calcul de dt
    temps0 = for_discret%TempsMoins
    temps1 = for_discret%TempsPlus
    dtemps = temps1-temps0
!   contrôle de rk5 : découpage successif, erreur maximale
    call jevech('PCARCRI', 'L', icarcr)
!   nombre d'itérations maxi (ITER_INTE_MAXI=-20 par défaut)
    nbdecp = abs(int(zr(icarcr)))
!   tolérance de convergence (RESI_INTE=1.0E-06 par défaut)
    errmax = zr(icarcr+2)
!
! --------------------------------------------------------------------------------------------------
!   Relation de comportement de choc
!
!   Coordonnees du discret dans le repère local
    xl(:) = 0.0
    if (for_discret%ndim .eq. 3) then
        call utpvgl(for_discret%nno, 3, for_discret%pgl, zr(igeom), xl)
    else if (for_discret%ndim .eq. 2) then
        call ut2vgl(for_discret%nno, 2, for_discret%pgl, zr(igeom), xl)
    end if
!
!   Caractéristiques du matériau
!    1          2          3          4          5         6        7        8     9
!   'RIGI_NOR','RIGI_TAN','AMOR_NOR','AMOR_TAN','COULOMB','DIST_1','DIST_2','JEU','CONTACT'
    valre1(:) = 0.0
    nbpar = 0; nompar = ' '; valpar = 0.0
!   Si mot_cle RIGI_NOR ==> rignor = valre1(1) sinon rignor = raide(1)
!   Si mot_cle RIGI_TAN ==> rigtan = valre1(2) sinon rigtan = 0.0
    rignor = raide(1)
    rigtan = 0.0
    call rcvala(zi(imat), ' ', 'DIS_CONTACT', nbpar, nompar, &
                [valpar], nbre1, nomre1, valre1, codre1, &
                0, nan='NON')
!
    if (nint(valre1(9)) .ne. 0) then
        messak(1) = 'DIS_CONTACT'
        messak(2) = 'DIS_CONTACT'
        messak(3) = '"1D"'
        call utmess('F', 'DISCRETS_35', nk=3, valk=messak)
    end if
!
    if (codre1(1) .eq. 0) rignor = valre1(1)
    if (codre1(2) .eq. 0) rigtan = valre1(2)
    coulom = valre1(5)
!   Paramètres de la loi de comportement
    ldcpar(:) = 0.0d0
    ldcpar(ikn) = rignor
    ldcpar(ikt) = rigtan
    ldcpar(imu) = coulom
    if (codre1(3) .eq. 0) ldcpar(icn) = valre1(3)
    if (codre1(4) .eq. 0) ldcpar(ict) = valre1(4)
    ldcpar(iky) = raide(2)
    ldcpar(ikz) = raide(3)
!   calcul du jeu final
    LgDiscret = 0.0
    if (for_discret%nno .eq. 2) then
!       Longueur du discret
        xd(1:3) = xl(1+for_discret%ndim:2*for_discret%ndim)-xl(1:for_discret%ndim)
        LgDiscret = xd(1)
        Dist12 = -valre1(6)-valre1(7)
    else
        Dist12 = valre1(8)-valre1(6)
    end if
    ldcpar(ijeu) = LgDiscret+Dist12
!   La loi complète
    ldcpai(1) = 2
!   si raideur (tangente = faible) ou (coulom=0) ==> pas de frottement, pas de seuil
!       juste du contact, avec/sans amortissement normal
    if ((abs(rigtan) <= r8prem()) .or. (abs(coulom) <= r8prem())) then
        ldcpai(1) = 1
    end if
!   Traitement de l'évolution du jeu
    jadre1 = 0; nbout = 0; jcodre1 = 0
    call rcadlv(' ', 1, 1, '+', zi(imat), &
                ' ', 'DIS_CONTACT', 'INST_COMP_INIT', 0, [' '], &
                [0.d0], jadre1, nbout, jcodre1, 0)
    if (jcodre1 .eq. 0 .and. nbout .eq. 2) then
        inst0 = zr(jadre1); inst1 = zr(jadre1+1)
        ASSERT(inst0 < inst1)
        evoljeu0 = min(max(0.0, (temps0-inst0)/(inst1-inst0)), 1.0)
        evoljeu1 = min(max(0.0, (temps1-inst0)/(inst1-inst0)), 1.0)
!
        ldcpai(2) = 1
!       Si le jeu évolue, alors ni frottement, ni amortissement, juste du contact
        if (abs(evoljeu1-evoljeu0) > r8prem()) then
            ldcpai(1) = 3
        end if
    else
!       Pas de fonction d'évolution du jeu
        ldcpai(2) = 0
        evoljeu0 = 1.0
        evoljeu1 = 1.0
    end if
!
!   Équations du système :
!              1   2   3   4   5   6   7   8   9   10    11    12   13   14
!       yy   : Ux, Uy, Uz, Fx, Fy, Fz, vx, vy, vz, Uyan, Uzan, Fcy, Fcz, jeu
!       vari :                         5   6   7   3     4     1    2    8
    Correspond(:) = [12, 13, 10, 11, 7, 8, 9, 14]
    y0(:) = 0.0
    dy0(:) = 0.0
    do ii = 1, nbcorr
        y0(Correspond(ii)) = varmo(ii)
    end do
!   Les dérivées
    if (for_discret%nno .eq. 1) then
        y0(1) = for_discret%ulm(1)+dpe(1)
        y0(2) = for_discret%ulm(2)+dpe(2)
        y0(4) = zr(icontm)
        y0(5) = zr(icontm+1)
        dy0(1) = for_discret%dul(1)/dtemps
        dy0(2) = for_discret%dul(2)/dtemps
        dy0(7) = (dvl(1)+dve(1)-y0(7))/dtemps
        dy0(8) = (dvl(2)+dve(2)-y0(8))/dtemps
        if (for_discret%ndim .eq. 3) then
            y0(3) = for_discret%ulm(3)+dpe(3)
            y0(6) = zr(icontm+2)
            dy0(3) = for_discret%dul(3)/dtemps
            dy0(9) = (dvl(3)+dve(3)-y0(9))/dtemps
        end if
    else
        y0(1) = ( &
                for_discret%ulm(1+for_discret%nc)-for_discret%ulm(1)+dpe(1+for_discret%nc)-dpe(1&
                &) &
                )
        y0(2) = ( &
                for_discret%ulm(2+for_discret%nc)-for_discret%ulm(2)+dpe(2+for_discret%nc)-dpe(2&
                &) &
                )
        y0(4) = zr(icontm)
        y0(5) = zr(icontm+1)
        dy0(1) = (for_discret%dul(1+for_discret%nc)-for_discret%dul(1))/dtemps
        dy0(2) = (for_discret%dul(2+for_discret%nc)-for_discret%dul(2))/dtemps
        dy0(7) = (dvl(1+for_discret%nc)-dvl(1)+dve(1+for_discret%nc)-dve(1)-y0(7))/dtemps
        dy0(8) = (dvl(2+for_discret%nc)-dvl(2)+dve(2+for_discret%nc)-dve(2)-y0(8))/dtemps
        if (for_discret%ndim .eq. 3) then
            y0(3) = ( &
                    for_discret%ulm(3+for_discret%nc)-for_discret%ulm(3)+dpe(3+for_discret%nc)-d&
                    &pe(3) &
                    )
            y0(6) = zr(icontm+2)
            dy0(3) = (for_discret%dul(3+for_discret%nc)-for_discret%dul(3))/dtemps
            dy0(9) = (dvl(3+for_discret%nc)-dvl(3)+dve(3+for_discret%nc)-dve(3)-y0(9))/dtemps
        end if
    end if
!
!    dy0(14) = lejeu*(evoljeu1-evoljeu0)/dtemps
!   calcul de la vitesse d'évolution du jeu
    if (nint(varmo(idebut)) .eq. 0) then
        y0(14) = LgDiscret
    end if
    dy0(14) = Dist12*(evoljeu1-evoljeu0)/dtemps
!
    force(:) = 0.0
!   Prédiction en dynamique, on retourne les efforts précédents
    Prediction = ((iterat .eq. 1) .and. (for_discret%option .eq. 'FULL_MECA'))
    Prediction = Prediction .or. ((iterat .eq. 1) .and. (for_discret%option .eq. 'RAPH_MECA'))
!
!   Soit on intègre le jeu soit on prend sa valeur
!       ldcpai(2) = 1 : intégration du jeu
!       ldcpai(2) = 0 : valeur finale
    xjeu = y0(14)*float(ldcpai(2))+ldcpar(ijeu)*(1.0-float(ldcpai(2)))
!
    if (Prediction .and. Dynamique .and. (ldcpai(2) .eq. 0)) then
        r8bid = y0(1)+dy0(1)*dtemps+xjeu
        raide(1) = 0.0
        if (r8bid <= 0.0) then
            raide(1) = rignor
            if (for_discret%option .eq. 'RAPH_MECA') then
                force(1:3) = y0(4:6)
                force(1) = raide(1)*r8bid
            end if
        end if
        goto 888
    end if
!
!   Si le discret est en contact        ==> raideur de contact
!                 n'est pas en contact  ==> raideur nulle
!       mise à jour de la raideur en sortie de la ldc
    if (y0(1)+xjeu <= 0.0) then
        raide(1) = rignor
    else
        raide(1) = 0.0
    end if
!   Si on a du contact initial : on retourne l'effort de contact
    if ((xjeu <= 0.0) .and. (nint(varmo(idebut)) .eq. 0)) then
        y0(4) = raide(1)*xjeu
    end if
!   Pas de normalisation des équations : on laisse l'algo faire, c'est prévu :)
!   On intègre
    iret = 0
    call rk5adp(nbequa, ldcpar, ldcpai, ldcpac, temps0, &
                dtemps, nbdecp, errmax, y0, dy0, &
                ldc_dis_contact_frot, resu, iret)
!   resu(1:nbequa)              : variables intégrées
!   resu(nbequa+1:2*nbequa)     : d(resu)/d(t) a t+dt
    if (iret .ne. 0) goto 999
!   Les efforts
    force(1:3) = resu(4:6)
!   Les variables internes
    do ii = 1, nbcorr
        varpl(ii) = resu(Correspond(ii))
    end do
    varpl(idebut) = 1.0
!   Les raideurs
    xjeu = resu(14)*float(ldcpai(2))+ldcpar(ijeu)*(1.0-float(ldcpai(2)))
    if (resu(1)+xjeu <= 0.0) then
        raide(1) = rignor
        if (for_discret%option .ne. 'RAPH_MECA') then
            deplac = resu(1)-y0(1)
            if (abs(deplac) > r8prem()) then
                raide(1) = abs((resu(4)-y0(4))/deplac)
            end if
            deplac = resu(2)-y0(2)
            if (abs(deplac) > r8prem()) then
                raide(2) = abs((resu(5)-y0(5))/deplac)
            end if
            deplac = resu(3)-y0(3)
            if (abs(deplac) > r8prem()) then
                raide(3) = abs((resu(6)-y0(6))/deplac)
            end if
        end if
    else
        raide(1) = 0.0
        force(1) = 0.0
    end if
!
888 continue
! --------------------------------------------------------------------------------------------------
!   Actualisation de la matrice tangente : klv(i,i) = raide(i)
    call diklvraid(for_discret%nomte, klv, raide)
    if (for_discret%lMatr) then
        call jevech('PMATUUR', 'E', imatsym)
        if (for_discret%ndim .eq. 3) then
            call utpslg(for_discret%nno, for_discret%nc, for_discret%pgl, klv, zr(imatsym))
        else
            call ut2mlg(for_discret%nno, for_discret%nc, for_discret%pgl, klv, zr(imatsym))
        end if
    end if
!
    if (for_discret%lVect .or. for_discret%lSigm) then
! Demi-matrice klv transformée en matrice pleine klc
        call vecma(klv, for_discret%nbt, klc, for_discret%neq)
! Calcul de fl = klc.dul (incrément d'effort)
        call pmavec('ZERO', for_discret%neq, klc, for_discret%dul, fl)
    end if
    !
! calcul des efforts généralisés et des forces nodales
    if (for_discret%lSigm) then
        call jevech('PCONTPR', 'E', icontp)
! Attention aux signes des efforts sur le premier noeud pour MECA_DIS_TR_L et MECA_DIS_T_L
        if (for_discret%nno .eq. 1) then
            do ii = 1, for_discret%neq
                zr(icontp-1+ii) = fl(ii)+zr(icontm-1+ii)
            end do
        else if (for_discret%nno .eq. 2) then
            do ii = 1, for_discret%nc
                zr(icontp-1+ii) = -fl(ii)+zr(icontm-1+ii)
                zr(icontp-1+ii+for_discret%nc) = fl(ii+for_discret%nc)+zr(icontm-1+ii+for_discre&
                                                 &t%nc)
            end do
        end if
        if (for_discret%nno .eq. 1) then
            zr(icontp-1+1) = force(1)
            zr(icontp-1+2) = force(2)
            if (for_discret%ndim .eq. 3) then
                zr(icontp-1+3) = force(3)
            end if
        else if (for_discret%nno .eq. 2) then
            zr(icontp-1+1) = force(1)
            zr(icontp-1+1+for_discret%nc) = force(1)
            zr(icontp-1+2) = force(2)
            zr(icontp-1+2+for_discret%nc) = force(2)
            if (for_discret%ndim .eq. 3) then
                zr(icontp-1+3) = force(3)
                zr(icontp-1+3+for_discret%nc) = force(3)
            end if
        end if
    end if
! calcul des forces nodales
    if (for_discret%lVect) then
        call jevech('PVECTUR', 'E', ifono)
! Attention aux signes des efforts sur le premier noeud pour MECA_DIS_TR_L et MECA_DIS_T_L
        if (for_discret%nno .eq. 1) then
            do ii = 1, for_discret%neq
                fl(ii) = fl(ii)+zr(icontm-1+ii)
            end do
        else if (for_discret%nno .eq. 2) then
            do ii = 1, for_discret%nc
                fl(ii) = fl(ii)-zr(icontm-1+ii)
                fl(ii+for_discret%nc) = fl(ii+for_discret%nc)+zr(icontm-1+ii+for_discret%nc)
            end do
        end if
        if (for_discret%nno .eq. 1) then
            fl(1) = force(1)
            fl(2) = force(2)
            if (for_discret%ndim .eq. 3) then
                fl(3) = force(3)
            end if
        else if (for_discret%nno .eq. 2) then
            fl(1) = -force(1)
            fl(1+for_discret%nc) = force(1)
            fl(2) = -force(2)
            fl(2+for_discret%nc) = force(2)
            if (for_discret%ndim .eq. 3) then
                fl(3) = -force(3)
                fl(3+for_discret%nc) = force(3)
            end if
        end if
!       Forces nodales aux noeuds 1 et 2 (repère global)
        if (for_discret%nc .ne. 2) then
            call utpvlg(for_discret%nno, for_discret%nc, for_discret%pgl, fl, zr(ifono))
        else
            call ut2vlg(for_discret%nno, for_discret%nc, for_discret%pgl, fl, zr(ifono))
        end if
    end if
!   mise à jour des variables internes
    if (for_discret%lVari) then
        call jevech('PVARIPR', 'E', ivarip)
        if (for_discret%nno .eq. 1) then
            do ii = 1, nbvari
                zr(ivarip+ii-1) = varpl(ii)
            end do
        else
            do ii = 1, nbvari
                zr(ivarip+ii-1) = varpl(ii)
                zr(ivarip+ii-1+nbvari) = varpl(ii)
            end do
        end if
    end if
!
999 continue
end subroutine
