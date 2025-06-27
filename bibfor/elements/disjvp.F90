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
subroutine disjvp(for_discret, iret)
!
! person_in_charge: jean-luc.flejou at edf.fr
! --------------------------------------------------------------------------------------------------
!
!        MODELE ELASTO PLASTIQUE ENDOMMAGEANT EN FLEXION (1DDL=RZ)
!
!                            Kp
!         Ke(Yp,Ym)  |-----|/|/|/|-----|
!      ----|/|/|/|---|                 |-----
!                    |-----|______-----|
!                          Xm,theta_p
!
! IN    for_discret : voir l'appel
!     rappel pour mémoire :
!     option   : option de calcul
!     nomte    : nom terme élémentaire
!     ndim     : dimension du problème
!     nbt      : nombre de terme dans la matrice de raideur
!     nno      : nombre de noeuds de l'élément
!     nc       : nombre de composante par noeud
!     ulm      : déplacement moins
!     dul      : incrément de déplacement
!     pgl      : matrice de passage de global à local
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
#include "asterfort/diraidklv.h"
#include "asterfort/diklvraid.h"
#include "asterfort/infdis.h"
#include "asterfort/jevech.h"
#include "asterfort/pmavec.h"
#include "asterfort/rcvala.h"
#include "asterfort/rk5adp.h"
#include "asterfort/tecach.h"
#include "asterfort/tecael.h"
#include "asterfort/ut2mgl.h"
#include "asterfort/ut2mlg.h"
#include "asterfort/ut2vlg.h"
#include "asterfort/ut2vgl.h"
#include "asterfort/utpslg.h"
#include "asterfort/utpsgl.h"
#include "asterfort/utpvlg.h"
#include "asterfort/utpvgl.h"
#include "asterfort/vecma.h"
#include "asterfort/ldc_disjvp.h"
#include "asterfort/utmess.h"
#include "blas/dcopy.h"
!
!
    type(te0047_dscr), intent(in) :: for_discret
    integer(kind=8), intent(out) :: iret
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: imat, ivarim, jdc, irep, jtp, jtm, idepen
    integer(kind=8) :: imatsym, ifono, iretlc, icontp, ivarip
    integer(kind=8) :: iadzi, iazk24, igeom
    integer(kind=8) :: icarcr
    integer(kind=8) :: icontm, ii
    character(len=24) :: messak(5)
!
!
    real(kind=8) :: klc(for_discret%neq*for_discret%neq), klv(for_discret%nbt)
    real(kind=8) :: dpe(for_discret%neq)
    real(kind=8) :: fl(for_discret%neq)
    real(kind=8) :: force(3), raide(6)
    real(kind=8) :: r8bid
    character(len=8) :: k8bid
!
! --------------------------------------------------------------------------------------------------
    integer(kind=8), parameter :: nbre1 = 8
    real(kind=8) :: valre1(nbre1)
    integer(kind=8) :: codre1(nbre1)
    character(len=8) :: nomre1(nbre1)
    integer(kind=8) :: nbpar
    real(kind=8) :: valpar
    character(len=8) :: nompar
    data nomre1/'KE', 'KP', 'KDP', 'KDM', 'RDP', 'RDM', 'MYP', 'MYM'/
!
! --------------------------------------------------------------------------------------------------
!   Pour l'intégration de la loi de comportement
    real(kind=8) :: temps0, temps1, dtemps
!   Paramètres de la loi :     KE      KP    KDP    KDM     RDP    RDM    MYP    MYM
    integer(kind=8), parameter :: ike = 1, ikdp = 3, ikdm = 4, irdp = 5, irdm = 6
!   integer, parameter      :: ike=1, ikp=2, ikdp=3, ikdm=4, irdp=5, irdm=6, imyp=7, imym=8
    integer(kind=8), parameter :: nbpara = 8
    real(kind=8) :: ldcpar(nbpara)
    integer(kind=8) :: ldcpai(1)
    character(len=8) :: ldcpac(1)
!   Équations du système
    integer(kind=8), parameter :: nbequa = 7
    real(kind=8) :: y0(nbequa), dy0(nbequa), resu(nbequa*2), errmax, ynorme(nbequa)
    integer(kind=8) :: nbdecp
!   Variables internes
    integer(kind=8), parameter :: nbvari = 9, nbcorr = 6, idebut = nbvari, iddp = 7, iddm = 8
    integer(kind=8) :: Correspond(nbcorr)
    real(kind=8) :: varmo(nbvari), varpl(nbvari)
!
!   système d'équations :
    integer(kind=8), parameter :: imoment = 1, itheta = 2, ithetap = 3, idp = 4, idm = 5, ixm = 6
    integer(kind=8), parameter :: idiss = 7
! --------------------------------------------------------------------------------------------------
    real(kind=8) :: xl(7), deplac, Dp, Dm
    blas_int :: b_incx, b_incy, b_n
! --------------------------------------------------------------------------------------------------
!   Seulement sur SEG2
    if (for_discret%nomte .ne. 'MECA_DIS_TR_L') then
        messak(1) = for_discret%nomte
        messak(2) = 'NON_LINEAR'
        messak(3) = for_discret%type_comp
        messak(4) = for_discret%rela_comp
        call tecael(iadzi, iazk24)
        messak(5) = zk24(iazk24-1+3)
        call utmess('F', 'DISCRETS_11', nk=5, valk=messak)
    end if
!
    iret = 0
!   Paramètres en entrée
    call jevech('PGEOMER', 'L', igeom)
!   récupération du matériau
    call jevech('PMATERC', 'L', imat)
!   variables a t-
    call jevech('PCONTMR', 'L', icontm)
!   récupération des caractéristiques élastique
    call jevech('PCADISK', 'L', jdc)
    call infdis('REPK', irep, r8bid, k8bid)
!   seulement en repère local : irep = 2
!   les caractéristiques sont toujours dans le repère local. on fait seulement une copie
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
!   les caractéristiques sont toujours dans le repère local. on fait seulement une copie
!   Récupération des termes diagonaux : raide = klv(i,i)
    call diraidklv(for_discret%nomte, raide, klv)
!
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
!
!   Variables internes
    call jevech('PVARIMR', 'L', ivarim)
    do ii = 1, nbvari
        varmo(ii) = zr(ivarim+ii-1)
        varpl(ii) = varmo(ii)
    end do
!
!   loi de comportement non-linéaire : récupération du temps + et - , calcul de dt
    call jevech('PINSTPR', 'L', jtp)
    call jevech('PINSTMR', 'L', jtm)
    temps0 = zr(jtm)
    temps1 = zr(jtp)
    dtemps = temps1-temps0
!   contrôle de rk5 : découpage successif, erreur maximale
    call jevech('PCARCRI', 'L', icarcr)
!   nombre d'itérations maxi (ITER_INTE_MAXI=-20 par defaut)
    nbdecp = abs(int(zr(icarcr)))
!   tolérance de convergence (RESI_INTE=1.0E-06 par défaut)
    errmax = zr(icarcr+2)
!
! --------------------------------------------------------------------------------------------------
!   Relation de comportement
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
!    1    2    3      4    5     6      7     8
!   'KE','KP','KDP','KDM','RDP','RDM','MYP','MYM'
!   récupère tous les paramètres
    valre1(:) = 0.0
    nbpar = 0
    nompar = ' '
    valpar = 0.d0
    call rcvala(zi(imat), ' ', 'JONC_ENDO_PLAS', nbpar, nompar, &
                [valpar], nbre1, nomre1, valre1, codre1, &
                0, nan='NON')
!   recuperation des parametres materiaux
    ldcpar(1:nbpara) = valre1(1:nbre1)
!
!   comportement non-linéaire suivant le RZ local
!   équations du système :
!              1      2       3    4   5   6   7
!       yy   : M, theta, theta_p,  dp, dm, Xm, Diss
!       vari :        1       2    3   4   5    6
    Correspond(:) = [2, 3, 4, 5, 6, 7]
    y0(:) = 0.0
    dy0(:) = 0.0
    do ii = 1, nbcorr
        y0(Correspond(ii)) = varmo(ii)
    end do
!
!   récupération du moment précédent, suivant l'axe z local
    y0(imoment) = zr(icontm+5)
!   récupération de la rotation précédente, suivant l'axe z local
    y0(itheta) = for_discret%ulm(6+for_discret%nc)-for_discret%ulm(6)+dpe(6+for_discret%nc)-dpe(6&
                 &)
!   initialisation de Yp et Ym
    if (nint(varmo(idebut)) .eq. 0) then
        y0(idp) = ldcpar(irdp)
        y0(idm) = abs(ldcpar(irdm))
    end if
!   récupération de l'increment de rotation, suivant l'axe z local
    dy0(itheta) = (for_discret%dul(6+for_discret%nc)-for_discret%dul(6))/dtemps
!   Normalisation des équations, par défaut 1
!   ynorme(:)=1.0
    ynorme(imoment) = ldcpar(ike)*min(ldcpar(irdp), abs(ldcpar(irdm)))*0.6
    ynorme(ixm) = ldcpar(ike)*min(ldcpar(irdp), abs(ldcpar(irdm)))*0.6
    ynorme(itheta) = min(ldcpar(irdp), abs(ldcpar(irdm)))*0.6
    ynorme(ithetap) = min(ldcpar(irdp), abs(ldcpar(irdm)))*0.6
    ynorme(idp) = min(ldcpar(irdp), abs(ldcpar(irdm)))*0.6
    ynorme(idm) = min(ldcpar(irdp), abs(ldcpar(irdm)))*0.6
    ynorme(idiss) = (ldcpar(ike)*min(ldcpar(irdp), abs(ldcpar(irdm)))**2)*0.6
!   Integration de la loi de comportement
    call rk5adp(nbequa, ldcpar, ldcpai, ldcpac, temps0, &
                dtemps, nbdecp, errmax, y0, dy0, &
                ldc_disjvp, resu, iret, ynorme=ynorme)
!   resu(1:nbeq)            : variables intégrées
!   resu(nbeq+1:2*nbeq)     : d(resu)/d(t) a t+dt
    if (iret .ne. 0) goto 999
!   Les efforts
    force(1) = resu(imoment)
!   Les variables internes
    do ii = 1, nbcorr
        varpl(ii) = resu(Correspond(ii))
    end do
    Dp = (1.0-ldcpar(ikdp)/ldcpar(ike))*(1.0-ldcpar(irdp)/resu(idp))
    Dm = (1.0-ldcpar(ikdm)/ldcpar(ike))*(1.0-abs(ldcpar(irdm))/resu(idm))
    varpl(iddp) = Dp
    varpl(iddm) = Dm
    varpl(idebut) = 1.0
!
!   Calcul des raideurs (raideur secante)
    deplac = resu(itheta)-y0(itheta)
    raide(6) = ldcpar(ike)
    if (abs(deplac) > r8prem()) then
        raide(6) = abs((resu(imoment)-y0(imoment))/deplac)
    end if
! --------------------------------------------------------------------------------------------------
!   Actualisation de la matrice tangente : klv(i,i) = raide(i)
    call diklvraid(for_discret%nomte, klv, raide)
    if (for_discret%lMatr) then
        call jevech('PMATUUR', 'E', imatsym)
        if (for_discret%ndim .eq. 3) then
            call utpslg(for_discret%nno, for_discret%nc, for_discret%pgl, klv, zr(imatsym))
        else if (for_discret%ndim .eq. 2) then
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
!       calcul des efforts généralisés, des forces nodales
        call jevech('PVECTUR', 'E', ifono)
        call jevech('PCONTPR', 'E', icontp)
!       efforts généralisés aux noeuds 1 et 2 (repère local)
!       on change le signe des efforts sur le premier noeud pour les MECA_DIS_TR_L et MECA_DIS_T_L
        do ii = 1, for_discret%nc
            zr(icontp-1+ii) = -fl(ii)+zr(icontm-1+ii)
            zr(icontp-1+ii+for_discret%nc) = fl(ii+for_discret%nc)+zr(icontm-1+ii+for_discret%nc)
            fl(ii) = fl(ii)-zr(icontm-1+ii)
            fl(ii+for_discret%nc) = fl(ii+for_discret%nc)+zr(icontm-1+ii+for_discret%nc)
        end do
        zr(icontp+5) = force(1)
        zr(icontp+for_discret%nc+5) = force(1)
        fl(6) = -force(1)
        fl(6+for_discret%nc) = force(1)
    end if
!   forces nodales aux noeuds 1 et 2 (repère global)
    if (for_discret%lVect) then
        if (for_discret%nc .ne. 2) then
            call utpvlg(for_discret%nno, for_discret%nc, for_discret%pgl, fl, zr(ifono))
        else
            call ut2vlg(for_discret%nno, for_discret%nc, for_discret%pgl, fl, zr(ifono))
        end if
    end if
!   Mise à jour des variables internes
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
999 continue
end subroutine
