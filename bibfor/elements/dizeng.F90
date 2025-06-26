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
subroutine dizeng(for_discret, iret)
!
! --------------------------------------------------------------------------------------------------
!
!        MODÈLE DE D'AMORTISSEUR DE ZENZER GÉNÉRALISÉ
!
!                         e2
!          e1     |-----======-----|
!      ---=====---|                |-----
!                 |--=====----=]---|
!                      e3    n3,a3
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
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8miem.h"
#include "asterfort/assert.h"
#include "asterfort/diraidklv.h"
#include "asterfort/diklvraid.h"
#include "asterfort/infdis.h"
#include "asterfort/jevech.h"
#include "asterfort/pmavec.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rk5adp.h"
#include "asterfort/tecael.h"
#include "asterfort/ut2mlg.h"
#include "asterfort/ut2vlg.h"
#include "asterfort/utmess.h"
#include "asterfort/utpslg.h"
#include "asterfort/utpvlg.h"
#include "asterfort/vecma.h"
#include "asterfort/zengen.h"
#include "blas/dcopy.h"
!
    type(te0047_dscr), intent(in) :: for_discret
    integer(kind=8), intent(out) :: iret
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: imat, ivarim, jdc, irep, jtp, jtm, ifono, icontp, ivarip, iadzi, iazk24
    integer(kind=8) :: icarcr
    integer(kind=8) :: icontm, ii, neq
    real(kind=8) :: r8bid, raide(6), fl(12), klv(78), klc(144), raideurDeno
    character(len=8) :: k8bid
    character(len=24) :: messak(5)
!   pour la loi de comportement
    integer(kind=8) :: nbpara
!   SOUPL_1  RAIDE_2  SOUPL_3  RAID_VISQ   PUIS_VISQ
    parameter(nbpara=5)
    real(kind=8) :: ldcpar(nbpara)
    integer(kind=8) :: ldcpai(1)
    character(len=8) :: ldcpac(1)
    real(kind=8) :: temps0, temps1, dtemps
!   équations du système : sigma, epsivis, epsi, puiss
    integer(kind=8) :: nbequa, nbdecp
    parameter(nbequa=4)
    real(kind=8) :: y0(nbequa), dy0(nbequa), resu(nbequa*2)
    real(kind=8) :: errmax
!
    real(kind=8) :: precis
    parameter(precis=1.0e-08)
!
!   paramètres issus de DEFI_MATERIAU
    integer(kind=8), parameter :: nbcar = 8, ie1 = 1, ie2 = 2, ie3 = 3, in3 = 4
    integer(kind=8), parameter :: ia3 = 5, is1 = 6, is2 = 7, is3 = 8
    character(len=16) :: nomcar(nbcar)
    real(kind=8) :: valcar(nbcar)
    integer(kind=8) :: codcar(nbcar)
    blas_int :: b_incx, b_incy, b_n
    data nomcar/'K1', 'K2', 'K3', 'C', 'PUIS_ALPHA', 'UNSUR_K1', 'UNSUR_K2', 'UNSUR_K3'/
! --------------------------------------------------------------------------------------------------
!
    iret = 0
    neq = for_discret%nno*for_discret%nc
!   récupération du matériau
    call jevech('PMATERC', 'L', imat)
!   variables a t-
    call jevech('PVARIMR', 'L', ivarim)
    call jevech('PCONTMR', 'L', icontm)
!   récupération des caractéristiques élastique
    call jevech('PCADISK', 'L', jdc)
    call infdis('REPK', irep, r8bid, k8bid)
!   seulement en repère local : irep = 2
    if (irep .ne. 2) then
        messak(1) = for_discret%nomte
        messak(2) = 'NON_LINEAR'
        messak(3) = for_discret%type_comp
        messak(4) = for_discret%rela_comp
        call tecael(iadzi, iazk24)
        messak(5) = zk24(iazk24-1+3)
        call utmess('I', 'DISCRETS_5', nk=5, valk=messak)
    end if
!   les caractéristiques sont toujours dans le repère local. on fait seulement une copie
    b_n = to_blas_int(for_discret%nbt)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, zr(jdc), b_incx, klv, b_incy)
!   Récupération des termes diagonaux : raide = klv(i,i)
    call diraidklv(for_discret%nomte, raide, klv)
!   les incréments de déplacement sont nuls
!       ==> récupération de la matrice tangente précédente, si possible
!       ==> si pas possible, calcul d'une tangente pas trop mauvaise,
!           après lecture des paramètres
    if (for_discret%lMatrPred) then
!       tangente précédente si elle existe
        if (abs(zr(ivarim+3)) .gt. r8miem()) then
            raide(1) = zr(ivarim+3)
            resu(1) = zr(ivarim)
            resu(2) = zr(ivarim+1)
            resu(4) = zr(ivarim+2)
            goto 800
        end if
    end if
!
!   récupère tous les paramètres
    valcar(:) = 0.0
    call rcvalb('FPG1', 1, 1, '+', zi(imat), &
                ' ', 'DIS_VISC', 0, ' ', [0.0d0], &
                nbcar, nomcar, valcar, codcar, 0)
!   examen des codcar. assert pas nécessaire ? car un_parmi(a,b) dans les catalogues
    ASSERT(codcar(ie1)+codcar(is1) .eq. 1)
    ASSERT(codcar(ie2)+codcar(is2) .eq. 1)
    ASSERT(codcar(ie3)+codcar(is3) .eq. 1)
    ASSERT(codcar(in3) .eq. 0)
    ASSERT(codcar(ia3) .eq. 0)
!
!   paramètres de la loi de comportement
!     souple1  raide2  souple3  raide_vi  puiss_vi
!
!   cara 1 : souplesse = 1/raideur
    if (codcar(ie1) .eq. 0) then
        ldcpar(1) = 1.0/valcar(ie1)
    else
        ldcpar(1) = valcar(is1)
    end if
!   cara 2 : raideur = 1/souplesse
    if (codcar(ie2) .eq. 0) then
        ldcpar(2) = valcar(ie2)
    else
        ldcpar(2) = 1.0/valcar(is2)
    end if
!   cara 3 : souplesse = 1/raideur
    if (codcar(ie3) .eq. 0) then
        ldcpar(3) = 1.0/valcar(ie3)
    else
        ldcpar(3) = valcar(is3)
    end if
!
    raideurDeno = (ldcpar(1)+ldcpar(3)+ldcpar(2)*ldcpar(1)*ldcpar(3))
    if (raideurDeno .le. r8miem()) then
        messak(1) = for_discret%nomte
        messak(2) = 'NON_LINEAR'
        messak(3) = for_discret%type_comp
        messak(4) = for_discret%rela_comp
        call tecael(iadzi, iazk24)
        messak(5) = zk24(iazk24-1+3)
        call utmess('F', 'DISCRETS_4', nk=5, valk=messak)
    end if
!
    ldcpar(4) = valcar(in3)
    ldcpar(5) = valcar(ia3)
!
!   les incréments de déplacement sont nuls
!       ==> la récupération de la matrice tangente précédente a échouée
!       ==> calcul d'une tangente pas trop mauvaise
    if (for_discret%lMatrPred) then
        raide(1) = (1.0+ldcpar(2)*ldcpar(3))/raideurDeno
        resu(1) = zr(ivarim)
        resu(2) = zr(ivarim+1)
        resu(4) = zr(ivarim+2)
        goto 800
    end if
!
!   loi de comportement non-linéaire : récupération du temps + et - , calcul de dt
    call jevech('PINSTPR', 'L', jtp)
    call jevech('PINSTMR', 'L', jtm)
    temps0 = zr(jtm)
    temps1 = zr(jtp)
    dtemps = temps1-temps0
!   contrôle de rk5 : découpage successif, erreur maximale
    call jevech('PCARCRI', 'L', icarcr)
!   nombre d'itérations maxi  (ITER_INTE_MAXI=-20 par défaut)
    nbdecp = abs(nint(zr(icarcr)))
!   tolérance de convergence (RESI_INTE)
    errmax = zr(icarcr+2)
!   comportement non-linéaire suivant le x local
!   équations du système :
!              1       2         3     4
!       yy   : sigma, epsivisq, epsi,  puiss
!       vari : sigma, epsivisq, puiss, tangente
    if (for_discret%nno .eq. 1) then
        y0(3) = for_discret%ulm(1)
        dy0(3) = for_discret%dul(1)/dtemps
    else
        y0(3) = for_discret%ulm(1+for_discret%nc)-for_discret%ulm(1)
        dy0(3) = (for_discret%dul(1+for_discret%nc)-for_discret%dul(1))/dtemps
    end if
!   récupération de l'effort précédent, suivant l'axe x local
    y0(1) = zr(icontm)
!   récupération des variables internes : epsivis  puiss  tangente
    y0(2) = zr(ivarim+1)
    y0(4) = zr(ivarim+2)
    call rk5adp(nbequa, ldcpar, ldcpai, ldcpac, temps0, &
                dtemps, nbdecp, errmax, y0, dy0, &
                zengen, resu, iret)
!   resu(1:nbeq)            : variables intégrées
!   resu(nbeq+1:2*nbeq)     : d(resu)/d(t) a t+dt
    if (iret .ne. 0) goto 999
!   calcul de la tangente au comportement
    if (abs(resu(nbequa+3)) .gt. precis) then
        raide(1) = resu(nbequa+1)/resu(nbequa+3)
    else
        raide(1) = resu(nbequa+1)
    end if
    if (abs(raide(1)) .lt. precis) then
        raide(1) = (1.0+ldcpar(2)*ldcpar(3))/raideurDeno
    end if
!
800 continue
!   Actualisation de la matrice tangente : klv(i,i) = raide(i)
    call diklvraid(for_discret%nomte, klv, raide)
    if (for_discret%lMatr) then
        call jevech('PMATUUR', 'E', imat)
        if (for_discret%ndim .eq. 3) then
            call utpslg(for_discret%nno, for_discret%nc, for_discret%pgl, klv, zr(imat))
        else if (for_discret%ndim .eq. 2) then
            call ut2mlg(for_discret%nno, for_discret%nc, for_discret%pgl, klv, zr(imat))
        end if
    end if
    !
    if (for_discret%lSigm .or. for_discret%lVect) then
! demi-matrice klv transformée en matrice pleine klc
        call vecma(klv, for_discret%nbt, klc, neq)
! calcul de fl = klc.dul (incrément d'effort)
        call pmavec('ZERO', neq, klc, for_discret%dul, fl)
    end if
! calcul des efforts généralisés
    if (for_discret%lSigm) then
        call jevech('PCONTPR', 'E', icontp)
! Attention aux signes des efforts sur le premier noeud pour MECA_DIS_TR_L et MECA_DIS_T_L
        if (for_discret%nno .eq. 1) then
            do ii = 1, for_discret%nc
                zr(icontp-1+ii) = fl(ii)+zr(icontm-1+ii)
            end do
            zr(icontp) = resu(1)
        else if (for_discret%nno .eq. 2) then
            do ii = 1, for_discret%nc
                zr(icontp-1+ii) = -fl(ii)+zr(icontm-1+ii)
                zr(icontp-1+ii+for_discret%nc) = fl(ii+for_discret%nc)+zr(icontm-1+ii+for_discre&
                                                 &t%nc)
            end do
            zr(icontp) = resu(1)
            zr(icontp+for_discret%nc) = resu(1)
        end if
    end if
! calcul des forces nodales
    if (for_discret%lVect) then
        call jevech('PVECTUR', 'E', ifono)
! Attention aux signes des efforts sur le premier noeud pour MECA_DIS_TR_L et MECA_DIS_T_L
        if (for_discret%nno .eq. 1) then
            do ii = 1, for_discret%nc
                fl(ii) = fl(ii)+zr(icontm-1+ii)
            end do
            fl(1) = resu(1)
        else if (for_discret%nno .eq. 2) then
            do ii = 1, for_discret%nc
                fl(ii) = fl(ii)-zr(icontm-1+ii)
                fl(ii+for_discret%nc) = fl(ii+for_discret%nc)+zr(icontm-1+ii+for_discret%nc)
            end do
            fl(1) = -resu(1)
            fl(1+for_discret%nc) = resu(1)
        end if
! forces nodales aux noeuds 1 et 2 (repère global)
        if (for_discret%nc .ne. 2) then
            call utpvlg(for_discret%nno, for_discret%nc, for_discret%pgl, fl, zr(ifono))
        else
            call ut2vlg(for_discret%nno, for_discret%nc, for_discret%pgl, fl, zr(ifono))
        end if
    end if
! mise à jour des variables internes : sigma  epsivis  puiss tangente
    if (for_discret%lVari) then
        call jevech('PVARIPR', 'E', ivarip)
        zr(ivarip) = resu(1)
        zr(ivarip+1) = resu(2)
        zr(ivarip+2) = resu(4)
        zr(ivarip+3) = raide(1)
        if (for_discret%nno .eq. 2) then
            zr(ivarip+4) = resu(1)
            zr(ivarip+4+1) = resu(2)
            zr(ivarip+4+2) = resu(4)
            zr(ivarip+4+3) = raide(1)
        end if
    end if
999 continue
end subroutine
