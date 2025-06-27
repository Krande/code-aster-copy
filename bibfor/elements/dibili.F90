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
subroutine dibili(for_discret, iret)
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
#include "asterfort/dinon4.h"
#include "asterfort/diklvraid.h"
#include "asterfort/dinonc.h"
#include "asterfort/infdis.h"
#include "asterfort/jevech.h"
#include "asterfort/pmavec.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/tecael.h"
#include "asterfort/ut2mlg.h"
#include "asterfort/ut2vlg.h"
#include "asterfort/utmess.h"
#include "asterfort/utpslg.h"
#include "asterfort/utpvlg.h"
#include "asterfort/vecma.h"
#include "blas/dcopy.h"
!
    type(te0047_dscr), intent(in) :: for_discret
    integer(kind=8), intent(out) :: iret
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: imat, ivarim, jdc, irep, iadzi, iazk24, neq, ivarip, ii, ifono, icontp
    integer(kind=8) :: icontm, iretvc
    real(kind=8) :: r8bid, klv(78), raide(6), ulp(12), temper, temp1, temp2, klc(144), fl(12)
    character(len=8) :: k8bid
    character(len=24) :: messak(5)
!
!   loi bi-linéaire sur 6 composantes
!       nbparc : 3 paramètres par composante
!       nbvint : 1 variables internes par composantes
    integer(kind=8) :: nbparc, nbvint
    parameter(nbparc=3, nbvint=1*6)
!       nbpart : nombre paramètres total
!       valpar : valeur paramètres de la loi
!       nompar : nom des paramètres de la loi
    integer(kind=8) :: nbpart
    parameter(nbpart=nbparc*6)
    real(kind=8) :: valpar(nbpart), coeflo(6, nbparc), vardnl(nbvint)
    integer(kind=8) :: codret(nbpart)
    character(len=8) :: nompar(nbpart)
    aster_logical :: okdire(6)
    blas_int :: b_incx, b_incy, b_n
!   nbparc paramètres par composante
    data nompar/'KDEB_DX', 'KFIN_DX', 'FPRE_DX', &
        'KDEB_DY', 'KFIN_DY', 'FPRE_DY', &
        'KDEB_DZ', 'KFIN_DZ', 'FPRE_DZ', &
        'KDEB_RX', 'KFIN_RX', 'FPRE_RX', &
        'KDEB_RY', 'KFIN_RY', 'FPRE_RY', &
        'KDEB_RZ', 'KFIN_RZ', 'FPRE_RZ'/
!
! --------------------------------------------------------------------------------------------------
!
    iret = 0
!   récupération du matériau
    call jevech('PMATERC', 'L', imat)
!   variables a t-
    call jevech('PVARIMR', 'L', ivarim)
    call jevech('PCONTMR', 'L', icontm)
!   récupération des caractéristiques
    call jevech('PCADISK', 'L', jdc)
!
    call infdis('REPK', irep, r8bid, k8bid)
!   seulement en repere local : irep = 2
    if (irep .ne. 2) then
        messak(1) = for_discret%nomte
        messak(2) = 'NON_LINEAR'
        messak(3) = for_discret%type_comp
        messak(4) = for_discret%rela_comp
        call tecael(iadzi, iazk24)
        messak(5) = zk24(iazk24-1+3)
        call utmess('F', 'DISCRETS_5', nk=5, valk=messak)
    end if
!   récupère tous les paramètres
!   température : si 2 noeuds ==> moyenne
    call rcvarc(' ', 'TEMP', '+', 'RIGI', 1, &
                1, temp1, iretvc)
    temper = temp1
    if (for_discret%nno .eq. 2) then
        call rcvarc(' ', 'TEMP', '+', 'RIGI', 2, &
                    1, temp2, iretvc)
        temper = (temp1+temp2)*0.5d0
    end if
    valpar(:) = 0.0d0
    call rcvalb('FPG1', 1, 1, '+', zi(imat), &
                ' ', 'DIS_BILI_ELAS', 1, 'TEMP', [temper], &
                nbpart, nompar, valpar, codret, 0)
!   les caractéristiques sont toujours dans le repère local on fait seulement une copie
    b_n = to_blas_int(for_discret%nbt)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, zr(jdc), b_incx, klv, b_incy)
!   si un ddl n'est pas affecte d'un comportement non-linéaire
!   il est donc élastique dans cette direction. ==> dinonc
    raide(:) = 0.0d0
    coeflo(:, :) = -1.0d0
!   examen des codret, valpar. on affecte raide, les paramètres
    call dinonc(for_discret%nomte, codret, valpar, klv, raide, &
                nbparc, coeflo, okdire)
!   loi de comportement non-linéaire
    neq = for_discret%nno*for_discret%nc
    ulp(1:12) = for_discret%ulm(1:12)+for_discret%dul(1:12)
    vardnl(:) = 0.0d0
    call dinon4(neq, for_discret%ulm, for_discret%dul, ulp, for_discret%nno, &
                for_discret%nc, zr(ivarim), raide, nbparc, coeflo, &
                okdire, vardnl)
!   actualisation de la matrice quasi-tangente
    call diklvraid(for_discret%nomte, klv, raide)
!   actualisation de la matrice quasi-tangente
    if (for_discret%lMatr) then
        call jevech('PMATUUR', 'E', imat)
        if (for_discret%ndim .eq. 3) then
            call utpslg(for_discret%nno, for_discret%nc, for_discret%pgl, klv, zr(imat))
        else if (for_discret%ndim .eq. 2) then
            call ut2mlg(for_discret%nno, for_discret%nc, for_discret%pgl, klv, zr(imat))
        end if
    end if
    !
    if (for_discret%lVect .or. for_discret%lSigm) then
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
            do ii = 1, neq
                zr(icontp-1+ii) = fl(ii)+zr(icontm-1+ii)
            end do
        else if (for_discret%nno .eq. 2) then
            do ii = 1, for_discret%nc
                zr(icontp-1+ii) = -fl(ii)+zr(icontm-1+ii)
                zr(icontp-1+ii+for_discret%nc) = fl(ii+for_discret%nc)+zr(icontm-1+ii+for_discre&
                                                 &t%nc)
            end do
        end if
    end if
! calcul des forces nodales
    if (for_discret%lVect) then
        call jevech('PVECTUR', 'E', ifono)
! Attention aux signes des efforts sur le premier noeud pour MECA_DIS_TR_L et MECA_DIS_T_L
        if (for_discret%nno .eq. 1) then
            do ii = 1, neq
                fl(ii) = fl(ii)+zr(icontm-1+ii)
            end do
        else if (for_discret%nno .eq. 2) then
            do ii = 1, for_discret%nc
                fl(ii) = fl(ii)-zr(icontm-1+ii)
                fl(ii+for_discret%nc) = fl(ii+for_discret%nc)+zr(icontm-1+ii+for_discret%nc)
            end do
        end if
! forces nodales aux noeuds 1 et 2 (repère global)
        if (for_discret%nc .ne. 2) then
            call utpvlg(for_discret%nno, for_discret%nc, for_discret%pgl, fl, zr(ifono))
        else
            call ut2vlg(for_discret%nno, for_discret%nc, for_discret%pgl, fl, zr(ifono))
        end if
    end if
! mise à jour des variables internes
    if (for_discret%lVari) then
        call jevech('PVARIPR', 'E', ivarip)
        do ii = 1, nbvint
            zr(ivarip+ii-1) = vardnl(ii)
        end do
    end if
end subroutine
