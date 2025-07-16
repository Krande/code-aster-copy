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
! person_in_charge:jean-luc.flejou@edf.fr
!
subroutine dichoc_galet_elasnl(for_discret, iret)
!
! --------------------------------------------------------------------------------------------------
!
!     RELATION DE COMPORTEMENT "CHOC_ELAS_TRAC" : COMPORTEMENT DISCRET CHOC NON-LINEAIRE
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
#include "asterfort/diraidklv.h"
#include "asterfort/diklvraid.h"
#include "asterfort/infdis.h"
#include "asterfort/jevech.h"
#include "asterfort/pmavec.h"
#include "asterfort/rcvala.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
#include "asterfort/utpsgl.h"
#include "asterfort/utpslg.h"
#include "asterfort/utpvgl.h"
#include "asterfort/utpvlg.h"
#include "asterfort/vecma.h"
#include "blas/dcopy.h"
#include "asterfort/Behaviour_type.h"
!
    type(te0047_dscr), intent(in) :: for_discret
    integer(kind=8), intent(out) :: iret
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: jdc, irep, imat, ivarim, ii, igeom, ivarip
    integer(kind=8) :: ifono, imatri
    integer(kind=8) :: icontp, icontm, iadzi, iazk24
!
    real(kind=8) :: klc(for_discret%neq*for_discret%neq)
    real(kind=8) :: klv(for_discret%nbt), fl(for_discret%neq), raide(6)
    real(kind=8) :: r8bid
    character(len=16), pointer :: compor(:) => null()
    character(len=8) :: k8bid
    character(len=24) :: messak(6)
!
! --------------------------------------------------------------------------------------------------
    integer(kind=8), parameter :: nbres = 3
    real(kind=8) :: valres(nbres)
    integer(kind=8) :: codres(nbres)
    character(len=8) :: nomres(nbres)
    integer(kind=8) :: nbpar
    real(kind=8) :: valpar
    character(len=8) :: nompar
!
! --------------------------------------------------------------------------------------------------
!   Variables internes
    integer(kind=8), parameter :: nbvari = 1
    real(kind=8) :: varmo(nbvari), varpl(nbvari)
! --------------------------------------------------------------------------------------------------
    real(kind=8) :: xl(6), xd(3), ld0, ldm, ldp, forcem, forcep, raidemp, ktang
    blas_int :: b_incx, b_incy, b_n
! --------------------------------------------------------------------------------------------------
!
    iret = 0
!   Seulement en 3D, sur un segment, avec seulement de la translation
    if ((for_discret%nomte(1:12) .ne. 'MECA_DIS_T_L') .or. (for_discret%ndim .ne. 3) .or. &
        (for_discret%nno .ne. 2) .or. (for_discret%nc .ne. 3)) then
        call jevech('PCOMPOR', 'L', vk16=compor)
        messak(1) = for_discret%nomte
        messak(2) = for_discret%option
        messak(3) = compor(INCRELAS)
        messak(4) = compor(RELA_NAME)
        call tecael(iadzi, iazk24)
        messak(5) = zk24(iazk24-1+3)
        call utmess('F', 'DISCRETS_22', nk=5, valk=messak)
    end if
! --------------------------------------------------------------------------------------------------
!   Paramètres en entrée
    call jevech('PCADISK', 'L', jdc)
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imat)
!
    call infdis('REPK', irep, r8bid, k8bid)
!   irep = 1 = matrice en repère global ==> passer en local
    if (irep .eq. 1) then
        call utpsgl(for_discret%nno, for_discret%nc, for_discret%pgl, zr(jdc), klv)
    else
        b_n = to_blas_int(for_discret%nbt)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(jdc), b_incx, klv, b_incy)
    end if
!   Récupération des termes diagonaux : raide = klv(i,i)
    call diraidklv(for_discret%nomte, raide, klv)
!   Variables internes
    call jevech('PVARIMR', 'L', ivarim)
    do ii = 1, nbvari
        varmo(ii) = zr(ivarim+ii-1)
    end do
    varpl(:) = varmo(:)
! --------------------------------------------------------------------------------------------------
!   Coordonnees du discret dans le repère local
    xl(:) = 0.0
    call utpvgl(for_discret%nno, 3, for_discret%pgl, zr(igeom), xl)
! --------------------------------------------------------------------------------------------------
!   Caractéristiques du matériau
    nbpar = 0; nompar = ' '; valpar = 0.d0
    valres(:) = 0.0; nomres(:) = ' '
    nomres(1) = 'DIST_1'; nomres(2) = 'DIST_2'
    call rcvala(zi(imat), ' ', 'DIS_CHOC_ELAS', nbpar, nompar, &
                [valpar], 2, nomres, valres, codres, &
                0, nan='NON')
!   longueur du discret
    xd(1:3) = xl(1+for_discret%ndim:2*for_discret%ndim)-xl(1:for_discret%ndim)
    ld0 = xd(1)-valres(1)-valres(2)
! --------------------------------------------------------------------------------------------------
!   Instant '-' : Longueur du discret
    ldm = ld0+for_discret%ulm(1+for_discret%nc)-for_discret%ulm(1)
!   Instant '+' : Longueur du discret
    ldp = ldm+for_discret%dul(1+for_discret%nc)-for_discret%dul(1)
!   calcul des efforts
    forcem = 0.0; forcep = 0.0; raidemp = raide(1)
    varpl(1) = 0.0; ktang = 0.0
!   Contact si longueur du discret <=0
!   Si pas contact force=0 et ktang=0
!   Instant '-' : Contact ou pas
    if (ldm <= 0.0) then
        nomres(1) = 'FX'; nomres(2) = 'RIGIP'
        call rcvala(zi(imat), ' ', 'DIS_CHOC_ELAS', 1, 'DX', &
                    [-ldm], 2, nomres, valres, codres, &
                    1)
        forcem = -valres(1)
        raidemp = valres(2)
    end if
!   Instant '+' : Contact ou pas
    if (ldp <= 0.0) then
        nomres(1) = 'FX'; nomres(2) = 'RIGIP'
        call rcvala(zi(imat), ' ', 'DIS_CHOC_ELAS', 1, 'DX', &
                    [-ldp], 2, nomres, valres, codres, &
                    1)
        forcep = -valres(1)
        raidemp = valres(2)
        varpl(1) = 1.0
    end if
!   Contact en '-' ou en '+' : Raideur 'corde'
    if ((ldp <= 0.0) .or. (ldm <= 0.0)) then
        if (abs(ldp-ldm) <= r8prem()) then
            ktang = raidemp
        else
            ktang = abs((forcep-forcem)/(ldp-ldm))
        end if
    end if
!   Actualisation de la raideur
    raide(1) = ktang
!
! --------------------------------------------------------------------------------------------------
!   Actualisation de la matrice tangente : klv(i,i) = raide(i)
    call diklvraid(for_discret%nomte, klv, raide)
!   Actualisation de la matrice quasi-tangente
    if (for_discret%rigi) then
        call jevech('PMATUUR', 'E', imatri)
        call utpslg(for_discret%nno, for_discret%nc, for_discret%pgl, klv, zr(imatri))
    end if
!
!   Efforts généralisés, Forces nodales
!       Sortie : efforts généralisés
!       Calcul : Forces nodales  (Mise à jour après)
    if (for_discret%resi) then
        call jevech('PCONTMR', 'L', icontm)
        call jevech('PCONTPR', 'E', icontp)
!       demi-matrice klv transformée en matrice pleine klc
        call vecma(klv, for_discret%nbt, klc, for_discret%neq)
!       calcul de fl = klc.dul (incrément d'effort)
        call pmavec('ZERO', for_discret%neq, klc, for_discret%dul, fl)
!       Efforts généralisés aux noeuds dans le repère local
!       Changement du signe des efforts sur le premier noeud pour les MECA_DIS_T_L
        do ii = 1, for_discret%nc
            zr(icontp-1+ii) = -fl(ii)+zr(icontm-1+ii)
            zr(icontp-1+ii+for_discret%nc) = fl(ii+for_discret%nc)+zr(icontm-1+ii+for_discret%nc&
                                             &)
            fl(ii) = fl(ii)-zr(icontm-1+ii)
            fl(ii+for_discret%nc) = fl(ii+for_discret%nc)+zr(icontm-1+ii+for_discret%nc)
        end do
!       Surcharge par l'effort 'réel' dans le discret
        zr(icontp-1+1) = forcep
        zr(icontp-1+1+for_discret%nc) = forcep
        fl(1) = -forcep
        fl(1+for_discret%nc) = forcep
!       Forces nodales aux noeuds dans le repère global
        call jevech('PVECTUR', 'E', ifono)
        call utpvlg(for_discret%nno, for_discret%nc, for_discret%pgl, fl, zr(ifono))
!       Variables internes : mise à jour
        call jevech('PVARIPR', 'E', ivarip)
        do ii = 1, nbvari
            zr(ivarip+ii-1) = varpl(ii)
            zr(ivarip+ii-1+nbvari) = varpl(ii)
        end do
    end if
!
end subroutine
