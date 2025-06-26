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
subroutine te0357(optioz, nomtez)
!
! --------------------------------------------------------------------------------------------------
!
!         COMPORTEMENT LINÉAIRE ET NON-LINÉAIRE POUR LES DISCRETS
!
! --------------------------------------------------------------------------------------------------
!
! Éléments concernés
!     MECA_DIS_TR_L     : sur une maille a 2 noeuds
!     MECA_DIS_T_L      : sur une maille a 2 noeuds
!     MECA_DIS_TR_N     : sur une maille a 1 noeud
!     MECA_DIS_T_N      : sur une maille a 1 noeud
!     MECA_2D_DIS_TR_L  : sur une maille a 2 noeuds
!     MECA_2D_DIS_T_L   : sur une maille a 2 noeuds
!     MECA_2D_DIS_TR_N  : sur une maille a 1 noeud
!     MECA_2D_DIS_T_N   : sur une maille a 1 noeud
!
! --------------------------------------------------------------------------------------------------
!
!   IN
!       optioz : nom de l'option a calculer
!       nomtez : nom du type_element
!
! --------------------------------------------------------------------------------------------------
!
    use Behaviour_module, only: behaviourOption
    use te0047_type
!
    implicit none
    character(len=*) :: optioz, nomtez
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/infdis.h"
#include "asterfort/jevech.h"
#include "asterfort/matrot.h"
#include "asterfort/tecach.h"
#include "asterfort/tecael.h"
#include "asterfort/ut2mgl.h"
#include "asterfort/ut2vlg.h"
#include "asterfort/utmess.h"
#include "asterfort/utpsgl.h"
#include "asterfort/utpvlg.h"
#include "asterfort/vecma.h"
#include "blas/dcopy.h"
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: lorien
    integer(kind=8) :: iadzi, iazk24, ibid, infodi, codret, ii, jj
    !
    integer(kind=8) :: neq, ifono, jdc, irep
    real(kind=8) :: r8bid, klv(78), klc(12, 12), fl(12), dulth(12)
    !
    character(len=8) :: k8bid
    character(len=24) :: messak(5)
    !
    type(te0047_dscr) :: DD
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!&<
    DD%option = optioz
    DD%nomte = nomtez
! on vérifie que les caractéristiques ont été affectées
! le code du discret
    call infdis('CODE', ibid, r8bid, DD%nomte)
! le code stocké dans la carte
    call infdis('TYDI', infodi, r8bid, k8bid)
    if (infodi .ne. ibid) then
        call utmess('F+', 'DISCRETS_25', sk=DD%nomte)
        call infdis('DUMP', ibid, r8bid, 'F+')
    end if
! discret de type raideur
    call infdis('DISK', infodi, r8bid, k8bid)
    if (infodi .eq. 0) then
        call utmess('A+', 'DISCRETS_27', sk=DD%nomte)
        call infdis('DUMP', ibid, r8bid, 'A+')
    end if
! Discrets symétriques ou pas
    call infdis('SYMK', DD%syme, r8bid, k8bid)
    if (DD%syme .ne. 1) then
        call utmess('A+', 'DISCRETS_32', 2, valk=[optioz, nomtez])
        call infdis('DUMP', ibid, r8bid, 'A+')
    end if
!
! Récupérations de plein d'informations sur les discrets
    call getDiscretInformations(DD)
    !
    ASSERT((DD%ndim .eq. 2) .or. (DD%ndim .eq. 3))
    !
! récupération des orientations (angles nautiques)
! orientation de l'élément et déplacements dans les repères g et l
    call tecach('ONO', 'PCAORIE', 'L', codret, iad=lorien)
    if (codret .ne. 0) then
        messak(1) = DD%nomte
        messak(2) = DD%option
        messak(3) = DD%type_comp
        messak(4) = DD%rela_comp
        call tecael(iadzi, iazk24)
        messak(5) = zk24(iazk24-1+3)
        call utmess('F', 'DISCRETS_6', nk=5, valk=messak)
    end if
    !
! matrice pgl de passage repère global -> repère local
    call matrot(zr(lorien), DD%pgl)
    call jevech('PCADISK', 'L', jdc)
    !
    call infdis('REPK', irep, r8bid, k8bid)
! absolu vers local ?
! irep = 1 . la matrice est dans le repère global ==> passer en local
    if (irep .eq. 1) then
        if (DD%ndim .eq. 3) then
            call utpsgl(DD%nno, DD%nc, DD%pgl, zr(jdc), klv)
        else if (DD%ndim .eq. 2) then
            call ut2mgl(DD%nno, DD%nc, DD%pgl, zr(jdc), klv)
        end if
    else
        b_n = to_blas_int(DD%nbt)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(jdc), b_incx, klv, b_incy)
    end if
! demi-matrice klv transformée en matrice pleine klc
    neq = DD%nno*DD%nc
    call vecma(klv, DD%nbt, klc, neq)
! calcul de fl = klc.dul (incrément d'effort)
!   Correction thermique : Dilatation = alpha*(T+ - Tref)*xl
    dulth(:) = 0.0
    dulth(1) = -DD%Dilatation
    dulth(DD%nc+1) = -dulth(1)
    fl(:) = 0.0
    do ii = 1, DD%nc
        do jj = 1, DD% nc
            fl(ii) = fl(ii) + klc(ii, jj )*dulth(jj)
            fl(ii+DD%nc) = fl(ii+DD%nc) + klc(ii+DD%nc,jj+DD%nc)*dulth(jj+DD%nc)
        end do
    end do
    !
    call jevech('PVECTUR', 'E', ifono)
! Passage du repere local au repere glocal
    if (DD%nc .ne. 2) then
        call utpvlg(DD%nno, DD%nc, DD%pgl, fl, zr(ifono))
    else
        call ut2vlg(DD%nno, DD%nc, DD%pgl, fl, zr(ifono))
    end if
!&>
end subroutine
