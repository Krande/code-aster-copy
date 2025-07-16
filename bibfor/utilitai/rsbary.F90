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
subroutine rsbary(lr8, nr8, tous, lexi, x, i1, i2, iposit, prec, crit)

    use searchlist_module, only: allUnique, almostEqual, getUnique
    implicit none
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/utmess.h"
!
!     ARGUMENTS:
!     ----------
    real(kind=8), intent(in) :: x, prec, lr8(*)
    integer(kind=8), intent(in) :: nr8
    aster_logical, intent(in) :: tous, lexi(*)
    integer(kind=8), intent(out) :: i1, i2, iposit
    character(len=8), intent(in) :: crit
!
! ----------------------------------------------------------------------
!     BUT:
!      TROUVER DANS UNE LISTE DE R8 QUELS SONT LES 2 REELS LES PLUS
!      PROCHES DU REEL X DONNE. (POUR FAIRE UN BARYCENTRE)
!     (ON PEUT NE PAS PRENDRE EN COMPTE TOUS LES REELS DE LA LISTE GRACE
!       A L'ARGUMENT LEXI)
!     IN:
!     NR8    : NOMBRE DE REELS DANS LA LISTE LR8.
!     LR8    : LISTE DE REELS (PAS FORCEMENT ORDONNEE).
!     TOUS   : INDIQUE QUE TOUS LES REELS DE LA LISTE SONT A CONSIDERER.
!     LEXI   : INDIQUE QUELS SONT LES REELS A CONSIDERER (SI TOUS=FALSE)
!              SI TOUS=.TRUE. CET ARGUMENT EST INUTILISE.
!      X     : REEL DONT ON CHERCHE LES COORDONEES BARYCENTRIQUES.
!
!     OUT:
!     I1,I2  : INDICES DES 2 REELS DE LA LISTE QUI "ENCADRENT" X
!              (LR8(I1) =< LR8(I2))
!              (EVENTUELLEMENT I1 PEUT ETRE EGAL A I2)
!     IPOSIT : CODE LA POSITION DE X PAR RAPPORT A LR8(I1) ET LR8(I2)
!         IPOSIT=0  -->  LR8(I1)  =<   X  =<   LR8(I2)
!         IPOSIT=1  -->  LR8(I1) =<  LR8(I2) =< X  (PROL_DR)
!         IPOSIT=-1 -->  X =< LR8(I1) =<  LR8(I2)  (PROL_GA)
!         IPOSIT=-2 -->  ERREUR : LA LISTE DE REELS EST VIDE.
!
!
! ----------------------------------------------------------------------
    integer(kind=8) :: imax, imin, ipp, ip, is, iss, ns
    real(kind=8) :: xmax, xmin
    real(kind=8), pointer :: v_diff(:) => null()
    aster_logical, pointer :: v_active(:) => null()

!-----------------------------------------------------------------------
    AS_ALLOCATE(vl=v_active, size=nr8)
!
    v_active = ASTER_FALSE
    if (tous) then
        v_active = ASTER_TRUE
    else
        v_active = lexi(1:nr8)
    end if

! DEB-------------------------------------------------------------------
!
!     --------XI-----XPP--XP-------X---XS-------XSS-------XJ-->
!
!     ON APPELLE : XP : LE REEL PRECEDENT X DANS LA LISTE
!                XPP: LE REEL PRECEDENT XP DANS LA LISTE
!                XS : LE REEL SUIVANT  X DANS LA LISTE
!                XSS: LE REEL SUIVANT  XS DANS LA LISTE
!               XMAX: LE REEL MAX DE LA LISTE
!               XMIN: LE REEL MIN DE LA LISTE
!
    i1 = 0
    i2 = 0
    ns = count(v_active)

    if (ns .eq. 0) then
        ! CAS LISTE VIDE
        iposit = -2
    else if (ns .eq. 1) then
        ! CAS EXCEPTION LISTE AVEC UN SEUL INSTANT
        i1 = MAXLOC(lr8(1:nr8), dim=1, mask=v_active)
        i2 = i1
        xmax = MAXVAL(lr8(1:nr8), dim=1, mask=v_active)
        xmin = xmax
        if (almostEqual(x, xmax, prec, crit)) then
            iposit = 0
        else if (x .gt. xmax) then
            iposit = 1
        else if (x .lt. xmin) then
            iposit = -1
        else
            ASSERT(.false.)
        end if
    else
        i1 = getUnique(x, lr8(1:nr8), prec, crit, v_active)
        if (i1 .ne. 0) then
            ! CAS EXCEPTION OU ON DEMANDE UN INSTANT EXISTANT DANS LA LISTE
            iposit = 0
            i2 = i1
        else
            ! INTERPOLATION
            !
            !
            ! VERIFICATION UNICITE DE CHAQUE INSTANT
            if (.not. allUnique(lr8(1:nr8), prec, crit, v_active)) then
                call utmess('F', 'UTILITAI_44', sk=crit, sr=prec)
            end if
            !
            ! RECHERCHE DE XMAX ET XMIN:
            xmax = MAXVAL(lr8(1:nr8), dim=1, mask=v_active)
            xmin = MINVAL(lr8(1:nr8), dim=1, mask=v_active)
            !
            if ((x .ge. xmin) .and. (x .le. xmax)) then
                ! 1ER CAS X EST INCLU DANS L'INTERVALLE DE LA LISTE:
                AS_ALLOCATE(vr=v_diff, size=nr8)
                v_diff = lr8(1:nr8)-x
                ip = MAXLOC(v_diff, dim=1, mask=(v_diff .lt. 0.d0 .and. v_active))
                is = MINLOC(v_diff, dim=1, mask=(v_diff .gt. 0.d0 .and. v_active))
                AS_DEALLOCATE(vr=v_diff)
                i1 = ip
                i2 = is
                iposit = 0
            else if (x .gt. xmax) then
                ! 2EME CAS X EST A DROITE DE L'INTERVALLE DE LA LISTE:
                imax = MAXLOC(lr8(1:nr8), dim=1, mask=v_active)
                ip = imax
                v_active(imax) = ASTER_FALSE
                ipp = MAXLOC(lr8(1:nr8), dim=1, mask=v_active)
                i1 = ipp
                i2 = ip
                iposit = 1
            else if (x .lt. xmin) then
                ! 3EME CAS X EST A GAUCHE DE L'INTERVALLE DE LA LISTE:
                imin = MINLOC(lr8(1:nr8), dim=1, mask=v_active)
                is = imin
                v_active(imin) = ASTER_FALSE
                iss = MINLOC(lr8(1:nr8), dim=1, mask=v_active)
                i1 = is
                i2 = iss
                iposit = -1
            else
                ASSERT(.false.)
            end if
        end if
    end if
!
!
    AS_DEALLOCATE(vl=v_active)

end subroutine
