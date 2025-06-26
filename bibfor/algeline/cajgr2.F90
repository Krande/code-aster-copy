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
subroutine cajgr2(igrap, vr, cocaj1, cocaj2)
    implicit none
! CALCUL DES COEFFICIENTS ADIMENSIONNELS DE FORCE D'AMORTISSEMENT
! GRAPPE2
!-----------------------------------------------------------------------
!  IN : IGRAP  : INDICE CARACTERISTIQUE DE LA CONFIGURATION
!                EXPERIMENTALE DE REFERENCE
!  IN : VR     : VITESSE REDUITE
! OUT : COCAJ1 : COEFFICIENT ADIMENSIONNEL DE FORCE D'AMORTISSEMENT
!                POUR UN MOUVEMENT DE TRANSLATION
! OUT : COCAJ2 : COEFFICIENT ADIMENSIONNEL DE FORCE D'AMORTISSEMENT
!                POUR UN MOUVEMENT DE ROTATION
!-----------------------------------------------------------------------
!     UN COMMON AJOUTE POUR RESORBER UNE GLUTE ANTIQUE (VOIR HISTOR):
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/ulisop.h"
#include "asterfort/ulopen.h"
#include "asterfort/wkvect.h"
    character(len=8) :: typflu
    common/kop144/typflu
!
    integer(kind=8) :: igrap
    real(kind=8) :: vr, cocaj1, cocaj2
!
    integer(kind=8) :: ncamax, nbomax, unit, nbloc, iflag
    real(kind=8) :: coeca1(10, 20, 11), coeca2(10, 20, 11)
    real(kind=8) :: boca1(10, 20), boca2(10, 20), borne1(10, 20)
    real(kind=8) :: zero, coef1(10, 20, 11), coef2(10, 20, 11), borne2(10, 20)
    character(len=16) :: k16nom
    character(len=24) :: nom1, nom2
    save borne1, coeca1, borne2, coeca2
! ----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iret, iunit, j, k, l, nb1
    integer(kind=8) :: nb2, nbmax
    real(kind=8) :: vr2, vr3
!-----------------------------------------------------------------------
    call jemarq()
    ncamax = 11
    nbomax = 20
    nbmax = 10
    zero = 0.0d0
!
    nom1 = '&&CAJGR2.FLAG'
    nom2 = typflu//'.UNIT_GRAPPES'
!
! --- ON TESTE L'EXISTENCE DU VECTEUR DES COEFFICIENTS
!     ===============================================
    call jeexin(nom1, iret)
    if (iret .eq. 0) then
!
! --- LECTURE DU FICHIER DE DONNEES
!     =============================
        call jeveuo(nom2, 'L', iunit)
        unit = zi(iunit-1+1)
        k16nom = ' '
        if (ulisop(unit, k16nom) .eq. 0) then
            call ulopen(unit, ' ', ' ', 'NEW', 'O')
        end if
!
! ---    BLOC D'INITIALISATION
        do i = 1, nbmax
            do j = 1, nbomax
                boca1(i, j) = zero
                boca2(i, j) = zero
                borne1(i, j) = zero
                borne2(i, j) = zero
                do k = 1, ncamax
                    coeca1(i, j, k) = zero
                    coeca2(i, j, k) = zero
                    coef1(i, j, k) = zero
                    coef2(i, j, k) = zero
                end do
            end do
        end do
!
        read (unit, *) nbloc
        do l = 1, nbloc
            read (unit, *) nb1
            if (nb1 .ne. 0) then
                read (unit, *) (boca1(l, i), i=1, nb1)
            end if
            do i = 1, nb1+1
                read (unit, *) (coef1(l, i, j), j=1, ncamax)
            end do
            read (unit, *) nb2
            if (nb2 .ne. 0) then
                read (unit, *) (boca2(l, i), i=1, nb2)
            end if
            do i = 1, nb2+1
                read (unit, *) (coef2(l, i, j), j=1, ncamax)
            end do
            read (unit, *)
            do i = 1, nbmax
                do j = 1, nbomax
                    borne1(i, j) = boca1(i, j)
                    borne2(i, j) = boca2(i, j)
                    do k = 1, ncamax
                        coeca1(i, j, k) = coef1(i, j, k)
                        coeca2(i, j, k) = coef2(i, j, k)
                    end do
                end do
            end do
            call jedetr(nom1)
            call wkvect(nom1, 'V V I', 1, iflag)
            zi(iflag-1+1) = 1
        end do
    end if
!
!-----1.CONFIG. ECOULEMENT ASCENDANT TIGE DE COMMANDE CENTREE
!
    if (igrap .eq. 1) then
!
        if (vr .lt. borne1(1, 1)) then
            cocaj1 = coeca1(1, 1, 8)+coeca1(1, 1, 9)*vr
        else
            cocaj1 = coeca1(1, 2, 8)+coeca1(1, 2, 9)*vr
        end if
!
        cocaj2 = coeca2(1, 1, 8)+coeca2(1, 1, 9)*vr
!
!-----2.CONFIG. ECOULEMENT ASCENDANT TIGE DE COMMANDE EXCENTREE
!
    else if (igrap .eq. 2) then
!
        cocaj1 = coeca1(2, 1, 8)+coeca1(2, 1, 9)*vr
        cocaj2 = coeca2(2, 1, 8)+coeca2(2, 1, 9)*vr
!
!-----3.CONFIG. ECOULEMENT DESCENDANT TIGE DE COMMANDE CENTREE
!
    else if (igrap .eq. 3) then
!
        if (vr .lt. borne1(3, 1)) then
            vr2 = vr*vr
            vr3 = vr2*vr
            cocaj1 = coeca1(3, 1, 8)+coeca1(3, 1, 9)*vr+coeca1(3, 1, 10)*vr2+coeca1(3, 1, 11)*vr3
        else if (vr .lt. borne1(3, 2)) then
            cocaj1 = coeca1(3, 2, 8)
        else
            cocaj1 = coeca1(3, 3, 8)+coeca1(3, 3, 9)*vr
        end if
!
        if (vr .lt. borne2(3, 1)) then
            cocaj2 = coeca2(3, 1, 8)+coeca2(3, 1, 9)*vr
        else if (vr .lt. borne2(3, 2)) then
            vr2 = vr*vr
            cocaj2 = coeca2(3, 2, 8)+coeca2(3, 2, 9)*vr+coeca2(3, 2, 10)*vr2
        else
            cocaj2 = coeca2(3, 3, 8)+coeca2(3, 3, 9)*vr
        end if
!
!-----4.CONFIG. ECOULEMENT DESCENDANT TIGE DE COMMANDE EXCENTREE
!
    else
!
        if (vr .lt. borne1(4, 1)) then
            vr2 = vr*vr
            vr3 = vr2*vr
            cocaj1 = coeca1(4, 1, 8)+coeca1(4, 1, 9)*vr+coeca1(4, 1, 10)*vr2+coeca1(4, 1, 11)*vr3
        else if (vr .lt. borne1(4, 2)) then
            vr2 = vr*vr
            vr3 = vr2*vr
            cocaj1 = coeca1(4, 2, 8)+coeca1(4, 2, 9)*vr+coeca1(4, 2, 10)*vr2+coeca1(4, 2, 11)*vr3
        else
            cocaj1 = coeca1(4, 3, 8)+coeca1(4, 3, 9)*vr
        end if
!
        if (vr .lt. borne2(4, 1)) then
            cocaj2 = coeca2(4, 1, 8)+coeca2(4, 1, 9)*vr
        else if (vr .lt. borne2(4, 2)) then
            vr2 = vr*vr
            cocaj2 = coeca2(4, 2, 8)+coeca2(4, 2, 9)*vr+coeca2(4, 2, 10)*vr2
        else
            cocaj2 = coeca2(4, 3, 8)+coeca2(4, 3, 9)*vr
        end if
!
    end if
!
!     FERMETURE DU FICHIER
    if (iret .eq. 0) call ulopen(-unit, ' ', ' ', ' ', ' ')
    call jedema()
end subroutine
