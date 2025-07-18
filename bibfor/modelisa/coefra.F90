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
subroutine coefra(ipas, ires, x, xsi0, ck)
!   CALCUL DU COEFFICIENT DE RAIDEUR AJOUTEE CK EN FONCTION DE LA
!   VITESSE REDUITE  (FAISCEAU DE TUBES SOUS ECOULEMENT TRANSVERSE)
!-----------------------------------------------------------------------
    implicit none
!  IN   : IPAS      : TYPE DE PAS
!  IN   : IRES      : TYPE DE RESEAU DU POINT COURANT
!  IN   : X         : VITESSE REDUITE
!  IN   : NBORCK :
!  IN   : BORNCK  :
!  IN   : COEFCK   :
!  OUT  : CK        : COEFFICIENT DE MASSE AJOUTEE
!-----------------------------------------------------------------------
!
#include "jeveux.h"
#include "asterfort/ckatrc.h"
#include "asterfort/coefrl.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: ipas, ires
    integer(kind=8) :: jborne, jcoeff, jvired
    real(kind=8) :: ck, xsi0
!
    integer(kind=8) :: nborck, nckmax, iret
    real(kind=8) :: zero, bornck(20), coefck(20, 11)
    character(len=24) :: nom1, nom2, nom3
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ipas1, ires1, j, k
    real(kind=8) :: x
!-----------------------------------------------------------------------
    call jemarq()
!
    x = dble(abs(x))
    nckmax = 11
    zero = 0.0d0
!
    nom1 = '&&COEFRA.CKI'
    nom2 = '&&COEFRA.CKR1'
    nom3 = '&&COEFRA.CKR2'
!
    if (ires .eq. 0) then
        ck = zero
        goto 1000
    end if
!
! --- ON TESTE L'EXISTENCE DU VECTEUR DE COEFFICIENTS
!     POUR LA CORRELATION RELATIVE A IPAS ET IRES
!     ===============================================
    call jeexin(nom2, iret)
    if (iret .eq. 0) then
!
! --- LECTURE DU FICHIER DE DONNEES
!     =============================
        call coefrl(nom1, nom2, nom3, nckmax, ipas, &
                    ires, bornck, nborck, coefck, ipas1, &
                    ires1)
    else
        call jeveuo(nom1, 'L', jborne)
        call jeveuo(nom2, 'L', jcoeff)
        call jeveuo(nom3, 'L', jvired)
        ipas1 = zi(jborne-1+1)
        ires1 = zi(jborne-1+2)
        nborck = zi(jborne-1+3)
        if (ipas1 .eq. ipas .and. ires1 .eq. ires) then
            k = 1
            do i = 1, nborck
                bornck(i) = zr(jcoeff+i-1)
                do j = 1, nckmax
                    coefck(i, j) = zr(jcoeff+nborck+k-1)
                    k = k+1
                end do
            end do
        else
            call jedetr(nom1)
            call jedetr(nom2)
            call jedetr(nom3)
            call coefrl(nom1, nom2, nom3, nckmax, ipas, &
                        ires, bornck, nborck, coefck, ipas1, &
                        ires1)
        end if
    end if
    if (ipas1 .ne. ipas .or. ires1 .ne. ires) then
        call utmess('F', 'MODELISA4_31')
    end if
!
! **********************************************************************
! ***                  FAISCEAU EN PAS CARRE LIGNE                   ***
! **********************************************************************
!
    if (ipas .eq. 1) then
!
        if (ires .ge. 1 .and. ires .le. 1000) then
!
            if (x .lt. bornck(1)) then
                ck = 0.d0
            else
                if (x .lt. bornck(nborck)) then
                    do i = 2, nborck
                        if (x .ge. bornck(i-1) .and. x .lt. bornck(i)) then
                            ck = coefck(i-1, 1)/(x*x*x*x*x*x*x)+coefck(i-1, 2)/(x*x*x*x*x*x)+co&
                                 &efck(i-1, 3)/(x*x*x*x*x)+coefck(i-1, 4)/(x*x*x*x)+coefck(i-1&
                                 &, 5)/(x*x*x)+coefck(i-1, 6)/(x*x)+coefck(i-1, 7)/(x)+coefck&
                                 &(i-1, 8)+coefck(i-1, 9)*(x)+coefck(i-1, 10)*(x*x)+coefck(i-1&
                                 &, 11)*(x*x*x)
                            goto 140
                        end if
                    end do
140                 continue
                else
                    ck = coefck(nborck, 1)/(x*x*x*x*x*x*x)+coefck(nborck, 2)/(x*x*x*x*x*x)+coe&
                         &fck(nborck, 3)/(x*x*x*x*x)+coefck(nborck, 4)/(x*x*x*x)+coefck(nborck&
                         &, 5)/(x*x*x)+coefck(nborck, 6)/(x*x)+coefck(nborck, 7)/(x)+coefck(&
                         &nborck, 8)+coefck(nborck, 9)*(x)+coefck(nborck, 10)*(x*x)+coefck(nb&
                         &orck, 11)*(x*x*x)
                end if
            end if
!
! --- CELLULE DE TUBES VIBRANTS EN DEBUT DE FAISCEAU VISCACHE1.
!
        else if (ires .eq. 1001) then
!
            if (x .lt. bornck(1)) then
                ck = 0.d0
            else if (x .lt. bornck(2)) then
                ck = coefck(1, 8)+coefck(1, 9)*x
            else if (x .lt. bornck(3)) then
                ck = coefck(2, 8)+coefck(2, 9)*x
            else if (x .lt. bornck(4)) then
                ck = coefck(3, 8)+coefck(3, 7)/x+coefck(3, 6)/(x*x)+coefck(3, 5)/(x*x*x)
            else
                ck = coefck(4, 8)
            end if
!
! --- CELLULE DE TUBES VIBRANTS EN MILIEU DE FAISCEAU CLOTAIRE.
!     (PROFIL DE VITESSE REEL)
!
        else if (ires .eq. 1002) then
!
            if (x .lt. bornck(1)) then
                ck = 0.d0
            else if (x .lt. bornck(2)) then
                ck = coefck(1, 8)+coefck(1, 9)*x
            else if (x .lt. bornck(3)) then
                ck = coefck(2, 8)+coefck(2, 7)/x+coefck(2, 6)/(x*x)+coefck(2, 5)/(x*x*x)
            else if (x .lt. bornck(4)) then
                ck = coefck(3, 8)+coefck(3, 7)/x+coefck(3, 6)/(x*x)+coefck(3, 5)/(x*x*x)
            else if (x .lt. bornck(5)) then
                ck = coefck(4, 8)+coefck(4, 9)*x
            else
                ck = coefck(5, 8)
            end if
!
! --- CELLULE DE TUBES VIBRANTS EN MILIEU DE FAISCEAU CLOTAIRE.
!     (PROFIL DE VITESSE UNIFORME)
!
        else if (ires .eq. 1003) then
!
            if (x .lt. bornck(1)) then
                ck = 0.d0
            else if (x .lt. bornck(2)) then
                ck = coefck(1, 8)+coefck(1, 9)*x
            else if (x .lt. bornck(3)) then
                ck = coefck(2, 8)+coefck(2, 7)/x+coefck(2, 6)/(x*x)+coefck(2, 5)/(x*x*x)
            else if (x .lt. bornck(4)) then
                ck = coefck(3, 8)+coefck(3, 7)/x+coefck(3, 6)/(x*x)+coefck(3, 5)/(x*x*x)
            else if (x .lt. bornck(5)) then
                ck = coefck(4, 8)+coefck(4, 9)*x
            else
                ck = coefck(5, 8)
            end if
!
! --- TUBE UNIQUE VIBRANT EN MILIEU DE FAISCEAU RIGIDE VISCACHE1.
!
        else if (ires .eq. 1004) then
!
            if (x .lt. bornck(1)) then
                ck = 0.d0
            else if (x .lt. bornck(2)) then
                ck = coefck(1, 8)+coefck(1, 9)*x
            else if (x .lt. bornck(3)) then
                ck = coefck(2, 8)+coefck(2, 7)/x+coefck(2, 6)/(x*x)
            else if (x .lt. bornck(4)) then
                ck = coefck(3, 8)+coefck(3, 9)*x
            else
                ck = coefck(4, 8)
            end if
!
! --- TUBE UNIQUE VIBRANT EN DEBUT DE FAISCEAU RIGIDE VISCACHE1.
!
        else if (ires .eq. 1005) then
!
            if (x .lt. bornck(1)) then
                ck = 0.d0
            else if (x .lt. bornck(2)) then
                ck = coefck(1, 8)+coefck(1, 9)*x
            else if (x .lt. bornck(3)) then
                ck = coefck(2, 8)+coefck(2, 7)/x+coefck(2, 6)/(x*x)+coefck(2, 5)/(x*x*x)
            else if (x .lt. bornck(4)) then
                ck = coefck(3, 8)+coefck(3, 7)/x+coefck(3, 6)/(x*x)+coefck(3, 5)/(x*x*x)
            else
                ck = coefck(4, 8)+coefck(4, 7)/x+coefck(4, 6)/(x*x)+coefck(4, 5)/(x*x*x)
            end if
!
! --- TUBE ROMPU
!
        else if (ires .eq. 1006) then
!
            if (x .le. bornck(1)) then
                ck = 0.0d0
            else
                call ckatrc(x, xsi0, coefck, ck)
            end if
!
! --- TUBE UNIQUE VIBRANT EN MILIEU DE FAISCEAU RIGIDE TANAKA.
!
        else if (ires .eq. 1007) then
!
            if (x .lt. bornck(1)) then
                ck = 0.d0
            else
                ck = coefck(1, 8)+coefck(1, 7)/x+coefck(1, 6)/(x*x)+coefck(1, 5)/(x*x*x)
            end if
!
! --- TUBE UNIQUE VIBRANT EN MILIEU DE FAISCEAU RIGIDE DIVA EAU.
!
        else if (ires .eq. 1008) then
!
            if (x .lt. bornck(1)) then
                ck = 0.d0
            else if (x .lt. bornck(2)) then
                ck = coefck(1, 8)+coefck(1, 9)*x
            else if (x .lt. bornck(3)) then
                ck = coefck(2, 8)+coefck(2, 7)/x
            else if (x .lt. bornck(4)) then
                ck = coefck(3, 8)+coefck(3, 7)/x+coefck(3, 6)/(x*x)+coefck(3, 5)/(x*x*x)
            else
                ck = coefck(4, 8)+coefck(4, 7)/x+coefck(4, 6)/(x*x)+coefck(4, 5)/(x*x*x)
            end if
!
!
! --- COEFFICIENT VISCACHE 2 CFD 90%
!
        else if (ires .eq. 1101) then
!
            if (x .lt. bornck(1)) then
                ck = 0.d0
            else if (x .lt. bornck(2)) then
                ck = coefck(1, 8)+coefck(1, 9)*x
            else if (x .lt. bornck(3)) then
                ck = coefck(2, 8)+coefck(2, 7)/x+coefck(2, 6)/(x*x)
            else
                ck = coefck(3, 8)+coefck(3, 7)/x+coefck(3, 6)/(x*x)+coefck(3, 5)/(x*x*x)
            end if
!
!
! --- COEFFICIENT VISCACHE 2 CFD 85 %
!
        else if (ires .eq. 1102) then
!
            if (x .lt. bornck(1)) then
                ck = 0.d0
            else if (x .lt. bornck(2)) then
                ck = coefck(1, 8)+coefck(1, 9)*x
            else
                ck = coefck(2, 8)+coefck(2, 7)/x+coefck(2, 6)/(x*x)+coefck(2, 5)/(x*x*x)+coef&
                     &ck(2, 4)/(x*x*x*x)+coefck(2, 3)/(x*x*x*x*x)
            end if
!
!
! --- COEFFICIENT VISCACHE 2 CFD 80 %
!
        else if (ires .eq. 1103) then
!
            if (x .lt. bornck(1)) then
                ck = 0.d0
            else if (x .lt. bornck(2)) then
                ck = coefck(1, 8)+coefck(1, 9)*x
            else if (x .lt. bornck(3)) then
                ck = coefck(2, 8)+coefck(2, 7)/x+coefck(2, 6)/(x*x)+coefck(2, 5)/(x*x*x)
            else
                ck = coefck(3, 8)+coefck(3, 7)/x+coefck(3, 6)/(x*x)+coefck(3, 5)/(x*x*x)+coef&
                     &ck(3, 4)/(x*x*x*x)
            end if
!
!
! --- COEFFICIENT VISCACHE 2 CFD 50 %
!
        else if (ires .eq. 1104) then
!
            if (x .lt. bornck(1)) then
                ck = 0.d0
            else if (x .lt. bornck(2)) then
                ck = coefck(1, 8)+coefck(1, 9)*x
            else
                ck = coefck(2, 8)+coefck(2, 7)/x+coefck(2, 6)/(x*x)+coefck(2, 5)/(x*x*x)
            end if
!
!
! --- COEFFICIENT VISCACHE 2 CFD 20 %
!
        else if (ires .eq. 1105) then
!
            if (x .lt. bornck(1)) then
                ck = 0.d0
            else if (x .lt. bornck(2)) then
                ck = coefck(1, 8)+coefck(1, 9)*x
            else
                ck = coefck(2, 8)+coefck(2, 7)/x+coefck(2, 6)/(x*x)+coefck(2, 5)/(x*x*x)+coef&
                     &ck(2, 4)/(x*x*x*x)
            end if
!
!
! --- COEFFICIENT VISCACHE 2 CFD 10 %
!
        else if (ires .eq. 1106) then
!
            if (x .lt. bornck(1)) then
                ck = 0.d0
            else if (x .lt. bornck(2)) then
                ck = coefck(1, 8)+coefck(1, 9)*x
            else
                ck = coefck(2, 8)+coefck(2, 7)/x+coefck(2, 6)/(x*x)+coefck(2, 5)/(x*x*x)
            end if
!
!
! --- COEFFICIENT VISCACHE 2 TUM 90 %
!
        else if (ires .eq. 1201) then
!
            if (x .lt. bornck(1)) then
                ck = 0.d0
            else if (x .lt. bornck(2)) then
                ck = coefck(1, 8)+coefck(1, 9)*x
            else
                ck = coefck(2, 8)+coefck(2, 7)/x+coefck(2, 6)/(x*x)+coefck(2, 5)/(x*x*x)+coef&
                     &ck(2, 4)/(x*x*x*x)
            end if
!
!
! --- COEFFICIENT VISCACHE 2 TUM 86 %
!
        else if (ires .eq. 1202) then
!
            if (x .lt. bornck(1)) then
                ck = 0.d0
            else if (x .lt. bornck(2)) then
                ck = coefck(1, 8)+coefck(1, 9)*x
            else
                ck = coefck(2, 8)+coefck(2, 7)/x+coefck(2, 6)/(x*x)
            end if
!
!
! --- COEFFICIENT VISCACHE 2 TUM 82 %
!
        else if (ires .eq. 1203) then
!
            if (x .lt. bornck(1)) then
                ck = 0.d0
            else if (x .lt. bornck(2)) then
                ck = coefck(1, 8)+coefck(1, 9)*x
            else
                ck = coefck(2, 8)+coefck(2, 7)/x+coefck(2, 6)/(x*x)+coefck(2, 5)/(x*x*x)
            end if
!
!
! --- COEFFICIENT VISCACHE 2 TUM 50 %
!
        else if (ires .eq. 1204) then
!
            if (x .lt. bornck(1)) then
                ck = 0.d0
            else if (x .lt. bornck(2)) then
                ck = coefck(1, 8)+coefck(1, 9)*x
            else
                ck = coefck(2, 8)+coefck(2, 7)/x+coefck(2, 6)/(x*x)
            end if
!
!
! --- COEFFICIENT VISCACHE 2 TUM 20 %
!
        else if (ires .eq. 1205) then
!
            if (x .lt. bornck(1)) then
                ck = 0.d0
            else if (x .lt. bornck(2)) then
                ck = coefck(1, 8)+coefck(1, 9)*x
            else
                ck = coefck(2, 8)+coefck(2, 7)/x+coefck(2, 6)/(x*x)+coefck(2, 5)/(x*x*x)
            end if
!
!
! --- COEFFICIENT VISCACHE 2 TUM 10 %
!
        else if (ires .eq. 1206) then
!
            if (x .lt. bornck(1)) then
                ck = 0.d0
            else if (x .lt. bornck(2)) then
                ck = coefck(1, 8)+coefck(1, 9)*x
            else
                ck = coefck(2, 8)+coefck(2, 7)/x+coefck(2, 6)/(x*x)+coefck(2, 5)/(x*x*x)
            end if
!
        else
!
            call utmess('F', 'MODELISA4_32')
!
        end if
!
! **********************************************************************
! ***               FAISCEAU EN PAS TRIANGULAIRE LIGNE               ***
! **********************************************************************
!
    else if (ipas .eq. 2) then
!
        if (ires .ge. 1 .and. ires .le. 1000) then
!
            if (x .lt. bornck(1)) then
                ck = 0.d0
            else
                if (x .lt. bornck(nborck)) then
                    do i = 2, nborck
                        if (x .ge. bornck(i-1) .and. x .lt. bornck(i)) then
                            ck = coefck(i-1, 1)/(x*x*x*x*x*x*x)+coefck(i-1, 2)/(x*x*x*x*x*x)+co&
                                 &efck(i-1, 3)/(x*x*x*x*x)+coefck(i-1, 4)/(x*x*x*x)+coefck(i-1&
                                 &, 5)/(x*x*x)+coefck(i-1, 6)/(x*x)+coefck(i-1, 7)/(x)+coefck&
                                 &(i-1, 8)+coefck(i-1, 9)*(x)+coefck(i-1, 10)*(x*x)+coefck(i-1&
                                 &, 11)*(x*x*x)
                            goto 160
                        end if
                    end do
160                 continue
                else
                    ck = coefck(nborck, 1)/(x*x*x*x*x*x*x)+coefck(nborck, 2)/(x*x*x*x*x*x)+coe&
                         &fck(nborck, 3)/(x*x*x*x*x)+coefck(nborck, 4)/(x*x*x*x)+coefck(nborck&
                         &, 5)/(x*x*x)+coefck(nborck, 6)/(x*x)+coefck(nborck, 7)/(x)+coefck(&
                         &nborck, 8)+coefck(nborck, 9)*(x)+coefck(nborck, 10)*(x*x)+coefck(nb&
                         &orck, 11)*(x*x*x)
                end if
            end if
!
! --- CELLULE DE TUBES VIBRANTS EN DEBUT DE FAISCEAU VISCACHE1.
!
        else if (ires .eq. 1001) then
            if (x .lt. bornck(1)) then
                ck = 0.d0
            else if (x .lt. bornck(2)) then
                ck = coefck(1, 8)+coefck(1, 9)*x
            else if (x .lt. bornck(3)) then
                ck = coefck(2, 8)+coefck(2, 9)*x+coefck(2, 10)*(x*x)+coefck(2, 11)*(x*x*x)
            else if (x .lt. bornck(4)) then
                ck = coefck(3, 8)+coefck(3, 7)/x+coefck(3, 6)/(x*x)
            else
                ck = coefck(4, 8)
            end if
!
! --- CELLULE DE TUBES VIBRANTS EN MILIEU DE FAISCEAU VISCACHE1.
!
        else if (ires .eq. 1002) then
            if (x .lt. bornck(1)) then
                ck = 0.d0
            else if (x .lt. bornck(2)) then
                ck = coefck(1, 8)+coefck(1, 9)*x
            else if (x .lt. bornck(3)) then
                ck = coefck(2, 8)+coefck(2, 7)/x+coefck(2, 6)/(x*x)+coefck(2, 5)/(x*x*x)
            else if (x .lt. bornck(4)) then
                ck = coefck(3, 8)+coefck(3, 7)/x+coefck(3, 6)/(x*x)
            else if (x .lt. bornck(5)) then
                ck = coefck(4, 8)+coefck(4, 7)/x+coefck(4, 6)/(x*x)+coefck(4, 5)/(x*x*x)
            else if (x .lt. bornck(6)) then
                ck = coefck(5, 8)+coefck(5, 9)*x
            else
                ck = coefck(6, 8)
            end if
!
! --- TUBE UNIQUE VIBRANT EN MILIEU DE FAISCEAU RIGIDE VISCACHE1.
!
        else if (ires .eq. 1003) then
            if (x .lt. bornck(1)) then
                ck = 0.d0
            else if (x .lt. bornck(2)) then
                ck = coefck(1, 8)+coefck(1, 9)*x
            else if (x .lt. bornck(3)) then
                ck = coefck(2, 8)+coefck(2, 7)/x
            else if (x .lt. bornck(4)) then
                ck = coefck(3, 8)+coefck(3, 7)/x+coefck(3, 6)/(x*x)+coefck(3, 5)/(x*x*x)
            else if (x .lt. bornck(5)) then
                ck = coefck(4, 8)+coefck(4, 7)/x+coefck(4, 6)/(x*x)
            else if (x .lt. bornck(6)) then
                ck = coefck(5, 8)+coefck(5, 7)/x
            else if (x .lt. bornck(7)) then
                ck = coefck(6, 8)+coefck(6, 9)*x
            else
                ck = coefck(7, 8)
            end if
!
! --- TUBE UNIQUE VIBRANT EN DEBUT DE FAISCEAU RIGIDE VISCACHE1.
!
        else if (ires .eq. 1004) then
            if (x .lt. bornck(1)) then
                ck = 0.d0
            else if (x .lt. bornck(2)) then
                ck = coefck(1, 8)+coefck(1, 9)*x
            else if (x .lt. bornck(3)) then
                ck = coefck(2, 8)+coefck(2, 7)/x+coefck(2, 6)/(x*x)
            else if (x .lt. bornck(4)) then
                ck = coefck(3, 8)+coefck(3, 9)*x
            else if (x .lt. bornck(5)) then
                ck = coefck(4, 8)+coefck(4, 7)/x+coefck(4, 6)/(x*x)+coefck(4, 5)/(x*x*x)
            else if (x .lt. bornck(6)) then
                ck = coefck(5, 8)+coefck(5, 9)*x
            else
                ck = coefck(6, 8)
            end if
!
        end if
!
    end if
!
! --- INVERSION DE CK POUR CONVENTION DE SIGNE ET FIN DE COEFR.
!
    ck = -ck
!
1000 continue
    call jedema()
end subroutine
