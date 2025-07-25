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
subroutine mtdete(option, method, lmat, mantis, expo, &
                  cmod)
    implicit none
#include "jeveux.h"
#include "asterc/r8miem.h"
#include "asterfort/almulr.h"
#include "asterfort/assert.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: lmat, expo, option
    real(kind=8) :: mantis
    complex(kind=8) :: cmod
    character(len=24) :: method
!     OPTION=1
!       CALCUL DU DETERMINANT D'UNE MATRICE ASSEMBLEE REELLE DECOMPOSEE
!       SOUS FORME L*D*LT
!       RESULTAT SOUS FORME MANTIS * (10**EXPO)
!              AVEC 10**(-30)<MANTISSE<10**(+30)
!     OPTION=2
!       CALCUL DU DETERMINANT NORMALISE D'UNE MATRICE COMPLEXE DECOM
!       POSEE SOUS FORME L*D*LT
!       RESULTAT SOUS FORME CMOD
!
!     SI METHODE='LDLT' OU 'MULT_FRONT': ON UTILISE L'OBJET .DIGS
!     SI METHODE='MUMPS': ON UTILISE L'OBJET '&&AMUMP.DETERMINANT'
!     ------------------------------------------------------------------
!
!     ------------------------------------------------------------------
!
!     ------------------------------------------------------------------
    integer(kind=8) :: i, neq, iret, ldiag, nbneg, ibid, ipiv, itrent, info34, ifm, niv
    integer(kind=8) :: ie
    real(kind=8) :: trent, trent1, rinf12, rinf13, rmin, rauxx, rauxy, rauxm
    complex(kind=8) :: cun, caux
    character(len=24) :: nomdia, kpiv
    data nomdia/'                   .DIGS'/
!     ------------------------------------------------------------------
!
!
    call jemarq()
    call infniv(ifm, niv)
    cun = dcmplx(1.d0, 0.d0)
    rmin = r8miem()*100
! --- TEST DES PARAMETRES D'ENTREES
    if (option .eq. 1) then
! --- LA MATRICE DOIT ETRE REELLE
        if (zi(lmat+3) .ne. 1) then
            ASSERT(.false.)
        end if
    else if (option .eq. 2) then
! --- LA MATRICE DOIT ETRE COMPLEXE
        if (zi(lmat+3) .ne. 2) then
            ASSERT(.false.)
        end if
    else
! --- MAUVAISE OPTION DE CALCUL
        ASSERT(.false.)
    end if
!
!
    if ((method(1:4) .eq. 'LDLT') .or. (method(1:10) .eq. 'MULT_FRONT')) then
! --- INITS. LDLT/MF
        nomdia(1:19) = zk24(zi(lmat+1))
        neq = zi(lmat+2)
        call jeexin(nomdia, iret)
        if (iret .eq. 0) then
            call utmess('F', 'MODELISA2_9', sk=nomdia)
        end if
        call jeveuo(nomdia, 'L', ldiag)
        ldiag = ldiag+neq
!
        if (option .eq. 1) then
! --- CALCUL DET AVEC LDLT/MF
            call almulr('ZERO', zr(ldiag), neq, mantis, expo)
            nbneg = 0
            do i = 0, neq-1
                if (zr(ldiag+i) .le. 0.d0) nbneg = nbneg+1
            end do
            call jedetr(nomdia)
            if (niv .ge. 2) write (ifm, *) '<MTDETE 1 LDLT/MF>  MANTIS/EXPO  ', mantis, expo
!
        else if (option .eq. 2) then
! --- CALCUL DET NORMALISE AVEC LDLT/MF
            cmod = cun
            do i = 1, neq
                caux = zc(ldiag+i-1)
                rauxx = dble(cmod)
                rauxy = dimag(cmod)
                rauxm = sqrt(rauxx*rauxx+rauxy*rauxy)
                if (rauxm .lt. rmin) rauxm = 1.d0
                caux = caux/rauxm
                cmod = cmod*caux
            end do
!
            if (niv .ge. 2) write (ifm, *) '<MTDETE 2 LDLT/MF>  CMOD  ', cmod
        end if
!
    else if (method(1:5) .eq. 'MUMPS') then
! --- INITS. MUMPS
        kpiv = '&&AMUMP.DETERMINANT'
        call jeexin(kpiv, ibid)
        if (ibid .ne. 0) then
            call jeveuo(kpiv, 'L', ipiv)
        else
            ASSERT(.false.)
        end if
! --- LE DETERMINANT ISSU DE MUMPS EST STOCKE SOUS LA FORME:
! ---                 MANTISSE * (2**EXP)   MANTISSE COMPLEXE
! --- CF. ROUTINE AMUMPU.F OPTION=4
! --- ON LE MET SOUS LA FORME ASTER:
! ---                 MANTISSE * (10**EXP) AVEC  MANTISSE REEL
! ---             10**(-30)<MANTISSE<10**(+30)
! --- CF. ROUTINES ALMULR
! --- INIT. ASTER DU PROCESSUS DES ROUTINES MTDETE/ALMULR
        trent = 1.d30
        trent1 = 1.d-30
        itrent = 30
!
! --- DONNEES ISSUES DE MUMPS
        rinf12 = zr(ipiv)
        rinf13 = zr(ipiv+1)
        info34 = nint(zr(ipiv+2))
!
        if (option .eq. 1) then
! --- CALCUL DET AVEC LDLT/MF
!
! --- ON S'ATTEND A UN DETERMINANT REEL, CETTE VALEUR EST SUSPECTE !
            if ((rinf12 .gt. trent1) .and. (rinf13 .gt. (0.05d0*rinf12))) then
                ASSERT(.false.)
            end if
            mantis = rinf12
            expo = 0
            do i = 1, info34
                mantis = mantis*2.d0
                if (abs(mantis) .ge. trent) then
                    mantis = mantis*trent1
                    expo = expo+itrent
                else if (abs(mantis) .le. trent1) then
                    mantis = mantis*trent
                    expo = expo-itrent
                end if
            end do
            if (abs(mantis) .gt. rmin) then
                ie = nint(log10(abs(mantis)))
                mantis = mantis/(10.d0**ie)
                expo = expo+ie
            else
                mantis = 0.d0
                expo = 1
            end if
!
            if (niv .ge. 2) then
                write (ifm, *) '<MTDETE 1 MUMPS> RINF ', rinf12, rinf13, &
                    info34
                write (ifm, *) '<MTDETE 1 MUMPS>  MANTIS/EXPO  ',&
     &                   mantis, expo
            end if
!
        else if (option .eq. 2) then
! --- CALCUL DET NORMALISE POUR LDLT/MF
            rauxm = sqrt(rinf12*rinf12+rinf13*rinf13)
            if (rauxm .lt. rmin) rauxm = 1.d0
            cmod = dcmplx(rinf12, rinf13)/rauxm
!
            if (niv .ge. 2) then
                write (ifm, *) '<MTDETE 2 MUMPS> RINF ', rinf12, rinf13, &
                    info34
                write (ifm, *) '<MTDETE 2 MUMPS>  CMOD  ', cmod
            end if
        end if
!
! --- FIN IF METHOD
    end if
!
    call jedema()
end subroutine
