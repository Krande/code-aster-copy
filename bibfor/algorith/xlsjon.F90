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

subroutine xlsjon(ino, jlsn, nfiss, jfisco, jonc_no)
    implicit none
!
#include "jeveux.h"
    integer(kind=8) :: ino, jlsn, nfiss, jfisco
    real(kind=8) :: jonc_no(nfiss)
!
!       CALCUL DE LA LSN D'UN NOEUD DE L'ELEMENT PARENT EN TENANT COMPTE
!       DES JONCTIONS POUR CHAQUE FISSURE
!
!     ENTREE
!       INO      : NUMERO DU NOEUD COURANT DE L'ELT PARENT
!       NFISS    : NOMBRE DE FISSURE DANS LE SUPPORT DU NOEUD
!     SORTIE
!       jonc_no    : FONCTION HEAVISIDE AU NOEUD POUR CHAQUE FISSURE EN TENANT COMPTE
!                  DES JONCTIONS
!
    integer(kind=8) :: nfimax
    parameter(nfimax=10)
    integer(kind=8) :: ifiss, i, nfisc, ifisc
    integer(kind=8) :: fisco(2*nfimax), fisc(2*nfimax)
    real(kind=8) :: ljonc(nfimax+1)
!
!.....................................................................
!
!
    if (nfiss .eq. 1) then
!
        if (zr(jlsn-1+ino) .lt. 0.d0) then
            jonc_no(1) = -1.d0
        else
            jonc_no(1) = +1.d0
        end if
!
    else
!
        fisc(1:2*nfimax) = 0
        fisco(1:2*nfimax) = 0
        if (nfiss .gt. 1) fisco(1:2*nfiss) = zi(jfisco:(jfisco+2*nfiss-1))
        jonc_no(1:nfiss) = 0.d0
        ljonc(1:nfimax) = 0.d0
        do ifiss = 1, nfiss
!   REMPLISSAGE DE FISC POUR CHAQUE FISSURE
            fisc(1:2*nfiss) = 0
            ifisc = ifiss
            nfisc = 0
80          continue
            if (fisco(2*ifisc-1) .gt. 0) then
                nfisc = nfisc+1
                fisc(2*(nfisc-1)+2) = fisco(2*ifisc)
                ifisc = fisco(2*ifisc-1)
                fisc(2*(nfisc-1)+1) = ifisc
                goto 80
            end if
!
            do i = 1, nfisc
                ljonc(i) = zr(jlsn-1+(ino-1)*nfiss+fisc(2*i-1))
            end do
            ljonc(nfisc+1) = zr(jlsn-1+(ino-1)*nfiss+ifiss)
!
!   MISE A ZERO POUR LA FONCTION JONCTION AU NIVEAU DU BRANCHEMENT
            do i = 1, nfisc
                if (fisc(2*i) .gt. 0.d0 .and. fisc(2*i)*ljonc(i) .ge. 0.d0) goto 300
                if (fisc(2*i) .lt. 0.d0 .and. fisc(2*i)*ljonc(i) .gt. 0.d0) goto 300
            end do
!
            if (ljonc(nfisc+1) .lt. 0.d0) then
                jonc_no(ifiss) = -1.d0
            else if (ljonc(nfisc+1) .ge. 0.d0) then
                jonc_no(ifiss) = +1.d0
            end if
300         continue
        end do
!
    end if
!
!
end subroutine
