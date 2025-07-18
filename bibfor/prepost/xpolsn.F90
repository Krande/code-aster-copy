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
subroutine xpolsn(elrefp, ino, n, jlsn, jlst, &
                  ima, iad, igeom, nfiss, ndime, &
                  ndim, jconx1, jconx2, fisco, co, &
                  lsn, lst)
! aslint: disable=W1306
    implicit none
!
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/xpoffo.h"
    integer(kind=8) :: n, jlsn, jlst, ndim, ndime, ino, nfiss
    integer(kind=8) :: iad, igeom, ima, jconx1, jconx2, fisco(*)
    character(len=8) :: elrefp
    real(kind=8) :: co(3), lsn(nfiss), lst(nfiss)
!
!              RECUPERATION DES COORDONNEES RELLES ET
!                 CALCUL DES LEVEL SET DE L'ENTITE
!
!   IN
!     ELREFP : ÉLÉMENT DE RÉFÉRENCE PARENT
!     INO   : NUMERO DU POINT D'INTER OU DU NOEUD
!     N      : NOMBRE DE NOEUDS DE L'ÉLÉMENT PARENT
!     JLSN   : ADRESSE DU CHAM_NO_S DE LA LEVEL NORMALE
!     JLST   : ADRESSE DU CHAM_NO_S DE LA LEVEL TANGENTE
!     IMA    : NUMERO DE MAILLE COURANTE PARENT
!     IAD    : POINTEUR DES COORDONÉES DE INO
!     IGEOM  : COORDONNÉES DES NOEUDS DE L'ÉLÉMENT PARENT
!     NFISS  : NOMBRE DE FISSURES "VUES" PAR L'ÉLÉMENT PARENT
!     NDIME  : DIMENSION TOPOLOGIQUE DE LA MAILLE PARENT
!     NDIM   : DIMENSION DU MAILLAGE
!     JCONX1 : ADRESSE DE LA CONNECTIVITE DU MAILLAGE SAIN
!              (CONNECTIVITE QUADRATIQUE SI LAGRANGES DE CONTACT
!              AUX ARETES)
!     JCONX2 : LONGUEUR CUMULEE DE LA CONNECTIVITE DU MAILLAGE SAIN
!              (CONNECTIVITE QUADRATIQUE SI LAGRANGES DE CONTACT
!              AUX ARETES)
!     FISCO  : CONNECTIVITE FISSURE/DDL
!   OUT
!     CO     : COORDONNEES DU NOEUD OU DU POINT
!     LSN    : LEVEL SET NORMALE DU NOEUD OU DU POINT
!     LST    : LEVEL SET TANGENTE DU NOEUD OU DU POINT
!
    real(kind=8) :: ff(n), somlsn(nfiss)
    integer(kind=8) :: i, j, ifiss, ifisc, nfisc, fisc(2*nfiss)
!
    call jemarq()
    co(:) = 0.d0
!
!     CAS D'UN NOEUD
    if (ino .lt. 1000) then
        i = zi(jconx1-1+zi(jconx2+ima-1)+ino-1)
        co(1) = zr(iad-1+3*(i-1)+1)
        co(2) = zr(iad-1+3*(i-1)+2)
        co(3) = zr(iad-1+3*(i-1)+3)
!
        do ifiss = 1, nfiss
!
            lsn(ifiss) = zr(jlsn-1+nfiss*(ino-1)+ifiss)
            lst(ifiss) = zr(jlst-1+nfiss*(ino-1)+ifiss)
!
!     TRAITEMENT SPECIAL POUR LES JONCTIONS
            do i = 1, 2*nfiss
                fisc(i) = 0
            end do
            ifisc = ifiss
            nfisc = 0
80          continue
            if (fisco(2*ifisc-1) .gt. 0) then
!     STOCKAGE DES FISSURES SUR LESQUELLES IFISS SE BRANCHE
                nfisc = nfisc+1
                fisc(2*(nfisc-1)+2) = fisco(2*ifisc)
                ifisc = fisco(2*ifisc-1)
                fisc(2*(nfisc-1)+1) = ifisc
                goto 80
            end if
            do i = 1, nfisc
                if (fisco(2*i)*zr(jlsn-1+(ino-1)*nfiss+fisco(2*i-1)) .gt. 0) lsn(ifiss) = &
                    1.d0
            end do
!
        end do
!
!     CAS D'UN POINT D'INTERSECTION OU D'UN POINT MILIEU
    else if (ino .gt. 1000) then
        do i = 1, ndim
            co(i) = zr(iad-1+i)
        end do
!
!       FF : FONCTIONS DE FORMES
        call xpoffo(ndim, ndime, elrefp, n, igeom, &
                    co, ff)
!
        do ifiss = 1, nfiss
!
            lsn(ifiss) = 0.d0
            lst(ifiss) = 0.d0
            do i = 1, n
                lsn(ifiss) = lsn(ifiss)+zr(jlsn-1+nfiss*(i-1)+ifiss)*ff( &
                             i)
                lst(ifiss) = lst(ifiss)+zr(jlst-1+nfiss*(i-1)+ifiss)*ff( &
                             i)
            end do
!
!     TRAITEMENT SPECIAL POUR LES JONCTIONS
            do i = 1, 2*nfiss
                fisc(i) = 0
            end do
            ifisc = ifiss
            nfisc = 0
90          continue
            if (fisco(2*ifisc-1) .gt. 0) then
!     STOCKAGE DES FISSURES SUR LESQUELLES IFISS SE BRANCHE
                nfisc = nfisc+1
                fisc(2*(nfisc-1)+2) = fisco(2*ifisc)
                ifisc = fisco(2*ifisc-1)
                fisc(2*(nfisc-1)+1) = ifisc
                goto 90
            end if
            somlsn(:) = 0.d0
            do i = 1, n
                do j = 1, nfisc
                    somlsn(j) = somlsn(j)+ff(i)*zr(jlsn-1+(i-1)*nfiss+fisco(2*j-1))
                end do
            end do
            do i = 1, nfisc
                if (fisco(2*i)*somlsn(i) .gt. 0) lsn(ifiss) = 1.d0
            end do
!
        end do
!
    end if
!
    call jedema()
!
end subroutine
