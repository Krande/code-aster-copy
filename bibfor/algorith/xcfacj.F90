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
subroutine xcfacj(ptint, ptmax, ipt, ainter, lsn, &
                  igeom, nno, ndim, nfiss, ifiss, &
                  fisco, nfisc, typma)
!
! aslint: disable=W1306
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/confac.h"
#include "asterfort/elref1.h"
#include "asterfort/elrfvf.h"
#include "asterfort/padist.h"
#include "asterfort/reereg.h"
#include "asterfort/xajpin.h"
!
    integer(kind=8) :: ptmax, ipt, igeom, nno, ndim
    real(kind=8) :: lsn(*), ptint(*), ainter(*)
    integer(kind=8) :: nfiss, ifiss, fisco(*), nfisc
    character(len=8) :: typma
!
!              TROUVER LES PTS D'INTER ENTRE LES JONCTIONS DE FISSURE
!                 ET LES FACES POUR LES ELEMENTS MULTI-HEAVISIDE
!
!     ENTREE
!       PTINT   : COORDONNEES DES POINTS D'INTERSECTION
!       PTMAX    : NOMBRE MAX DE POINTS D'INTERSECTION
!       IPT      : COMPTEUR DE NOMBRE DE POINTS D'INTERSECTION
!       AINTER   : INFOS SUR LES ARETES ASSOCIEES
!       LSN      : VALEURS DE LA LEVEL SET NORMALE
!       IGEOM    : ADRESSE DES COORDONNEES DES NOEUDS DE L'ELT PARENT
!       NNO      : NOMBRE DE NOEUDS DE L'ELEMENT
!       NDIM     : DIMENSION DE L'ESPACE
!       NFISS    : NOMBRE DE FISSURES VUES DANS L'ÉLÉMENT
!       IFISS    : NUMÉRO DE LA FISSURE EN COURS
!       FISCO    : NUM ET COEF DES FISS SUR LESQUELLES IFISS SE BRANCHE
!       NFISC    : NOMBRE DE FISSURES SUR LESQUELLES IFISS SE BRANCHE
!       TYPMA    : TYPE DE LA MAILLE ASSOCIEE A L'ELEMENT
!
!     SORTIE
!       PTINT    : COORDONNEES DES POINTS D'INTERSECTION
!       IPT      : COMPTEUR DE NOMBRE DE POINTS D'INTERSECTION
!       AINTER   : INFOS SUR LES ARETES ASSOCIEES
!
!     ------------------------------------------------------------------
!
    character(len=8) :: elref, alias
    real(kind=8) :: lsna, lsnb, lsja, lsjb, lsj
    real(kind=8) :: a(3), b(3), c(3), mp(2), prec, ff(nno)
    real(kind=8) :: loncar, dst
    real(kind=8) :: m(3), somlsn, epsi(2), coorma(8)
    integer(kind=8) :: i, nbf, ibid, ifq, j
    integer(kind=8) :: fa(6, 8), ibid3(12, 3), ifisc, jfisc, ino
    integer(kind=8) :: nnof, na, nb, iret
    aster_logical :: chgsgn, lajpf, ajout, c1, c2
! ----------------------------------------------------------------------
!
    call elref1(elref)
    prec = 100.d0*r8prem()
    call confac(typma, ibid3, ibid, fa, nbf)
!
!     BOUCLE SUR LES FACES
    do ifq = 1, nbf
!
        lajpf = ASTER_FALSE
        if (fa(ifq, 4) .eq. 0) then
            nnof = 3
            alias = 'TR3'
        else
            nnof = 4
            alias = 'QU4'
        end if
        ASSERT(nnof .le. 4)
!
!       RECHERCHE DES INTERSECTION ENTRE LSN ET LES LSJ SUR LA FACE
!
        somlsn = 0.d0
        do i = 1, nnof
            somlsn = somlsn+abs(lsn(fa(ifq, i)))
        end do
        if (somlsn .eq. 0.d0) goto 200
!
        do ifisc = 1, nfisc
            chgsgn = ASTER_FALSE
            do i = 1, nnof
                na = fa(ifq, i)
                if (i .eq. 1) then
                    j = nnof
                else
                    j = i-1
                end if
                nb = fa(ifq, j)
                lsna = lsn((na-1)*nfiss+ifiss)
                lsnb = lsn((nb-1)*nfiss+ifiss)
                lsja = lsn((na-1)*nfiss+fisco(2*ifisc-1))*fisco(2*ifisc)
                lsjb = lsn((nb-1)*nfiss+fisco(2*ifisc-1))*fisco(2*ifisc)
                coorma(2*i-1) = lsja
                coorma(2*i) = lsna
!           SI LE FOND COINCIDE AVEC UN COTE DE LA FACE, ON SORT
                if (lsna .eq. 0.d0 .and. lsnb .eq. 0.d0 .and. lsja .eq. 0.d0 .and. lsjb .eq. &
                    0.d0) goto 110
!           ON ACCEPTE SI LE FRONT EST SUR UN NOEUD OU UN PT DE L'ARETE
                if (lsna .eq. 0.d0 .and. lsja .eq. 0.d0 .or. lsna .eq. 0.d0 .and. lsnb .eq. &
                    0.d0 .and. (lsja*lsjb) .lt. r8prem()) chgsgn = ASTER_TRUE
!           ON ACCEPTE SI UNE ARETE DE LA FACETTE EST COUPÉE
                c1 = abs(lsna-lsnb) .gt. r8prem()
                if (c1) then
                    c1 = lsjb-lsnb*(lsja-lsjb)/(lsna-lsnb) .lt. prec
                end if
                c2 = abs(lsna-lsnb) .le. r8prem() .and. (lsja*lsjb) .lt. r8prem()
                if (c1 .or. c2) then
                    chgsgn = ASTER_TRUE
                end if
            end do
            if (.not. chgsgn) goto 110
!
!         ON CHERCHE SUR LA MAILLE LE POINT CORRESPONDANT À LSN=LSJ=0
            mp(1) = 0.d0
            mp(2) = 0.d0
            call reereg('C', alias, nnof, coorma, mp, &
                        2, epsi, iret)
            if (iret .eq. 1) goto 110
!         ON NE PREND PAS EN COMPTE LES POINTS QUI SORTENT DU DOMAINE
            if (alias .eq. 'QU4') then
                if (abs(epsi(1)) .gt. 1.d0) goto 110
                if (abs(epsi(2)) .gt. 1.d0) goto 110
            else if (alias .eq. 'TR3') then
                if (epsi(1) .lt. 0.d0) goto 110
                if (epsi(2) .lt. 0.d0) goto 110
                if (epsi(1)+epsi(2) .gt. 1.d0) goto 110
            end if
            mp(1) = epsi(1)
            mp(2) = epsi(2)
            call elrfvf(alias, mp, ff)
            do jfisc = ifisc+1, nfisc
                lsj = 0
                do j = 1, nnof
                    ino = fa(ifq, j)
                    lsj = lsj+lsn((ino-1)*nfiss+fisco(2*jfisc-1))*fisco(2*jfisc)*ff(j)
                end do
                if (lsj .gt. 0) goto 110
            end do
            lajpf = ASTER_TRUE
            do i = 1, ndim
                m(i) = 0
                do j = 1, nnof
                    ino = fa(ifq, j)
                    m(i) = m(i)+zr(igeom-1+ndim*(ino-1)+i)*ff(j)
                end do
            end do
!
110         continue
        end do
!
        if (lajpf) then
!       POUR IGNORER LES POINTS CONFONDUS AVEC CEUX
!       DETECTES DANS XCFACE LORSQUE LE PT EST EXACT SUR UNE ARETE
            do j = 1, ipt
                dst = padist(ndim, m, ptint(ndim*(j-1)+1))
                if (dst .le. r8prem()) lajpf = ASTER_FALSE
            end do
        end if
!
        if (lajpf) then
!       ON AJOUTE A LA LISTE LE POINT M
            do i = 1, ndim
                a(i) = zr(igeom-1+ndim*(fa(ifq, 1)-1)+i)
                b(i) = zr(igeom-1+ndim*(fa(ifq, 2)-1)+i)
                c(i) = zr(igeom-1+ndim*(fa(ifq, 3)-1)+i)
            end do
            loncar = (padist(ndim, a, b)+padist(ndim, a, c))/2.d0
            call xajpin(ndim, ptint, ptmax, ipt, ibid, &
                        m, loncar, ainter, 0, 0, &
                        0.d0, ajout)
        end if
200     continue
    end do
!
!
end subroutine
