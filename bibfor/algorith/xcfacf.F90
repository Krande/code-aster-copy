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
subroutine xcfacf(ptint, ptmax, ipt, ainter, lsn, &
                  lst, igeom, nno, ndim, typma, &
                  noma, nmaabs)
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8maem.h"
#include "asterc/r8prem.h"
#include "asterfort/confac.h"
#include "asterfort/elref1.h"
#include "asterfort/intfac.h"
#include "asterfort/padist.h"
#include "asterfort/xajpin.h"
    integer(kind=8) :: ptmax, ipt, igeom, nno, ndim, nmaabs
    real(kind=8) :: lsn(nno), lst(nno), ptint(*), ainter(*)
    character(len=8) :: typma, noma
!
! person_in_charge: samuel.geniaut at edf.fr
!              TROUVER LES PTS D'INTERSECTION ENTRE LE FOND DE FISSURE
!                 ET LES FACES POUR LES ELEMENTS EN FOND DE FISSURE
!
!     ENTREE
!       PTINT    : COORDONNEES DES POINTS D'INTERSECTION
!       PTMAX    : NOMBRE MAX DE POINTS D'INTERSECTION
!       IPT      : COMPTEUR DE NOMBRE DE POINTS D'INTERSECTION
!       AINTER   : INFOS SUR LES ARETES ASSOCIEES
!       LSN      : VALEURS DE LA LEVEL SET NORMALE
!       LST      : VALEURS DE LA LEVEL SET TANGENTE
!       IGEOM    : ADRESSE DES COORDONNEES DES NOEUDS DE L'ELT PARENT
!       NNO      : NOMBRE DE NOEUDS DE L'ELEMENT
!       NDIM     : DIMENSION DE L'ESPACE
!       TYPMA    : TYPE DE LA MAILLE ASSOCIEE A L'ELEMENT
!       NOMA     : NOM DU MAILLAGE
!       NMAABS   : INDICE DE LA MAILLE
!
!     SORTIE
!       PTINT    : COORDONNEES DES POINTS D'INTERSECTION
!       IPT      : COMPTEUR DE NOMBRE DE POINTS D'INTERSECTION
!       AINTER   : INFOS SUR LES ARETES ASSOCIEES
!
!     ------------------------------------------------------------------
!
    character(len=8) :: elref
    real(kind=8) :: rbid, maxlsn, minlsn, maxlst, minlst
    real(kind=8) :: a(3), b(3), c(3)
    real(kind=8) :: loncar, dst
    real(kind=8) :: m(3)
    integer(kind=8) :: i, nbf, ibid, ifq, j, codret
    integer(kind=8) :: fa(6, 8), ibid3(12, 3), indptf(3)
    aster_logical :: ajout
! ----------------------------------------------------------------------
!
    call elref1(elref)
!
!     INITIALISATION DES MIN ET MAX
    maxlsn = -1.d0*r8maem()
    minlsn = r8maem()
    maxlst = -1.d0*r8maem()
    minlst = r8maem()
!
!     RECHERCHE DU MIN ET MAX DES LEVEL SETS SUR LES NOEUDS
    do i = 1, nno
        maxlsn = max(lsn(i), maxlsn)
        maxlst = max(lst(i), maxlst)
        minlsn = min(lsn(i), minlsn)
        minlst = min(lst(i), minlst)
    end do
!
!     SI CE N'EST PAS UN ELEMENT EN FOND DE FISSURE, ON SORT
!     EN FAIT, CE TEST NE PERMET PAS DE DETECTER TOUS LES CAS
!     IL SE PEUT DONC QUE L'ON NE SORTE PAS MAIS QUE L'ON NE
!     SOIT PAS SUR UNE ELEMENT EN FOND DE FISSURE
    if (minlsn*maxlsn .ge. 0.d0 .or. minlst*maxlst .ge. 0.d0) goto 999
!
    call confac(typma, ibid3, ibid, fa, nbf)
!
!     BOUCLE SUR LES FACES
    rbid = 0.d0
    do ifq = 1, nbf
!
!       RECHERCHE DES INTERSECTION ENTRE LE FOND DE FISSURE ET LA FACE
        call intfac(noma, nmaabs, ifq, fa, nno, &
                    lst, lsn, ndim, 'NON', ibid, &
                    ibid, igeom, m, indptf, [rbid], &
                    [rbid], codret)
        if (codret .eq. 0) goto 200
!
!       POUR IGNORER LES POINTS CONFONDUS AVEC CEUX
!       DETECTES DANS XCFACE LORSQUE LE PT EST EXACT SUR UNE ARETE
        do j = 1, ipt
            dst = padist(ndim, m, ptint(ndim*(j-1)+1))
            if (dst .le. r8prem()) goto 200
        end do
!
!       LONGUEUR CARACTERISTIQUE
        do i = 1, ndim
            a(i) = zr(igeom-1+ndim*(fa(ifq, 1)-1)+i)
            b(i) = zr(igeom-1+ndim*(fa(ifq, 2)-1)+i)
            c(i) = zr(igeom-1+ndim*(fa(ifq, 3)-1)+i)
        end do
        loncar = (padist(ndim, a, b)+padist(ndim, a, c))/2.d0
!
!       ON AJOUTE A LA LISTE LE POINT M
!
        call xajpin(ndim, ptint, ptmax, ipt, ibid, &
                    m, loncar, ainter, 0, 0, &
                    0.d0, ajout)
!
200     continue
    end do
!
999 continue
end subroutine
