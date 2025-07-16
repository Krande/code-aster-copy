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
! person_in_charge: samuel.geniaut at edf.fr
! aslint: disable=W1306
!
subroutine intfac(noma, nmaabs, ifq, fa, nno, &
                  lst, lsn, ndim, grad, jglsn, &
                  jglst, igeom, m, indptf, gln, &
                  glt, codret)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/elrfvf.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/reereg.h"
!
    integer(kind=8) :: ifq, fa(6, 8), nno, ndim, jglsn, jglst, igeom, codret
    integer(kind=8) :: indptf(3), nmaabs
    real(kind=8) :: lsn(nno), lst(nno), m(ndim), gln(ndim), glt(ndim)
    character(len=3) :: grad
    character(len=8) :: noma
!
!              TROUVER LES PTS D'INTERSECTION ENTRE LE FOND DE FISSURE
!                 ET UNE FACE POUR LES ELEMENTS EN FOND DE FISSURE
!
!     ENTREE
!  NOMA   : NOM DU MAILLAGE
!  NMAABS : NUMERO DE LA MAILLE
!  IFQ    : NUMERO LOCAL DE LA FACE DE LA MAILLE
!  FA     : COONECTIVITE DES FACES DE LA MAILLE
!  NNO    : NOMBRE DE NOEUDS DE LA MAILLE
!  LSN    : VECTEUR LOCAL DES VALEURS NODALES DE LA MAILLE POUR LSN
!  LST    : VECTEUR LOCAL DES VALEURS NODALES DE LA MAILLE POUR LST
!  NDIM   : DIMENSION DE L'ESPACE
!  GRAD   : SI 'OUI' : ON CALCULE AUSSI LES GRADIENTS AU POINT TROUVE
!  JGLSN  : ADRESSE DU VECTEUR LOCAL DES VALEURS NODALES DE GRAD DE LSN
!  JGLST  : ADRESSE DU VECTEUR LOCAL DES VALEURS NODALES DE GRAD DE LST
!  IGEOM  : ADRESSE DU VECTEUR LOCAL DES COORDONNEES DES NOEUDS
!
!     SORTIE
!  M      : POINT TROUVE
!  INDPTF : VECTEUR INDICE DU POINT M TROUVE
!  GLN    : GRAD DE LSN EN M (SI DEMANDE)
!  GLT    : GRAD DE LST EN M (SI DEMANDE)
!  CODRET : CODE RETOUR = 1 SI ON A BIEN TROUVE UN POINT M
!                       = 0 SI ON N'A PAS PU TROUVE UN UNIQUE POINT M
!
!     ------------------------------------------------------------------
!
    integer(kind=8) :: nnof, i, j, k, ino, iret, jconx2, numnoa, numnob
    real(kind=8) :: coorma(8), prec, mp(2), epsi(2), ff(nno), lsta, lsna, lstb
    real(kind=8) :: lsnb, solsn, a(ndim), b(ndim), mem(3), memo, normab, coeffk
    real(kind=8) :: prec2, length(12)
    character(len=8) :: alias
    aster_logical :: chgsgn, c1, c2
    integer(kind=8), pointer :: connex(:) => null()
! ----------------------------------------------------------------------
!
    call jemarq()
!
    m(1:ndim) = 0.d0
    if (grad .eq. 'OUI') then
        gln(1:ndim) = 0.d0
        glt(1:ndim) = 0.d0
    end if

    indptf = 0
!
    prec = r8prem()
    prec2 = 1.d-4
    codret = 0
!
!     INITIALISATION DES COORDONNéES (LS) DES NOEUDS DE LA FACE
    coorma = 0.d0
    lsta = 0.d0
    lsna = 0.d0
    lstb = 0.d0
    lsnb = 0.d0
!
!     RECUPERATION DES DONNEES SUR LE MAILLAGE
    call jeveuo(noma//'.CONNEX', 'L', vi=connex)
    call jeveuo(jexatr(noma//'.CONNEX', 'LONCUM'), 'L', jconx2)
!
!     NOMBRE DE SOMMETS DE LA FACE
    if (fa(ifq, 4) .eq. 0) then
        nnof = 3
        alias = 'TR3'
    else
        nnof = 4
        alias = 'QU4'
    end if
    ASSERT(nnof .le. 4)
!
!     NOEUDS SOMMETS DE LA FACE : FA(IFQ,1) ... FA(IFQ,NNOF)
!
    chgsgn = .false.
!
!     SI LA FACE COINCIDE AVEC LA SURFACE DE LA FISSURE, ON SORT
!     (CAD SI LES LSN DES SOMMETS DE LA FACE SONT TOUS NULS)
    solsn = 0.d0
    do i = 1, nnof
        solsn = solsn+abs(lsn(fa(ifq, i)))
    end do
    if (solsn .eq. 0.d0) goto 999
!
    do i = 1, nnof
        if (i .eq. 1) then
            j = nnof
        else
            j = i-1
        end if
        lsta = lst(fa(ifq, i))
        lsna = lsn(fa(ifq, i))
        lstb = lst(fa(ifq, j))
        lsnb = lsn(fa(ifq, j))
        coorma(2*i-1) = lsta
        coorma(2*i) = lsna
!
!       SI LE FOND COINCIDE AVEC UN COTE DE LA FACE, ON SORT
        if (lsna .eq. 0.d0 .and. lsnb .eq. 0.d0 .and. lsta .eq. 0.d0 .and. lstb .eq. 0.d0) then
            goto 999
        end if
!
!       ON ACCEPTE TOUT DE SUITE LA FACE SI LE FRONT COINCIDE
!       AVEC L'UN DES NOEUDS DE LA FACE
        if (lsna .eq. 0.d0 .and. lsta .eq. 0.d0) then
            chgsgn = .true.
            indptf(1) = 1
            indptf(2) = connex(zi(jconx2+nmaabs-1)+fa(ifq, i)-1)
            cycle
        end if
!
!       ON ACCEPTE TOUT DE SUITE LA FACE SI LE FRONT COINCIDE
!       AVEC UN POINT D'UNE ARETE DE LA FACE
        if (lsna .eq. 0.d0 .and. lsnb .eq. 0.d0 .and. lsta*lstb .lt. prec) then
            chgsgn = .true.
            indptf(1) = 2
            indptf(2) = connex(zi(jconx2+nmaabs-1)+fa(ifq, i)-1)
            indptf(3) = connex(zi(jconx2+nmaabs-1)+fa(ifq, j)-1)
            cycle
        end if
!
!       ON NE CONSERVE QUE LES FACES COUPEES
!       SUR UNE SEULE ARETE PAR LE DEMI-PLAN (LSN=0/LST<0)
!       POUR CELA, ON CONTROLE LE SIGNE DE LST(I) OU I EST LE POINT
!       D'INTERSECTION DE L'ARETE AVEC LSN=0
!       (ON A LSN(I)=LST(B)-LSN(B)*(LST(A)-LST(B))/(LSN(A)-LSN(B)))
        if (lsna*lsnb .lt. prec) then
            c1 = abs((lsna-lsnb)) .gt. prec
            if (c1) then
                c1 = lstb-(lsnb*(lsta-lstb)/(lsna-lsnb)) .lt. prec
            end if
            c2 = abs((lsna-lsnb)) .le. r8prem() .and. (lsta*lstb) .lt. r8prem()
            if (c1 .or. c2) then
                chgsgn = .true.
                indptf(1) = 3
                indptf(2) = 0
                indptf(3) = 0
            end if
        end if
    end do
!
    if (.not. chgsgn) goto 999
!
!     ON CHERCHE SUR LA MAILLE LE POINT CORRESPONDANT à LSN=LST=0
    mp(1) = 0.d0
    mp(2) = 0.d0
    call reereg('C', alias, nnof, coorma, mp, &
                2, epsi, iret, toler=prec*1.d4)
!
    if (iret .eq. 1) goto 999
!
!     ON NE PREND PAS EN COMPTE LES POINTS QUI SORTENT DU DOMAINE,
!     AVEC UNE TOLERANCE IDENTIQUE A CELLE UTILISEE POUR LA RECHERCHE
!     DU POINT SOLUTION FOURNI PAR REEREG
    if (alias .eq. 'QU4') then
        if (abs(epsi(1)) .gt. (1.d0+prec*1.d4)) goto 999
        if (abs(epsi(2)) .gt. (1.d0+prec*1.d4)) goto 999
    else if (alias .eq. 'TR3') then
        if (epsi(1) .lt. (0.d0-prec*1.d4)) goto 999
        if (epsi(2) .lt. (0.d0-prec*1.d4)) goto 999
        if (epsi(1)+epsi(2) .gt. (1.d0+prec*1.d4)) goto 999
    end if
!
    mp(1) = epsi(1)
    mp(2) = epsi(2)
!     ON DOIT MAINTENANT MULTIPLIER LES COORD. PARAM. DE M PAR CHACUNE
!     DES FF DES NOEUDS DE L'éLéMENT POUR OBTENIR LES COORD. CART.
    call elrfvf(alias, mp, ff)
    do i = 1, ndim
        do j = 1, nnof
            ino = fa(ifq, j)
            m(i) = m(i)+zr(igeom-1+ndim*(ino-1)+i)*ff(j)
            if (grad .eq. 'OUI') then
                glt(i) = glt(i)+zr(jglst-1+ndim*(ino-1)+i)*ff(j)
                gln(i) = gln(i)+zr(jglsn-1+ndim*(ino-1)+i)*ff(j)
            end if
        end do
    end do
!
!     TRAITEMENT DES POINTS M PROCHES DES SOMMETS (FIT TO VERTEX)
    if ((indptf(1) .eq. 2) .or. (indptf(1) .eq. 3)) then
        do i = 1, nnof
            memo = 0.d0
            do j = 1, ndim
                a(j) = zr(igeom-1+ndim*(fa(ifq, i)-1)+j)
                memo = memo+(a(j)-m(j))**2
            end do
            length(3*(i-1)+1) = sqrt(memo)
            length(3*(i-1)+2) = connex(zi(jconx2+nmaabs-1)+fa(ifq, &
                                                              i)-1)
            length(3*(i-1)+3) = 0
        end do
!       ON TRIE LE VECTEUR LENGTH
        do i = 1, nnof-1
            do j = i+1, nnof
                if (length(3*(j-1)+1) .lt. length(3*(i-1)+1)) then
                    do k = 1, 3
                        mem(k) = length(3*(i-1)+k)
                    end do
                    do k = 1, 3
                        length(3*(i-1)+k) = length(3*(j-1)+k)
                    end do
                    do k = 1, 3
                        length(3*(j-1)+k) = mem(k)
                    end do
                end if
            end do
        end do
!       M EST PROCHE D'UN SOMMET ? SI OUI, ON LE REPLACE SUR LE SOMMET
        if (length(1) .lt. (prec2*length(4))) then
            indptf(1) = 1
            indptf(2) = int(length(2))
            indptf(3) = 0
            goto 222
        end if
    end if
!
!     TRAITEMENT DES POINTS M PROCHES DES ARETES (FIT TO VERTEX)
    if (indptf(1) .eq. 3) then
        do i = 1, nnof
            do j = 1, ndim
                a(j) = zr(igeom-1+ndim*(fa(ifq, i)-1)+j)
            end do
            numnoa = connex(zi(jconx2+nmaabs-1)+fa(ifq, i)-1)
            if (i .eq. nnof) then
                do j = 1, ndim
                    b(j) = zr(igeom-1+ndim*(fa(ifq, 1)-1)+j)
                end do
                numnob = connex(zi(jconx2+nmaabs-1)+fa(ifq, 1)-1)
            else
                do j = 1, ndim
                    b(j) = zr(igeom-1+ndim*(fa(ifq, i+1)-1)+j)
                end do
                numnob = connex(zi(jconx2+nmaabs-1)+fa(ifq, i+1)-1)
            end if
            normab = 0.d0
            coeffk = 0.d0
            memo = 0.d0
            do k = 1, ndim
                normab = normab+(b(k)-a(k))**2
                coeffk = coeffk+(b(k)-a(k))*(m(k)-a(k))
            end do
            do k = 1, ndim
                memo = memo+(a(k)-m(k)+(coeffk/normab)*(b(k)-a(k)))**2
            end do
            length(3*(i-1)+1) = memo
            length(3*(i-1)+2) = numnoa
            length(3*(i-1)+3) = numnob
        end do
!       ON TRIE LE VECTEUR LENGTH
        do i = 1, nnof-1
            do j = i+1, nnof
                if (length(3*(j-1)+1) .lt. length(3*(i-1)+1)) then
                    do k = 1, 3
                        mem(k) = length(3*(i-1)+k)
                    end do
                    do k = 1, 3
                        length(3*(i-1)+k) = length(3*(j-1)+k)
                    end do
                    do k = 1, 3
                        length(3*(j-1)+k) = mem(k)
                    end do
                end if
            end do
        end do
!       M EST PROCHE D'UNE ARETE ? SI OUI, ON LE REPLACE SUR L'ARETE
        if (length(1) .lt. (prec2*length(4))) then
            indptf(1) = 2
            indptf(2) = int(length(2))
            indptf(3) = int(length(3))
            goto 222
        end if
    end if
!
222 continue
!
!     TOUT S'EST BIEN PASSE
    codret = 1
!
999 continue
!
    call jedema()
end subroutine
