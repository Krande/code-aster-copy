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
subroutine xdecou(ndim, elp, nnop, nnose, it, &
                  pintt, cnset, lsn, fisco, igeom, &
                  nfiss, ifiss, pinter, ninter, npts, &
                  ainter, lonref, nfisc)
! aslint: disable=W1306
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/conare.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/lteatt.h"
#include "asterfort/reeref.h"
#include "asterfort/xajpin.h"
#include "asterfort/xinter.h"
#include "asterfort/xxmmvd.h"
    real(kind=8) :: lsn(*), pintt(*), lonref, pinter(*), ainter(*)
    integer(kind=8) :: ndim, nnop, nnose, it, cnset(*), ninter, igeom, npts
    integer(kind=8) :: nfiss, ifiss, fisco(*), nfisc
    character(len=8) :: elp
! person_in_charge: samuel.geniaut at edf.fr
!                      TROUVER LES PTS D'INTERSECTION ENTRE LES ARETES
!                      ET LE PLAN DE FISSURE
!
!     ENTREE
!       NDIM     : DIMENSION DE L'ESPACE
!       ELP      : ELEMENT DE REFERENCE PARENT
!       NNOP     : NOMBRE DE NOEUDS DE L'ELEMENT PARENT
!       NNOSE    : NOMBRE DE NOEUDS DU SOUS TETRA
!       IT       : INDICE DU TETRA EN COURS
!       PINTT    :
!       CNSET    : CONNECTIVITEE DES NOEUDS DU TETRA
!       LSN      : VALEURS DE LA LEVEL SET NORMALE
!       FISCO    :
!       IGEOM    : ADRESSE DES COORDONNEES DES NOEUDS DE L'ELT PARENT
!       NFIS     : NOMBRE DE FISSURES "VUES" PAR L'ÉLÉMENT
!       NFISS
!       IFISS
!       LONREF   :
!       NFISC    :
!
!     SORTIE
!       PINTER   : COORDONNEES DES POINTS D'INTERSECTION
!       NINTER   : NB DE POINTS D'INTERSECTION
!       NPTS     : NB DE PTS D'INTERSECTION COINCIDANT AVEC UN NOEUD
!       AINTER   : INFOS ARETE ASSOCIEES AU POINTS D'INTERSECTION
!     ------------------------------------------------------------------
!
    real(kind=8) :: a(3), b(3), c(3), lsna, lsnb, lsnm, tampor(4)
    real(kind=8) :: somlsn(nfisc+1), ff(nnop), lsnelp(8)
    real(kind=8) :: rbid2(ndim), cref(3)
    integer(kind=8) :: ar(12, 3), nbar, nta, ntb, na, nb, nm, ins
    integer(kind=8) :: ia, i, j, ipt, ibid, pp, pd, k, ptmax
    integer(kind=8) :: ndime, a1, a2
    integer(kind=8) :: mxstac
    character(len=8) :: typma
    integer(kind=8) :: zxain
    aster_logical :: axi, papillon, ajout
    parameter(mxstac=1000)
!
! ----------------------------------------------------------------------
!
!
!     VERIF QUE LES TABLEAUX LOCAUX DYNAMIQUES NE SONT PAS TROP GRANDS
!     (VOIR CRS 1404)
    ASSERT(nnop .le. mxstac)
    ASSERT(nfisc .le. mxstac)
    ASSERT(ndim .le. mxstac)
!
    zxain = xxmmvd('ZXAIN')
    call elrefe_info(fami='RIGI', ndim=ndime)
!
    axi = lteatt('AXIS', 'OUI')
!
    if (ndim .eq. 3) then
!
        if (ndime .eq. 3) then
            typma = 'TETRA4'
            ptmax = 4
        else if (ndime .eq. 2) then
            typma = 'TRIA3'
            ptmax = 3
        else if (ndime .eq. 1) then
            typma = 'SEG2'
            ptmax = 2
        end if
!
    else if (ndim .eq. 2) then
!
        if (ndime .eq. 2) then
            typma = 'TRIA3'
            ptmax = 2
        else if (ndime .eq. 1) then
            typma = 'SEG2'
            ptmax = 2
        end if
!
    else if (ndim .eq. 1) then
!
        ASSERT(.false.)
!
    end if
!
!     VECTEUR REEL A ZXAIN COMPOSANTES, POUR CHAQUE PT D'INTER :
!     - NUMERO ARETE CORRESPONDANTE (0 SI C'EST UN NOEUD SOMMET)
!     - VRAI NUMERO NOEUD CORRESPONDANT (SERT QUE POUR NOEUD SOMMET)
!     - LONGUEUR DE L'ARETE
!     - POSITION DU PT SUR L'ARETE
!     - ARETE VITALE (NE SERT A RIEN ICI)
!
!     COMPTEUR DE POINT INTERSECTION ET POINT D'INTERSECTION SOMMENT
    ipt = 0
    ins = 0
!
!     SOMME DES LSN SUR LES NOEUDS DU SE
    somlsn(:) = 0.d0
    do k = 1, nnose
        na = cnset(nnose*(it-1)+k)
        if (na .lt. 1000) then
            do i = 1, nfisc
                somlsn(i) = somlsn(i)+lsn((na-1)*nfiss+fisco(2*i-1))
            end do
        else
!         RECUP COOR GLOBALES
            a(:) = 0.d0
            do i = 1, ndim
                a(i) = pintt(ndim*(na-1001)+i)
            end do
!           CALCUL DES FF
            call reeref(elp, nnop, zr(igeom), a, ndim, &
                        rbid2, ff)
!           INTERPOLATION LSN
            do j = 1, nnop
                do i = 1, nfisc
                    somlsn(i) = somlsn(i)+ff(j)*lsn((j-1)*nfiss+fisco(2* &
                                                                      i-1))
                end do
            end do
        end if
    end do
!  SI ON EST PAS DU COTÉ INTERSECTÉ, ON SORT
    do i = 1, nfisc
        if (fisco(2*i)*somlsn(i) .gt. 0) goto 999
    end do
!
    call conare(typma, ar, nbar)
!
!     BOUCLE SUR LES ARETES POUR DETERMINER LES POINTS D'INTERSECTION
!
    do ia = 1, nbar
!       NUM NO DU SOUS-ELEMENT
        nta = ar(ia, 1)
        ntb = ar(ia, 2)
!       NUM NO OU POINT D'INTER DE L'ELEMENT PARENT
        na = cnset(nnose*(it-1)+nta)
        nb = cnset(nnose*(it-1)+ntb)
!
        a(:) = 0.d0
        b(:) = 0.d0
        do i = 1, ndim
            if (na .lt. 1000) then
                a(i) = zr(igeom-1+ndim*(na-1)+i)
            else
                a(i) = pintt(ndim*(na-1001)+i)
            end if
            if (nb .lt. 1000) then
                b(i) = zr(igeom-1+ndim*(nb-1)+i)
            else
                b(i) = pintt(ndim*(nb-1001)+i)
            end if
        end do
!        LONGAR=PADIST(NDIM,A,B)
!
        if (na .lt. 1000) then
            lsna = lsn((na-1)*nfiss+ifiss)
        else
!         CALCUL DES FF
            call reeref(elp, nnop, zr(igeom), a, ndim, &
                        rbid2, ff)
!         INTERPOLATION LSN
            lsna = 0
            do i = 1, nnop
                lsna = lsna+ff(i)*lsn((i-1)*nfiss+ifiss)
            end do
            if (abs(lsna) .lt. lonref*1.d-4) lsna = 0
        end if
        if (nb .lt. 1000) then
            lsnb = lsn((nb-1)*nfiss+ifiss)
        else
!         CALCUL DES FF
            call reeref(elp, nnop, zr(igeom), b, ndim, &
                        rbid2, ff)
!         INTERPOLATION LSN
            lsnb = 0
            do i = 1, nnop
                lsnb = lsnb+ff(i)*lsn((i-1)*nfiss+ifiss)
            end do
            if (abs(lsnb) .lt. lonref*1.d-4) lsnb = 0
        end if
!
        if ((lsna*lsnb) .le. 0) then
            if (lsna .eq. 0) then
!           ON AJOUTE A LA LISTE LE POINT A
                call xajpin(ndim, pinter, ptmax, ipt, ins, &
                            a, lonref, ainter, 0, na, &
                            0.d0, ajout)
            end if
            if (lsnb .eq. 0) then
!           ON AJOUTE A LA LISTE LE POINT B
                call xajpin(ndim, pinter, ptmax, ipt, ins, &
                            b, lonref, ainter, 0, nb, &
                            0.d0, ajout)
            end if
            if (lsna .ne. 0 .and. lsnb .ne. 0) then
!           INTERPOLATION DES COORDONNEES DE C
                do i = 1, ndim
                    c(i) = a(i)-lsna/(lsnb-lsna)*(b(i)-a(i))
                end do
                nm = 0
                lsnm = (lsna+lsnb)/2.d0
                do i = 1, nnop
                    lsnelp(i) = lsn((i-1)*nfiss+ifiss)
                end do
                call xinter(ndim, ndime, elp, zr(igeom), lsnelp, &
                            na, nb, nm, pintt, pintt, &
                            lsna, lsnb, lsnm, cref, c)
!           POSITION DU PT D'INTERSECTION SUR L'ARETE
!            ALPHA=PADIST(NDIM,A,C)
!           ON AJOUTE A LA LISTE LE POINT C
                call xajpin(ndim, pinter, ptmax, ipt, ibid, &
                            c, lonref, ainter, ia, 0, &
                            0.d0, ajout)
            end if
        end if
    end do
!
999 continue
    ninter = ipt
    npts = ins
    ASSERT(ninter .ge. npts .and. ninter .le. ptmax)
!
!     TRI DES POINTS D'INTERSECTION PAR ORDRE CROISSANT DES ARETES
    do pd = 1, ninter-1
        pp = pd
        do i = pp, ninter
            if (ainter(zxain*(i-1)+1) .lt. ainter(zxain*(pp-1)+1)) pp = i
        end do
        do k = 1, 4
            tampor(k) = ainter(zxain*(pp-1)+k)
            ainter(zxain*(pp-1)+k) = ainter(zxain*(pd-1)+k)
            ainter(zxain*(pd-1)+k) = tampor(k)
        end do
        do k = 1, ndim
            tampor(k) = pinter(ndim*(pp-1)+k)
            pinter(ndim*(pp-1)+k) = pinter(ndim*(pd-1)+k)
            pinter(ndim*(pd-1)+k) = tampor(k)
        end do
    end do
!
!      TRI DES POINTS POUR QUE LE POLYGONE IP1,IP2,IP3,IP4 SOIT CONVEXE
!      IP1 IP2 ET IP3 ONT UN SOMMET EN COMMUN
!      IP1 ET IP4 N ONT PAS DE SOMMET COMMUN
    if (ninter .eq. 4 .and. npts .eq. 0) then
        a1 = nint(ainter(1))
        do ia = 2, 3
            a2 = nint(ainter(zxain*(ia-1)+1))
            papillon = .true.
            do i = 1, 2
                do j = 1, 2
                    if (ar(a1, i) .eq. ar(a2, j)) papillon = .false.
                end do
            end do
            if (papillon) then
!        CONFIGURATION RENCONTREE PAR EXEMPLE DANS SSNV510C
                do k = 1, (zxain-1)
                    tampor(k) = ainter(zxain*(ia-1)+k)
                    ainter(zxain*(ia-1)+k) = ainter(zxain*(4-1)+k)
                    ainter(zxain*(4-1)+k) = tampor(k)
                end do
                do k = 1, ndim
                    tampor(k) = pinter(ndim*(ia-1)+k)
                    pinter(ndim*(ia-1)+k) = pinter(ndim*(4-1)+k)
                    pinter(ndim*(4-1)+k) = tampor(k)
                end do
            end if
        end do
    end if
end subroutine
