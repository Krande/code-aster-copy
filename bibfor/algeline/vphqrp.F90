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
subroutine vphqrp(mat, neq, mxeq, icode, w, &
                  z, iz, wk, mxiter, ier, &
                  nitqr)
    implicit none
#include "asterfort/vpzbal.h"
#include "asterfort/vpzech.h"
#include "asterfort/vpzhes.h"
#include "asterfort/vpzqrh.h"
#include "asterfort/vpzrbk.h"
    integer(kind=8) :: neq, mxeq, icode, iz, ier, nitqr
    real(kind=8) :: mat(mxeq, 1), wk(neq, 1), w(1), z(1)
!     CALCUL DE TOUTES LES VALEURS PROPRES D'UNE MATRICE COMPLETE REELLE
!     MISE SOUS FORME DE HESSENBERG PUIS RESOLUTION PAR LA METHODE QR
!     ------------------------------------------------------------------
! VAR MAT   : R8 : MATRICE REEL D'ORDRE NEQ DONT ON CALCULE LES VALEURS
!                  TE LES VECTEURS PROPRES
!              *** EN SORTIE MAT EST DETRUITE ***
! IN  NEQ   : IS : ORDRE DE LA MATRICE
! IN  MXEQ  : IS : DIMENSION EXACT DE LA MATRICE ( IE  MAT(MXEQ,1) )
! IN  ICODE : IS : CODE DE CALCUL
!          = 0, CALCUL DES VALEURS PROPRES SEULEMENT
!          = 1, CALCUL DES VALEURS ET VECTEURS PROPRES
! OUT W     : C8 : VECTEUR (COMPLEXE) DES VALEURS PROPRES DE LA MATRICE
! OUT Z     : C8 : MATRICE (COMPLEXE) DES VECTEURS PROPRES DE LA MATRICE
!                  LA J-IEME COLONNE CONTIENT LE VECTEUR ASSOCIE A LA
!                  LA J-IEME VALEUR PROPRE DE W
!                  IF ICODE = 0, Z N'EST PAS UTILISE
! IN  IZ    : IS : 1-ERE DIMENSION (EXACTE) DE LA MATRICE Z
! LOC WK    : R8 : ZONE DE TRAVAIL DE LONGUEUR 2*NEQ
! IN  MXITER: IS : NOMBRE MAX D'ITERARION POUR LE QR
!                    (30 EST UN BON NOMBRE)
! OUT IER   : IS : PARAMETRE  D'ERREUR
!             IER = 0 OK
!             IER = J >0 , NON CONVERGENCE POUR LA J-IEME VALEUR PROPRE
!                LES J PREMIERES VALEURS PROPRES NE SONT PAS CALCULEES
! OUT NITQR : NOMBRE D'ITERATIONS QR POUR ATTEINDRE LA CONVERGENCE
!     ------------------------------------------------------------------
!     REFERENCE: F.L. BAUER - J.H. WILKINSON - C. REINSCH
!        HANDBOOK FOR AUTOMATIC COMPUTATION - LINEAR ALGEBRA - VOL.2
!        PAGE XXX
!     ------------------------------------------------------------------
    integer(kind=8) :: jer, iz2, k, l, i, n2, iiz, npi, jw, j, is, ig, igz
    real(kind=8) :: z11
!
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: mxiter
!-----------------------------------------------------------------------
    ier = 0
    jer = 0
    iz2 = iz+iz
    if (icode .lt. 0 .or. icode .gt. 1) then
        icode = 0
    end if
    if (icode .ne. 0 .and. iz .lt. neq) then
        icode = 0
    end if
!
    k = neq
    l = neq
    n2 = 2
    if (icode .eq. 0) n2 = 1
!
!     --- EQUILIBRAGE DE LA MATRICE INITIALE ---
!
    call vpzbal(mat, neq, mxeq, wk, k, &
                l)
!
!     --- MISE SOUS FORME DE HESSENBERG ---
!
    if (icode .eq. 0) then
        iiz = 1
    else
        iiz = neq
    end if
    if (l .ne. 0) then
!        IF L <> 0, MAT EST DEJA SOUS FORME DE HESSENBERG
        call vpzhes(mat, k, l, neq, neq, &
                    wk(1, n2))
    end if
!
!     --- TRANSFORMATION INVERSE HESSENBERG - EQUILIBRAGE
!     --- POUR LE CALCUL DES VECTEURS PROPRES
!
    if (icode .eq. 1) then
!        --- FORMATION DE LA MATRICE IDENTITE DANS Z
        do i = 1, neq*neq
            z(i) = 0.d0
        end do
        do i = 1, neq*neq, neq+1
            z(i) = 1.d0
        end do
        call vpzrbk(z, mat, wk(1, n2), neq, neq, &
                    k, l)
    end if
!
!     --- CALCUL DES VALEURS PROPRES (ET DES VECTEURS PROPRES)
!
    if (icode .eq. 0 .and. neq .eq. 1) then
        z11 = z(1)
    end if
    nitqr = 0
    call vpzqrh(mat, neq, neq, k, l, &
                w(1), w(neq+1), z, iiz, mxiter, &
                jer, nitqr)
    if (icode .eq. 0 .and. neq .eq. 1) then
        z(1) = z11
    end if
!
!     --- TRANSFORMATION INVERSE EQUILIBRAGE - MATRICE INITIALE
!     --- POUR LE CALCUL DES VECTEURS PROPRES
!
    if (jer .eq. 0 .and. icode .eq. 1) then
        call vpzech(wk, z, k, l, neq, &
                    neq, neq)
    end if
!
!     --- CONVERSION DES VALEURS PROPRES (W) EN FORMAT COMPLEXE
!
    do i = 1, neq
        npi = neq+i
        wk(i, 1) = w(npi)
    end do
    jw = neq+neq
    j = neq
    do i = 1, neq
        w(jw-1) = w(j)
        w(jw) = wk(j, 1)
        jw = jw-2
        j = j-1
    end do
!
!        TRAITEMENT DES VECTEURS PROPRES
!        IE. : LES CONVERTIR  EN FORMAT COMPLEXE DANS Z(IZ,NEQ)
!
    if (icode .eq. 1) then
!
        j = neq
60      continue
        if (j .lt. 1) goto 999
        if (w(j+j) .ne. 0.d0) then
!           TRANSLATER LA PAIRE DE VECTEURS COMPLEXES CONJUGUES
            is = iz2*(j-1)+1
            ig = neq*(j-2)+1
            igz = ig+neq
!           TRANSLATER LE VECTEUR COMPLEXE CONJUGE
            do i = 1, neq
                z(is) = z(ig)
                z(is+1) = -z(igz)
                is = is+2
                ig = ig+1
                igz = igz+1
            end do
!           TRANSLATER LE VECTEUR COMPLEXE
            is = iz2*(j-2)+1
            ig = is+iz2
            do i = 1, neq
                z(is) = z(ig)
                z(is+1) = -z(ig+1)
                is = is+2
                ig = ig+2
            end do
            j = j-2
            goto 60
        end if
!        TRANSLATER LE VECTEUR REEL
        is = iz2*(j-1)+neq+neq
        ig = neq*j
        do i = 1, neq
            z(is-1) = z(ig)
            z(is) = 0.d0
            is = is-2
            ig = ig-1
        end do
        j = j-1
        goto 60
    end if
!
999 continue
    ier = jer
end subroutine
