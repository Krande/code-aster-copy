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
subroutine aveppr(nbordr, vwork, tdisp, kwork, sommw, &
                  tspaq, i, jveppr, jvsipn)
!
    implicit none
#include "jeveux.h"
#include "asterfort/jacobi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
!
    integer(kind=8) :: nbordr, tdisp, kwork, sommw, tspaq, i
    real(kind=8) :: vwork(tdisp)
    integer(kind=8) :: jveppr, jvsipn
!, veppr(nbordr), vsipn(nbordr)
! ----------------------------------------------------------------------
! BUT: CALCULER LA DEF PLAS PRINCIPALE ET CONT NOMRMALE ASSOCIE
! ----------------------------------------------------------------------
! ARGUMENTS :
!  NBORDR : IN   I  : NOMBRE DE NUMEROS D'ORDRE.
!  VWORK  : IN   R  : VECTEUR DE TRAVAIL CONTENANT
!                     L'HISTORIQUE DES TENSEURS DES CONTRAINTES
!                     ATTACHES A CHAQUE POINT DE GAUSS DES MAILLES
!                     DU <<PAQUET>> DE MAILLES.
!  TDISP  : IN   I  : TAILLE DU VECTEUR DE TRAVAIL.
!  KWORK  : IN   I  : KWORK = 0 ON TRAITE LA 1ERE MAILLE DU PAQUET
!                               MAILLES OU LE 1ER NOEUD DU PAQUET DE
!                               NOEUDS;
!                     KWORK = 1 ON TRAITE LA IEME (I>1) MAILLE DU PAQUET
!                               MAILLES OU LE IEME NOEUD DU PAQUET
!                               DE NOEUDS.
!  SOMMW  : IN   I  : SOMME DES POINTS DE GAUSS OU DES NOEUDS DES N
!                     MAILLES PRECEDANT LA MAILLE COURANTE.
!  TSPAQ  : IN   I  : TAILLE DU SOUS-PAQUET DU <<PAQUET>> DE MAILLES
!                     OU DE NOEUDS COURANT.
!  I      : IN   I  : IEME POINT DE GAUSS OU IEME NOEUD.
! JEPPR : OUT  R  : VECTEUR CONTENANT LA  A
!                     TOUS LES INSTANTS.
! ----------------------------------------------------------------------
    integer(kind=8) :: iordr, adrs, nvp, nperm, nitjac, j, iordre, itype
    integer(kind=8) :: decal
    real(kind=8) :: sig(6), tol, toldyn, ar(6), br(6)
    real(kind=8) :: vecpro(3, 3), valpro(3), nm1x, nm1y, nm1z, sxm, sym, szm
    real(kind=8) :: jacaux(3), epsp(6)
!     ------------------------------------------------------------------
!
!234567                                                              012
!
    call jemarq()
!
! ---------------------------------------------------------------
! RECUPER LES CONTRAINTES ET DEFORMATION
!
    decal = 18
    do iordr = 1, nbordr
        adrs = (iordr-1)*tspaq+kwork*sommw*decal+(i-1)*decal
!
        do j = 1, 6
            sig(j) = vwork(adrs+j)
            epsp(j) = vwork(adrs+j+6+6)
        end do
!
        nvp = 3
        nperm = 12
        tol = 1.d-10
        toldyn = 1.d-2
        itype = 0
        iordre = 1
        ar(1) = epsp(1)
        ar(2) = epsp(4)
        ar(3) = epsp(5)
        ar(4) = epsp(2)
        ar(5) = epsp(6)
        ar(6) = epsp(3)
        br(1) = 1.d0
        br(2) = 0.d0
        br(3) = 0.d0
        br(4) = 1.d0
        br(5) = 0.d0
        br(6) = 1.d0
!
        call jacobi(nvp, nperm, tol, toldyn, ar, &
                    br, vecpro, valpro, jacaux, nitjac, &
                    itype, iordre)
!
        zr(jveppr+iordr) = valpro(1)
!
        nm1x = vecpro(1, 1)
        nm1y = vecpro(2, 1)
        nm1z = vecpro(3, 1)
! CALCvect_F = [SIG].vect_n
!
        sxm = sig(1)*nm1x+sig(4)*nm1y+sig(5)*nm1z
        sym = sig(4)*nm1x+sig(2)*nm1y+sig(6)*nm1z
        szm = sig(5)*nm1x+sig(6)*nm1y+sig(3)*nm1z
!
! CALCNORM = vect_F.vect_n
!
        zr(jvsipn+iordr) = abs(sxm*nm1x+sym*nm1y+szm*nm1z)
!
    end do
!
    call jedema()
!
end subroutine
