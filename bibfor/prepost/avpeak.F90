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
subroutine avpeak(jvalax, nbvec, nbordr, pseuil, iflag, &
                  npoin, jvalpo, jvalor)
! person_in_charge: van-xuan.tran at edf.fr
    implicit none
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
    integer(kind=8) :: nbvec, nbordr, npoin(nbvec), jvalor
    integer(kind=8) :: iflag(nbvec), jvalax, jvalpo
    real(kind=8) :: pseuil
! ----------------------------------------------------------------------
! BUT: EXTRAIRE LES PICS D'UNE FONCTION.
! ----------------------------------------------------------------------
! ARGUMENTS:
! VALAXE    IN   R  : VECTEUR CONTENANT L'HISTORIQUE DES PROJECTIONS
!                     POUR TOUS LES VECTEURS NORMAUX (n) ET TOUS LES
!                     NUMEROS D'ORDRE.
! NBVEC     IN   I  : NOMBRE DE VECTEURS NORMAUX.
! NBORDR    IN   I  : NOMBRE DE NUMERO D'ORDRE.
! PSEUIL    IN   R  : SEUIL DE DTECTION DES PICS
! IFLAG     IN   I  : VECTEUR DE DRAPEAUX QUI INDIQUE :
!                      - IFLAG(i) = 0 --> CAS GENERAL ;
!                      - IFLAG(i) = 1 --> CAS OU LES POINTS DANS LE
!                                         PLAN DE CISAILLEMENT SONT
!                                         ALIGNES VERTICALEMENT ;
!                      - IFLAG(i) = 2 --> CAS OU LES POINTS DANS LE
!                                         PLAN DE CISAILLEMENT SONT
!                                         ALIGNES HORIZONTALEMENT ;
!                      - IFLAG(i) = 3 --> CAS OU LES POINTS DANS LE
!                                         PLAN DE CISAILLEMENT SONT
!                                         CONTENUS DANS UN CADRE DE
!                                         COTES INFERIEURS A EPSILO.
! NPOIN     OUT  I  : NOMBRE DE PICS DETECTES POUR TOUS LES VECTEURS
!                     NORMAUX.
! VALPOI    OUT  R  : VALEUR DES PICS DETECTES POUR TOUS LES VECTEURS
!                     NORMAUX.
! VALORD    OUT  I  : NUMEROS D'ORDRE ASSOCIES AUX PICS DETECTES POUR
!                     TOUS LES VECTEURS NORMAUX.
!
!-----------------------------------------------------------------------
!     ------------------------------------------------------------------
    integer(kind=8) :: ivect, iordr, pass, sortie, ordmax, ordmin
!
    real(kind=8) :: vmin, vmax, valeur, epsilo
!
!-----------------------------------------------------------------------
!234567                                                              012
!-----------------------------------------------------------------------
    epsilo = 1.0d-7
!-----------------------------------------------------------------------
!
    call jemarq()
!
    do ivect = 1, nbvec
!
        if (iflag(ivect) .eq. 3) then
            goto 10
        end if
!
! ----- LE PREMIER POINT EST UN PIC -----
!
        npoin(ivect) = 1
        zr(jvalpo+(ivect-1)*nbordr+1) = zr(jvalax+(ivect-1)*nbordr+1)
        zi(jvalor+(ivect-1)*nbordr+1) = 1
        vmax = zr(jvalpo+(ivect-1)*nbordr+1)
        vmin = zr(jvalpo+(ivect-1)*nbordr+1)
        pass = 0
        sortie = 2
!
! ----- RECHERCHE DES PICS INTERMEDIAIRES -----
!
        do iordr = 2, nbordr
            valeur = zr(jvalax+(ivect-1)*nbordr+iordr)
            if (vmax .lt. valeur) then
                vmax = valeur
                ordmax = iordr
            end if
            if (vmin .gt. valeur) then
                vmin = valeur
                ordmin = iordr
            end if
            if (pass .eq. 0) then
                if ((valeur-vmin) .gt. pseuil) then
                    sortie = 1
                    pass = 1
                end if
                if ((vmax-valeur) .gt. pseuil) then
                    sortie = 0
                    pass = 1
                end if
            end if
            if ((sortie .eq. 1) .and. ((vmax-valeur)-pseuil .gt. epsilo)) then
                npoin(ivect) = npoin(ivect)+1
                zr(jvalpo+(ivect-1)*nbordr+npoin(ivect)) = vmax
                zi(jvalor+(ivect-1)*nbordr+npoin(ivect)) = ordmax
                vmin = valeur
                ordmin = iordr
                sortie = 0
            end if
            if ((sortie .eq. 0) .and. ((valeur-vmin)-pseuil .gt. epsilo)) then
                npoin(ivect) = npoin(ivect)+1
                zr(jvalpo+(ivect-1)*nbordr+npoin(ivect)) = vmin
                zi(jvalor+(ivect-1)*nbordr+npoin(ivect)) = ordmin
                vmax = valeur
                ordmax = iordr
                sortie = 1
            end if
        end do
!
        if (sortie .eq. 0) then
            npoin(ivect) = npoin(ivect)+1
            zr(jvalpo+(ivect-1)*nbordr+npoin(ivect)) = vmin
            zi(jvalor+(ivect-1)*nbordr+npoin(ivect)) = ordmin
        end if
        if (sortie .eq. 1) then
            npoin(ivect) = npoin(ivect)+1
            zr(jvalpo+(ivect-1)*nbordr+npoin(ivect)) = vmax
            zi(jvalor+(ivect-1)*nbordr+npoin(ivect)) = ordmax
        end if
10      continue
    end do
!
!
    call jedema()
end subroutine
