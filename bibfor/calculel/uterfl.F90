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
subroutine uterfl(ndim, iflup, iflum, ino, mno, &
                  jno, nsomm, jac, term22, aux, &
                  ltheta, valthe, valunt, niv, ifm, &
                  xn, yn, zn, valfp, valfm, &
                  ityp, noe)
! person_in_charge: olivier.boiteau at edf.fr
!-----------------------------------------------------------------------
!    - FONCTION REALISEE: UTILITAIRE DE CALCUL DE L'ERREUR DUE A LA
!                         CONDITION DE FLUX. POUR AERER TE0003.
!
! IN NDIM : DIMENSION DU CALCUL
! IN IFLUP/M : ADRESSE JEVEUX DU FLUX +/-
! IN INO/JNO/MNO : NUMERO DES NOEUDS DE L'ARETE (INO NUMERO FACE EN 3D)
! IN NSOMM  : NBRE DE SOMMETS DE L'ARETE OU DE LA FACE.
! IN JAC  : JACOBIEN
! IN LTHETA/VALTHE/UNT : PARAMETRES THETA METHODE
! IN NIV/IFM : PARAMETRES D'IMPRESSION
! IN ITYP : TYPE DE FACE (EN 3D)
! IN XN/YN/ZN : COMPOSANTES DE LA NORMALE AUX POINTS D'INTEGRATION
! IN VALFP/M  : VALEURS DU FLUX AUX INSTANTS + ET -
! IN NOE : TABLEAU NUMEROS DE NOEUDS PAR FACE ET PAR TYPE D'ELEMENT 3D
! OUT TERM22 : CONTRIBUTION DE L'ERREUR
! OUT AUX : TERME DE NORMALISATION
!   -------------------------------------------------------------------
!     SUBROUTINES APPELLEES:
!       AUCUNE.
!
!     FONCTIONS INTRINSEQUES:
!       AUCUNE.
!   -------------------------------------------------------------------
!     ASTER INFORMATIONS:
!       24/09/01 (OB): CREATION POUR SIMPLIFIER TE0003.F.
!----------------------------------------------------------------------
! CORPS DU PROGRAMME
! aslint: disable=W1504
    implicit none
!
! DECLARATION PARAMETRES D'APPELS
#include "asterf_types.h"
#include "jeveux.h"
    integer(kind=8) :: iflup, iflum, ndim, ino, mno, jno, nsomm, ityp, noe(9, 6, 3), ifm
    integer(kind=8) :: niv
    real(kind=8) :: jac(9), term22, aux, valthe, valunt, xn(9), yn(9), zn(9)
    real(kind=8) :: valfp(9), valfm(9)
    aster_logical :: ltheta
!
!
! DECLARATION VARIABLES LOCALES
    integer(kind=8) :: ij, in, iino
    real(kind=8) :: aux1, aux2, aux3, aux4, term23
!
! INIT.
    term22 = 0.d0
    aux = 0.d0
!
    if (ndim .eq. 2) then
!
! CAS 2D
! POINT DE DEPART: INO
        ij = iflup+(ino-1)*ndim
        aux1 = valthe*valfp(1)
        aux2 = valthe*(zr(ij)*xn(1)+zr(ij+1)*yn(1))
        term23 = aux1+aux2
        if (niv .eq. 2) write (ifm, *) ' ZR INO P', aux1, aux2
        if (ltheta) then
            ij = ij+iflum-iflup
            aux3 = valunt*valfm(1)
            aux4 = valunt*(zr(ij)*xn(1)+zr(ij+1)*yn(1))
            aux1 = aux1+aux3
            term23 = term23+aux3+aux4
            if (niv .eq. 2) write (ifm, *) ' ZR INO M', aux3, aux4
        end if
        term22 = jac(1)*term23*term23
        aux = jac(1)*aux1*aux1
!
! POINT EXTREME: JNO
        ij = iflup+(jno-1)*ndim
        aux1 = valthe*valfp(2)
        aux2 = valthe*(zr(ij)*xn(2)+zr(ij+1)*yn(2))
        term23 = aux1+aux2
        if (niv .eq. 2) write (ifm, *) ' ZR JNO P', aux1, aux2
        if (ltheta) then
            ij = ij+iflum-iflup
            aux3 = valunt*valfm(2)
            aux4 = valunt*(zr(ij)*xn(2)+zr(ij+1)*yn(2))
            aux1 = aux1+aux3
            term23 = term23+aux3+aux4
            if (niv .eq. 2) write (ifm, *) ' ZR JNO M', aux3, aux4
        end if
        term22 = term22+jac(2)*term23*term23
        aux = aux+jac(2)*aux1*aux1
!
! POINT MILIEU SI NECESSAIRE: MNO
        if (nsomm .eq. 3) then
            ij = iflup+(mno-1)*ndim
            aux1 = valthe*valfp(3)
            aux2 = valthe*(zr(ij)*xn(3)+zr(ij+1)*yn(3))
            term23 = aux1+aux2
            if (niv .eq. 2) write (ifm, *) ' ZR MNO P', aux1, aux2
            if (ltheta) then
                ij = ij+iflum-iflup
                aux3 = valunt*valfm(3)
                aux4 = valunt*(zr(ij)*xn(3)+zr(ij+1)*yn(3))
                aux1 = aux1+aux3
                term23 = term23+aux3+aux4
                if (niv .eq. 2) write (ifm, *) ' ZR MNO M', aux3, aux4
            end if
            term22 = term22+jac(3)*term23*term23
            aux = aux+jac(3)*aux1*aux1
        end if
    else
!
! CAS 3D
        do in = 1, nsomm
            iino = noe(in, ino, ityp)
            ij = iflup+(iino-1)*ndim
            aux1 = valthe*valfp(in)
            aux2 = valthe*(zr(ij)*xn(in)+zr(ij+1)*yn(in)+zr(ij+2)*zn(in))
            term23 = aux1+aux2
            if (niv .eq. 2) write (ifm, *) ' ZR IINO P/IN ', aux1, aux2, in
            if (ltheta) then
                ij = ij+iflum-iflup
                aux3 = valunt*valfm(in)
                aux4 = valunt*(zr(ij)*xn(in)+zr(ij+1)*yn(in)+zr(ij+2)*zn(in))
                aux1 = aux1+aux3
                term23 = term23+aux3+aux4
                if (niv .eq. 2) write (ifm, *) ' ZR JNO M    ', aux3, aux4
            end if
            term22 = term22+term23*term23*jac(in)
            aux = aux+aux1*aux1*jac(in)
        end do
!
    end if
!
end subroutine
