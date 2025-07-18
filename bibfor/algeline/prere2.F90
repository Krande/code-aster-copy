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

subroutine prere2(solveu, base, iret, matpre, matass, &
                  npvneg, istop)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elg_kellag.h"
#include "asterfort/elg_preres.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/prere1.h"
!-----------------------------------------------------------------------
    integer(kind=8) :: npvneg, istop, iret
    character(len=1) :: base
    character(len=*) :: matass, matpre, solveu
!-----------------------------------------------------------------------
! person_in_charge: jacques.pellet at edf.fr
! BUT : FACTORISER UNE MATR_ASSE (LDLT/MULT_FRONT/MUMPS)
!       OU FABRIQUER UNE MATRICE DE PRECONDITIONNEMENT (GCPC,PETSC)
!
! SOLVEU (K19) IN : OBJET SOLVEUR (OU ' ')
! BASE (K1)    IN : BASE SUR LAQUELLE ON CREE LA MATRICE FACTORISEE
!                  (OU LA MATRICE DE PRECONDITIONNEMENT)
! IRET (I)     OUT : CODE_RETOUR :
!             /0 -> OK (PAR DEFAUT AVEC GCPC/PETSC)
!             /2 -> LA FACTORISATION N'A PAS PU SE FAIRE
!                   JUSQU'AU BOUT.
!             /1 -> LA FACTORISATION EST ALLEE AU BOUT
!                   MAIS ON A PERDU BEAUCOUP DE DECIMALES
!             /3 -> LA FACTORISATION EST ALLEE AU BOUT
!                   MAIS ON NE SAIT PAS DIRE SI ON A PERDU DES DECIMALES
!
! MATPRE(K19) IN/JXVAR : MATRICE DE PRECONDITIONNEMENT (SI GCPC)
! MATASS(K19) IN/JXVAR : MATRICE A FACTORISER OU A PRECONDITIONNER
! NPVNEG (I) OUT : NBRE DE TERMES DIAGONAUX NEGATIFS DE LA FACTORISEE
!          CE NBRE N'EST LICITE QUE SI LA MATRICE EST REELLE SYMETRIQUE
!          ET N'EST FOURNI QUE PAR UN SOLVEUR DIRECT: LDLT, MF OU MUMPS
! ISTOP (I)  IN: COMPORTEMENT EN CAS DE DETECTION DE SINGULARITE. CE
!                PARAMETRE N'A D'UTILITE QU'AVEC UN SOLVEUR DIRECT
!                  /0 -> SI IRET>0 : ERREUR <F>
!                  /1 -> SI IRET=1 : ALARME <A>
!                        SI IRET=2 : ERREUR <F>
!                  /2 -> LE PROGRAMME NE S'ARRETE PAS
!                        ET N'IMPRIME AUCUN MESSAGE.
!                 /-9999 -> ON PREND LA VALEUR PREVUE DS LA SD_SOLVEUR
!                        POUR STOP_SINGULIER (VALEUR 0 OU 1 SEULEMENT)
!                 /AUTRE --> ASSERT
!-----------------------------------------------------------------------
! cette routine est une surcouche de la routine prere1.
! elle est necessaire pour traiter le cas elim_lagr='oui'
!----------------------------------------------------------------------
    character(len=3) :: kellag
    character(len=19) :: solve1, matas1
!----------------------------------------------------------------------
    call jemarq()
!
    matas1 = matass
    solve1 = solveu
!
!   1. CALCUL DE KELLAG :
!   -------------------------------------
    call elg_kellag(matas1, solve1, kellag)
!
!
!   2. SI ELIM_LAGR /= 'OUI', ON APPELLE SIMPLEMENT PRERE1 :
!   --------------------------------------------------------
    if (.not. (kellag .eq. 'OUI')) then
        call prere1(solve1, base, iret, matpre, matas1, &
                    npvneg, istop)
    else
        call elg_preres(solve1, base, iret, matpre, matas1, &
                        npvneg, istop)
!
    end if
!
    call jedema()
end subroutine
