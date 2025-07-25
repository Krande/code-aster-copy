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

subroutine resou2(matass, matpre, solveu, chcine, nsecm, &
                  chsecm, chsolu, base, rsolu, csolu, &
                  criter, prepos, istop, iret)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/elg_kellag.h"
#include "asterfort/elg_resoud.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/resou1.h"
!-----------------------------------------------------------------------
!
    character(len=*) :: matass, matpre, solveu, chcine
    integer(kind=8) :: nsecm
    character(len=*) :: chsecm, chsolu, base
    real(kind=8) :: rsolu(*)
    complex(kind=8) :: csolu(*)
    character(len=*) :: criter
    aster_logical :: prepos
    integer(kind=8) :: istop, iret
!
!-----------------------------------------------------------------------
! BUT : RESOUDRE UN SYSTEME LINEAIRE D'EQUATIONS (REEL OU COMPLEXE)
!-----------------------------------------------------------------------
!
! ARGUMENTS :
!------------
!
! REMARQUES : ON PEUT APPELER RESOUD DE 2 FACONS
!   1) AVEC NSECM = 0 + CHSECM, CHSOLU, BASE
!   2) AVEC NSECM > 0 + RSOLU (OU CSOLU) + (CHSECM=CHSOLU=' ')
!
! IN/JXIN  K19 MATASS : MATR_ASSE PREMIER MEMBRE DU SYSTEME LINEAIRE
! IN/JXIN  K19 MATPRE : MATR_ASSE DE PRECONDITIONNEMENT
!                       POUR SOLVEUR ITERATIF GCPC (OU ' ' SINON)
! IN/JXIN  K19 SOLVEU : SD_SOLVEUR (OU ' ')
!                       SI SOLVEU=' ' ON PREND LE SOLVEUR DE MATASS
! IN/JXIN  K*  CHCINE : CHAMP ASSOCIE AUX CHARGES CINEMATIQUES (OU ' ')
! IN       I   NSECM  : / 0 => ON UTILISE CHSECM, CHSOLU, BASE
!                       / N => ON UTILISE RSOLU (OU CSOLU)
!                         N : NOMBRE DE SECONDS MEMBRES
! IN/JXIN  K*  CHSECM : CHAMP SECOND MEMBRE DU SYSTEME LINEAIRE
! IN/JXOUT K*  CHSOLU : CHAMP SOLUTION DU SYSTEME LINEAIRE
! IN       K*  BASE   : BASE SUR LAQUELLE ON CREE CHSOLU
! IN/OUT   R   RSOLU  : TABLEAU (*,NSECM)
!           EN ENTREE : VECTEUR DE REELS CONTENANT LES SECONDS MEMBRES
!           EN SORTIE : VECTEUR DE REELS CONTENANT LES SOLUTIONS
! IN/OUT   C   CSOLU  : TABLEAU (*,NSECM)
!                       IDEM RSOLU POUR LES COMPLEXES
! IN/JXOUT K*  CRITER : SD_CRITER (CRITERES DE CONVERGENCE)
!                       POUR SOLVEUR ITERATIF GCPC (OU ' ' SINON)
! IN       L   PREPOS : / .TRUE.  => ON FAIT LES PRE ET POST-TRAITEMENTS
!                                    DU SMB ET DE LA SOLUTION
!                       / .FALSE. => ON NE FAIT AUCUN TRAITEMENT
!                                    (EN MODAL PAR EXEMPLE)
! IN       I   ISTOP  : COMPORTEMENT EN CAS D'ERREUR (CE PARAMETRE N'A
!                       D'UTILITE QUE POUR UN SOLVEUR ITERATIF)
!                       / 0     : ON S'ARRETE EN <F>
!                       / 2     : ON CONTINUE SANS MESSAGE D'ERREUR
!                       / -9999 : ON PREND LA VALEUR DEFINIE DANS LA
!                                 SD_SOLVEUR POUR STOP_SINGULIER
! OUT      I   IRET   : CODE RETOUR
!                       / 0 : OK (PAR DEFAUT POUR SOLVEURS DIRECTS)
!                       / 1 : ECHEC (NOMBRE MAX. D'ITERATIONS ATTEINT)
!-----------------------------------------------------------------------
! Cette routine est une surcouche de la routine resou1.
! Elle est necessaire pour traiter le cas ELIM_LAGR='OUI'
! ----------------------------------------------------------------------
    character(len=19) :: matas1, solve1
    character(len=3) :: kellag
! ----------------------------------------------------------------------
!
    call jemarq()
    matas1 = matass
    solve1 = solveu
!
!
!   1. CALCUL DE KELLAG :
!   -------------------------------------
    call elg_kellag(matas1, solve1, kellag)
!
!
!   2. SI ELIM_LAGR /= 'OUI', ON APPELLE SIMPLEMENT RESOU1 :
!   --------------------------------------------------------
    if (.not. (kellag .eq. 'OUI')) then
        call resou1(matas1, matpre, solve1, chcine, nsecm, &
                    chsecm, chsolu, base, rsolu, csolu, &
                    criter, prepos, istop, iret)
    else
        call elg_resoud(matas1, matpre, nsecm, chsecm, chsolu, &
                        base, rsolu, csolu, criter, prepos, &
                        istop, iret)
    end if
!
!
    call jedema()
end subroutine
