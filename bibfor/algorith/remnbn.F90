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
subroutine remnbn(basmod, nbmod, nbddr, nbdax, flexdr, &
                  flexga, flexax, tetgd, tetax, cmode, &
                  vecmod, neq, beta)
    implicit none
#include "jeveux.h"
#include "asterfort/dcapno.h"
#include "asterfort/jedema.h"
#include "asterfort/jelibe.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
    integer(kind=8) :: nbmod, nbddr, nbdax, neq
    real(kind=8) :: beta
    complex(kind=8) :: cmode(nbmod+nbddr+nbdax), vecmod(neq)
    character(len=8) :: basmod
    character(len=24) :: flexdr, flexga, flexax, tetgd, tetax
!-----------------------------------------------------------------------
!
!  BUT:     <RESTITUTION MAC-NEAL DE BAS NIVEAU>
!
!  CALCUL DU VECTEUR COMPLEXE EN DDL PHYSIQUES A PARTIR DU MODE
!   COMPLEXE ISSU DU CALCUL CYCLIQUE ET DES VECTEURS DES NUMERO ORDRE
!  DES DEFORMEES RELATIVES AUX INTERFACES DROITE GAUCHE ET AXE DE TYPE
!  MAC-NEAL
!
!-----------------------------------------------------------------------
!
! BASMOD   /I/: NOM UT DE LA BASE MODALE EN AMONT
! NBMOD    /I/: NOMBRE DE MODES PROPRES UTILISES POUR CALCUL CYCLIQUE
! NBDDR    /I/: NOMBRE DE DEFORMEES INTERFACE DROITE (ET GAUCHE)
! NBDAX    /I/: NOMBRE DE DEFORMEES INTERFACE AXE
! FLEXDR   /I/: NOM K24 FLEXIBILITE TOUTES-LIGNES/COLONNES-DROITES
! FELXGA   /I/: NOM K24 FELXIBILITE TOUTES-LIGNES/COLONNES-GAUCHES
! FLEXAX   /I/: NOM K24 FELXIBILITE TOUTES-LIGNES/COLONNES-AXES
! TETGD    /I/: NOM K24 DE LA MATRICE DE PASSAGE GAUCHE-DROITE
! TETAX    /I/: NOM K24 DE LA MATRICE DE PASSAGE AXE-AXE
! IORD     /I/: VECTEUR DES NUMEROS ORDRE DEFORMEES DE DROITE
! IORG     /I/: VECTEUR DES NUMEROS ORDRE DES DEFORMEES DE GAUCHE
! IORA     /I/: VECTEUR DES NUMEROS ORDRE DES DEFORMEES  AXE
! CMODE    /I/: MODE COMPLEXES ISSU DU CALCUL CYCLIQUE
! VECMOD   /I/: VECTEUR MODAL COMPLEXE EN DDL PHYSIQUE
! NEQ      /I/: NOMBRE DE DDL PHYSIQUES ASSEMBLES
! BETA     /I/: DEPHASAGE INTER-SECTEUR
!
!
!
!
    integer(kind=8) :: i, j, iad, llcham, llfdr, llfga, lltgd, lltax, llfax
    real(kind=8) :: abeta, bbeta
    complex(kind=8) :: dephc, cfact, cmult
    character(len=24) :: chaval
!
!-----------------------------------------------------------------------
!
    call jemarq()
!
! --- MISE A ZERO DU MODE PROPRES RESULTAT
!
    do i = 1, neq
        vecmod(i) = dcmplx(0.d0, 0.d0)
    end do
!
    abeta = cos(beta)
    bbeta = sin(beta)
    dephc = dcmplx(abeta, bbeta)
!
! --- CONTRIBUTION DES MODES PROPRES
!
    do i = 1, nbmod
        call dcapno(basmod, 'DEPL    ', i, chaval)
        call jeveuo(chaval, 'L', llcham)
        do j = 1, neq
            cfact = dcmplx(zr(llcham+j-1), 0.d0)
            vecmod(j) = vecmod(j)+cmode(i)*cfact
        end do
        call jelibe(chaval)
    end do
!
! --- CONTRIBUTION DES DEFORMEES DE DROITE
!
    call jeveuo(flexdr, 'L', llfdr)
    do i = 1, nbddr
        do j = 1, neq
            iad = llfdr+((i-1)*neq)+j-1
            cfact = dcmplx(zr(iad), 0.d0)
            vecmod(j) = vecmod(j)+cmode(i+nbmod)*cfact
        end do
    end do
    call jelibe(flexdr)
!
! --- CONTRIBUTION DES DEFORMEES DE GAUCHE
!
    call jeveuo(tetgd, 'L', lltgd)
    call jeveuo(flexga, 'L', llfga)
    do i = 1, nbddr
!
        cmult = dcmplx(0.d0, 0.d0)
        do j = 1, nbddr
            iad = lltgd+((j-1)*nbddr)+i-1
            cfact = dcmplx(zr(iad), 0.d0)*cmode(j+nbmod)
            cmult = cmult-cfact
        end do
!
        iad = llfga+((i-1)*neq)
        do j = 1, neq
            cfact = dcmplx(zr(iad+j-1), 0.d0)
            vecmod(j) = vecmod(j)+dephc*cfact*cmult
        end do
!
    end do
    call jelibe(flexga)
    call jelibe(tetgd)
!
! --- EVENTUELLE CONTRIBUTION DE L'AXE
!
    if (nbdax .gt. 0) then
        call jeveuo(tetax, 'L', lltax)
        call jeveuo(flexax, 'L', llfax)
        do i = 1, nbdax
!
            cmult = dcmplx(1.d0, 0.d0)
            do j = 1, nbdax
                iad = lltax+((j-1)*nbdax)+i-1
                cfact = dcmplx(zr(iad), 0.d0)*cmode(nbmod+nbddr+j)
                cmult = cmult-(dephc*cfact)
            end do
!
            iad = llfax+((i-1)*neq)
            do j = 1, neq
                cfact = dcmplx(2*zr(iad+j-1), 0.d0)
                vecmod(j) = vecmod(j)+cfact*cmult
            end do
!
        end do
        call jelibe(flexax)
        call jelibe(tetax)
    end if
!
    call jedema()
end subroutine
