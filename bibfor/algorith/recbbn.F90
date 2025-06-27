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
subroutine recbbn(basmod, nbmod, nbddr, nbdax, tetgd, &
                  iord, iorg, iora, cmode, vecmod, &
                  neq, beta)
    implicit none
#include "jeveux.h"
#include "asterfort/dcapno.h"
#include "asterfort/jedema.h"
#include "asterfort/jelibe.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
    integer(kind=8) :: nbmod, nbddr, nbdax, iord(nbddr), iorg(nbddr), iora(nbdax), neq
    real(kind=8) :: beta
    complex(kind=8) :: cmode(nbmod+nbddr+nbdax), vecmod(neq)
    character(len=8) :: basmod
    character(len=24) :: tetgd
!-----------------------------------------------------------------------
!
!  BUT:        < RESTITUTION CRAIG_BAMPTON BAS NIVEAU >
!
!    CALCUL DU VECTEUR COMPLEXE EN DDL PHYSIQUES A PARTIR DU MODE
!   COMPLEXE ISSU DU CALCUL CYCLIQUE ET DES VECTEURS DES NUMERO ORDRE
!  DES DEFORMEES RELATIVES AUX INTERFACES DROITE GAUCHE ET AXE DE TYPE
!   CRAIG-BAMPTON
!  AINSI QUE DE LA MATRICE TETGD
!
!-----------------------------------------------------------------------
!
! BASMOD   /I/: NOM UT DE LA BASE MODALE EN AMONT
! NBMOD    /I/: NOMBRE DE MODES PROPRES UTILISES POUR CALCUL CYCLIQUE
! NBDDR    /I/: NOMBRE DE DEFORMEES INTERFACE DROITE (ET GAUCHE)
! NBDAX    /I/: NOMBRE DE DEFORMEES INTERFACE AXE
! TETGD    /I/: NOM K24 DE MATRICE DE PASSAGE GAUCHE-DROITE
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
    integer(kind=8) :: i, j, iad, llcham, lltgd
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
    do i = 1, nbddr
        call dcapno(basmod, 'DEPL    ', iord(i), chaval)
        call jeveuo(chaval, 'L', llcham)
        do j = 1, neq
            cfact = dcmplx(zr(llcham+j-1), 0.d0)
            vecmod(j) = vecmod(j)+cmode(i+nbmod)*cfact
        end do
        call jelibe(chaval)
    end do
!
! --- CONTRIBUTION DES DEFORMEES DE GAUCHE
!
    call jeveuo(tetgd, 'L', lltgd)
!
    do i = 1, nbddr
        call dcapno(basmod, 'DEPL    ', iorg(i), chaval)
        call jeveuo(chaval, 'L', llcham)
!
        cmult = dcmplx(0.d0, 0.d0)
        do j = 1, nbddr
            iad = lltgd+((j-1)*nbddr)+i-1
            cfact = dcmplx(zr(iad), 0.d0)*cmode(j+nbmod)
            cmult = cmult+cfact
        end do
!
        do j = 1, neq
            cfact = dcmplx(zr(llcham+j-1), 0.d0)
            vecmod(j) = vecmod(j)+dephc*cfact*cmult
        end do
!
        call jelibe(chaval)
    end do
    call jelibe(tetgd)
!
! --- EVENTUELLE CONTRIBUTION DE L'AXE
!
    if (nbdax .gt. 0) then
        do i = 1, nbdax
            call dcapno(basmod, 'DEPL    ', iora(i), chaval)
            call jeveuo(chaval, 'L', llcham)
            do j = 1, neq
                cfact = dcmplx(2*zr(llcham+j-1), 0.d0)
                vecmod(j) = vecmod(j)+cfact*cmode(i+nbmod+nbddr)
            end do
            call jelibe(chaval)
        end do
    end if
!
    call jedema()
end subroutine
