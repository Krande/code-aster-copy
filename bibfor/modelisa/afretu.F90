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

subroutine afretu(iprno, lonlis, klisno, noepou, noma, &
                  vale1, nbcoef, idec, coef, nomddl, &
                  lisrel)
    implicit none
#include "jeveux.h"
#include "asterfort/afrela.h"
#include "asterfort/dismoi.h"
#include "asterfort/imprel.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/char8_to_int.h"
!
    integer(kind=8) :: lonlis, iprno(*), idec, nbcoef
    real(kind=8) :: coef(nbcoef)
    character(len=8) :: klisno(lonlis), noepou, noma, nomddl(nbcoef)
    character(len=24) :: vale1
    character(len=19) :: lisrel
!     RACCORD (COQUE OU 3D)_TUYAU : UNE RELATION LINEAIRE
!
    integer(kind=8) :: ino, ival, idch1, nbterm, i, nbec

    character(len=8) :: betaf
    character(len=16) :: motfac
    real(kind=8) :: beta
    complex(kind=8) :: betac
    complex(kind=8), pointer :: coec(:) => null()
    real(kind=8), pointer :: coer(:) => null()
    integer(kind=8), pointer :: dime(:) => null()
    real(kind=8), pointer :: direct(:) => null()
    character(len=8), pointer :: lisddl(:) => null()
    character(len=8), pointer :: lisno(:) => null()
!
    call jemarq()
!
    motfac = 'LIAISON_ELEM'
    betaf = '&FOZERO'
    beta = 0.0d0
    betac = (0.0d0, 0.0d0)
    motfac = 'LIAISON_ELEM'
    call jeveuo(vale1, 'L', idch1)
    call dismoi('NB_EC', 'DEPL_R', 'GRANDEUR', repi=nbec)
    if (nbec .gt. 10) then
        call utmess('F', 'MODELISA_94')
    end if
!
! --- CREATION DES TABLEAUX DE TRAVAIL NECESSAIRES A L'AFFECTATION
! --- NBTERM : MAJORANT DU NOMBRE DE TERMES DANS UNE RELATION
!
    nbterm = 3*lonlis+nbcoef
! ---     VECTEUR DU NOM DES NOEUDS
    AS_ALLOCATE(vk8=lisno, size=nbterm)
! ---     VECTEUR DU NOM DES DDLS
    AS_ALLOCATE(vk8=lisddl, size=nbterm)
! ---     VECTEUR DES COEFFICIENTS REELS
    AS_ALLOCATE(vr=coer, size=nbterm)
! ---     VECTEUR DES COEFFICIENTS COMPLEXES
    AS_ALLOCATE(vc=coec, size=nbterm)
! ---     VECTEUR DES DIRECTIONS DES DDLS A CONTRAINDRE
    AS_ALLOCATE(vr=direct, size=3*nbterm)
! ---     VECTEUR DES DIMENSIONS DE CES DIRECTIONS
    AS_ALLOCATE(vi=dime, size=nbterm)
!
!        RELATIONS ENTRE LES NOEUDS DE COQUE ET LE NOEUD NOEPOU
!        DE TUYAU SUR LE DDL NOMDDL
!
    do i = 1, lonlis
        ino = char8_to_int(klisno(i))
!           ADRESSE DE LA PREMIERE CMP DU NOEUD INO DANS LES CHAMNO
        ival = iprno((ino-1)*(nbec+2)+1)
!
        lisno(1+3*(i-1)+1-1) = klisno(i)
        lisno(1+3*(i-1)+2-1) = klisno(i)
        lisno(1+3*(i-1)+3-1) = klisno(i)
!
        lisddl(1+3*(i-1)+1-1) = 'DX'
        lisddl(1+3*(i-1)+2-1) = 'DY'
        lisddl(1+3*(i-1)+3-1) = 'DZ'
!
! RACCORD  3D_TUYAU : IDEC=0 DANS TOUS LES APPELS A AFRETU
! RACCORD COQ_TUYAU : IDEC=0 OU 3 DANS LES APPELS A AFRETU
!
        if (idec .ge. 0) then
            coer(1+3*(i-1)+1-1) = zr(idch1+ival-1+idec+0)
            coer(1+3*(i-1)+2-1) = zr(idch1+ival-1+idec+1)
            coer(1+3*(i-1)+3-1) = zr(idch1+ival-1+idec+2)
        end if
!
!   Gestion des cas WI1 et WO1, où les conditions d'orthogonalité
!   imposent wi1 = vo1 et wo1 = -vi1 comme coeff. devant les cos(phi) et sin(phi)
!   L'entier "idec" est détourné de son utilisation
!
!       Cas WI1 = VO1
        if (idec .eq. -1) then
            ival = 6*(i-1)
            coer(1+3*(i-1)+1-1) = zr(idch1+ival+0+0)
            coer(1+3*(i-1)+2-1) = zr(idch1+ival+0+1)
            coer(1+3*(i-1)+3-1) = zr(idch1+ival+0+2)
        end if
!
!       Cas WO1 = -VI1
        if (idec .eq. -2) then
            ival = 6*(i-1)
            coer(1+3*(i-1)+1-1) = zr(idch1+ival+3+0)
            coer(1+3*(i-1)+2-1) = zr(idch1+ival+3+1)
            coer(1+3*(i-1)+3-1) = zr(idch1+ival+3+2)
        end if
    end do
!
    do i = 1, nbcoef
        lisno(1+3*lonlis+i-1) = noepou
        lisddl(1+3*lonlis+i-1) = nomddl(i)
        coer(1+3*lonlis+i-1) = coef(i)
    end do
!
    call afrela(coer, coec, lisddl, lisno, dime, &
                direct, nbterm, beta, betac, betaf, &
                'REEL', 'REEL', 0.d0, lisrel)
    call imprel(motfac, nbterm, coer, lisddl, lisno, &
                beta)
!
!
! --- MENAGE
!
    AS_DEALLOCATE(vk8=lisno)
    AS_DEALLOCATE(vk8=lisddl)
    AS_DEALLOCATE(vr=coer)
    AS_DEALLOCATE(vc=coec)
    AS_DEALLOCATE(vr=direct)
    AS_DEALLOCATE(vi=dime)
!
    call jedema()
end subroutine
