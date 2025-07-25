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

subroutine simono()
    implicit none
!
!     OPERATEUR :   CALC_CHAR_SEISME
!
!     CREE LE VECTEUR SECOND MEMBRE DANS LE CAS D'UN CALCUL SISMIQUE
!     STRUCTURE : MONO-APPUI
!
!     ------------------------------------------------------------------
!
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mrmult.h"
#include "asterfort/mtdscr.h"
#include "asterfort/pteddl.h"
#include "asterfort/utmess.h"
#include "asterfort/vtcrem.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
    integer(kind=8) :: lmat, neq
    real(kind=8) :: xnorm, depl(6)
    character(len=8) :: tabcmp(6), masse
    character(len=14) :: nume
    character(len=16) :: type, nomcmd
    character(len=19) :: resu
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, in, jvec, nbd
    integer(kind=8) :: nbdir, nbv
    integer(kind=8), pointer :: ddl(:) => null()
    real(kind=8), pointer :: vale(:) => null()
!-----------------------------------------------------------------------
    data tabcmp/'DX', 'DY', 'DZ', 'DRX', 'DRY', 'DRZ'/
!     ------------------------------------------------------------------
!
! --- RECUPERATION DES ARGUMENTS DE LA COMMANDE
!
    call jemarq()
    resu = ' '
    call getres(resu, type, nomcmd)
!
! --- MATRICE DE MASSE
!
    call getvid(' ', 'MATR_MASS', scal=masse, nbret=nbv)
    call mtdscr(masse)
    call jeveuo(masse//'           .&INT', 'E', lmat)
    call dismoi('NOM_NUME_DDL', masse, 'MATR_ASSE', repk=nume)
    call dismoi('NB_EQUA', masse, 'MATR_ASSE', repi=neq)
!
! --- QUELLE EST LA DIRECTION ?
!
    call getvr8(' ', 'DIRECTION', nbval=0, nbret=nbd)
    nbdir = -nbd
    call getvr8(' ', 'DIRECTION', nbval=nbdir, vect=depl, nbret=nbd)
!
!     --- ON NORMALISE LE VECTEUR ---
    xnorm = 0.d0
    do i = 1, nbdir
        xnorm = xnorm+depl(i)*depl(i)
    end do
    xnorm = sqrt(xnorm)
    if (xnorm .lt. 0.d0) then
        call utmess('F', 'ALGORITH9_81')
    end if
    do i = 1, nbdir
        depl(i) = depl(i)/xnorm
    end do
!
    call wkvect('&&SIMONO.VECTEUR', 'V V R', neq, jvec)
    AS_ALLOCATE(vi=ddl, size=neq*nbdir)
    call pteddl('NUME_DDL', nume, nbdir, tabcmp, neq, &
                tabl_equa=ddl)
    do i = 1, nbdir
        do in = 0, neq-1
            zr(jvec+in) = zr(jvec+in)-ddl(1+(i-1)*neq+in)*depl(i)
        end do
    end do
!
!     --- CREATION DU CHAMNO ---
!
    call vtcrem(resu, masse, 'G', 'R')
    call jeveuo(resu//'.VALE', 'E', vr=vale)
!
    call mrmult('ZERO', lmat, zr(jvec), vale, 1, &
                .true._1)
!
!      CALL WKVECT('&&SIMONO.DDL.BLOQUE','V V I',NEQ,IDDL)
!      CALL TYPDDL('BLOQ',NUME,NEQ,ZI(IDDL),NBACT,NBBLO,NBLAG,NBLIAI)
!      DO 40 IN = 0,NEQ-1
!         ZR(IDCHM+IN) = ( 1 - ZI(IDDL+IN) ) * ZR(IDCHM+IN)
! 40   CONTINUE
!
! --- MENAGE
    call jedetr('&&SIMONO.VECTEUR')
    AS_DEALLOCATE(vi=ddl)
!
    call jedema()
end subroutine
