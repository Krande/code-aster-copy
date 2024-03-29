! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

subroutine sinoz1(modele, sigma, signo)
    implicit none
!
!     ARGUMENTS:
!     ----------
! ......................................................................
!     BUT:
!           CALCUL DES CONTRAINTES AUX NOEUDS PAR LA METHODE ZZ1
!     ENTREES:
!        MODELE : NOM DU MODELE
!        SIGMA  : NOM DU CHAMP DE CONTRAINTES AUX POINTS DE GAUSS
!        SIGNO  : NOM DU CHAMP DE CONTRAINTES AUX NOEUDS
!
!----------------------------------------------------------------------
!
! ----------------------- DECLARATIONS --------------------------------
!
#include "jeveux.h"
#include "asterfort/asasve.h"
#include "asterfort/asmatr.h"
#include "asterfort/assert.h"
#include "asterfort/crcnct.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/me2zme.h"
#include "asterfort/memzme.h"
#include "asterfort/numero.h"
#include "asterfort/preres.h"
#include "asterfort/resoud.h"
#include "asterfort/utmess.h"
!
    character(len=1) :: typres, k1bid
    character(len=8) :: modele
    character(len=14) :: nupgm
    character(len=8) :: licmp(6), ma
    character(len=19) :: infcha
    character(len=19) :: solveu, vecele, matpre, k19bid, criter, masselK19
    character(len=24) :: signo, sigma, massel
    character(len=24) :: nume, vecass, vect(6)
    real(kind=8) :: rcmp(6)
    integer :: ibid, jvect, nbcmp, repdim
    complex(kind=8) :: cbid
    integer :: iret
!
!
!-----------------------------------------------------------------------
    integer :: i, ieq, ier, indeq, jprno
    integer :: jvecas, nbno
    real(kind=8), pointer :: sig(:) => null()
    real(kind=8), pointer :: sixx(:) => null()
    real(kind=8), pointer :: sixy(:) => null()
    real(kind=8), pointer :: sixz(:) => null()
    real(kind=8), pointer :: siyy(:) => null()
    real(kind=8), pointer :: siyz(:) => null()
    real(kind=8), pointer :: sizz(:) => null()
    integer, pointer :: slvi(:) => null()
    integer, pointer :: nueq(:) => null()
    cbid = dcmplx(0.d0, 0.d0)
!-----------------------------------------------------------------------
    call jemarq()
!
    call dismoi('DIM_GEOM', modele, 'MODELE', repi=repdim)
!
    if ((repdim .ne. 2) .or. (repdim .ne. 3)) then
        if (repdim .eq. 1) then
            call utmess('F', 'CALCULEL_75')
        end if
        if (repdim .gt. 3) then
            call utmess('F', 'CALCULEL_76')
        end if
    end if
!
    if (repdim .eq. 2) then
        nbcmp = 4
    else if (repdim .eq. 3) then
        nbcmp = 6
    else
        ASSERT(.false.)
    end if
!
    massel = '&&MASSEL'
!
!     CALCUL DE LA MATRICE DE MASSE ZZ1 (A 1 CMP)
    call memzme(modele, massel(1:19))
!
    typres = 'R'
!
!
!     --  APPEL A NUMER2 POUR CONSTRUIRE UN NUME_DDL
!         SUR LA GRANDEUR SIZZ_R (1 CMP)
    nupgm = '&&NUME'
    infcha = '&&SINOZ1.INFCHA'
!
!     -- CREATION DU SOLVEUR :
    solveu = '&&OP0042.SOLVEUR'

    call numero(nupgm, 'VV', &
                modelocz='DDL_NOZ1', modelz=modele, &
                nb_matr_elem=1, list_matr_elem=massel)

    masselK19 = massel(1:19)
!
    call asmatr(1, masselK19, ' ', nupgm, &
                infcha, 'ZERO', 'V', 1, '&&MASSAS')
!
!     CALCUL DES SECONDS MEMBRES
    vecele = '&&VECELE'
    call me2zme(modele, sigma(1:19), vecele)
!
!
!     ASSEMBLAGE DES SECONDS MEMBRES
    nume = '&&NUME'
    vecass = '??????'
    call asasve(vecele, nume, typres, vecass)
    call jeveuo(vecass, 'L', jvecas)
    do i = 1, nbcmp
        vect(i) = zk24(jvecas-1+i)
    end do
!
!
!     RESOLUTIONS SANS DIRICHLET
!     -- ON FORCE STOP_SINGULIER='NON' MAIS POURQUOI ??
    call jeveuo(solveu//'.SLVI', 'E', vi=slvi)
    slvi(3) = 1
    matpre = '&&SINOZ1.MATPRE'
    call preres(solveu, 'V', ier, matpre, '&&MASSAS', &
                ibid, -9999)
!
    k19bid = ' '
    k1bid = ' '
    criter = ' '
    do i = 1, nbcmp
        call jeveuo(vect(i) (1:19)//'.VALE', 'E', jvect)
        call resoud('&&MASSAS', matpre, solveu, k19bid, 1, &
                    k19bid, k19bid, k1bid, zr(jvect), [cbid], &
                    criter, .true._1, 0, iret)
    end do
!
!   CREATION DU CHAM_NO_SIEF_R A PARTIR DES 4 CHAM_NO_SIZZ_R (A 1 CMP)
!
    do i = 1, nbcmp
        rcmp(i) = 0.d0
    end do
    licmp(1) = 'SIXX'
    licmp(2) = 'SIYY'
    licmp(3) = 'SIZZ'
    licmp(4) = 'SIXY'
    licmp(5) = 'SIXZ'
    licmp(6) = 'SIYZ'
    call dismoi('NOM_MAILLA', sigma(1:19), 'CHAM_ELEM', repk=ma)
    call crcnct('G', signo, ma, 'SIEF_R', nbcmp, &
                licmp, rcmp)
    call jeveuo(signo(1:19)//'.VALE', 'E', vr=sig)
    call jeveuo(vect(1) (1:19)//'.VALE', 'E', vr=sixx)
    call jeveuo(vect(2) (1:19)//'.VALE', 'E', vr=siyy)
    call jeveuo(vect(3) (1:19)//'.VALE', 'E', vr=sizz)
    call jeveuo(vect(4) (1:19)//'.VALE', 'E', vr=sixy)
    if (nbcmp .eq. 6) then
        call jeveuo(vect(5) (1:19)//'.VALE', 'E', vr=sixz)
        call jeveuo(vect(6) (1:19)//'.VALE', 'E', vr=siyz)
    end if
    call jeveuo(jexnum(nume(1:14)//'.NUME.PRNO', 1), 'L', jprno)
    call jeveuo(nume(1:14)//'.NUME.NUEQ', 'L', vi=nueq)
!
    call dismoi('NB_NO_MAILLA', ma, 'MAILLAGE', repi=nbno)
    do i = 1, nbno
        indeq = zi(jprno-1+3*(i-1)+1)
        ieq = nueq(indeq)
        sig(nbcmp*(i-1)+1) = sixx(ieq)
        sig(nbcmp*(i-1)+2) = siyy(ieq)
        sig(nbcmp*(i-1)+3) = sizz(ieq)
        sig(nbcmp*(i-1)+4) = sixy(ieq)
        if (nbcmp .eq. 6) then
            sig(nbcmp*(i-1)+5) = sixz(ieq)
            sig(nbcmp*(i-1)+6) = siyz(ieq)
        end if
    end do
!
    call detrsd('MATR_ASSE', '&&MASSAS')
    call jedetr(nupgm//'.&LMODCHAR')
    call detrsd('NUME_DDL', nupgm)
!
    do i = 1, nbcmp
        call detrsd('CHAMP_GD', zk24(jvecas-1+i))
    end do
    call jedetr(vecass)
!
!
    call jedema()
end subroutine
