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
subroutine sinoz1(model, sigmElga, sigmNoeu)
!
    implicit none
!
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
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
#include "asterfort/numoch.h"
#include "asterfort/preres.h"
#include "asterfort/resoud.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    character(len=8), intent(in) :: model
    character(len=24), intent(in) :: sigmElga, sigmNoeu
!
! --------------------------------------------------------------------------------------------------
!
!     BUT:
!           CALCUL DES CONTRAINTES AUX NOEUDS PAR LA METHODE ZZ1
!     ENTREES:
!        MODELE : NOM DU MODELE
!        SIGMA  : NOM DU CHAMP DE CONTRAINTES AUX POINTS DE GAUSS
!        SIGNO  : NOM DU CHAMP DE CONTRAINTES AUX NOEUDS
!
! --------------------------------------------------------------------------------------------------
!
    character(len=1) :: typres, k1bid
    character(len=14) :: numeDof
    character(len=8) :: licmp(6), mesh
    character(len=19), parameter :: listLoad = '&&SINOZ1.INFCHA'
    character(len=19), parameter :: vectElem = '&&VECELE'
    character(len=19) :: solveu, matpre, k19bid, criter, matrElem19
    character(len=24), parameter :: matrElem = '&&MASSEL'
    character(len=24) :: vecass, vect(6)
    integer(kind=8), parameter :: nbMatrElem = 1
    character(len=24), pointer :: listMatrElem(:) => null()
    integer(kind=8) :: nbLigr
    character(len=24), pointer :: listLigr(:) => null()
    real(kind=8) :: rcmp(6)
    integer(kind=8) :: ibid, jvect, nbcmp, repdim
    complex(kind=8) :: cbid
    integer(kind=8) :: iret
    integer(kind=8) :: i, ieq, ier, indeq, jprno
    integer(kind=8) :: jvecas, nbno
    real(kind=8), pointer :: sig(:) => null()
    real(kind=8), pointer :: sixx(:) => null()
    real(kind=8), pointer :: sixy(:) => null()
    real(kind=8), pointer :: sixz(:) => null()
    real(kind=8), pointer :: siyy(:) => null()
    real(kind=8), pointer :: siyz(:) => null()
    real(kind=8), pointer :: sizz(:) => null()
    integer(kind=8), pointer :: slvi(:) => null()
    integer(kind=8), pointer :: nueq(:) => null()
    cbid = dcmplx(0.d0, 0.d0)
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    call dismoi('DIM_GEOM', model, 'MODELE', repi=repdim)
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

!     CALCUL DE LA MATRICE DE MASSE ZZ1 (A 1 CMP)
    call memzme(model, matrElem(1:19))
!
    typres = 'R'
!
!     -- CREATION DU SOLVEUR :
    solveu = '&&OP0042.SOLVEUR'

! - Create list of LIGREL from matr_elem
    AS_ALLOCATE(vk24=listMatrElem, size=nbMatrElem)
    listMatrElem(1) = matrElem
    call numoch(listMatrElem, nbMatrElem, listLigr, nbLigr)
    AS_DEALLOCATE(vk24=listMatrElem)

! - Numbering
    numeDof = '&&NUME'
    call numero(numeDof, 'VV', &
                nbLigr, listLigr, &
                modeLocZ_="DDL_NOZ1")
    AS_DEALLOCATE(vk24=listLigr)

! - Assemblying
    matrElem19 = matrElem(1:19)
    call asmatr(1, matrElem19, ' ', numeDof, &
                listLoad, 'ZERO', 'V', 1, '&&MASSAS')

! - CALCUL DES SECONDS MEMBRES
    call me2zme(model, sigmElga(1:19), vectElem)

!  ASSEMBLAGE DES SECONDS MEMBRES
    vecass = '??????'
    call asasve(vectElem, numeDof, typres, vecass)
    call jeveuo(vecass, 'L', jvecas)
    do i = 1, nbcmp
        vect(i) = zk24(jvecas-1+i)
    end do

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
    rcmp = 0.d0
    licmp(1) = 'SIXX'
    licmp(2) = 'SIYY'
    licmp(3) = 'SIZZ'
    licmp(4) = 'SIXY'
    licmp(5) = 'SIXZ'
    licmp(6) = 'SIYZ'
    call dismoi('NOM_MAILLA', sigmElga(1:19), 'CHAM_ELEM', repk=mesh)
    call crcnct('G', sigmNoeu, mesh, 'SIEF_R', nbcmp, &
                licmp, rcmp)
    call jeveuo(sigmNoeu(1:19)//'.VALE', 'E', vr=sig)
    call jeveuo(vect(1) (1:19)//'.VALE', 'E', vr=sixx)
    call jeveuo(vect(2) (1:19)//'.VALE', 'E', vr=siyy)
    call jeveuo(vect(3) (1:19)//'.VALE', 'E', vr=sizz)
    call jeveuo(vect(4) (1:19)//'.VALE', 'E', vr=sixy)
    if (nbcmp .eq. 6) then
        call jeveuo(vect(5) (1:19)//'.VALE', 'E', vr=sixz)
        call jeveuo(vect(6) (1:19)//'.VALE', 'E', vr=siyz)
    end if
    call jeveuo(jexnum(numeDof(1:14)//'.NUME.PRNO', 1), 'L', jprno)
    call jeveuo(numeDof(1:14)//'.NUME.NUEQ', 'L', vi=nueq)
!
    call dismoi('NB_NO_MAILLA', mesh, 'MAILLAGE', repi=nbno)
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
    call jedetr(numeDof//'.&LMODCHAR')
    call detrsd('NUME_DDL', numeDof)
!
    do i = 1, nbcmp
        call detrsd('CHAMP_GD', zk24(jvecas-1+i))
    end do
    call jedetr(vecass)
!
!
    call jedema()
end subroutine
