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
subroutine mearcc(option, mo, chin, chout)
    implicit none
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/celces.h"
#include "asterfort/cescel.h"
#include "asterfort/cescre.h"
#include "asterfort/cesexi.h"
#include "asterfort/cesred.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/srlima.h"
#include "asterfort/utmess.h"
#include "asterfort/int_to_char8.h"
!
    character(len=8) :: mo
    character(len=16) :: option
    character(len=24) :: chin, chout
!
!     BUT: REDUIRE LE CHAMP DE CONTRAINTES (ELNO) DES ELEMENTS 3D
!          A LEURS FACES POUR L'OPTION SIRO_ELEM
!
!     IN  MO     : NOM DU MODELE
!     IN  OPTION : NOM DE L'OPTION
!     IN  CHIN   : CHAMP DE CONTRAINTE ELNO DES ELEMENTS 3D
!     OUT CHOUT  : CHAMP DE CONTRAINTES ELNO REDUIT AUX MAILLES DE PEAU
!
!
    integer(kind=8) :: nbcmp, nbnomx
    parameter(nbcmp=6, nbnomx=9)
!
    integer(kind=8) :: nbma, ibid, jma2d, jma3d, ndim
    integer(kind=8) :: jcesd3, jcesk3, jcesl3, jcesd2, ima
    integer(kind=8) :: jcesk2, jcesl2, jcesc2, jlcnx, jcnx, ipt, icp, ino2, ino3
    integer(kind=8) :: jco3, jco2, npt3, npt2, ipt2, ipt3, k, npt
    integer(kind=8) :: iad3, iad2, nucmp, numasu, numavo, nbcmp2
!
    character(len=8) :: ma, comp(nbcmp), k8b, nomasu, nomavo, valk(2)
    character(len=19) :: chous, chins
    character(len=24) :: mail2d, mail3d, mailto, ligrmo
    integer(kind=8), pointer :: pt3d(:) => null()
    real(kind=8), pointer :: cesv2(:) => null()
    real(kind=8), pointer :: cesv3(:) => null()
    character(len=8), pointer :: cesc3(:) => null()
!
    data comp/'SIXX', 'SIYY', 'SIZZ', 'SIXY', 'SIXZ', 'SIYZ'/
!--------------------------------------------------------------------------
!
    call jemarq()
!
    mail2d = '&&MEARCC.MAILLE_FACE'
    mail3d = '&&MEARCC.MAILLE_3D_SUPP'
    mailto = '&&MEARCC.MAILLE_2D_3D'
!
    call dismoi('NOM_MAILLA', mo, 'MODELE', repk=ma)
    call dismoi('DIM_GEOM', ma, 'MAILLAGE', repi=ndim)
    nbcmp2 = nbcmp
    if (ndim .eq. 2) nbcmp2 = 4
!    ASSERT(ndim.eq.3)
!
!     RECUPERATION DES MAILLES DE FACES ET DES MAILLES 3D SUPPORT
    call srlima(mo, mail2d, mail3d, mailto, nbma)
    call jeveuo(mail2d, 'L', jma2d)
    call jeveuo(mail3d, 'L', jma3d)
!
!     =============================================================
! --- POUR CHAQUE NOEUD DE LA MAILLE 3D SUPPORT, ON RECOPIE LES
!     CONTRAINTES AUX NOEUDS DE LA MAILLE DE PEAU
!     =============================================================
!
!     TRANSFORMATION DU CHAMP 3D (IN) EN CHAMP SIMPLE
    chins = '&&MEARCC.CHIN_S'
    call celces(chin, 'V', chins)
    call jeveuo(chins//'.CESV', 'L', vr=cesv3)
    call jeveuo(chins//'.CESD', 'L', jcesd3)
    call jeveuo(chins//'.CESK', 'L', jcesk3)
    call jeveuo(chins//'.CESL', 'L', jcesl3)
    call jeveuo(chins//'.CESC', 'L', vk8=cesc3)
!
!     CREATION DU CHAMP 2D (OUT) SIMPLE
    chous = '&&MEARCC.CHOUT_S'
    call cescre('V', chous, 'ELNO', ma, 'SIEF_R', &
                nbcmp2, comp, [-1], [-1], [-nbcmp2])
!
    call jeveuo(chous//'.CESV', 'E', vr=cesv2)
    call jeveuo(chous//'.CESD', 'E', jcesd2)
    call jeveuo(chous//'.CESK', 'E', jcesk2)
    call jeveuo(chous//'.CESL', 'E', jcesl2)
    call jeveuo(chous//'.CESC', 'E', jcesc2)
!
!
!   -- correspondance maille 2d / maille 3d: zi(jpt3d)
!      pour chaque point de la maille 2d, on cherche le
!      point de la maille 3d correspondant
!   ---------------------------------------------------------
    AS_ALLOCATE(vi=pt3d, size=nbma*nbnomx)
    call jeveuo(ma//'.CONNEX', 'L', jcnx)
    call jeveuo(jexatr(ma//'.CONNEX', 'LONCUM'), 'L', jlcnx)
    do ima = 1, nbma
        if (zi(jma3d+ima-1) .eq. 0) goto 2
!
        jco3 = jcnx+zi(jlcnx-1+zi(jma3d+ima-1))-1
        jco2 = jcnx+zi(jlcnx-1+zi(jma2d+ima-1))-1
        npt3 = zi(jcesd3-1+5+4*(zi(jma3d+ima-1)-1)+1)
        npt2 = zi(jcesd2-1+5+4*(zi(jma2d+ima-1)-1)+1)
        k = 0
        do ipt2 = 1, npt2
            ino2 = zi(jco2+ipt2-1)
            do ipt3 = 1, npt3
                ino3 = zi(jco3+ipt3-1)
                if (ino3 .eq. ino2) then
                    k = k+1
                    pt3d(1+nbnomx*(ima-1)+k-1) = ipt3
                    goto 110
                end if
            end do
110         continue
        end do
2       continue
    end do
!
!
!   -- remplissage du champ simple 3d
!   ----------------------------------
    do ima = 1, nbma
        if (zi(jma3d+ima-1) .eq. 0) goto 1
!
        npt = zi(jcesd2-1+5+4*(zi(jma2d+ima-1)-1)+1)
        k8b = int_to_char8(zi(jma2d+ima-1))
        k8b = int_to_char8(zi(jma3d+ima-1))
        do ipt = 1, npt
            do icp = 1, nbcmp2
                nucmp = indik8(cesc3, comp(icp), 1, zi(jcesd3+1) &
                               )
!
                call cesexi('C', jcesd3, jcesl3, zi(jma3d+ima-1), pt3d(1+nbnomx*(ima-1)+ipt-1), &
                            1, nucmp, iad3)
                if (iad3 .eq. 0) then
                    numasu = zi(jma2d+ima-1)
                    numavo = zi(jma3d+ima-1)
                    nomasu = int_to_char8(numasu)
                    nomavo = int_to_char8(numavo)
                    valk(1) = nomavo
                    valk(2) = nomasu
                    call utmess('F', 'CALCULEL5_52', nk=2, valk=valk)
                end if
                call cesexi('S', jcesd2, jcesl2, zi(jma2d+ima-1), ipt, &
                            1, nucmp, iad2)
                cesv2(1-iad2-1) = cesv3(iad3)
                zl(jcesl2-iad2-1) = .true.
            end do
        end do
1       continue
    end do
!
    call cesred(chous, nbma, zi(jma2d), 0, [k8b], &
                'V', chous)
!
    call dismoi('NOM_LIGREL', mo, 'MODELE', repk=ligrmo)
!
    call cescel(chous, ligrmo, option, 'PSIG3D', 'OUI', &
                ibid, 'V', chout, 'F', ibid)
!
    AS_DEALLOCATE(vi=pt3d)
    call jedetr(mail2d)
    call jedetr(mail3d)
    call jedetr(mailto)
    call detrsd('CHAM_ELEM_S', chous)
    call detrsd('CHAM_ELEM_S', chins)
!
    call jedema()
!
end subroutine
