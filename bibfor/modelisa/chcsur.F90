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

subroutine chcsur(chcinez, chamnosz, type, modelz, gran_name)
    implicit none
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisdg.h"
#include "asterfort/gettco.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
!
    character(len=1) :: type
    character(len=8) :: gran_name
    character(len=*) :: chcinez, chamnosz, modelz
! OBJET : CREATION D"UNE CHARGE CINEMATIQUE.
!        1) LE .REFE DE LA CHARGE DOIT DEJA EXISTER
!        2) MISE A JOUR DE : .AFCI ET .AFCV
!-----------------------------------------------------------------------
! OUT  CHCINE  K*19    : NOM DE LA CHARGE CINEMATIQUE
! IN   CNS     K*19    : NOM D'UN CHAM_NO_S CONTENANT LES DEGRES IMPOSES
! IN   TYPE    K*1     : 'R','C' OU 'F' TYPE DE LA CHARGE
! IN   MO      K*      : NOM DU MODELE
! IN   NOMGD   K*      : NOM DE LA GRANDEUR
!-----------------------------------------------------------------------
!
    integer(kind=8) :: nb_affe_cine, i_affe_cine, i_cmp_chmx, nb_cmp_chmx
    integer(kind=8) ::  i_node, nb_node, nbec, i_cmp
    integer(kind=8) :: jcnsd, jcnsv, jcnsl, iaprnm, jcnsc, i_cmp_mx
    integer(kind=8) :: nb_cmp_mx
    character(len=8) :: model, cmp_name, nommai
    character(len=16) :: sdtyp
    character(len=19) :: chcine, chamnos
    character(len=24) :: cafci, cafcv
    integer(kind=8), pointer :: corres(:) => null()
    integer(kind=8), pointer :: afci(:) => null()
    real(kind=8), pointer :: afcv_r(:) => null()
    complex(kind=8), pointer :: afcv_c(:) => null()
    character(len=8), pointer :: afcv_f(:) => null()
    character(len=8), pointer :: cata_gd_nomcmp(:) => null()
!
! --- DEBUT -----------------------------------------------------------
!
    call jemarq()
!
    model = modelz
    call dismoi('NB_EC', gran_name, 'GRANDEUR', repi=nbec)
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=nommai)
    call gettco(nommai, sdtyp)
    call jeveuo(model//'.MODELE    .PRNM', 'L', iaprnm)
!
    chcine = chcinez
    chamnos = chamnosz
    cafci = chcine(1:19)//'.AFCI'
    cafcv = chcine(1:19)//'.AFCV'
!
! - Access to CHAM_NO_S
!
    call jeveuo(chamnos//'.CNSD', 'L', jcnsd)
    call jeveuo(chamnos//'.CNSC', 'L', jcnsc)
    call jeveuo(chamnos//'.CNSV', 'L', jcnsv)
    call jeveuo(chamnos//'.CNSL', 'L', jcnsl)
!
    nb_node = zi(jcnsd)
    nb_cmp_chmx = zi(jcnsd+1)
!
! - Set bijection between local component index in CHAM_NO_S and global component index in GRANDEUR
!
    call jeveuo(jexnom('&CATA.GD.NOMCMP', gran_name), 'L', vk8=cata_gd_nomcmp)
    call jelira(jexnom('&CATA.GD.NOMCMP', gran_name), 'LONMAX', nb_cmp_mx)
    AS_ALLOCATE(vi=corres, size=nb_cmp_mx)
    do i_cmp_mx = 1, nb_cmp_mx
        cmp_name = cata_gd_nomcmp(i_cmp_mx)
        i_cmp_chmx = indik8(zk8(jcnsc), cmp_name, 1, nb_cmp_chmx)
        corres(i_cmp_mx) = i_cmp_chmx
    end do
!
! - Number of affected values
!
    nb_affe_cine = 0
    do i_cmp_chmx = 1, nb_cmp_chmx
        do i_node = 1, nb_node
            if (zl(jcnsl+(i_node-1)*nb_cmp_chmx+i_cmp_chmx-1)) then
                nb_affe_cine = nb_affe_cine+1
            end if
        end do
    end do
!
! - Create datastructure
!
    call wkvect(cafci, 'G V I', (3*nb_affe_cine+1), vi=afci)
    if (type .eq. 'R') then
        call wkvect(cafcv, 'G V R', max(nb_affe_cine, 1), vr=afcv_r)
    else if (type .eq. 'C') then
        call wkvect(cafcv, 'G V C', max(nb_affe_cine, 1), vc=afcv_c)
    else if (type .eq. 'F') then
        call wkvect(cafcv, 'G V K8', max(nb_affe_cine, 1), vk8=afcv_f)
    end if
!
! - Set datastructure
!
    i_affe_cine = 0
    do i_node = 1, nb_node
        i_cmp = 0
        do i_cmp_mx = 1, nb_cmp_mx
            if (exisdg(zi(iaprnm-1+nbec*(i_node-1)+1), i_cmp_mx)) then
                i_cmp = i_cmp+1
                i_cmp_chmx = corres(i_cmp_mx)
                if (i_cmp_chmx .ne. 0) then
                    if (zl(jcnsl+(i_node-1)*nb_cmp_chmx+i_cmp_chmx-1)) then
                        i_affe_cine = i_affe_cine+1
                        afci(3*(i_affe_cine-1)+2) = i_node
                        afci(3*(i_affe_cine-1)+3) = i_cmp
                        if (type .eq. 'R') then
                            afcv_r(i_affe_cine) = zr(jcnsv+(i_node-1)*nb_cmp_chmx+i_cmp_chmx-1)
                        else if (type .eq. 'C') then
                            afcv_c(i_affe_cine) = zc(jcnsv+(i_node-1)*nb_cmp_chmx+i_cmp_chmx-1)
                        else if (type .eq. 'F') then
                            afcv_f(i_affe_cine) = zk8(jcnsv+(i_node-1)*nb_cmp_chmx+i_cmp_chmx-1)
                        else
                            ASSERT(.false.)
                        end if
                    end if
                end if
            end if
        end do
    end do
!
    if (i_affe_cine .eq. 0 .and. sdtyp .ne. 'MAILLAGE_P') then
        call utmess('F', 'CALCULEL_9')
    end if
    afci(1) = i_affe_cine
!
    AS_DEALLOCATE(vi=corres)
    call jedema()
end subroutine
