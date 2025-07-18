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

subroutine cmpcha(fieldz, cmp_name, cata_to_field, field_to_cata, nb_cmpz, &
                  nb_cmp_mxz)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisdg.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
!
! person_in_charge: nicolas.sellenet at edf.fr
!
    character(len=*), intent(in) :: fieldz
    character(len=8), pointer :: cmp_name(:)
    integer(kind=8), pointer :: cata_to_field(:)
    integer(kind=8), pointer :: field_to_cata(:)
    integer(kind=8), optional, intent(out) :: nb_cmpz
    integer(kind=8), optional, intent(out) :: nb_cmp_mxz
!
! --------------------------------------------------------------------------------------------------
!
! Create objects for global components (catalog) <=> local components (field)
!
! --------------------------------------------------------------------------------------------------
!
! In  field         : name of field
! Out nb_cmp_mx     : number of components in GRANDEUR (catalog)
! Out nb_cmp        : number of components in GRANDEUR (field)
! Out cmp_name      : pointer to name of components in field
! Out cata_to_field : pointer to converter from global components (catalog) to local (field)
! Out field_to_cata : pointer to converter from local components (field) to global (catalog)
!
!
! --------------------------------------------------------------------------------------------------
!
!  EXEMPLE :
!  SI cmp_name EST UN CHAMP DE DEPL_R NE CONTENANT QUE 'DY' ET 'DRZ' :
!    NCMP=2
!    NOMCP=('DY','DRZ')
!    CORR1=(0,1,0,0,2,0,...)
!    CORR2=(2,5)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nb_grel, jceld, nb_ec, jcmpgd, nb_cmp, nb_cmp_mx
    integer(kind=8) :: igr, imolo, jmolo, idx_gd, nb_pt, i_pt, k, iadg, i_cmp
    integer(kind=8) :: jdesc, long, jprno, nb_node, i_node, ncmpp
    integer(kind=8) :: ngrmx, nbedit, igd, ient, debgd, dg(200), ior, kpt, kcmp
    aster_logical :: diff, l_pmesh
    character(len=8) :: gran_name, mesh
    character(len=16) :: typsd
    character(len=19) :: field, nume_equa
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    field = fieldz
!
! - Informations about field
!
    call dismoi('TYPE_CHAMP', fieldz, 'CHAMP', repk=typsd)
    if (typsd .eq. 'NOEU') then
        call dismoi('NOM_GD', field, 'CHAM_NO', repk=gran_name)
    else if (typsd(1:2) .eq. 'EL') then
        call dismoi('NOM_GD', field, 'CHAM_ELEM', repk=gran_name)
    else if (typsd .eq. 'CART') then
        call dismoi('NOM_GD', field, 'CARTE', repk=gran_name)
    else if (typsd .eq. 'GEOM') then
        gran_name = "GEOM_R"
    else
        ASSERT(.false.)
    end if
    l_pmesh = .false.
!
    call dismoi('NB_EC', gran_name, 'GRANDEUR', repi=nb_ec)
    call dismoi('NB_CMP_MAX', gran_name, 'GRANDEUR', repi=nb_cmp_mx)
    call dismoi('NUM_GD', gran_name, 'GRANDEUR', repi=idx_gd)
    call jeveuo(jexnum('&CATA.GD.NOMCMP', idx_gd), 'L', jcmpgd)
!
!
!
!     -- 1. POUR ECONOMISER LES APPELS A EXISDG, ON VA CALCULER
!           UN DESCRIPTEUR_GRANDEUR (DG) "ENVELOPPE" DE TOUS LES
!           POINTS DU CHAMP.
!     ----------------------------------------------------------------
    ASSERT(nb_ec .le. 200)
    dg(1:200) = 0
!
!
!     -- 1.1 CAS DES CHAM_NO
!     ----------------------------------------------------------------
    if (typsd .eq. 'NOEU') then
        call dismoi('NOM_MAILLA', field, 'CHAM_NO', repk=mesh)
        l_pmesh = isParallelMesh(mesh)
        call dismoi('NB_NO_MAILLA', mesh, 'MAILLAGE', repi=nb_node)
!
!        -- 1.1.2 CAS DES CHAM_NO A NUME_EQUA:
        call dismoi('NUME_EQUA', field, 'CHAM_NO', repk=nume_equa)
        call jeveuo(jexnum(nume_equa//'.PRNO', 1), 'L', jprno)
        do i_node = 1, nb_node
            ncmpp = zi(jprno-1+(i_node-1)*(nb_ec+2)+2)
            if (ncmpp .ne. 0) then
                iadg = jprno-1+(i_node-1)*(nb_ec+2)+3
                do k = 1, nb_ec
                    dg(k) = ior(dg(k), zi(iadg-1+k))
                end do
            end if
        end do
!
!     -- 1.1 CAS DES CHAM_GEOM
!     ----------------------------------------------------------------
    elseif (typsd .eq. 'GEOM') then
        call jeveuo(field//'.DESC', 'L', jdesc)
        ASSERT(zi(jdesc-1+2) == -3)
        call jelira(field//'.DESC', 'LONMAX', long)
        ASSERT(long .eq. (2+nb_ec))
        iadg = jdesc-1+3
        do k = 1, nb_ec
            dg(k) = zi(iadg-1+k)
        end do
        !
!
!
!     -- 1.2 CAS DES CHAM_ELEM
!     ----------------------------------------------------------------
    else if (typsd(1:2) .eq. 'EL') then
        call dismoi('NOM_MAILLA', field, 'CHAMP', repk=mesh)
        l_pmesh = isParallelMesh(mesh)
        call jeveuo(field//'.CELD', 'L', jceld)
        nb_grel = zi(jceld-1+2)
!
        do igr = 1, nb_grel
            imolo = zi(jceld-1+zi(jceld-1+4+igr)+2)
            if (imolo .ne. 0) then
                call jeveuo(jexnum('&CATA.TE.MODELOC', imolo), 'L', jmolo)
                ASSERT(zi(jmolo-1+1) .le. 3)
                ASSERT(zi(jmolo-1+2) .eq. idx_gd)
                diff = (zi(jmolo-1+4) .gt. 10000)
                nb_pt = mod(zi(jmolo-1+4), 10000)
                do i_pt = 1, nb_pt
                    kpt = 1
                    if (diff) kpt = i_pt
                    iadg = jmolo-1+4+(kpt-1)*nb_ec+1
                    do k = 1, nb_ec
                        dg(k) = ior(dg(k), zi(iadg-1+k))
                    end do
                end do
            end if
        end do
!
!
!     -- 1.3 CAS DES CARTES
!     ----------------------------------------------------------------
    else if (typsd .eq. 'CART') then
        call dismoi('NOM_MAILLA', field, 'CARTE', repk=mesh)
        l_pmesh = isParallelMesh(mesh)
        call jeveuo(field//'.DESC', 'L', jdesc)
        ngrmx = zi(jdesc-1+2)
        nbedit = zi(jdesc-1+3)
        do igd = 1, nbedit
            ient = zi(jdesc-1+3+2*igd)
            if (ient .ne. 0) then
                debgd = 3+2*ngrmx+(igd-1)*nb_ec+1
                do k = 1, nb_ec
                    dg(k) = ior(dg(k), zi(jdesc-1+debgd-1+k))
                end do
            end if
        end do
!
    else
        ASSERT(.false.)
    end if
!
! - Create cata_to_field object
!
    AS_ALLOCATE(vi=cata_to_field, size=nb_cmp_mx)
    nb_cmp = 0
    do i_cmp = 1, nb_cmp_mx
        if (exisdg(dg(1), i_cmp)) then
            nb_cmp = nb_cmp+1
            cata_to_field(i_cmp) = nb_cmp
        end if
    end do
    if (.not. l_pmesh) then
        ASSERT(nb_cmp .ne. 0)
    end if
!
! - Create field_to_cata and cmp_name objects
!
    if (nb_cmp .ne. 0) then
        AS_ALLOCATE(vi=field_to_cata, size=nb_cmp)
        AS_ALLOCATE(vk8=cmp_name, size=nb_cmp)
    end if
    kcmp = 0
    do i_cmp = 1, nb_cmp_mx
        if (cata_to_field(i_cmp) .gt. 0) then
            kcmp = kcmp+1
            field_to_cata(kcmp) = i_cmp
            cmp_name(kcmp) = zk8(jcmpgd-1+i_cmp)
        end if
    end do
    ASSERT(kcmp .eq. nb_cmp)
!
    if (present(nb_cmpz)) then
        nb_cmpz = nb_cmp
    end if
    if (present(nb_cmp_mxz)) then
        nb_cmp_mxz = nb_cmp_mx
    end if
!
    call jedema()
end subroutine
