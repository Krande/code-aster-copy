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

subroutine xmele3(mesh, model, ligrel, nfiss, chelem, &
                  param, option, list_func_acti)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/cescel.h"
#include "asterfort/cescre.h"
#include "asterfort/cesexi.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/infdbg.h"
#include "asterfort/isfonc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=8), intent(in) :: mesh
    character(len=8), intent(in) :: model
    character(len=*), intent(in) :: param
    character(len=*), intent(in) :: option
    integer, intent(in) :: nfiss
    character(len=19), intent(in) :: chelem
    character(len=19), intent(in) :: ligrel
    integer, intent(in) :: list_func_acti(*)
!
! ----------------------------------------------------------------------
!
! ROUTINE XFEM (METHODE XFEM - CREATION CHAM_ELEM)
!
! CREATION CHAM_ELEM RELATIFS AU CONTACT TYPE "MORTAR"
!
! ----------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  model            : name of model
! In  list_func_acti   : list of active functionnalities
! IN  NFISS  : NOMBRE TOTAL DE FISSURES
! IN  LIGREL : NOM DU LIGREL DES MAILLES TARDIVES
! IN  CHELEM : NOM DU CHAM_ELEM A CREER
! IN  PARAM  : NOM DE PARAMETRE
!
!
!
!
    integer :: ifm, niv
    character(len=8) :: k8bid
    integer :: ibid, iad, i, ima, ifis
    integer :: ino, itypma, nno
    integer :: ndim
    integer :: nbma, nmaenr
    character(len=8) :: nomfis, nomgd, typma, licmp3(3), licmp5(5)
    integer :: jcesl, jcesd, ncmp, icmp, jnbsp
    character(len=24) :: grp
    integer ::  jgrp, iret, ib1, ipt
    character(len=19) :: chelsi, chnbsp
    real(kind=8) :: valr
    character(len=8), pointer :: fiss(:) => null()
    integer, pointer :: typmail(:) => null()
    real(kind=8), pointer :: cesv(:) => null()
    aster_logical :: lxthm
!
    data licmp3/'X1', 'X2', 'X3'/
    data licmp5/'X1', 'X2', 'X3', 'X4', 'X5'/
!
! ----------------------------------------------------------------------
!
    call jemarq()
    call infdbg('XFEM', ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<XFEM  > CREATION DU CHAM_ELEM PINDCOI '
    end if
!
! --- INITIALISATIONS CHAMPS SIMPLES DE TRAVAIL
!
    chelsi = '&&XMELE3.CES'
!
! --- RECUPERATION DES INFOS SUR LE MAILLAGE ET LE MODELE
!
    call jeveuo(model//'.FISS', 'L', vk8=fiss)
    call dismoi('NB_MA_MAILLA', mesh, 'MAILLAGE', repi=nbma)
    call dismoi('DIM_GEOM', model, 'MODELE', repi=ndim)
    call jeveuo(mesh//'.TYPMAIL', 'L', vi=typmail)
!
    lxthm = isfonc(list_func_acti, 'THM')
    chnbsp = '&&XMELE3.NBSP'
    call wkvect(chnbsp, 'V V I', nbma, jnbsp)
!
    ASSERT(param .eq. 'PCOHES')
    nomgd = 'NEUT_R'
!
    if (lxthm) then
        ncmp = 5
    else
        ncmp = 3
    end if
    if (lxthm) then
! --- REMPLISSAGE DES SOUS POINTS POUR LA MULTU-FISSURATION
        do ifis = 1, nfiss
            nomfis = fiss(ifis)
            grp = nomfis(1:8)//'.MAILFISS.CONT'
            call jeexin(grp, iret)
            if (iret .ne. 0) then
                call jeveuo(grp, 'L', jgrp)
                call jelira(grp, 'LONMAX', nmaenr, k8bid)
                do i = 1, nmaenr
                    ima = zi(jgrp-1+i)
                    ASSERT(ima .le. nbma)
                    zi(jnbsp-1+ima) = 3
                end do
            end if
        end do
    end if
!
! --- TEST EXISTENCE DU CHAM_ELEM OU NON
!
    call exisd('CHAM_ELEM', chelem, iret)
    if (iret .eq. 0) then
        if (lxthm) then
            call cescre('V', chelsi, 'ELNO', mesh, nomgd, &
                        ncmp, licmp5, [-1], zi(jnbsp), [-ncmp])
        else
            call cescre('V', chelsi, 'ELNO', mesh, nomgd, &
                        ncmp, licmp3, [-1], [-1], [-ncmp])
        end if
!
!
! --- ACCES AU CHAM_ELEM_S
!
        call jeveuo(chelsi//'.CESD', 'L', jcesd)
        call jeveuo(chelsi//'.CESL', 'E', jcesl)
        call jeveuo(chelsi//'.CESV', 'E', vr=cesv)
!
! --- ENRICHISSEMENT DU CHAM_ELEM_S POUR LA MULTIFISSURATION
!
        do ifis = 1, nfiss
!
! --- ACCES FISSURE COURANTE
!
            nomfis = fiss(ifis)
            grp = nomfis(1:8)//'.MAILFISS.CONT'
            call jeexin(grp, iret)
            valr = 0.d0
!
! --- ON COPIE LES CHAMPS CORRESP. AUX ELEM. DE CONTACT
!
            if (iret .ne. 0) then
                call jeveuo(grp, 'L', jgrp)
                call jelira(grp, 'LONMAX', nmaenr, k8bid)
                do i = 1, nmaenr
                    ima = zi(jgrp-1+i)
                    itypma = typmail(ima)
                    call jenuno(jexnum('&CATA.TM.NOMTM', itypma), typma)
                    call dismoi('NBNO_TYPMAIL', typma, 'TYPE_MAILLE', repi=nno)
!
! --- RECOPIE EFFECTIVE DES CHAMPS
!
                    do ino = 1, nno
                        do icmp = 1, ncmp
                            if (lxthm) then
                                do ipt = 1, 3
                                    call cesexi('S', jcesd, jcesl, ima, ino, &
                                                ipt, icmp, iad)
                                    zl(jcesl-1+abs(iad)) = .true.
                                    cesv(abs(iad)) = valr
                                end do
                            else
                                call cesexi('S', jcesd, jcesl, ima, ino, &
                                            1, icmp, iad)
                                zl(jcesl-1+abs(iad)) = .true.
                                cesv(abs(iad)) = valr
                            end if
                        end do
                    end do
                end do
            end if
        end do
!
! --- CONVERSION CHAM_ELEM_S -> CHAM_ELEM
!
! on autorise un prolongement par zero
! sinon, il faudrait mettre NON a la place de OUI
        call cescel(chelsi, ligrel, option, param, 'NON', &
                    ib1, 'V', chelem, 'F', ibid)
!
! --- MENAGE
!
        call detrsd('CHAM_ELEM_S', chelsi)
    end if
!
    call jedetr(chnbsp)
!
    call jedema()
!
end subroutine
