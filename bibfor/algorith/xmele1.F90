! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

subroutine xmele1(mesh  , model, ds_contact, ligrel, nfiss,&
                  chelem, param, option,   list_func_acti)
!
use NonLin_Datastructure_type
!
implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/celces.h"
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
#include "asterfort/jeveuo.h"
#include "asterfort/mminfi.h"
#include "asterfort/mminfl.h"
#include "asterfort/mminfr.h"
#include "asterfort/wkvect.h"
#include "asterfort/xmelin.h"
#include "asterfort/xxconi.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=8), intent(in) :: mesh
    character(len=8), intent(in) :: model
    type(NL_DS_Contact), intent(in) :: ds_contact
    integer, intent(in) :: nfiss
    character(len=19), intent(in) :: chelem
    character(len=19), intent(in) :: ligrel
    character(len=*), intent(in) :: param
    character(len=*), intent(in) :: option
    integer, intent(in) :: list_func_acti(*)
!
! ----------------------------------------------------------------------
!
! ROUTINE XFEM (METHODE XFEM - CREATION CHAM_ELEM)
!
! CREATION CHAM_ELEM RELATIFS AU CONTACT
!
! ----------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  model            : name of model
! In  ds_contact       : datastructure for contact management
! In  list_func_acti   : list of active functionnalities
! IN  NFISS  : NOMBRE TOTAL DE FISSURES
! IN  LIGREL : NOM DU LIGREL DES MAILLES TARDIVES
! IN  CHELEM : NOM DU CHAM_ELEM A CREER
! IN  PARAM  : NOM DE PARAMETRE
!
! ----------------------------------------------------------------------
!
    integer :: ifm, niv
    integer :: ibid, iad, ispt, i, ima, ifis, izone
    integer :: ifima, jcesd1, jcesl1
    integer :: ndim, nface, nfisc, nnint, npg, typint, nfisc2
    integer :: nbma, nmaenr, jcesd2, jcesl2
    character(len=8) :: nomfis, nomgd, elc, nomfi2, licmp3(3), licmp5(5)
    integer :: jcesl, jcesv, jcesd, ncmp, icmp
    character(len=24) :: grp
    integer :: jgrp, iret, jnbsp, ifiss
    aster_logical :: vall, isint, lxthm
    character(len=19) :: chelsi, cmafis, faclon, chnbsp
    real(kind=8) :: valr
    character(len=8), pointer :: fiss(:) => null()
    character(len=8), pointer :: cesv1(:) => null()
    integer, pointer :: cesv2(:) => null()
    integer, pointer :: xfem_cont(:) => null()
!
    data licmp3    / 'X1', 'X2', 'X3'/
    data licmp5    / 'X1', 'X2', 'X3', 'X4', 'X5'/
!
! ----------------------------------------------------------------------
!
    call jemarq()
    call infdbg('XFEM', ifm, niv)
    if (niv .ge. 2) then
        write (ifm,*) '<XFEM  > CREATION DU CHAM_ELEM PINDCOI '
    endif
!
! --- INITIALISATIONS CHAMPS SIMPLES DE TRAVAIL
!
    chelsi = '&&XMELE1.CES'
    faclon = '&&XMELE1.FACLO'
    chnbsp = '&&XMELE1.NBSP'
    cmafis = '&&XMELE1.XMAFIS'
!
! --- RECUPERATION DES INFOS SUR LE MAILLAGE ET LE MODELE
!
    call jeveuo(model//'.FISS', 'L', vk8=fiss)
    call dismoi('NB_MA_MAILLA', mesh, 'MAILLAGE', repi=nbma)
    call dismoi('DIM_GEOM', model, 'MODELE', repi=ndim)
    call jeveuo(model//'.XFEM_CONT', 'L', vi=xfem_cont)
!
! --- ELEMENT DE REFERENCE ASSOCIE A UNE FACETTE DE CONTACT
!
    lxthm=isfonc(list_func_acti,'THM')
!
    if (ndim .eq. 3) then
        if (lxthm) then
            elc='TR6'
        else
            elc='TR3'
        endif
    else if (ndim.eq.2) then
        if (xfem_cont(1) .le. 2 .and. .not.lxthm) then
            elc='SE2'
        else
            elc='SE3'
        endif
    endif
!
! --- INFOS SUR LE TYPE DE LA VALEUR A INITIALISER
!
    isint = (param.eq.'PINDCOI'.or. param.eq.'PMEMCON'.or. param.eq.'PGLISS')
    if (isint) then
        nomgd = 'NEUT_I'
    else
        nomgd = 'NEUT_R'
    endif
!
    if (param .eq. 'PCOHES') then
       if (lxthm) then
          ncmp = 5
       else
          ncmp = 3
       endif
    else
        ncmp = 1
    endif
!
! --- CREATION DES OBJETS DE TRAVAIL
!
    call celces(model//'.TOPOFAC.LO', 'V', faclon)
    call jeveuo(faclon(1:19)//'.CESD', 'L', jcesd2)
    call jeveuo(faclon(1:19)//'.CESL', 'L', jcesl2)
    call jeveuo(faclon(1:19)//'.CESV', 'L', vi=cesv2)
    call celces(model//'.XMAFIS', 'V', cmafis)
    call jeveuo(cmafis(1:19)//'.CESD', 'L', jcesd1)
    call jeveuo(cmafis(1:19)//'.CESL', 'L', jcesl1)
    call jeveuo(cmafis(1:19)//'.CESV', 'L', vk8=cesv1)
    call wkvect(chnbsp, 'V V I', nbma, jnbsp)
!
! --- TEST EXISTENCE DU CHAM_ELEM OU NON
! --- SI LE CHAMP EXISTE, ON NE FAIT RIEN
!
    call exisd('CHAM_ELEM', chelem, iret)
    if (iret .eq. 0) then
!
! ---  BOUCLE SUR LES FISSURES
!
        do ifiss = 1, nfiss
!
! --- RECUPERATION NOMBRE DE POINTS DE GAUSS PAR FACETTE
!
            nomfis = fiss(ifiss)
            izone = xxconi(ds_contact%sdcont_defi,nomfis,'MAIT')
            typint = mminfi(ds_contact%sdcont_defi,'INTEGRATION',izone )
            call xmelin(elc, typint, nnint)
!
! --- RECUP LISTE DES MAILLES DE CONTACT POUR LA FISSURE
!
            grp = nomfis(1:8)//'.MAILFISS.CONT'
            call jeexin(grp, iret)
            nmaenr = 0
            if (iret .ne. 0) then
                call jeveuo(grp, 'L', jgrp)
                call jelira(grp, 'LONMAX', nmaenr)
            endif
!
! --- BOUCLE SUR LES MAILLES DE CONTACT DE LA FISSURE
!
            do i = 1, nmaenr
                ima = zi(jgrp-1+i)
                ASSERT(ima.le.nbma)
! --- NOMBRE DE FISSURES VUE PAR LA MAILLE
                nfisc = zi(jcesd2-1+5+4*(ima-1)+2)
                nfisc2 = zi(jcesd1-1+5+4*(ima-1)+2)
                ASSERT(nfisc.eq.nfisc2)
                nface = 0
! --- BOUCLE SUR CES FISSURES
                do ifima = 1, nfisc
                    call cesexi('S', jcesd1, jcesl1, ima, 1,&
                                ifima, 1, iad)
                    nomfi2 = cesv1(iad)
!
! --- S IL S AGIT DE LA FISSURE COURANTE
! --- ON RECUPERE LE NOMBRE DE FACETTES ET ON INCREMENTE
! --- LE NOMBRE DE SOUS-POINTS
!
                    if (nomfis .eq. nomfi2) then
                        ASSERT(nface.eq.0)
                        call cesexi('S', jcesd2, jcesl2, ima, 1,&
                                    ifima, 2, iad)
                        nface = cesv2(iad)
                        ASSERT(nface.le.30)
                        zi(jnbsp-1+ima) = zi(jnbsp-1+ima)+nface*nnint
                    endif
                end do
            end do
        end do
!
        if (param .eq. 'PCOHES') then
           if (lxthm) then
              call cescre('V', chelsi, 'ELEM', mesh, nomgd,&
                          ncmp, licmp5, [-1], zi( jnbsp), [-ncmp])
           else
              call cescre('V', chelsi, 'ELEM', mesh, nomgd,&
                          ncmp, licmp3, [-1], zi( jnbsp), [-ncmp])
           endif
        else
            call cescre('V', chelsi, 'ELEM', mesh, nomgd,&
                        1, 'X1', [-1], zi( jnbsp), [-ncmp])
        endif
!
! --- RAZ VECTEUR DE DIMENSIONNEMENT
!
    call jedetr(chnbsp)
    call wkvect(chnbsp, 'V V I', nbma, jnbsp)
!
! --- ACCES AU CHAM_ELEM_S
!
    call jeveuo(chelsi//'.CESD', 'L', jcesd)
    call jeveuo(chelsi//'.CESL', 'E', jcesl)
    call jeveuo(chelsi//'.CESV', 'E', jcesv)
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
        if (iret .ne. 0) call jeveuo(grp, 'L', jgrp)
!
! --- ZONE DE CONTACT IZONE CORRESPONDANTE
!
        izone = xxconi(ds_contact%sdcont_defi,nomfis,'MAIT')
        typint = mminfi(ds_contact%sdcont_defi,'INTEGRATION',izone )
        call xmelin(elc, typint, nnint)
!
! --- CONTACT INIT
!
        if (isint) then
            if (param .eq. 'PINDCOI' .or. param .eq. 'PMEMCON') then
                vall = mminfi(ds_contact%sdcont_defi,'CONTACT_INIT',izone ).eq.1
            else if (param.eq.'PGLISS') then
                vall = mminfl(ds_contact%sdcont_defi,'GLISSIERE_ZONE',izone )
            else
                ASSERT(.false.)
            endif
        else
            if (param .eq. 'PSEUIL') then
                valr = mminfr(ds_contact%sdcont_defi,'SEUIL_INIT',izone )
            else if (param(1:5).eq.'PCOHE') then
                valr = 0.d0
            else
                ASSERT(.false.)
            endif
            vall = valr.ne.0.d0
        endif
!
        if (vall) then
!
! --- ON COPIE LES CHAMPS CORRESP. AUX ELEM. DE CONTACT
!
            call jeexin(grp, iret)
            if (iret .ne. 0) then
                call jeveuo(grp, 'L', jgrp)
                call jelira(grp, 'LONMAX', nmaenr)
                do i = 1, nmaenr
                    ima = zi(jgrp-1+i)
!
! --- INDICE LOCAL DE LA FISSURE COURANTE
! --- ET DECALAGE CORRESPONDANT
!
                    nfisc = zi(jcesd1-1+5+4*(ima-1)+2)
                    do ifima = 1, nfisc
                        call cesexi('S', jcesd1, jcesl1, ima, 1,&
                                    ifima, 1, iad)
                        nomfi2 = cesv1(iad)
                        if (nomfis .eq. nomfi2) then
                            call cesexi('S', jcesd2, jcesl2, ima, 1,&
                                        ifima, 2, iad)
                            nface = cesv2(iad)
                            npg = nface*nnint
                        endif
                    end do
!
! --- RECOPIE EFFECTIVE DES CHAMPS
!
                    do ispt = 1, npg
                        do icmp = 1, ncmp
                            call cesexi('S', jcesd, jcesl, ima, 1,&
                                        zi(jnbsp-1+ima)+ispt, icmp, iad)
                            zl(jcesl-1+abs(iad)) = .true.
                            if (isint) zi(jcesv-1+abs(iad)) = 1
                            if (.not.isint) zr(jcesv-1+abs(iad)) = valr
                        end do
                    end do
!
! --- INCREMENTATION REPERAGE POUR LES FISSURES SUIVANTES
!
                    zi(jnbsp-1+ima) = zi(jnbsp-1+ima) + npg
!
                end do
            endif
        endif
    end do
!
! --- CONVERSION CHAM_ELEM_S -> CHAM_ELEM
!
    call cescel(chelsi, ligrel, option, param, 'OUI',&
                ibid, 'V', chelem, 'F', ibid)
!
! --- MENAGE
!
    call detrsd('CHAM_ELEM_S', chelsi)
    endif
    call detrsd('CHAM_ELEM_S', cmafis)
    call detrsd('CHAM_ELEM_S', faclon)
    call jedetr(chnbsp)
!
    call jedema()
!
end subroutine
