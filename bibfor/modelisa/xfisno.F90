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

subroutine xfisno(noma, modelx)
! person_in_charge: patrick.massin at edf.fr
    implicit none
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/celces.h"
#include "asterfort/cescel.h"
#include "asterfort/cescre.h"
#include "asterfort/cesexi.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
!
    character(len=8) :: noma, modelx
!
!----------------------------------------------------------------------
!  BUT: CREATION D'UN CHAMPS ELNO QUI ASSOCIE POUR CHAQUE NOEUD LE
!       NUMÉRO DE FISSURE LOCALE AU DDL HEAVISIDE
!
!----------------------------------------------------------------------
!
!     ARGUMENTS/
!  NOMA       IN       K8 : MAILLAGE
!  MODELX     IN/OUT   K8 : MODELE XFEM
!
!
!
!
!
    integer(kind=8) :: jlcnx, jcesfd, jcesfl, jcesd, jcesl
    integer(kind=8) :: jcesd2, jcesl2
    integer(kind=8) :: nbma, ima, nbno, ino, nheav, iheav, nfiss, ifiss
    integer(kind=8) :: ibid, iad, nncp
    integer(kind=8), pointer :: p_mail_affe(:) => null()
    character(len=19) :: fissno, ces, cesf, ligrel, ces2, heavno
    aster_logical :: lcont
    integer(kind=8), pointer :: xfem_cont(:) => null()
    integer(kind=8), pointer :: cesfv(:) => null()
    integer(kind=8), pointer :: cesv2(:) => null()
    integer(kind=8), pointer :: cesv(:) => null()
    integer(kind=8), pointer :: nbsp2(:) => null()
    integer(kind=8), pointer :: nbsp(:) => null()
!     ------------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATIONS
!
    call dismoi('NOM_LIGREL', modelx, 'MODELE', repk=ligrel)
    fissno = modelx(1:8)//'.FISSNO'
    ces = '&&XFISNO.FISSNO'
    cesf = '&&XFISNO.STNO'
!
! --- LE CONTACT EST-IL DÉCLARÉ
!
    call jeveuo(modelx(1:8)//'.XFEM_CONT', 'L', vi=xfem_cont)
    ASSERT(xfem_cont(1) .le. 1 .or. xfem_cont(1) .eq. 2 .or. xfem_cont(1) .eq. 3)
    lcont = (xfem_cont(1) .eq. 1 .or. xfem_cont(1) .eq. 2 .or. xfem_cont(1) .eq. 3)
    if (lcont) then
        heavno = modelx(1:8)//'.HEAVNO'
        ces2 = '&&XFISNO.HEAVNO'
    end if
!
! --- TRANSFO CHAM_ELEM -> CHAM_ELEM_S DE STANO
!
    call celces(modelx(1:8)//'.STNO', 'V', cesf)
!
    call jeveuo(cesf//'.CESD', 'L', jcesfd)
    call jeveuo(cesf//'.CESV', 'L', vi=cesfv)
    call jeveuo(cesf//'.CESL', 'L', jcesfl)
!
! --- RECUPERATION DE LA LISTE DES MAILLES AFFECTEES PAR DES EF
!
    call jeveuo(ligrel//'.TYFE', 'L', vi=p_mail_affe)
!
! --- RECUPERATION DU NOMBRE DE FISSURES VUES
!
    call jeveuo('&&XTYELE.NBSP', 'L', vi=nbsp)
!
!
! --- RECUPERATION DU NOMBRE DE FONCTIONS HEAVISIDES
!
    call jeveuo('&&XTYELE.NBSP2', 'L', vi=nbsp2)
!
! --- CREATION DE LA SD ELNO FISSNO
!
    call cescre('V', ces, 'ELNO', noma, 'NEUT_I', &
                1, 'X1', [ibid], nbsp2, [-1])
!
    call jeveuo(ces//'.CESD', 'L', jcesd)
    call jeveuo(ces//'.CESV', 'E', vi=cesv)
    call jeveuo(ces//'.CESL', 'E', jcesl)
!
! --- SI CONTACT, CREATION DE LA SD ELNO HEAVNO
!
    if (lcont) then
        call cescre('V', ces2, 'ELNO', noma, 'NEUT_I', &
                    1, 'X1', [ibid], nbsp, [-1])
!
        call jeveuo(ces2//'.CESD', 'L', jcesd2)
        call jeveuo(ces2//'.CESV', 'E', vi=cesv2)
        call jeveuo(ces2//'.CESL', 'E', jcesl2)
    end if
!
! --- INFOS SUR LE MAILLAGE
!
    call dismoi('NB_MA_MAILLA', noma, 'MAILLAGE', repi=nbma)
    call jeveuo(jexatr(noma//'.CONNEX', 'LONCUM'), 'L', jlcnx)
!
    do ima = 1, nbma
!       la maille est-elle affectee par un ef ?
        if (p_mail_affe(ima) .eq. 0) cycle
        nfiss = nbsp(ima)
        nheav = nbsp2(ima)
        if (nfiss .ge. 2) then
            nbno = zi(jlcnx+ima)-zi(jlcnx-1+ima)
            do ino = 1, nbno
!
! --- PREMIERE PASSE, ON REMPLIT AVEC LES HEAVISIDES ACTIFS
!
                do ifiss = 1, nfiss
                    call cesexi('S', jcesfd, jcesfl, ima, ino, &
                                ifiss, 1, iad)
                    ASSERT(iad .gt. 0)
                    if (cesfv(iad) .eq. 1) then
                        do iheav = 1, nheav
                            call cesexi('S', jcesd, jcesl, ima, ino, &
                                        iheav, 1, iad)
                            if (iad .lt. 0) then
                                zl(jcesl-1-iad) = .true.
                                cesv(1-1-iad) = ifiss
                                if (lcont) then
                                    call cesexi('S', jcesd2, jcesl2, ima, ino, &
                                                ifiss, 1, iad)
                                    ASSERT(iad .lt. 0)
                                    zl(jcesl2-1-iad) = .true.
                                    cesv2(1-1-iad) = iheav
                                end if
                                goto 30
                            end if
                        end do
                    end if
30                  continue
                end do
!
! --- DEUXIEME PASSE, ON REMPLIT AVEC LES HEAVISIDES INACTIFS
!
                do ifiss = 1, nfiss
                    call cesexi('S', jcesfd, jcesfl, ima, ino, &
                                ifiss, 1, iad)
                    ASSERT(iad .gt. 0)
                    if (cesfv(iad) .eq. 0) then
                        do iheav = 1, nheav
                            call cesexi('S', jcesd, jcesl, ima, ino, &
                                        iheav, 1, iad)
                            if (iad .lt. 0) then
                                zl(jcesl-1-iad) = .true.
                                cesv(1-1-iad) = ifiss
                                if (lcont) then
                                    call cesexi('S', jcesd2, jcesl2, ima, ino, &
                                                ifiss, 1, iad)
                                    ASSERT(iad .lt. 0)
                                    zl(jcesl2-1-iad) = .true.
                                    cesv2(1-1-iad) = iheav
                                end if
                                goto 50
                            end if
                        end do
                    end if
50                  continue
                end do
!
! --- FIN DES DEUX PASSES, FISSNO EST DEFINI ENTIEREMENT POUR LE NOEUD
! ---                      HEAVNO N'EST PAS COMPLET
!
            end do
        end if
    end do
!
! --- CONVERSION CHAM_ELEM_S -> CHAM_ELEM
!
    call cescel(ces, ligrel, 'INI_XFEM_ELNO', 'PFISNO', 'OUI', &
                nncp, 'G', fissno, 'F', ibid)
    call detrsd('CHAM_ELEM_S', ces)
    if (lcont) then
        call cescel(ces2, ligrel, 'FULL_MECA', 'PHEAVNO', 'NAN', &
                    nncp, 'G', heavno, 'F', ibid)
        call detrsd('CHAM_ELEM_S', ces2)
    end if
    call detrsd('CHAM_ELEM_S', cesf)
    call jedema()
end subroutine
