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

subroutine xtmafi(ndim, fiss, nfiss, lismai, &
                  mesmai, nbma, mesh, model, typ_enr)
!
! person_in_charge: samuel.geniaut at edf.fr
!
! aslint: disable=W1306
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/int_to_char8.h"
    integer(kind=8) :: nfiss, nbma, ndim
    character(len=8) :: fiss(nfiss)
    character(len=24) :: lismai, mesmai
    character(len=8), optional, intent(in) :: mesh, model
    character(len=*), optional, intent(in) :: typ_enr
!
! ----------------------------------------------------------------------
!
! ROUTINE XFEM
!
!  CREATION DE LA LISTE TOTALE DES NOMS DES MAILLES FISSUREES
!  (MAILLES HEAV + CTIP + HECT) POUR LES FISSURES CONTENUES
!  DANS LA LISTE FISS
!  SI NDIM = 0 : ON LES PREND TOUTES
!  SI NDIM NON NUL : ON NE PREND QUE LES MAILLES DE DIMENSION NDIM
!
! ----------------------------------------------------------------------
!
! IN     NDIM   : DIMENSION DES MAILLES A LISTER
! IN     FISS   : LISTE DES NOMS DES SD FISS_XFEM
! IN     NFISS  : LONGUEUR DE FISS
! IN/OUT LISMAI : NOM DE LA LISTE CREEE CONTENANT LES NUMEROS DE MAILLES
! IN/OUT MESMAI : NOM DE LA LISTE CREEE CONTENANT LES NOMS DES MAILLES
! OUT    NBMA   : LONGUEUR DE MESMAI
! IN (o) mesh    : optionnel / nom du maillage
! IN (o) model   : optionnel / nom du modele
! IN (o) typ_enr : optionnel / parmi 'HEAV', 'CTIP', 'HECT'
!
! regles sur les arguments optionnels :
! -------------------------------------
! - mesh ou model doit etre present (ou exclusif)
! - si mesh est present  -> on prend toutes les mailles sachant ndim
! - si model est present -> on prend toutes les mailles affectees dans
!                           model sachant ndim
! - si typ_enr present (parmi 'HEAV', 'CTIP', 'HECT'), on ne garde que
!   les mailles de types typ_enr
!
    integer(kind=8) :: ifiss, kk, jgrp, nmaenr, i, ima, cpt, iret
    integer(kind=8) ::   ndime, jmad, mxstac
    character(len=8) :: noma, nomafi, nomail, k8_typ_enr, vk8_typ_enr(3)
    character(len=8) :: k8_test
    character(len=24) :: grp(nfiss, 3)
    integer(kind=8), pointer :: temi(:) => null()
    character(len=8), pointer :: temp(:) => null()
    integer(kind=8), pointer :: tmdim(:) => null()
    integer(kind=8), pointer :: typmail(:) => null()
    integer(kind=8), pointer :: p_mail_affe(:) => null()
    aster_logical :: lmesh, lmodel, l_mail_affe, l_group_ok
!
    parameter(mxstac=100)
!
    data vk8_typ_enr/'HEAV', 'CTIP', 'HECT'/
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
!   VERIF QUE LES TABLEAUX LOCAUX DYNAMIQUES NE SONT PAS TROP GRANDS
!   (VOIR CRS 1404)
    ASSERT(nfiss .le. mxstac)
!
    lmesh = .false.
    lmodel = .false.
!
! - Verification ou exclusif pour les arguments optionnels mesh / model
!
    ASSERT(present(mesh) .or. present(model))
    if (present(mesh)) then
        lmesh = .true.
        ASSERT(.not. present(model))
    end if
    if (present(model)) then
        lmodel = .true.
        ASSERT(.not. present(mesh))
    end if
!
! - Verification sur argument optionnel typ_enr (3 valeurs aurorisees)
!
    k8_typ_enr = ''
    if (present(typ_enr)) then
        k8_typ_enr = typ_enr
        ASSERT(indik8(vk8_typ_enr, k8_typ_enr, 1, 3) .gt. 0)
    end if
!
! - Recuperation de l'objet '.TYPMAIL' pour filtrer sur ndim
!
    if (present(mesh)) then
        noma = mesh
    else
        call dismoi('NOM_MAILLA', model, 'MODELE', repk=noma)
    end if
    call jeveuo('&CATA.TM.TMDIM', 'L', vi=tmdim)
    call jeveuo(noma//'.TYPMAIL', 'L', vi=typmail)
!
! - Si model present, recuperation de l'objet '.MAILLE' pour filtrer
!   sur les mailles affectees
!
    if (present(model)) then
        call jeveuo(model//'.MAILLE', 'L', vi=p_mail_affe)
    end if
!
! - Dimensionnement grossier de la liste
!
    cpt = 0
    do ifiss = 1, nfiss

!       verif coherence maillage in <-> maillage de definition de la fissure
        call dismoi('NOM_MAILLA', fiss(ifiss), 'FISS_XFEM', repk=nomafi)
        ASSERT(nomafi .eq. noma)
!
        grp(ifiss, 1) = fiss(ifiss)//'.MAILFISS.HEAV'
        grp(ifiss, 2) = fiss(ifiss)//'.MAILFISS.CTIP'
        grp(ifiss, 3) = fiss(ifiss)//'.MAILFISS.HECT'
!
        do kk = 1, 3
            k8_test = grp(ifiss, kk) (19:22)
            call jeexin(grp(ifiss, kk), iret)
            l_group_ok = .not. (present(typ_enr))
            l_group_ok = l_group_ok .or. (k8_typ_enr .eq. k8_test)
            l_group_ok = l_group_ok .and. (iret .ne. 0)
            if (l_group_ok) then
                call jelira(grp(ifiss, kk), 'LONMAX', nmaenr)
                cpt = cpt+nmaenr
            end if
        end do
!
    end do
    ASSERT(cpt .gt. 0)
!
! - Creation des listes temporaires
!
    AS_ALLOCATE(vk8=temp, size=cpt)
    AS_ALLOCATE(vi=temi, size=cpt)
!
! - Remplissage des listes temporaires
!
    nbma = 0
!
!   boucle sur les fissures
    do ifiss = 1, nfiss
!
!       boucle sur les 3 groupes HEAV, CTIP et HECT
        do kk = 1, 3
!
            k8_test = grp(ifiss, kk) (19:22)
            call jeexin(grp(ifiss, kk), iret)
            l_group_ok = .not. (present(typ_enr))
            l_group_ok = l_group_ok .or. (k8_typ_enr .eq. k8_test)
            l_group_ok = l_group_ok .and. (iret .ne. 0)
!
!           si les conditions pour scruter le groupe sont remplies
            if (l_group_ok) then
!
                call jeveuo(grp(ifiss, kk), 'L', jgrp)
                call jelira(grp(ifiss, kk), 'LONMAX', nmaenr)
!
!               boucle sur les mailles de chaque groupe
                do i = 1, nmaenr
!
                    ima = zi(jgrp-1+i)
!                   ndime : dimension topologique de la maille
                    ndime = tmdim(typmail(ima))
                    if ((ndim .eq. ndime) .or. (ndim .eq. 0)) then
!
!                       on retient la maille si mesh present
!                       ou si model present et maille affectee
                        l_mail_affe = .false.
                        if (lmodel) then
                            l_mail_affe = p_mail_affe(ima) .ne. 0
                        end if
                        if (lmesh .or. l_mail_affe) then
                            nomail = int_to_char8(ima)
                            nbma = nbma+1
                            temp(nbma) = nomail
                            temi(nbma) = ima
                        end if
!
                    end if
!
!               fin boucle sur les mailles de chaque groupe
                end do
!
            end if
!
!       fin boucle sur les 3 groupes HEAV, CTIP et HECT
        end do

!   fin boucle sur les fissures
    end do
!
! - Verification
!
    ASSERT(nbma .le. cpt)
    ASSERT(nbma .ge. 1)
!
! - Creation des listes definitives
!
    call wkvect(mesmai, 'V V K8', nbma, jmad)
    do i = 1, nbma
        zk8(jmad-1+i) = temp(i)
    end do
    call wkvect(lismai, 'V V I', nbma, jmad)
    do i = 1, nbma
        zi(jmad-1+i) = temi(i)
    end do
!
! - Menage
!
    AS_DEALLOCATE(vk8=temp)
    AS_DEALLOCATE(vi=temi)
!
    call jedema()
end subroutine
