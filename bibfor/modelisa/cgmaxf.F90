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

subroutine cgmaxf(mofaz, iocc, nomaz, lismaz, nbma)
    implicit none
#include "jeveux.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/gmgnre.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/wkvect.h"
#include "asterfort/xtmafi.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/int_to_char8.h"
!
    integer(kind=8) :: iocc, nbma
    character(len=*) :: mofaz, nomaz, lismaz
!
!       CGMAXF -- TRAITEMENT DE L'OPTION FISS_XFEM
!                 DU MOT FACTEUR CREA_GROUP_MA DE
!                 LA COMMANDE DEFI_GROUP
!
! -------------------------------------------------------
!  MOFAZ         - IN    - K16  - : MOT FACTEUR 'CREA_GROUP_MA'
!  IOCC          - IN    - I    - : NUMERO D'OCCURENCE DU MOT-FACTEUR
!  NOMAZ         - IN    - K8   - : NOM DU MAILLAGE
!  LISMAZ        - JXVAR - K24  - : NOM DE LA LISTE DE MAILLES
!                                   DU TYPE XFEM DEMANDE PAR
!                                   L'UTILISATEUR
!  NBMA          - OUT   -  I   - : LONGUEUR DE CETTE LISTE
! -------------------------------------------------------
!
!
    integer(kind=8) :: i, ii, ima, ifiss, ino, n
    integer(kind=8) :: nbno, nbnot, nfiss, nmax, nbmalo, nbmala
    integer(kind=8) :: jlmas, idlist, jtem3, jtem4
    integer(kind=8) :: ibid, test, valeno
    character(len=8) :: noma, nomail, fiss
    character(len=16) :: motfac, typgrp
    character(len=19) :: stno
    character(len=24) :: lismai, lismar, lisman, maifis
    character(len=8), pointer :: vfiss(:) => null()
    integer(kind=8), pointer :: tem5(:) => null()
    integer(kind=8), pointer :: vale(:) => null()
!
!     -----------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATIONS :
!     ================
    nbmala = 0
    nbma = 0
    motfac = mofaz
    noma = nomaz
    lismai = lismaz
!
! --  RECUPERATION DU TYPE DE MAILLE XFEM :
    call getvtx(motfac, 'TYPE_GROUP', iocc=iocc, scal=typgrp, nbret=ibid)
    call getvid(motfac, 'FISSURE', iocc=iocc, nbval=0, nbret=nfiss)
    nfiss = -nfiss
    AS_ALLOCATE(vk8=vfiss, size=nfiss)
    call getvid(motfac, 'FISSURE', iocc=iocc, nbval=nfiss, vect=vfiss, &
                nbret=ibid)
!
!
! --- TYPE DE MAILLE = 'HEAVISIDE', 'CRACKTIP' OU  'MIXTE'
!     ====================================================
    if ((typgrp .eq. 'HEAVISIDE') .or. (typgrp .eq. 'CRACKTIP') .or. (typgrp .eq. 'MIXTE')) then
!
        call wkvect('&&CGMAXF.TEM3', 'V V I', nfiss, jtem3)
        call wkvect('&&CGMAXF.TEM4', 'V V I', nfiss, jtem4)
!
! ---   TYPE DE MAILLE = 'HEAVISIDE'
        if (typgrp .eq. 'HEAVISIDE') then
            maifis = '.MAILFISS.HEAV'
!
! ---   TYPE DE MAILLE = 'CRACKTIP'
        else if (typgrp .eq. 'CRACKTIP') then
            maifis = '.MAILFISS.CTIP'
!
! ---   TYPE DE MAILLE =
        else if (typgrp .eq. 'MIXTE') then
            maifis = '.MAILFISS.HECT'
        end if
! ---   BOUCLE SUR TOUTES LES FISSURES
        do ifiss = 1, nfiss
            fiss = vfiss(ifiss)
            call jeveuo(fiss//maifis, 'L', zi(jtem3-1+ifiss))
            call jelira(fiss//maifis, 'LONMAX', n)
            zi(jtem4-1+ifiss) = n
            nbma = nbma+n
        end do
! ---   CREATION ET REMPLISSAGE DU VECTEUR DE SORTIE
        if (nbma .gt. 0) then
            call wkvect(lismai, 'V V I', nbma, idlist)
            do ifiss = 1, nfiss
                jlmas = zi(jtem3-1+ifiss)
                do i = 1, zi(jtem4-1+ifiss)
                    nbmala = nbmala+1
                    zi(idlist+nbmala-1) = zi(jlmas+i-1)
                    nomail = int_to_char8(zi(jlmas+i-1))
                end do
            end do
        end if
        call jedetr('&&CGMAXF.TEM3')
        call jedetr('&&CGMAXF.TEM4')
!
! --- TYPE DE MAILLE = 'XFEM'
!     ============================
    else if (typgrp .eq. 'XFEM') then
!
        lisman = '&&CGMAXF.TEM1'
        call xtmafi(0, vfiss, nfiss, lismai, &
                    lisman, nbma, mesh=noma)
        call jedetr(lisman)
!
! --- TYPE DE MAILLE = 'FISSUREE'
!     ============================
    else if (typgrp .eq. 'FISSUREE') then
!
        call dismoi('NB_NO_MAILLA', noma, 'MAILLAGE', repi=nbnot)
        call dismoi('NB_NO_MAX', '&CATA', 'CATALOGUE', repi=nmax)
!
!       POUR DIMENSIONNER GROSSIEREMENT LA LISTE DES MAILLES
        lisman = '&&CGMAXF.TEM1'
        lismar = '&&CGMAXF.TEM2'
        call xtmafi(0, vfiss, nfiss, lismar, &
                    lisman, nbma, mesh=noma)
        call jedetr(lismar)
        call jedetr(lisman)
!
        call wkvect('&&CGMAXF.TEM3', 'V V I', nbnot, jtem3)
        call wkvect('&&CGMAXF.TEM4', 'V V I', nmax, jtem4)
        AS_ALLOCATE(vi=tem5, size=nbma)
!
        nbmalo = 0
!
!       POUR CHAQUE FISSURE
        do ifiss = 1, nfiss
            fiss = vfiss(ifiss)
            stno = fiss//'.STNO'
            call jeveuo(stno//'.VALE', 'L', vi=vale)
!
!         RECUPERATION DE TOUTES MAILLES XFEM DE LA FISSURE COURANTE
            lisman = '&&CGMAXF.TEM1'
            call xtmafi(0, fiss, 1, lismar, &
                        lisman, nbmala, mesh=noma)
            call jeveuo(lismar, 'L', jlmas)
!
!         POUR CHAQUE MAILLE XFEM DE LA FISSURE COURANTE
            do ii = 1, nbmala
!           RECUPERATION DES NOEUDS
                call gmgnre(noma, nbnot, zi(jtem3), zi(jlmas+ii-1), 1, &
                            zi(jtem4), nbno, 'TOUS')
!           TRI DES NOEUDS SELON LEUR STATUS XFEM
                test = 1
                do ino = 1, nbno
                    valeno = vale(1+zi(jtem4+ino-1)-1)
                    if (valeno .eq. 0) then
                        test = 0
                    end if
                end do
!           MAILLES QUI REPOSENT SUR LES NOEUDS AU STATUS <> 0
                if (test .eq. 1) then
                    nbmalo = nbmalo+1
                    tem5(nbmalo) = zi(jlmas+ii-1)
                    nomail = int_to_char8(tem5(nbmalo))
                end if
            end do
            call jedetr(lismar)
            call jedetr(lisman)
            call jedetr(stno)
        end do
        call wkvect(lismai, 'V V I', nbmalo, idlist)
        do ima = 1, nbmalo
            zi(idlist-1+ima) = tem5(ima)
        end do
        call jedetr('&&CGMAXF.TEM3')
        call jedetr('&&CGMAXF.TEM4')
        AS_DEALLOCATE(vi=tem5)
        nbma = nbmalo
!
    end if
!
! --- FIN
!     ===
!
! --- MENAGE
!
    AS_DEALLOCATE(vk8=vfiss)
!
    call jedema()
!
end subroutine
