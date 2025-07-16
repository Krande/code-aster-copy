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
subroutine caechp(load, loadLigrel, mesh, model, geomDime, valeType)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/alcart.h"
#include "asterfort/assert.h"
#include "asterfort/char_nb_ligf.h"
#include "asterfort/char_crea_ligf.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/nocart.h"
#include "asterfort/paligi.h"
#include "asterfort/palima.h"
#include "asterfort/patrma.h"
#include "asterfort/utmess.h"
#include "asterfort/xtempc.h"
#include "asterfort/xtmafi.h"
#include "asterfort/xvelfm.h"
!
    character(len=8), intent(in) :: load
    character(len=19), intent(in) :: loadLigrel
    character(len=8), intent(in) :: mesh, model
    integer(kind=8), intent(in) :: geomDime
    character(len=4), intent(in) :: valeType
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Treatment of load ECHANGE_PAROI
!
! --------------------------------------------------------------------------------------------------
!
! In  load             : load
! In  loadLigrel
! In  mesh             : mesh
! In  model            : model
! In  geomDime         : dimension of space
! In  valeType         : affected value type (real, complex or function)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: keywordFact = 'ECHANGE_PAROI'
    integer(kind=8), parameter :: nfismx = 100, nbtymx = 7
    integer(kind=8) :: nechp, ibid, jvalv, iocc, nh, nt, j
    integer(kind=8) :: nbtyp, jlistt, nbm, nfiss, jma, ntcon
    aster_logical :: ltcon, lcoefh
    integer(kind=8) :: igrel, inema
    integer(kind=8) :: jligr, ncmp
    real(kind=8) :: t(3), cechpr
    character(len=8) :: k8b, cechpf, fiss(nfismx)
    character(len=24) :: liel, modelisa, llist1, llist2, llistt
    character(len=19) :: carte
    integer(kind=8) :: nb_elem_late, nb_noel_maxi
    character(len=24) :: mesmai, lismai
    character(len=8), pointer :: vncmp(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    call getfac(keywordFact, nechp)
    if (nechp .eq. 0) goto 999
!
    liel = loadLigrel(1:19)//'.LIEL'
    igrel = 0
    inema = 0
    call dismoi('MODELISATION', model, 'MODELE', repk=modelisa)

! - Count number of late elements
    call char_nb_ligf(mesh, keywordFact, 'Elem', nb_elem_late, nb_noel_maxi, &
                      suffix='_1')

! - Create <LIGREL> on late elements
    if (nb_elem_late .ne. 0) then
        call char_crea_ligf(mesh, loadLigrel, nb_elem_late, nb_noel_maxi)
    end if
!
!     LE MOT-CLE COEF_H EST-IL PRESENT ?
    lcoefh = .false.
    do iocc = 1, nechp
        if (valeType .eq. 'REEL') then
            call getvr8(keywordFact, 'COEF_H', iocc=iocc, scal=cechpr, nbret=nh)
        else if (valeType .eq. 'FONC') then
            call getvid(keywordFact, 'COEF_H', iocc=iocc, scal=cechpf, nbret=nh)
        end if
        if (nh .ne. 0) then
            lcoefh = .true.
            exit
        end if
    end do
!
!     SI LE MOT-CLE COEF_H EST PRESENT, ON ALLOUE ET PREAPRE LA CARTE
    if (lcoefh) then
        carte = load//'.CHTH.HECHP'
        if (valeType .eq. 'REEL') then
            call alcart('G', carte, mesh, 'COEH_R')
        else if (valeType .eq. 'FONC') then
            call alcart('G', carte, mesh, 'COEH_F')
        else
            call utmess('F', 'MODELISA2_37', sk=valeType)
        end if
!       NOM DE LA CMP DU COEFFICIENT D'ECHANGE DANS LA CARTE
        call jeveuo(carte//'.NCMP', 'E', vk8=vncmp)
        call jeveuo(carte//'.VALV', 'E', jvalv)
        ncmp = 1
        vncmp(1) = 'H'
    end if
!
! ----------------------------------------------------------------------
! --- BOUCLE SUR LES OCCURENCES DU MCF
! ----------------------------------------------------------------------
    do iocc = 1, nechp
!
!       RECUPERATION DU COEFFICIENT D'ECHANGE
        if (valeType .eq. 'REEL') then
            call getvr8(keywordFact, 'COEF_H', iocc=iocc, scal=cechpr, nbret=nh)
        else if (valeType .eq. 'FONC') then
            call getvid(keywordFact, 'COEF_H', iocc=iocc, scal=cechpf, nbret=nh)
        end if
!
!       RECUPERATION DU VECTEUR DE TRANSLATION POUR PATRMA
        t = 0.d0
        call getvr8(keywordFact, 'TRAN', iocc=iocc, nbval=3, vect=t, nbret=nt)
        call getvid(keywordFact, 'FISSURE', iocc=iocc, nbval=0, nbret=nfiss)
!
! ----------------------------------------------------------------------
! ----- CAS MOT-CLEF FISSURE (X-FEM)
! ----------------------------------------------------------------------
        if (nfiss .ne. 0) then
!
!         RECUPERATION DU NOM DES FISSURES
            nfiss = -nfiss
            call getvid(keywordFact, 'FISSURE', iocc=iocc, nbval=nfiss, vect=fiss, &
                        nbret=ibid)
!         VERIFICATION DE LA COHERENCE ENTRE LES FISSURES ET LE MODELE
            call xvelfm(nfiss, fiss, model)
!
!         ON SCRUTE LE MC TEMP_CONTINUE
            ltcon = .false.
            call getvtx(keywordFact, 'TEMP_CONTINUE', iocc=iocc, scal=k8b, nbret=ntcon)
!         VERIF DE COHERENCE AVEC LE MC COEF_H
            if (ntcon .eq. 1) then
                ASSERT(k8b(1:3) .eq. 'OUI' .and. nh .eq. 0)
                ltcon = .true.
            else
                ASSERT(nh .eq. 1 .and. ntcon .eq. 0)
            end if
!
! ----------------------------------------------------------------------
! ------- CAS TEMP_CONTINUE (X-FEM / PAS D'ECHANGE)
! ----------------------------------------------------------------------
            if (ltcon) then
!
                call xtempc(nfiss, fiss, valeType, load)
!
! ----------------------------------------------------------------------
! ------- CAS COEF_H (X-FEM / ECHANGE)
! ----------------------------------------------------------------------
            else
!
!           ON NOTE 0. OU '&FOZERO' DANS LA CARTE POUR TOUT LE MAILLAGE
                if (valeType .eq. 'REEL') then
                    zr(jvalv) = 0.d0
                else if (valeType .eq. 'FONC') then
                    zk8(jvalv) = '&FOZERO'
                end if
                call nocart(carte, 1, ncmp)
!
!           RECUPERATION DES MAILLES PRINCIPALES XFEM POUR FISS(1:NFISS)
                mesmai = '&&CAECHP.MES_MAILLES'
                lismai = '&&CAECHP.NUM_MAILLES'
                call xtmafi(geomDime, fiss, nfiss, lismai, &
                            mesmai, nbm, model=model)
                call jeveuo(mesmai, 'L', jma)
!
!           STOCKAGE DANS LA CARTE SUR CES MAILLES
                if (valeType .eq. 'REEL') then
                    zr(jvalv) = cechpr
                else if (valeType .eq. 'FONC') then
                    zk8(jvalv) = cechpf
                end if
                call nocart(carte, 3, ncmp, mode='NOM', nma=nbm, &
                            limano=zk8(jma))
!
!           MENAGE
                call jedetr(mesmai)
                call jedetr(lismai)
!
            end if
!
! ----------------------------------------------------------------------
! ----- CAS MOTS-CLEFS GROUP_MA_1... (PAROI MAILLEE)
! ----------------------------------------------------------------------
        else
!
            llist1 = '&&CAECHP.LLIST1'
            llist2 = '&&CAECHP.LLIST2'
            llistt = '&&CAECHP.LLIST.TRIE'
!
            call palima(mesh, keywordFact, 'GROUP_MA_1', 'MAILLE_1', iocc, &
                        llist1)
            call palima(mesh, keywordFact, 'GROUP_MA_2', 'MAILLE_2', iocc, &
                        llist2)
!
            call patrma(llist1, llist2, t, nbtymx, mesh, &
                        llistt, nbtyp)
!
!         MISE A JOUR DE LIGRCH ET STOCKAGE DANS LA CARTE
            do j = 1, nbtyp
                igrel = igrel+1
                call jeveuo(jexnum(llistt, j), 'L', jlistt)
                call paligi(modelisa, loadLigrel, igrel, inema, zi(jlistt))
!
!           STOCKAGE DANS LA CARTE
                call jeveuo(jexnum(liel, igrel), 'E', jligr)
                call jelira(jexnum(liel, igrel), 'LONMAX', nbm)
                nbm = nbm-1
                if (valeType .eq. 'REEL') then
                    zr(jvalv) = cechpr
                else if (valeType .eq. 'FONC') then
                    zk8(jvalv) = cechpf
                end if
                call nocart(carte, -3, ncmp, ligrel=loadLigrel, nma=nbm, &
                            limanu=zi(jligr))
            end do
!
!         MENAGE
            call jedetr(llist1)
            call jedetr(llist2)
            call jedetr(llistt)
! ------
        end if
!
    end do
!
999 continue
    call jedema()
end subroutine
