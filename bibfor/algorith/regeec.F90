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

subroutine regeec(nomres, resgen, nomsst)
    implicit none
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/dcapno.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvis.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/mgutdm.h"
#include "asterfort/nueq_chck.h"
#include "asterfort/refdcp.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rscrsd.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnoch.h"
#include "asterfort/rsorac.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"
#include "asterfort/vtcrea.h"
#include "asterfort/wkvect.h"
!
    character(len=8) :: nomres, resgen, nomsst
!
!  BUT : < RESTITUTION GENERALISEE ECLATEE >
!
!  RESTITUER EN BASE PHYSIQUE SUR UNE SOUS-STRUCTURE LES RESULTATS
!  ISSUS DE LA SOUS-STRUCTURATION GENERALE
!  LE CONCEPT RESULTAT EST UN RESULTAT COMPOSE "MODE_MECA"
!
!  /!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
!
!--   LES ROUTINES REGEEC ET REGE2C FONT LA MEME CHOSE, UNE EN REEL,
!--   L'AUTRE EN COMPLEXE. EN CAS DE MODIFICATION D'UNE DES ROUTINES,
!--   NE PAS OUBLIER DE REPORTER LE CHANGEMENT DANS L'AUTRE.
!
!  /!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
!-----------------------------------------------------------------------
!
! NOMRES /I/ : NOM K8 DU CONCEPT MODE MECA RESULTAT
! RESGEN /I/ : NOM K8 DU MODE_GENE AMONT
! NOMSST /I/ : NOM K8 DE LA SOUS-STRUCTURE SUR LAQUELLE ON RESTITUE
!
!
!
!
!
    integer(kind=8) :: i, iad, ibid, ieq, ier, iord, j, jbid, k, llchab, i_ligr_ss
    integer(kind=8) :: llchol, llors, llprs, vali(2), nbbas, nbddg, nbmod(1), nbsst
    integer(kind=8) :: neq, nno, numo, nusst, nutars, iadpar(8)
    integer(kind=8) :: elim, neqet, neqred, lmapro, lsilia, lsst, lmoet, i1, k1
    real(kind=8) :: freq, genek, genem, omeg2, rbid, genec, amor
    character(len=8) :: kbid, basmod, mailla, lint, model_gene
    character(len=16) :: depl, nompar(8), typres, quamod
    character(len=19) :: raid, numddl, chamne, nume_equa_gene
    character(len=14) :: nume_gene
    character(len=24) :: crefe(2), chamol, chamba
    character(len=24) :: valk(2), seliai, sizlia, sst
    character(len=3) :: typesca
    complex(kind=8) :: cbid
    character(len=24), pointer :: refa(:) => null()
    integer(kind=8), pointer :: nueq(:) => null()
    character(len=24), pointer :: refn(:) => null()
    real(kind=8), pointer :: vale(:) => null()
!
!-----------------------------------------------------------------------
    data depl/'DEPL            '/
    data nompar/'FREQ', 'RIGI_GENE', 'MASS_GENE', 'AMOR_GENE', 'OMEGA2',&
     &            'NUME_MODE', 'AMOR_REDUIT', 'TYPE_MODE'/
!-----------------------------------------------------------------------
!
    call jemarq()
    call titre()
!
! --- RECUPERATION DU MODELE GENERALISE
!
    call dismoi('REF_RIGI_PREM', resgen, 'RESU_DYNA', repk=raid)
!
    call jeveuo(raid//'.REFA', 'L', vk24=refa)
    nume_gene = refa(2) (1:14)
    nume_equa_gene = nume_gene(1:14)//'.NUME'
    call nueq_chck(nume_equa_gene, neqred)
    call jelibe(raid//'.REFA')
!
    call jeveuo(nume_equa_gene//'.REFN', 'L', vk24=refn)
    model_gene = refn(1) (1:8)
    call jelibe(nume_equa_gene//'.REFN')
!
! --- RECUPERATION NUMERO DE SOUS-STRUCTURE
!     ET DU NOEUD TARDIF CORRESPONDANT
!
    call jenonu(jexnom(model_gene//'      .MODG.SSNO', nomsst), nusst)
    if (nusst .eq. 0) then
        valk(1) = model_gene
        valk(2) = nomsst
        call utmess('F', 'ALGORITH14_25', nk=2, valk=valk)
    end if
!
!
!
!-- ON TESTE SI ON A EU RECOURS A L'ELIMINATION
!
    seliai = nume_gene(1:14)//'.ELIM.BASE'
    sizlia = nume_gene(1:14)//'.ELIM.TAIL'
    sst = nume_gene(1:14)//'.ELIM.NOMS'
!
    call jeexin(seliai, elim)
!
    if (elim .eq. 0) then
!
        call jenonu(jexnom(nume_equa_gene//'.LILI', '&SOUSSTR'), i_ligr_ss)
        call jeveuo(jexnum(nume_equa_gene//'.ORIG', i_ligr_ss), 'L', llors)
        call jelira(jexnum(nume_equa_gene//'.ORIG', i_ligr_ss), 'LONMAX', nbsst)
!
        nutars = 0
        do i = 1, nbsst
            if (zi(llors+i-1) .eq. nusst) nutars = i
        end do
!
!
        call jeveuo(jexnum(nume_equa_gene//'.PRNO', i_ligr_ss), 'L', llprs)
        nbddg = zi(llprs+(nutars-1)*2+1)
        ieq = zi(llprs+(nutars-1)*2)
!
    else
!
        call jelira(model_gene//'      .MODG.SSNO', 'NOMMAX', nbsst)
        call jeveuo(sst, 'L', ibid)
        do i1 = 1, nbsst
            if (nomsst .eq. zk8(ibid+i1-1)) then
                nusst = i1
            end if
        end do
        neqet = 0
        ieq = 0
!
        call jeveuo(seliai, 'L', lmapro)
        call jeveuo(sizlia, 'L', lsilia)
        call jeveuo(sst, 'L', lsst)
        ibid = 1
        do i = 1, nbsst
            neqet = neqet+zi(lsilia+i-1)
        end do
!
        ieq = 0
        do i1 = 1, nusst-1
            ieq = ieq+zi(lsilia+i1-1)
        end do
        call wkvect('&&MODE_ETENDU_REST_ELIM', 'V V R', neqet, lmoet)
!
    end if
!
! --- RECUPERATION DE LA BASE MODALE
!
    call mgutdm(model_gene, nomsst, ibid, 'NOM_BASE_MODALE', ibid, &
                basmod)
!
    call refdcp(basmod, nomres)
!
    call dismoi('NB_MODES_TOT', basmod, 'RESULTAT', repi=nbbas)
!
    if (elim .eq. 0) then
        if (nbbas .ne. nbddg) then
            valk(1) = basmod
            vali(1) = nbbas
            vali(2) = nbddg
            call utmess('F', 'ALGORITH14_26', sk=valk(1), ni=2, vali=vali)
        end if
    end if
!
    call dismoi('REF_INTD_PREM', basmod, 'RESU_DYNA', repk=lint)
!
    call dismoi('NOM_MAILLA', lint, 'INTERF_DYNA', repk=mailla)
    call dismoi('NOM_NUME_DDL', lint, 'INTERF_DYNA', repk=numddl)
    call dismoi('NB_EQUA', numddl, 'NUME_DDL', repi=neq)
!
    crefe(1) = mailla
    crefe(2) = numddl
!
! --- RECUPERATION NOMBRE DE MODES PROPRES CALCULES
!
    call rsorac(resgen, 'LONUTI', 0, rbid, kbid, &
                cbid, rbid, kbid, nbmod, 1, &
                ibid)
!
! --- ON RESTITUE SUR TOUS LES MODES OU SUR QUELQUES MODES:
!
!
    call getres(kbid, typres, quamod)
    if (quamod .ne. 'CALC_CORR_SSD') then
        call getvis(' ', 'NUME_ORDRE', nbval=0, nbret=nno)
    else
!-- SI ON APPELLE DEPUIS QUAL_MODL, ON RESTITUE TOUS LES MODES
        nno = 0
    end if
!
    if (nno .ne. 0) then
        nbmod(1) = -nno
        call wkvect('&&REGEEC.NUME', 'V V I', nbmod(1), jbid)
        call getvis(' ', 'NUME_ORDRE', nbval=nbmod(1), vect=zi(jbid), nbret=nno)
    else
        call wkvect('&&REGEEC.NUME', 'V V I', nbmod(1), jbid)
        do i = 1, nbmod(1)
            zi(jbid+i-1) = i
        end do
    end if
!
! --- ALLOCATION STRUCTURE DE DONNEES RESULTAT
!
    call rscrsd('G', nomres, 'MODE_MECA', nbmod(1))
!
! --- RESTITUTION PROPREMENT DITE
!
    call jeveuo(nume_equa_gene//'.NUEQ', 'L', vi=nueq)
!
! --- BOUCLE SUR LES MODES A RESTITUER
    do i = 1, nbmod(1)
        iord = zi(jbid+i-1)
!
! ----- REQUETTE NOM ET ADRESSE CHAMNO GENERALISE
        call dcapno(resgen, depl, iord, chamol)
        call dismoi('TYPE_SCA', chamol(1:19), 'CHAMP', repk=typesca)
        if (typesca .ne. "R") then
            call utmess('F', 'SOUSTRUC_84')
        end if
        call jeveuo(chamol, 'L', llchol)
!-- SI ELIMINATION, ON RESTITUE D'ABORD LES MODES GENERALISES
        if (elim .ne. 0) then
            do i1 = 1, neqet
                zr(lmoet+i1-1) = 0.d0
                do k1 = 1, neqred
                    zr(lmoet+i1-1) = zr(lmoet+i1-1)+zr(lmapro+(k1-1)* &
                                                       neqet+i1-1)*zr(llchol+k1-1)
                end do
            end do
            llchol = lmoet
        end if
!
! ----- REQUETTE NOM ET ADRESSE NOUVEAU CHAMNO
        call rsexch(' ', nomres, depl, i, chamne, &
                    ier)
        call vtcrea(chamne, crefe, 'G', 'R', neq)
        call jeveuo(chamne//'.VALE', 'E', vr=vale)
!
        call rsadpa(resgen, 'L', 7, nompar, iord, &
                    0, tjv=iadpar, styp=kbid)
        freq = zr(iadpar(1))
        genek = zr(iadpar(2))
        genem = zr(iadpar(3))
        genec = zr(iadpar(4))
        omeg2 = zr(iadpar(5))
        numo = zi(iadpar(6))
        amor = zr(iadpar(7))
!
! ----- BOUCLE SUR LES MODES PROPRES DE LA BASE
        if (elim .ne. 0) then
            ibid = nbbas
        else
            ibid = nbddg
        end if
        do j = 1, ibid
            call dcapno(basmod, depl, j, chamba)
            call jeveuo(chamba, 'L', llchab)
!
! ------- BOUCLE SUR LES EQUATIONS PHYSIQUES
            do k = 1, neq
                if (elim .ne. 0) then
                    iad = llchol+ieq+j-1
                else
                    iad = llchol+nueq(1+ieq+j-2)-1
                end if
                vale(k) = vale(k)+zr(llchab+k-1)*zr(iad)
            end do
            call jelibe(chamba)
        end do
        call rsnoch(nomres, depl, i)
        call rsadpa(nomres, 'E', 8, nompar, i, &
                    0, tjv=iadpar, styp=kbid)
        zr(iadpar(1)) = freq
        zr(iadpar(2)) = genek
        zr(iadpar(3)) = genem
        zr(iadpar(4)) = geneC
        zr(iadpar(5)) = omeg2
        zi(iadpar(6)) = numo
        zr(iadpar(7)) = amor
        zk16(iadpar(8)) = 'MODE_DYN'
!
        call jelibe(chamol)
    end do
!
    call jelibe(nume_equa_gene//'.NUEQ')
    call jedetr('&&REGEEC.NUME')
!
    call jedema()
end subroutine
