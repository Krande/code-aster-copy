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

subroutine rehaec(nomres, resgen, nomsst)
!    T. KERBER     DATE 14/05/93
!-----------------------------------------------------------------------
!  BUT:      < RESTITUTION HARMONIQUE ECLATEE >
    implicit none
!
!      RESTITUER EN BASE PHYSIQUE SUR UNE SOUS-STRUCTURE LES RESULTATS
!                ISSU DE LA SOUS-STRUCTURATION GENERALE
!  LE CONCEPT RESULTAT EST UN RESULTAT COMPOSE "DYNA_HARMO"
!
!-----------------------------------------------------------------------
!
! NOMRES   /I/: NOM K8 DU CONCEPT DYNA HARMO RESULTAT
! RESGEN   /I/: NOM K8 DU HARM_GENE AMONT
! NOMSST   /I/: NOM K8 DE LA SOUS-STRUCTURE SUR LAQUELLE ON RESTITUE
!
!
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dcapno.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
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
#include "asterfort/rstran.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"
#include "asterfort/vtcrea.h"
#include "asterfort/wkvect.h"
!
!
!
    real(kind=8) :: epsi
    character(len=4) :: champ(8)
    character(len=8) :: chmp(3), crit, interp, k8b, nomres, basmod, mailla
    character(len=8) :: lint, nomsst, model_gene, resgen, k8rep
    character(len=14) :: nume_gene
    character(len=19) :: numddl, nume_equa_gene, knume, kfreq, harmge
    character(len=24) :: crefe(2), chamba, chamno, seliai, sizlia, sst
    character(len=24) :: valk(2)
    integer(kind=8) :: itresu(3), elim, iret
    integer(kind=8) :: i, i1, iad, iarchi, ibid, ich, i_ligr_ss
    integer(kind=8) :: idresu, ieq, ire1, ire2, ire3
    integer(kind=8) :: iretou, j, jfreq, jnume, k, k1, ldnew
    integer(kind=8) :: lfreq, llchab, llors, llprs
    integer(kind=8) :: lmapro, lmoet, lsilia, lsst
    integer(kind=8) :: n1, nbcham, nbddg, nbfreq, nbsst, neq
    integer(kind=8) :: neqet, neqgen, neqred, nusst, nutars
    character(len=24), pointer :: refn(:) => null()
    integer(kind=8), pointer :: nueq(:) => null()
!-----------------------------------------------------------------------
!
! --- ECRITURE DU TITRE
!
    call jemarq()
    call titre()
!
! --- DETERMINATION DES CHAMPS A RESTITUER, PARMI DEPL, VITE ET ACCE
!
    harmge = resgen
! --- CALCUL DU NOMBRE DE CHAMPS A RESTITUER ET LEURS ADDRESSES
!
    call jeexin(resgen//'           .DEPL', ire1)
    call jeexin(resgen//'           .VITE', ire2)
    call jeexin(resgen//'           .ACCE', ire3)
!
    if (ire1 .eq. 0 .and. ire2 .eq. 0 .and. ire3 .eq. 0) then
        valk(1) = resgen
        call utmess('F', 'ALGORITH14_35', sk=valk(1))
    end if
!
    call getvtx(' ', 'TOUT_CHAM', nbval=0, nbret=n1)
    if (n1 .ne. 0) then
        call getvtx(' ', 'TOUT_CHAM', scal=k8rep, nbret=n1)
    else
        k8rep = ' '
    end if
    if (k8rep(1:3) .eq. 'OUI') then
        if (ire1 .eq. 0) then
            call utmess('F', 'ALGORITH10_44')
        end if
        if (ire2 .eq. 0) then
            call utmess('F', 'ALGORITH10_45')
        end if
        if (ire3 .eq. 0) then
            call utmess('F', 'ALGORITH10_46')
        end if
        nbcham = 3
        chmp(1) = 'DEPL'
        chmp(2) = 'VITE'
        chmp(3) = 'ACCE'
        call jeveuo(harmge//'.DEPL', 'L', itresu(1))
        call jeveuo(harmge//'.VITE', 'L', itresu(2))
        call jeveuo(harmge//'.ACCE', 'L', itresu(3))
    else
! ----  ON RECHERCHE LES CHAMPS QU'IL FAUT RESTITUER
        call getvtx(' ', 'NOM_CHAM', nbval=0, nbret=n1)
        nbcham = -n1
        call getvtx(' ', 'NOM_CHAM', nbval=nbcham, vect=champ, nbret=n1)
! ----   BOUCLE SUR LES CHAMPS DEMANDES
        do i = 1, nbcham
!
            if (champ(i) .eq. 'DEPL') then
                chmp(i) = 'DEPL'
                call jeexin(harmge//'.DEPL', iret)
                if (iret .eq. 0) then
                    call utmess('F', 'ALGORITH10_11')
                else
                    call jeveuo(harmge//'.DEPL', 'L', itresu(i))
                end if
            else if (champ(i) .eq. 'VITE') then
                chmp(i) = 'VITE'
                call jeexin(harmge//'.VITE', iret)
                if (iret .eq. 0) then
                    call utmess('F', 'ALGORITH10_12')
                else
                    call jeveuo(harmge//'.VITE', 'L', itresu(i))
                end if
            else if (champ(i) .eq. 'ACCE') then
                chmp(i) = 'ACCE'
                call jeexin(harmge//'.ACCE', iret)
                if (iret .eq. 0) then
                    call utmess('F', 'ALGORITH10_13')
                else
                    call jeveuo(harmge//'.ACCE', 'L', itresu(i))
                end if
            else
! ----        SI LE CHAMP N'EST PAS DEPL,VITE OU ACCE ON PLANTE
                call utmess('F', 'ALGORITH10_16')
            end if
        end do
    end if
!
!     NOMBRE DE CHAMPS SYMBOLIQUES CALCULES.
!     ON S'ASSURE QUE LEUR NOMBRE EST NON NUL.
!
    if (nbcham .eq. 0) then
        valk(1) = resgen
        call utmess('F', 'ALGORITH14_35', sk=valk(1))
    end if
!
! --- RECUPERATION DE LA NUMEROTATION ET DU MODELE GENERALISE
!
    call dismoi('NUME_DDL', harmge, 'RESU_DYNA', repk=nume_gene)
    nume_equa_gene = nume_gene//'.NUME'
    call jeveuo(nume_equa_gene//'.REFN', 'L', vk24=refn)
    model_gene = refn(1) (1:8)
    call jelibe(nume_equa_gene//'.REFN')
    call nueq_chck(nume_equa_gene, neqgen)
!
! --- RECUPERATION NUMERO DE LA SOUS-STRUCTURE
!
    call jenonu(jexnom(model_gene//'      .MODG.SSNO', nomsst), nusst)
    if (nusst .eq. 0) then
        valk(1) = model_gene
        valk(2) = nomsst
        call utmess('F', 'ALGORITH14_25', nk=2, valk=valk)
    end if
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
        call jeveuo(jexnum(nume_equa_gene//'.PRNO', i_ligr_ss), 'L', llprs)
        call jelira(jexnum(nume_equa_gene//'.ORIG', i_ligr_ss), 'LONMAX', nbsst)
!
        nutars = 0
        do i = 1, nbsst
            if (zi(llors+i-1) .eq. nusst) nutars = i
        end do
!
! --- NOMBRE DE MODES ET NUMERO DU PREMIER DDL DE LA SOUS-STRUCTURE
!
        nbddg = zi(llprs+(nutars-1)*2+1)
        ieq = zi(llprs+(nutars-1)*2)
!
    else
!
        neqet = 0
        ieq = 0
        call jelira(model_gene//'      .MODG.SSNO', 'NOMMAX', nbsst)
        neqred = neqgen
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
!
        call wkvect('&&MODE_ETENDU_REST_ELIM', 'V V C', neqet, lmoet)
    end if
!
! --- RECUPERATION D'INFORMATIONS SUR LA SOUS-STRUCTURE
!
    call mgutdm(model_gene, nomsst, ibid, 'NOM_BASE_MODALE', ibid, &
                basmod)
    if (elim .ne. 0) then
        call dismoi('NB_MODES_TOT', basmod, 'RESULTAT', repi=nbddg)
    end if
!
! -->AAC-->NORMALEMENT CE .REFD EST INCOHERENT AVEC CELUI DE DYNA_GENE
    call dismoi('REF_INTD_PREM', basmod, 'RESU_DYNA', repk=lint)
    call dismoi('NOM_MAILLA', lint, 'INTERF_DYNA', repk=mailla)
    call dismoi('NOM_NUME_DDL', lint, 'INTERF_DYNA', repk=numddl)
    call dismoi('NB_EQUA', numddl, 'NUME_DDL', repi=neq)
!
    crefe(1) = mailla
    crefe(2) = numddl
!
! --- RECUPERATION DES FREQUENCES
!
    call getvtx(' ', 'CRITERE', scal=crit, nbret=n1)
    call getvr8(' ', 'PRECISION', scal=epsi, nbret=n1)
    call getvtx(' ', 'INTERPOL', scal=interp, nbret=n1)
!
    knume = '&&RETREC.NUM_RANG'
    kfreq = '&&RETREC.FREQ'
    call rstran(interp, harmge, ' ', 1, kfreq, &
                knume, nbfreq, iretou)
    if (iretou .ne. 0) then
        call utmess('F', 'ALGORITH10_47')
    end if
    call jeexin(kfreq, iret)
    if (iret .gt. 0) then
        call jeveuo(kfreq, 'E', jfreq)
        call jeveuo(knume, 'E', jnume)
    end if
!
! --- ALLOCATION DE LA STRUCTURE DE DONNEES RESULTAT-COMPOSE
!
    call rscrsd('G', nomres, 'DYNA_HARMO', nbfreq)
!
! -------------------------------------
! --- RESTITUTION SUR BASE PHYSIQUE ---
! -------------------------------------
!
    call jeveuo(nume_equa_gene//'.NUEQ', 'L', vi=nueq)
!
    iarchi = 0
    if (interp(1:3) .ne. 'NON') then
        call utmess('F', 'ALGORITH3_86')
!
    else
        call jeexin(harmge//'.ORDR', iret)
!
        do i = 0, nbfreq-1
            iarchi = iarchi+1
!
            do ich = 1, nbcham
                idresu = itresu(ich)
!
!-- SI ELIMINATION, ON RESTITUE D'ABORD LES MODES GENERALISES
                if (elim .ne. 0) then
                    do i1 = 1, neqet
                        zc(lmoet+i1-1) = dcmplx(0.d0, 0.d0)
                        do k1 = 1, neqred
                            zc(lmoet+i1-1) = &
                                zc(lmoet+i1-1)+ &
                                zr(lmapro+(k1-1)*neqet+i1-1)* &
                                zc(idresu+k1-1+(zi(jnume+i)-1)*neqred)
                        end do
                    end do
                end if
!
                call rsexch(' ', nomres, chmp(ich), iarchi, chamno, &
                            iret)
                if (iret .eq. 0) then
                    call utmess('A', 'ALGORITH2_64', sk=chamno)
                else if (iret .eq. 100) then
                    call vtcrea(chamno, crefe, 'G', 'C', neq)
                else
                    ASSERT(.false.)
                end if
                chamno(20:24) = '.VALE'
                call jeveuo(chamno, 'E', ldnew)
!
! --- BOUCLE SUR LES MODES PROPRES DE LA BASE
!
                do j = 1, nbddg
                    call dcapno(basmod, 'DEPL', j, chamba)
                    call jeveuo(chamba, 'L', llchab)
!
                    if (elim .ne. 0) then
                        iad = lmoet+ieq+j-1
                    else
                        iad = idresu+(zi(jnume+i)-1)*neqgen+nueq(1+ &
                                                                 ieq+j-2)-1
                    end if
!
! --- BOUCLE SUR LES EQUATIONS PHYSIQUES
!
                    do k = 1, neq
                        zc(ldnew+k-1) = zc(ldnew+k-1)+zr(llchab+k-1)*zc(iad)
                    end do
                    call jelibe(chamba)
                end do
                call jelibe(chamno)
                call rsnoch(nomres, chmp(ich), iarchi)
            end do
            call rsadpa(nomres, 'E', 1, 'FREQ', iarchi, &
                        0, sjv=lfreq, styp=k8b)
            zr(lfreq) = zr(jfreq+i)
        end do
!
    end if
!
! --> AAC-->NORMALEMENT CE .REFD EST INCOHERENT AVEC CELUI DE DYNA_GENE
    call refdcp(basmod, nomres)
    call jelibe(nume_equa_gene//'.NUEQ')
    call jedetr('&&RETREC.NUM_RANG')
    call jedetr('&&RETREC.FREQ')
!
    goto 999
!
999 continue
    call jedema()
end subroutine
