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

subroutine regegl(nomres, resgen, mailsk, profno)
    implicit none
#include "jeveux.h"
#include "asterfort/dcapno.h"
#include "asterfort/dismoi.h"
#include "asterfort/genugl.h"
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
#include "asterfort/matrot.h"
#include "asterfort/mgutdm.h"
#include "asterfort/nueq_chck.h"
#include "asterfort/rotchm.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rscrsd.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnoch.h"
#include "asterfort/rsorac.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"
#include "asterfort/vtcrea.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
!
    character(len=8) :: nomres, resgen, mailsk
    character(len=19) :: profno
!
!  BUT : < RESTITUTION GENERALISEE GLOBALE >
!
!  RESTITUER EN BASE PHYSIQUE SUR UN MAILLAGE SQUELETTE LES RESULTATS
!  ISSUS DE LA SOUS-STRUCTURATION GENERALE.
!  LE CONCEPT RESULTAT EST UN RESULTAT COMPOSE "MODE_MECA"
!
!
!  /!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
!
!--   LES ROUTINES REGEGL ET REGEGC FONT LA MEME CHOSE, UNE EN REEL,
!--   L'AUTRE EN COMPLEXE. EN CAS DE MODIFICATION D'UNE DES ROUTINES,
!--   NE PAS OUBLIER DE REPORTER LE CHANGEMENT DANS L'AUTRE.
!
!  /!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
!
!-----------------------------------------------------------------------
!
! NOMRES  /I/ : NOM K8 DU CONCEPT MODE_MECA RESULTAT
! RESGEN  /I/ : NOM K8 DU RESULTAT GENERALISE AMONT
! MAILSK  /I/ : NOM K8 DU MAILLAGE SQUELETTE SUPPORT
! PROFNO  /I/ : NOM K19 DU NUME_EQUA DU SQUELETTE
!
!
!
!
!
    integer(kind=8) :: i, iad, iar, ibid, idep, ieq, ier, iord, iret, j, jbid, k, l
    integer(kind=8) ::  llchab, llchol, llind, llmass, i1, k1
    integer(kind=8) ::  llors, llprs, llrot
    integer(kind=8) ::  nbbas, nbcmp, nbcou, nbmas, nbmax, nbmod(1), nbnot, nbsst
    integer(kind=8) :: neq, neqs, nno, numo, nutars, elim, lmoet
    integer(kind=8) :: neqet, lmapro, neqred, lsilia, numsst, lsst, lrevind
    integer(kind=8) :: iadpar(12)
    integer(kind=8) :: vali(2), i_ligr_ss
    real(kind=8) :: compx, compy, compz, efmasx, efmasy, efmasz, freq, genek, fpartx, fparty
    real(kind=8) :: fpartz
    real(kind=8) :: genem, mat(3, 3), omeg2, rbid
    character(len=8) :: basmod, macrel, model_gene, kbid
    character(len=16) :: depl, nompar(12)
    character(len=14) :: nume_gene
    character(len=19) :: raid, numddl, chamne, nume_equa_gene
    character(len=24) :: crefe(2), chamol, chamba, indirf, seliai, sizlia, sst
    character(len=24) :: valk, nomsst, intf
    complex(kind=8) :: cbid
    real(kind=8), pointer :: rotx(:) => null()
    real(kind=8), pointer :: roty(:) => null()
    real(kind=8), pointer :: rotz(:) => null()
    real(kind=8), pointer :: trav(:) => null()
    integer(kind=8), pointer :: skeleton(:) => null()
    character(len=24), pointer :: refa(:) => null()
    real(kind=8), pointer :: vale(:) => null()
    character(len=8), pointer :: idc_type(:) => null()
    integer(kind=8), pointer :: nueq(:) => null()
    real(kind=8), pointer :: mael_iner_vale(:) => null()
    character(len=24), pointer :: refn(:) => null()
!
!-----------------------------------------------------------------------
    data depl/'DEPL            '/
    data nompar/'FREQ', 'RIGI_GENE', 'MASS_GENE', 'OMEGA2', 'NUME_MODE',&
     &            'MASS_EFFE_DX', 'MASS_EFFE_DY', 'MASS_EFFE_DZ',&
     &            'FACT_PARTICI_DX', 'FACT_PARTICI_DY', 'FACT_PARTICI_DZ',&
     &            'TYPE_MODE'/
!-----------------------------------------------------------------------
!
    call jemarq()
!
    indirf = '&&REGEGL.INDIR.SST'
!
!-----ECRITURE DU TITRE-------------------------------------------------
!
    call titre()
!
!-----VERIF SQUELETTE---------------------------------------------------
!
    call jeexin(mailsk//'.INV.SKELETON', iret)
    if (iret .eq. 0) then
        valk = mailsk
        call utmess('F', 'ALGORITH14_27', sk=valk)
    end if
    call jeveuo(mailsk//'.INV.SKELETON', 'L', vi=skeleton)
!
!-----RECUPERATION DU MODELE GENERALISE--------------------------------
!
    call dismoi('REF_RIGI_PREM', resgen, 'RESU_DYNA', repk=raid)
!
    call jeveuo(raid//'.REFA', 'L', vk24=refa)
    nume_gene = refa(2) (1:14)
    nume_equa_gene = nume_gene//'.NUME'
    call nueq_chck(nume_equa_gene, neqred)
    call jelibe(raid//'.REFA')
!
    call jeveuo(nume_equa_gene//'.REFN', 'L', vk24=refn)
    model_gene = refn(1) (1:8)
    call jelibe(nume_equa_gene//'.REFN')
!
    call jelira(model_gene//'      .MODG.SSNO', 'NOMMAX', nbsst)
    kbid = '  '
    call mgutdm(model_gene, kbid, 1, 'NB_CMP_MAX', nbcmp, &
                kbid)
!
!-----RECUPERATION DES ROTATIONS----------------------------------------
!
    AS_ALLOCATE(vr=rotx, size=nbsst)
    AS_ALLOCATE(vr=roty, size=nbsst)
    AS_ALLOCATE(vr=rotz, size=nbsst)
    do i = 1, nbsst
        call jeveuo(jexnum(model_gene//'      .MODG.SSOR', i), 'L', llrot)
        rotz(i) = zr(llrot)
        roty(i) = zr(llrot+1)
        rotx(i) = zr(llrot+2)
    end do
!
!-----CREATION DU PROF-CHAMNO-------------------------------------------
!
    call genugl(profno, indirf, model_gene, mailsk)
    call nueq_chck(profno, neq)
!
!-----RECUPERATION DU NOMBRE DE NOEUDS----------------------------------
!
    call dismoi('NB_NO_MAILLA', mailsk, 'MAILLAGE', repi=nbnot)
!
!-----RECUPERATION DE LA BASE MODALE------------------------------------
!
    crefe(1) = mailsk
    crefe(2) = profno
!
!-----RECUPERATION NOMBRE DE MODES PROPRES CALCULES---------------------
!
    call rsorac(resgen, 'LONUTI', 0, rbid, kbid, &
                cbid, rbid, kbid, nbmod, 1, &
                ibid)
!
!
! --- ON RESTITUE SUR TOUS LES MODES OU SUR QUELQUES MODES:
!
    call getvis(' ', 'NUME_ORDRE', nbval=0, nbret=nno)
    if (nno .ne. 0) then
        nbmod(1) = -nno
        call wkvect('&&REGEGL.NUME', 'V V I', nbmod(1), jbid)
        call getvis(' ', 'NUME_ORDRE', nbval=nbmod(1), vect=zi(jbid), nbret=nno)
    else
        call wkvect('&&REGEGL.NUME', 'V V I', nbmod(1), jbid)
        do i = 1, nbmod(1)
            zi(jbid+i-1) = i
        end do
    end if
!
!-----ALLOCATION STRUCTURE DE DONNEES RESULTAT--------------------------
!
    call rscrsd('G', nomres, 'MODE_MECA', nbmod(1))
!
!-- ON TESTE SI ON A EU RECOURS A L'ELIMINATION
!
    seliai = nume_gene(1:14)//'.ELIM.BASE'
    sizlia = nume_gene(1:14)//'.ELIM.TAIL'
    sst = nume_gene(1:14)//'.ELIM.NOMS'
!
    call jeexin(seliai, elim)
    if (elim .ne. 0) then
        neqet = 0
        nomsst = model_gene//'      .MODG.SSNO'
        call jeveuo(seliai, 'L', lmapro)
        call jeveuo(sizlia, 'L', lsilia)
        call jeveuo(sst, 'L', lsst)
        do i = 1, nbsst
            neqet = neqet+zi(lsilia+i-1)
        end do
        call wkvect('&&MODE_ETENDU_REST_ELIM', 'V V R', neqet, lmoet)
        !-- Recherche de la position dans MODGEN des SST de NUMGEN
        call wkvect('&&REGEGL.REVERSE_INDEX', 'V V I', nbsst, lrevind)
        do k = 1, nbsst
            call jenonu(jexnom(nomsst, zk8(lsst+k-1)), numsst)
            zi(lrevind+numsst-1) = k
        end do
    end if
!
!
!---RESTITUTION PROPREMENT DITE---------------------------------------
!
!
    call jeveuo(nume_equa_gene//'.NUEQ', 'L', vi=nueq)
    call jenonu(jexnom(nume_equa_gene//'.LILI', '&SOUSSTR'), i_ligr_ss)
    call jeveuo(jexnum(nume_equa_gene//'.ORIG', i_ligr_ss), 'L', llors)
    call jeveuo(jexnum(nume_equa_gene//'.PRNO', i_ligr_ss), 'L', llprs)
!
!  BOUCLE SUR LES MODES A RESTITUER
!
    do i = 1, nbmod(1)
        iord = zi(jbid+i-1)
!
!  REQUETTE NOM ET ADRESSE CHAMNO GENERALISE
!
        call dcapno(resgen, depl, iord, chamol)
        call jeveuo(chamol, 'L', llchol)
!
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
!  REQUETTE NOM ET ADRESSE NOUVEAU CHAMNO
!
        call rsexch(' ', nomres, depl, i, chamne, &
                    ier)
        call vtcrea(chamne, crefe, 'G', 'R', neq)
        call jeveuo(chamne//'.VALE', 'E', vr=vale)
!
        call rsadpa(resgen, 'L', 8, nompar, iord, &
                    0, tjv=iadpar, styp=kbid)
        freq = zr(iadpar(1))
        genek = zr(iadpar(2))
        genem = zr(iadpar(3))
        omeg2 = zr(iadpar(4))
        numo = zi(iadpar(5))
        efmasx = 0.d0
        efmasy = 0.d0
        efmasz = 0.d0
!
!  BOUCLE SUR LES SOUS-STRUCTURES
!
        do k = 1, nbsst
            call jeexin(jexnum(indirf, k), iret)
!
!
!  TEST SI LA SST GENERE DES DDL GLOBAUX
!
            if (iret .ne. 0) then
                kbid = '  '
                if (elim .ne. 0) then
                    numsst = zi(lrevind+k-1)
                    ieq = 0
                    do i1 = 1, numsst-1
                        ieq = ieq+zi(lsilia+i1-1)
                    end do
                else
                    numsst = k
!  RECUPERATION DU NUMERO TARDIF DE LA SST
                    do j = 1, nbsst
                        if (zi(llors+j-1) .eq. numsst) nutars = j
                    end do
                    ieq = zi(llprs+(nutars-1)*2)
                end if
!
                call mgutdm(model_gene, kbid, k, 'NOM_BASE_MODALE', ibid, &
                            basmod)
!
                call dismoi('REF_INTD_PREM', basmod, 'RESU_DYNA', repk=intf)
                call jeveuo(intf(1:8)//'.IDC_TYPE', 'L', vk8=idc_type)
                if (idc_type(1) .eq. 'AUCUN') then
                    vali(1) = k
                    vali(2) = k
                    call utmess('A', 'ALGORITH14_28', ni=2, vali=vali)
                end if
                call mgutdm(model_gene, kbid, k, 'NOM_MACR_ELEM', ibid, &
                            macrel)
                call jeveuo(jexnum(macrel//'.MAEL_MASS_VALE', 1), 'L', llmass)
                call jelira(macrel//'.MAEL_MASS_VALE', 'LONMAX', nbmax)
                call jelira(macrel//'.MAEL_MASS_VALE', 'LONUTI', nbmas)
                call jeveuo(macrel//'.MAEL_INER_VALE', 'L', vr=mael_iner_vale)
!
!           --- CALCUL DE LA MATRICE DE ROTAION ---
                call jeveuo(jexnum(model_gene//'      .MODG.SSOR', k), 'L', llrot)
                call matrot(zr(llrot), mat)
!
                call dismoi('NB_MODES_TOT', basmod, 'RESULTAT', repi=nbbas)
                kbid = '  '
                call mgutdm(model_gene, kbid, k, 'NOM_NUME_DDL', ibid, &
                            numddl)
                call dismoi('NB_EQUA', numddl, 'NUME_DDL', repi=neqs)
                AS_ALLOCATE(vr=trav, size=neqs)
!
!  BOUCLE SUR LES MODES PROPRES DE LA BASE
!
                do j = 1, nbbas
                    call dcapno(basmod, depl, j, chamba)
                    call jeveuo(chamba, 'L', llchab)
!
!  BOUCLE SUR LES EQUATIONS PHYSIQUES
!
                    if (elim .ne. 0) then
                        iad = llchol+ieq+j-1
                    else
                        iad = llchol+nueq(1+ieq+j-2)-1
                    end if
!
!-- DANS LE CAS ELIM, CHANGER LE IAD, le ZI(LLNUEQ EST PAS BON)
!
!           --- CALCUL DES MASSES EFFECTIVES ---
                    compx = mael_iner_vale(j)
                    compy = mael_iner_vale(1+nbbas+j-1)
                    compz = mael_iner_vale(1+2*nbbas+j-1)
!             --- UTILISATION DE MAT TRANSPOSEE (TRANSFORMATION INVERSE)
                    efmasx = efmasx+zr(iad)*(compx*mat(1, 1)+compy*mat(2, 1)+compz*mat(3, 1))
                    efmasy = efmasy+zr(iad)*(compx*mat(1, 2)+compy*mat(2, 2)+compz*mat(3, 2))
                    efmasz = efmasz+zr(iad)*(compx*mat(1, 3)+compy*mat(2, 3)+compz*mat(3, 3))
!
!  BOUCLE SUR LES DDL DE LA BASE
!
                    do l = 1, neqs
                        trav(l) = trav(l)+zr(llchab+l-1)* &
                                  zr(iad)
                    end do
!
                    call jelibe(chamba)
                end do
                call jeveuo(jexnum(indirf, k), 'L', llind)
                call jelira(jexnum(indirf, k), 'LONMAX', nbcou)
                nbcou = nbcou/2
                do l = 1, nbcou
                    idep = zi(llind+(l-1)*2)
                    iar = zi(llind+(l-1)*2+1)
                    vale(iar) = trav(idep)
                end do
                call jelibe(jexnum(indirf, k))
                AS_DEALLOCATE(vr=trav)
            end if
        end do
!
        fpartx = efmasx/genem
        fparty = efmasy/genem
        fpartz = efmasz/genem
        efmasx = efmasx*efmasx/genem
        efmasy = efmasy*efmasy/genem
        efmasz = efmasz*efmasz/genem
        call rsnoch(nomres, depl, i)
        call rsadpa(nomres, 'E', 12, nompar, i, &
                    0, tjv=iadpar, styp=kbid)
        zr(iadpar(1)) = freq
        zr(iadpar(2)) = genek
        zr(iadpar(3)) = genem
        zr(iadpar(4)) = omeg2
        zi(iadpar(5)) = numo
        zr(iadpar(6)) = efmasx
        zr(iadpar(7)) = efmasy
        zr(iadpar(8)) = efmasz
        zr(iadpar(9)) = fpartx
        zr(iadpar(10)) = fparty
        zr(iadpar(11)) = fpartz
        zk16(iadpar(12)) = 'MODE_DYN'
!
        call jelibe(chamol)
!
!  ROTATION DU CHAMPS AUX NOEUDS
!
        call rotchm(profno, vale, rotx, nbsst, skeleton, &
                    nbnot, nbcmp, 1)
        call rotchm(profno, vale, roty, nbsst, skeleton, &
                    nbnot, nbcmp, 2)
        call rotchm(profno, vale, rotz, nbsst, skeleton, &
                    nbnot, nbcmp, 3)

    end do
!
! --- MENAGE
    AS_DEALLOCATE(vr=rotx)
    AS_DEALLOCATE(vr=roty)
    AS_DEALLOCATE(vr=rotz)
!
    call jedetr('&&MODE_ETENDU_REST_ELIM')
    call jedetr('&&REGEGL.REVERSE_INDEX')
    call jedetr(indirf)
    call jelibe(nume_equa_gene//'.NUEQ')
    call jedetr('&&REGEGL.NUME')
!
    call jedema()
end subroutine
