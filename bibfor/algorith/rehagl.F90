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

subroutine rehagl(nomres, resgen, mailsk, profno)
!    T. KERBER     DATE 16/05/93
!-----------------------------------------------------------------------
!  BUT:      < RESTITUTION HARMONIQUE GLOBALE >
    implicit none
!
!      RESTITUER EN BASE PHYSIQUE SUR UN MAILLAGE SQUELETTE
!  LES RESULTATS ISSUS DE LA SOUS-STRUCTURATION GENERALE
!  LE CONCEPT RESULTAT EST UN RESULTAT COMPOSE "DYNA_HARMO"
!
!-----------------------------------------------------------------------
!
! NOMRES   /I/: NOM K8 DU CONCEPT DYNA_HARMO RESULTAT
! RESGEN   /I/: NOM K8 DU RESULTAT GENERALISE AMONT
! MAILSK   /I/: NOM K8 DU MAILLAGE SQUELETTE SUPPORT
! PROFNO   /I/: NOM K19 DU NUME_EQUA DU SQUELETTE
!
!
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dcapno.h"
#include "asterfort/dismoi.h"
#include "asterfort/genugl.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/mgutdm.h"
#include "asterfort/nueq_chck.h"
#include "asterfort/refdcp.h"
#include "asterfort/rotchc.h"
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
    character(len=6) :: pgc
    character(len=8) :: chmp(3), crit, interp, k8b, nomres, basmod, mailsk
    character(len=8) :: model_gene, resgen, k8bid, k8rep
    character(len=14) :: nume_gene
    character(len=19) :: numddl, nume_equa_gene, knume, kfreq, harmge, profno
    character(len=24) :: crefe(2), chamba, indirf, chamno, seliai, sizlia, sst
    character(len=24) :: valk, nomsst
    integer(kind=8) :: itresu(3), elim, neqet, neqred, lmapro, lsilia, lsst, lmoet
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, i1, iad, iar, iarchi, ibid, ich, i_ligr_ss
    integer(kind=8) :: idep, idresu, ieq, ire1
    integer(kind=8) :: ire2, ire3, iret, iretou, j, jfreq, jnume
    integer(kind=8) :: k, k1, l, ldnew, lfreq, llchab, llind
    integer(kind=8) :: llinsk, llors, llprs
    integer(kind=8) :: llrot, ltrotx, ltroty, ltrotz, ltvec, n1
    integer(kind=8) :: nbbas, nbcham, nbcmp, nbcou, nbfreq, nbnot
    integer(kind=8) :: nbsst, neq, neqgen, neqs, numsst, nutars
    integer(kind=8), pointer :: nueq(:) => null()
    character(len=24), pointer :: refn(:) => null()
!-----------------------------------------------------------------------
    data pgc/'REHAGL'/
!-----------------------------------------------------------------------
!
    call jemarq()
!
    indirf = '&&'//pgc//'.INDIR.SST'
!
! --- ECRITURE DU TITRE
    call titre()
!
! --- VERIFICATION SQUELETTE
    call jeexin(mailsk//'.INV.SKELETON', iret)
    if (iret .eq. 0) then
        valk = mailsk
        call utmess('F', 'ALGORITH14_27', sk=valk)
    end if
    call jeveuo(mailsk//'.INV.SKELETON', 'L', llinsk)
!
! --- DETERMINATION DES CHAMPS A RESTITUER, PARMI DEPL, VITE ET ACCE
!
    harmge = resgen
!
! --- CALCUL DU NOMBRE DE CHAMPS A RESTITUER ET LEURS ADDRESSES
!
    call jeexin(resgen//'           .DEPL', ire1)
    call jeexin(resgen//'           .VITE', ire2)
    call jeexin(resgen//'           .ACCE', ire3)
!
    if (ire1 .eq. 0 .and. ire2 .eq. 0 .and. ire3 .eq. 0) then
        valk = resgen
        call utmess('F', 'ALGORITH14_35', sk=valk)
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
!     NOMBRE DE CHAMPS SYMBOLIQUES CALCULES.
!     ON S'ASSURE QUE LEUR NOMBRE EST NON NUL.
!
    if (nbcham .eq. 0) then
        valk = resgen
        call utmess('F', 'ALGORITH14_35', sk=valk)
    end if
!
! --- RECUPERATION DE LA NUMEROTATION ET DU MODELE GENERALISE
    call dismoi('NUME_DDL', harmge, 'RESU_DYNA', repk=nume_gene)
    nume_equa_gene = nume_gene//'.NUME'
    call jeveuo(nume_equa_gene//'.REFN', 'L', vk24=refn)
    model_gene = refn(1) (1:8)
    call jelira(model_gene//'      .MODG.SSNO', 'NOMMAX', nbsst)
    k8bid = '  '
    call mgutdm(model_gene, k8bid, 1, 'NB_CMP_MAX', nbcmp, &
                k8bid)
    call nueq_chck(nume_equa_gene, neqgen)
!
! --- RECUPERATION DES ROTATIONS
    call wkvect('&&'//pgc//'ROTX', 'V V R', nbsst, ltrotx)
    call wkvect('&&'//pgc//'ROTY', 'V V R', nbsst, ltroty)
    call wkvect('&&'//pgc//'ROTZ', 'V V R', nbsst, ltrotz)
    do i = 1, nbsst
        call jeveuo(jexnum(model_gene//'      .MODG.SSOR', i), 'L', llrot)
        zr(ltrotz+i-1) = zr(llrot)
        zr(ltroty+i-1) = zr(llrot+1)
        zr(ltrotx+i-1) = zr(llrot+2)
    end do
!
! --- CREATION DU PROF-CHAMNO
    call genugl(profno, indirf, model_gene, mailsk)
    call nueq_chck(profno, neq)
!
! --- RECUPERATION DU NOMBRE DE NOEUDS
    call dismoi('NB_NO_MAILLA', mailsk, 'MAILLAGE', repi=nbnot)
!
! --- INFORMATIONS POUR CREATION DES CHAMNO A PARTIR DES .REFE
    crefe(2) = profno
!
! --- RECUPERATION DES FREQUENCES
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
!-- ON TESTE SI ON A EU RECOURS A L'ELIMINATION
!
    seliai = nume_gene(1:14)//'.ELIM.BASE'
    sizlia = nume_gene(1:14)//'.ELIM.TAIL'
    sst = nume_gene(1:14)//'.ELIM.NOMS'
!
    call jeexin(seliai, elim)
    if (elim .ne. 0) then
        neqet = 0
        neqred = neqgen
        nomsst = model_gene//'      .MODG.SSNO'
        call jeveuo(seliai, 'L', lmapro)
        call jeveuo(sizlia, 'L', lsilia)
        call jeveuo(sst, 'L', lsst)
        do i = 1, nbsst
            neqet = neqet+zi(lsilia+i-1)
        end do
        call wkvect('&&MODE_ETENDU_REST_ELIM', 'V V C', neqet, lmoet)
    end if
!
! -------------------------------------
! --- RESTITUTION SUR BASE PHYSIQUE ---
! -------------------------------------
!
    call jeveuo(nume_equa_gene//'.NUEQ', 'L', vi=nueq)
    call jenonu(jexnom(nume_equa_gene//'.LILI', '&SOUSSTR'), i_ligr_ss)
    call jeveuo(jexnum(nume_equa_gene//'.ORIG', i_ligr_ss), 'L', llors)
    call jeveuo(jexnum(nume_equa_gene//'.PRNO', i_ligr_ss), 'L', llprs)
!
    iarchi = 0
!
    if (interp(1:3) .ne. 'NON') then
        call utmess('F', 'ALGORITH3_86')
    else
!
!        CALL JEEXIN(HARMGE//'.ORDR',IRET)
!        IF (IRET.NE.0.AND.ZI(JNUME).EQ.1) IARCHI=0
!
        do i = 0, nbfreq-1
            iarchi = iarchi+1
!
            do ich = 1, nbcham
                idresu = itresu(ich)
!-- SI ELIMINATION, ON RESTITUE D'ABORD LES MODES GENERALISES
                if (elim .ne. 0) then
                    do i1 = 1, neqet
                        zc(lmoet+i1-1) = dcmplx(0.d0, 0.d0)
                        do k1 = 1, neqred
                            zc(lmoet+i1-1) = zc(lmoet+i1-1)+ &
                                             zr(lmapro+(k1-1)*neqet+i1-1)* &
                                             zc(idresu+k1-1+(zi(jnume+i)-1)*neqred)
                        end do
                    end do
                end if
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
! --- BOUCLE SUR LES SOUS-STRUCTURES
!
                do k = 1, nbsst
                    call jeexin(jexnum(indirf, k), iret)
!
! --- TEST SI LA SST GENERE DES DDL GLOBAUX
!
                    if (iret .ne. 0) then
!
! --- RECUPERATION DU NUMERO TARDIF DE LA SST
!
                        if (elim .ne. 0) then
                            call jenonu(jexnom(nomsst, zk8(lsst+k-1)), numsst)
                            ieq = 0
                            do i1 = 1, k-1
                                ieq = ieq+zi(lsilia+i1-1)
                            end do
                        else
                            numsst = k
!  RECUPERATION DU NUMERO TARDIF DE LA SST
                            nutars = 0
                            do j = 1, nbsst
                                if (zi(llors+j-1) .eq. numsst) nutars = j
                            end do
                            ieq = zi(llprs+(nutars-1)*2)
                        end if
                        k8bid = '  '
                        call mgutdm(model_gene, k8bid, numsst, 'NOM_BASE_MODALE', ibid, &
                                    basmod)
                        call dismoi('NB_MODES_TOT', basmod, 'RESULTAT', repi=nbbas)
                        k8bid = '  '
                        call mgutdm(model_gene, k8bid, numsst, 'NOM_NUME_DDL', ibid, &
                                    numddl)
                        call dismoi('NB_EQUA', numddl, 'NUME_DDL', repi=neqs)
                        call wkvect('&&'//pgc//'.TRAV', 'V V C', neqs, ltvec)
!
! --- BOUCLE SUR LES MODES PROPRES DE LA BASE
!
                        do j = 1, nbbas
                            call dcapno(basmod, 'DEPL', j, chamba)
                            call jeveuo(chamba, 'L', llchab)
!
                            if (elim .ne. 0) then
                                iad = lmoet+ieq+j-1
                            else
                                iad = idresu+(zi(jnume+i)-1)*neqgen+ &
                                      nueq(1+ieq+j-2)-1
                            end if
!
! --- BOUCLE SUR LES EQUATIONS PHYSIQUES
!
                            do l = 1, neqs
                                zc(ltvec+l-1) = zc(ltvec+l-1)+zr(llchab+ &
                                                                 l-1)*zc(iad)
                            end do
                        end do
                        call jeveuo(jexnum(indirf, numsst), 'L', llind)
                        call jelira(jexnum(indirf, numsst), 'LONMAX', nbcou)
                        nbcou = nbcou/2
                        do l = 1, nbcou
                            idep = zi(llind+(l-1)*2)
                            iar = zi(llind+(l-1)*2+1)
                            zc(ldnew+iar-1) = zc(ltvec+idep-1)
                        end do
                        call jedetr('&&'//pgc//'.TRAV')
                    end if
!
                end do
                call rsnoch(nomres, chmp(ich), iarchi)
!
! --- ROTATION DU CHAMP AUX NOEUDS
!
                call rotchc(profno, zc(ldnew), zr(ltrotz), nbsst, zi(llinsk), &
                            nbnot, nbcmp, 3)
                call rotchc(profno, zc(ldnew), zr(ltroty), nbsst, zi(llinsk), &
                            nbnot, nbcmp, 2)
                call rotchc(profno, zc(ldnew), zr(ltrotx), nbsst, zi(llinsk), &
                            nbnot, nbcmp, 1)
!
            end do
            call rsadpa(nomres, 'E', 1, 'FREQ', iarchi, &
                        0, sjv=lfreq, styp=k8b)
            zr(lfreq) = zr(jfreq+i)
        end do
!
    end if
!
    call refdcp(basmod, nomres)
!
! --- MENAGE
    call jedetr('&&'//pgc//'ROTX')
    call jedetr('&&'//pgc//'ROTY')
    call jedetr('&&'//pgc//'ROTZ')
    call jedetr('&&'//pgc//'.INDIR.SST')
!
    call jedema()
end subroutine
