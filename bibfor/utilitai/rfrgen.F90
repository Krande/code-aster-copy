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

subroutine rfrgen(trange)

    use DynaGene_module

    implicit none
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/gettco.h"
#include "asterfort/assert.h"
#include "asterfort/copmod.h"
#include "asterfort/dismoi.h"
#include "asterfort/extrac.h"
#include "asterfort/foattr.h"
#include "asterfort/foimpr.h"
#include "asterfort/fointe.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lxlgut.h"
#include "asterfort/mdgep2.h"
#include "asterfort/mdgep4.h"
#include "asterfort/ordonn.h"
#include "asterfort/posddl.h"
#include "asterfort/rfhge2.h"
#include "asterfort/rfmge1.h"
#include "asterfort/rstran.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"
#include "asterfort/utnono.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
!
    character(len=19) :: trange
!
!     OPERATEUR "RECU_FONCTION"  MOT CLE "RESU_GENE"
!     ------------------------------------------------------------------
    integer(kind=8) :: ibid
    integer(kind=8) :: ifm, niv
    character(len=24) :: valk(2), nogno
    character(len=4) :: interp(2), intres
    character(len=8) :: k8b, crit, noeud, cmp, noma, nomacc, basemo
    character(len=8) :: monmot(2), nonmot, nomno
    character(len=14) :: nume
    character(len=16) :: nomcmd, typcon, nomcha, tysd
    character(len=19) :: nomfon, knume, kinst, nomres, fonct
!     ------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, idbase, iddl, ie, shift
    integer(kind=8) :: ier, ierd, ii, inoeud, iordr
    integer(kind=8) :: iret, jfon, jinst
    integer(kind=8) ::  lfon, lg1, lg2, lordr, lpro
    integer(kind=8) :: lvar, n1, n2, n3, nbexci, nbinsg
    integer(kind=8) :: nbmode, nbordr, neq
    integer(kind=8) :: nfonct, ngn, numcmp, i_cham, i_bloc
    real(kind=8) :: alpha, epsi, rep, rep1(1)
    complex(kind=8) :: cbid
    real(kind=8), pointer :: vectgene(:) => null()
    integer(kind=8), pointer :: desc(:) => null()
    real(kind=8), pointer :: ipsd(:) => null()
    real(kind=8), pointer :: disc(:) => null()
    real(kind=8), pointer :: resu(:) => null()
    type(DynaGene) :: dyna_gene
    cbid = dcmplx(0.d0, 0.d0)
!-----------------------------------------------------------------------
    call jemarq()
!
    call infmaj()
    call infniv(ifm, niv)
!
    call getres(nomfon, typcon, nomcmd)
!
    call gettco(trange, tysd)
! TRAITEMENT DU MODE_GENE
    if (tysd .eq. 'MODE_GENE') then
        call getvtx(' ', 'NOM_PARA_RESU', scal=k8b, nbret=n1)
        call getvis(' ', 'NUME_CMP_GENE', scal=ibid, nbret=n2)
        if ((n1+n2) .ne. 0) then
            call rfmge1(trange)
        else
!CC  FONCTIONNALITE NON DEVELOPPEE
            ASSERT(.false.)
        end if
        goto 999
! TRAITEMENT DU HARM_GENE
    else if (tysd .eq. 'HARM_GENE') then
        call rfhge2(trange)
        goto 999
    end if
! TRAITEMENT DU TRAN_GENE
    nomres = trange
    interp(1) = 'LIN '
    interp(2) = 'LIN '
    intres = 'NON '
!
    call getvtx(' ', 'CRITERE', scal=crit, nbret=n1)
    call getvr8(' ', 'PRECISION', scal=epsi, nbret=n1)
    call getvtx(' ', 'INTERP_NUME', scal=intres, nbret=n1)
    call getvtx(' ', 'INTERPOL', nbval=2, vect=interp, nbret=n1)
    if (n1 .eq. 1) interp(2) = interp(1)
!
    noeud = ' '
    cmp = ' '
    call getvtx(' ', 'NOEUD', scal=noeud, nbret=n1)
    call getvtx(' ', 'NOM_CMP', scal=cmp, nbret=n2)
    call getvtx(' ', 'NOM_CHAM', scal=nomcha, nbret=n3)
!
    nomacc = 'INST'
    knume = '&&RFRGEN.NUME_ORDR'
    kinst = '&&RFRGEN.INSTANT'
    call rstran(intres, nomres, ' ', 1, kinst, &
                knume, nbordr, ie)
    if (ie .ne. 0) then
        call utmess('F', 'UTILITAI4_24')
    end if
    call jeexin(kinst, iret)
    if (iret .gt. 0) then
        call jeveuo(kinst, 'L', jinst)
        call jeveuo(knume, 'L', lordr)
    end if
!
!     --- REMPLISSAGE DU .PROL ---
!
    ASSERT(lxlgut(nomfon) .le. 24)
    call wkvect(nomfon//'.PROL', 'G V K24', 6, lpro)
    zk24(lpro) = 'FONCTION'
    zk24(lpro+1) = interp(1)//interp(2)
    zk24(lpro+2) = nomacc(1:8)
    zk24(lpro+3) = nomcha(1:4)
    zk24(lpro+4) = 'EE      '
    zk24(lpro+5) = nomfon

    call dyna_gene%init(trange(1:8))
!
!----------------------------------------------------------------------
!                            P T E M
!----------------------------------------------------------------------
!
    if (nomcha(1:4) .eq. 'PTEM') then
!
        call wkvect(nomfon//'.VALE', 'G V R', 2*nbordr, lvar)
        lfon = lvar+nbordr
        if (intres(1:3) .ne. 'NON') then
            do iordr = 0, nbordr-1
                call dyna_gene%get_values_by_disc(dyna_gene%ptem, zr(jinst+iordr), length=nbinsg, &
                                                  vr=resu)
                call dyna_gene%get_current_bloc(dyna_gene%ptem, i_bloc)
                call dyna_gene%get_values(dyna_gene%disc, i_bloc, vr=disc)
                call extrac(intres, epsi, crit, nbinsg-2, disc, &
                            zr(jinst+iordr), resu, 1, rep1, ierd)
                zr(lvar+iordr) = zr(jinst+iordr)
                zr(lfon+iordr) = rep1(1)
            end do
        else
            do iordr = 0, nbordr-1
                ii = zi(lordr+iordr)
                call dyna_gene%get_values_by_index(dyna_gene%ptem, ii, shift, vr=resu)
                zr(lvar+iordr) = zr(jinst+iordr)
                zr(lfon+iordr) = resu(ii-shift)
            end do
        end if
!
!----------------------------------------------------------------------
!                 D E P L   ---   V I T E   ---   A C C E
!----------------------------------------------------------------------
!
    else
        if (nomcha(1:4) .eq. 'DEPL') then
            i_cham = dyna_gene%depl
        else if (nomcha(1:4) .eq. 'VITE') then
            i_cham = dyna_gene%vite
        else if (nomcha(1:4) .eq. 'ACCE') then
            i_cham = dyna_gene%acce
        else
            call utmess('F', 'UTILITAI4_23', sk=nomcha)
        end if
!
        call jeveuo(nomres//'.DESC', 'L', vi=desc)
        nbmode = desc(2)
        call getvis(' ', 'NUME_CMP_GENE', scal=numcmp, nbret=n1)
!
        if (n1 .ne. 0) then
            if (numcmp .gt. nbmode) then
                call utmess('F', 'UTILITAI4_14')
            end if
            call wkvect(nomfon//'.VALE', 'G V R', 2*nbordr, lvar)
            lfon = lvar+nbordr
            if (intres(1:3) .ne. 'NON') then
                AS_ALLOCATE(vr=vectgene, size=nbmode)
                do iordr = 0, nbordr-1
                    call dyna_gene%get_values_by_disc(i_cham, zr(jinst+iordr), length=nbinsg, &
                                                      vr=resu)
                    call dyna_gene%get_current_bloc(i_cham, i_bloc)
                    call dyna_gene%get_values(dyna_gene%disc, i_bloc, vr=disc)

                    call extrac(intres, epsi, crit, nbinsg, disc, &
                                zr(jinst+iordr), resu, nbmode, vectgene, ierd)
                    zr(lvar+iordr) = zr(jinst+iordr)
                    zr(lfon+iordr) = vectgene(numcmp)
                end do
                AS_DEALLOCATE(vr=vectgene)
            else
                do iordr = 0, nbordr-1
                    ii = zi(lordr+iordr)
                    call dyna_gene%get_values_by_index(i_cham, ii, shift, vr=resu)
                    zr(lvar+iordr) = zr(jinst+iordr)
                    zr(lfon+iordr) = resu(nbmode*(ii-1-shift)+numcmp)
                end do
            end if
        else
            call dismoi('BASE_MODALE', nomres, 'RESU_DYNA', repk=basemo)
            call dismoi('NUME_DDL', basemo, 'RESU_DYNA', repk=nume)
            call dismoi('NOM_MAILLA', nume, 'NUME_DDL', repk=noma)
!
            call dismoi('NB_EQUA', nume, 'NUME_DDL', repi=neq)
            call wkvect('&&RFRGEN.VECT.PROPRE', 'V V R', neq*nbmode, idbase)
            call copmod(basemo, numer=nume, bmodr=zr(idbase), nbmodes=nbmode, nequa=neq)
!
            call getvtx(' ', 'GROUP_NO', scal=nogno, nbret=ngn)
            if (ngn .ne. 0) then
                call utnono(' ', noma, 'NOEUD', nogno, nomno, &
                            iret)
                if (iret .eq. 10) then
                    call utmess('F', 'ELEMENTS_67', sk=nogno)
                else if (iret .eq. 1) then
                    valk(1) = nogno
                    valk(2) = nomno
                    call utmess('A', 'SOUSTRUC_87', nk=2, valk=valk)
                end if
                noeud = nomno
            end if
            call posddl('NUME_DDL', nume, noeud, cmp, inoeud, &
                        iddl)
            if (inoeud .eq. 0) then
                lg1 = lxlgut(noeud)
                call utmess('F', 'UTILITAI_92', sk=noeud(1:lg1))
            else if (iddl .eq. 0) then
                lg1 = lxlgut(noeud)
                lg2 = lxlgut(cmp)
                valk(1) = cmp(1:lg2)
                valk(2) = noeud(1:lg1)
                call utmess('F', 'UTILITAI_93', nk=2, valk=valk)
            end if
!
!        --- RECHERCHE SI UNE ACCELERATION D'ENTRAINEMENT EXISTE ---
            nfonct = 0
            call getvid(' ', 'ACCE_MONO_APPUI', scal=fonct, nbret=nfonct)
            if (nfonct .ne. 0) then
                if (nomcha(1:4) .ne. 'ACCE') then
!           --- ACCE_MONO_APPUI COMPATIBLE UNIQUEMENT AVEC ACCELERATION
                    call utmess('F', 'UTILITAI4_26')
                    goto 999
                end if
                zk24(lpro+3) (5:8) = '_ABS'
            end if
!        --------------------------------------------------------------
            call wkvect(nomfon//'.VALE', 'G V R', 2*nbordr, lvar)
            lfon = lvar+nbordr
            if (intres(1:3) .ne. 'NON') then
                AS_ALLOCATE(vr=vectgene, size=nbmode)
                do iordr = 0, nbordr-1
                    call dyna_gene%get_values_by_disc(i_cham, zr(jinst+iordr), length=nbinsg, &
                                                      vr=resu)
                    call dyna_gene%get_current_bloc(i_cham, i_bloc)
                    call dyna_gene%get_values(dyna_gene%disc, i_bloc, vr=disc)
                    call extrac(intres, epsi, crit, nbinsg, disc, &
                                zr(jinst+iordr), resu, nbmode, vectgene, ierd)
                    call mdgep2(neq, nbmode, zr(idbase), vectgene, iddl, rep)
                    zr(lvar+iordr) = zr(jinst+iordr)
                    zr(lfon+iordr) = rep
                end do
                AS_DEALLOCATE(vr=vectgene)
            else
                do iordr = 0, nbordr-1
                    ii = zi(lordr+iordr)
                    call dyna_gene%get_values_by_index(i_cham, ii, shift, vr=resu)
                    vectgene => resu(nbmode*(ii-1-shift)+1:nbmode*(ii-shift))
                    call mdgep2(neq, nbmode, zr(idbase), vectgene, iddl, rep)
                    zr(lvar+iordr) = zr(jinst+iordr)
                    zr(lfon+iordr) = rep
                end do
            end if
            monmot(1) = 'NON'
            monmot(2) = 'NON'
            nonmot = 'NON'
            call getvtx(' ', 'MULT_APPUI', scal=monmot(1), nbret=n1)
            call getvtx(' ', 'CORR_STAT', scal=monmot(2), nbret=n2)
            if (monmot(1) .eq. 'OUI' .or. monmot(2) .eq. 'OUI') nonmot = 'OUI'
            if (nonmot(1:3) .eq. 'OUI') then
                call jeexin(nomres//'.F'//nomcha(1:3), iret)
                if (iret .eq. 0) then
                    call utmess('F', 'SEISME_45', sk=nomcha)
                end if
                call jeveuo(nomres//'.F'//nomcha(1:3), 'L', jfon)
                call jeveuo(nomres//'.IPSD', 'L', vr=ipsd)
                call jelira(nomres//'.F'//nomcha(1:3), 'LONMAX', nbexci)
                nbexci = nbexci/2
                do iordr = 0, nbordr-1
                    call mdgep4(neq, nbexci, ipsd, zr(lvar+iordr), zk8(jfon), &
                                iddl, rep)
                    zr(lfon+iordr) = zr(lfon+iordr)+rep
                end do
            end if
            call jedetr('&&RFRGEN.VECT.PROPRE')
!
!        --- PRISE EN COMPTE D'UNE ACCELERATION D'ENTRAINEMENT ---
            if (nfonct .ne. 0) then
                do i = 0, nbordr-1
                    iret = 0
                    call fointe('F', fonct, 1, 'INST', zr(jinst+i), &
                                alpha, ier)
!              --- ACCELERATION ABSOLUE = RELATIVE + ENTRAINEMENT ---
                    zr(lfon+i) = zr(lfon+i)+alpha
                end do
            end if
        end if
!     ---------------------------------------------------------------
    end if
    call jedetr(knume)
    call jedetr(kinst)

    call dyna_gene%free()

999 continue
!
    call foattr(' ', 1, nomfon)
!
!     --- VERIFICATION QU'ON A BIEN CREE UNE FONCTION ---
!         ET REMISE DES ABSCISSES EN ORDRE CROISSANT
    call ordonn(nomfon, 0)
!
    call titre()
    if (niv .gt. 1) call foimpr(nomfon, niv, ifm, 0, k8b)
!
    call jedema()

end subroutine
