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

subroutine ace_masse_repartie(nbocc, infdonn, grplmax, lmax, infcarte, nbdisc, zjdlm)
!
!
! --------------------------------------------------------------------------------------------------
!     AFFE_CARA_ELEM
!
!     MASSES RÉPARTIES
!
! --------------------------------------------------------------------------------------------------
! person_in_charge: jean-luc.flejou at edf.fr
!
    use cara_elem_parameter_module
    use cara_elem_info_type
    use cara_elem_carte_type
    implicit none
    integer(kind=8) :: nbocc
    type(cara_elem_info) :: infdonn
    character(len=24) :: grplmax(*)
    integer(kind=8) :: lmax
    type(cara_elem_carte) :: infcarte(*)
    integer(kind=8) :: nbdisc
    integer(kind=8) :: zjdlm(*)
!
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/calcul_cara_maille.h"
#include "asterfort/dismoi.h"
#include "asterfort/fointe.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/in_liste_entier.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/nocart.h"
#include "asterfort/affdis.h"
#include "asterfort/int_to_char8.h"
!
! --------------------------------------------------------------------------------------------------
 integer(kind=8) :: iocc, ii, jj, kk, iret, nbgrp, nb, ldgm, nm, nb_mail_grp, nb_noeu_grp, ifm, irep
    integer(kind=8) :: ndim, appui, ltypmail, imail, ntopo, isym, iv
    integer(kind=8) :: nfct, compte_maille, nb_noeud_uniq, ll, ncmp, nbsurchpoi1, nbsurchpoi2
    integer(kind=8) :: ivr(4)
    real(kind=8) :: lamasse, valfongro, surfacetotale, zero(6)
    real(kind=8) :: eta, maillesurf(2), maillecdg(3)
    character(len=8) :: typm, fongro, discretm, discretk, nommaille
    character(len=24) :: magrma, connex

    integer(kind=8) :: jdc(3), jdv(3), dimcar
    integer(kind=8) :: jdcinf, jdvinf
    character(len=19) :: cart(3), cartdi, masse_type
    logical :: repartition, ok
!
    integer(kind=8) :: nbnoeu, nbmail
    character(len=8) :: noma
!
! --------------------------------------------------------------------------------------------------
    integer(kind=8), pointer :: noeuds(:) => null()
    real(kind=8), pointer :: coord(:) => null()
!
    integer(kind=8), pointer :: nummaisur(:) => null()
    integer(kind=8), pointer :: lstnumnoe(:) => null()
    integer(kind=8), pointer :: lstnummaipoi1(:) => null()
    integer(kind=8), pointer :: lstpoi1surch(:) => null()
    real(kind=8), pointer :: lstcoenoe(:) => null()
! --------------------------------------------------------------------------------------------------
    integer(kind=8)           :: vmessi(4)
    character(len=24) :: vmessk(6)
! --------------------------------------------------------------------------------------------------
!
    if (nbocc .eq. 0) goto 999
!
    call jemarq()
!
    noma = infdonn%maillage
    nbnoeu = infdonn%nbnoeu
    ndim = infdonn%dimmod
    nbmail = infdonn%nbmail
    ivr(:) = infdonn%ivr(:)
!
!   Pour les discrets c'est obligatoirement du 2d ou 3d
    ASSERT((ndim .eq. 2) .or. (ndim .eq. 3))
!
!   Les cartes sont déjà construites : ace_crea_carte
    cartdi = infcarte(ACE_CAR_DINFO)%nom_carte
    jdcinf = infcarte(ACE_CAR_DINFO)%adr_cmp
    jdvinf = infcarte(ACE_CAR_DINFO)%adr_val
    dimcar = infcarte(ACE_CAR_DINFO)%nbr_cmp
!
    cart(1) = infcarte(ACE_CAR_DISCK)%nom_carte
    jdc(1) = infcarte(ACE_CAR_DISCK)%adr_cmp
    jdv(1) = infcarte(ACE_CAR_DISCK)%adr_val
!
    cart(2) = infcarte(ACE_CAR_DISCM)%nom_carte
    jdc(2) = infcarte(ACE_CAR_DISCM)%adr_cmp
    jdv(2) = infcarte(ACE_CAR_DISCM)%adr_val
!
    cart(3) = infcarte(ACE_CAR_DISCA)%nom_carte
    jdc(3) = infcarte(ACE_CAR_DISCA)%adr_cmp
    jdv(3) = infcarte(ACE_CAR_DISCA)%adr_val
!
    ifm = ivr(4)
!
    magrma = noma//'.GROUPEMA'
    connex = noma//'.CONNEX'
    call jeveuo(noma//'.TYPMAIL', 'L', ltypmail)
    call jeveuo(noma//'.COORDO    .VALE', 'L', vr=coord)
!
!   Comptage des mailles POI1 communes entre les différentes occurrences
    AS_ALLOCATE(vi=lstnummaipoi1, size=nbdisc)
    compte_maille = 0
!   Les mailles de la 1ère occurrence
    call getvtx('MASS_REP', 'GROUP_MA_POI1', iocc=1, nbval=lmax, vect=grplmax, nbret=nbgrp)
!   on éclate les GROUP_MA en mailles
    do ii = 1, nbgrp
        call jelira(jexnom(magrma, grplmax(ii)), 'LONUTI', nb)
        call jeveuo(jexnom(magrma, grplmax(ii)), 'L', ldgm)
        if (compte_maille+nb .gt. nbdisc) then
            vmessi(1) = 1
            vmessi(2) = nb
            vmessi(3) = nbdisc
            vmessi(4) = compte_maille
            vmessk(1) = grplmax(ii)
            call utmess('F', 'AFFECARAELEM_21', nk=1, valk=vmessk, ni=4, vali=vmessi)
        end if
        do jj = ldgm, ldgm+nb-1
            compte_maille = compte_maille+1
            lstnummaipoi1(compte_maille) = zi(jj)
        end do
    end do
!   Les mailles des autres occurrences
    if (nbocc .ge. 2) then
        AS_ALLOCATE(vi=lstpoi1surch, size=nbdisc)
        nbsurchpoi1 = 0
        nbsurchpoi2 = 0
        do iocc = 2, nbocc
            call getvtx('MASS_REP', 'GROUP_MA_POI1', iocc=iocc, nbval=lmax, vect=grplmax, &
                        nbret=nbgrp)
!           on éclate les GROUP_MA en mailles
            do ii = 1, nbgrp
                kk = compte_maille
                call jelira(jexnom(magrma, grplmax(ii)), 'LONUTI', nb)
                call jeveuo(jexnom(magrma, grplmax(ii)), 'L', ldgm)
                do jj = ldgm, ldgm+nb-1
                    imail = zi(jj)
                    if (in_liste_entier(imail, lstnummaipoi1(1:compte_maille))) then
                        if (.not. in_liste_entier(imail, lstpoi1surch(1:nbsurchpoi2))) then
                            if (nbsurchpoi2 .lt. nbdisc) then
                                nbsurchpoi2 = nbsurchpoi2+1
                                lstpoi1surch(nbsurchpoi2) = imail
                            end if
                        end if
                        nbsurchpoi1 = nbsurchpoi1+1
                    else
                        compte_maille = compte_maille+1
                        if (compte_maille .gt. nbdisc) then
                            vmessi(1) = iocc
                            vmessi(2) = nb
                            vmessi(3) = nbdisc
                            vmessi(4) = kk
                            vmessk(1) = grplmax(ii)
                            call utmess('F', 'AFFECARAELEM_21', nk=1, valk=vmessk, ni=4, &
                                        vali=vmessi)
                        end if
                        lstnummaipoi1(compte_maille) = imail
                    end if
                end do
            end do
        end do
        if (nbsurchpoi1 .ne. 0) then
            vmessk(1) = 'MASS_REP'
            vmessi(1) = nbocc
            vmessi(2) = nbsurchpoi1
            call utmess('I', 'AFFECARAELEM_20', nk=1, valk=vmessk, ni=2, vali=vmessi)
            if (ivr(3) .eq. 2) then
                write (ifm, '(A,I3,A)') 'Liste des ', nbsurchpoi2, ' mailles concernées :'
                do ii = 1, nbsurchpoi2
                    nommaille = int_to_char8(lstpoi1surch(ii))
                    if (mod(ii, 8) .eq. 0) then
                        write (ifm, '(A8)') nommaille
                    else
                        write (ifm, '(A8,2X)', advance="no") nommaille
                    end if
                end do
                write (ifm, *)
            end if
        end if
        AS_DEALLOCATE(vi=lstpoi1surch)
    end if
    AS_DEALLOCATE(vi=lstnummaipoi1)
!
    do iocc = 1, nbocc
!       Dimension de l'appui pour cette occurrence
        if (ndim .eq. 2) then
            appui = 1
        else
!           En 3D la dimension n'est pas encore déterminée
            appui = -1
        end if
!       Pour les messages
        vmessi(1) = iocc
        vmessk(1) = 'MASS_REP'
!
        nb_mail_grp = 0
        nb_noeu_grp = 0
        call getvtx('MASS_REP', 'GROUP_MA', iocc=iocc, nbval=lmax, vect=grplmax, nbret=nbgrp)
!       on éclate les GROUP_MA en mailles pour déterminer les noeuds concernés
        do ii = 1, nbgrp
            call jelira(jexnom(magrma, grplmax(ii)), 'LONUTI', nb)
            call jeveuo(jexnom(magrma, grplmax(ii)), 'L', ldgm)
            nb_mail_grp = nb_mail_grp+nb
            do jj = ldgm, ldgm+nb-1
                imail = zi(jj)
                call jenuno(jexnum('&CATA.TM.NOMTM', zi(ltypmail-1+imail)), typm)
                call dismoi('DIM_TOPO', typm, 'TYPE_MAILLE', repi=ntopo)
!               La dimension de la première maille définit l'appui si pas encore déterminé
                if (appui .eq. -1) then
                    appui = ntopo
                    if ((appui .ne. 1) .and. (appui .ne. 2)) then
                        vmessk(2) = grplmax(ii)
                        vmessk(3) = int_to_char8(imail)
                        vmessk(4) = typm
                        vmessi(2) = appui
                        call utmess('F', 'AFFECARAELEM_11', nk=4, valk=vmessk, ni=2, vali=vmessi)
                    end if
                end if
                if (ntopo .ne. appui) then
                    vmessk(2) = grplmax(ii)
!                   Maille courante qui va pas bien
                    vmessk(3) = int_to_char8(imail)
                    vmessk(4) = typm
                    vmessi(2) = appui
                    if (ndim .eq. 3) then
!                       Maille de définition de la topologie du groupe
                        vmessk(5) = int_to_char8(zi(ldgm))
                        call jenuno(jexnum('&CATA.TM.NOMTM', zi(ltypmail-1+zi(ldgm))), vmessk(6))
                        vmessi(3) = ntopo
!                       vmessi(1) N°occurrence
!                       vmessi(2) Topologie de la 1ère maille du groupe
!                       vmessi(3) Topologie de la maille courante
!                       vmessk(1) 'MASS_REP'
!                       vmessk(2) Nom du groupe
!                       vmessk(3) Nom de la maille courante
!                       vmessk(4) Type de la maille courante
!                       vmessk(5) Nom de la 1ère maille du groupe
!                       vmessk(6) Type de la 1ère maille du groupe
                        call utmess('F', 'AFFECARAELEM_23', nk=6, valk=vmessk, ni=3, vali=vmessi)
                    else
!                       vmessi(1) N°occurrence
!                       vmessi(2) Topologie =1 on est en 2D
!                       vmessi(3) Topologie de la maille courante
!                       vmessk(1) 'MASS_REP'
!                       vmessk(2) Nom du groupe
!                       vmessk(3) Nom de la maille courante
!                       vmessk(4) Type de la maille courante
                        call utmess('F', 'AFFECARAELEM_10', nk=4, valk=vmessk, ni=3, vali=vmessi)
                    end if
                end if
!               Noeuds de la maille
                call jelira(jexnum(connex, imail), 'LONMAX', nm)
                call jeveuo(jexnum(connex, imail), 'L', vi=noeuds)
!               Si pas le bon nombre de noeuds
                if (.not. in_liste_entier(nm, [2, 3, 4, 6, 7, 8, 9])) then
                    vmessk(2) = grplmax(ii)
                    vmessk(3) = int_to_char8(imail)
                    vmessk(4) = typm
                    vmessk(5) = '2 3 4 6 7 8 9'
                    vmessi(2) = nm
                    call utmess('F', 'AFFECARAELEM_13', nk=5, valk=vmessk, ni=2, vali=vmessi)
                end if
!               Tout semble ok
                nb_noeu_grp = nb_noeu_grp+nm
            end do
        end do
!
!       La masse à répartir
        call getvr8('MASS_REP', 'VALE', iocc=iocc, scal=lamasse)
!       Le type de masse : TOTALE LINEIQUE SURFACIQUE
        call getvtx('MASS_REP', 'TYPE', iocc=iocc, scal=masse_type)
!
        if (appui .eq. 2 .and. masse_type .eq. 'LINEIQUE') then
            call utmess('F', 'AFFECARAELEM_17', nk=1, valk=vmessk, ni=1, vali=vmessi)
        end if
        if (appui .eq. 1 .and. masse_type .eq. 'SURFACIQUE') then
            call utmess('F', 'AFFECARAELEM_18', nk=1, valk=vmessk, ni=1, vali=vmessi)
        end if
!
        if (masse_type .eq. 'TOTALE') then
            repartition = .true.
        else
            repartition = .false.
        end if
        call getvid('MASS_REP', 'FONC_MULT', iocc=iocc, scal=fongro, nbret=nfct)
        ASSERT(nfct .eq. 0 .or. nfct .eq. 1)
!
!       Numéro des mailles de surface, pour vérifier qu'il n'y a pas de doublon
        AS_ALLOCATE(vi=nummaisur, size=nb_mail_grp)
!       Numéro des noeuds de la surface, maille POI1, la masse pondérée
        AS_ALLOCATE(vi=lstnumnoe, size=nb_noeu_grp)
        AS_ALLOCATE(vi=lstnummaipoi1, size=nb_noeu_grp)
        AS_ALLOCATE(vr=lstcoenoe, size=nb_noeu_grp)
!
        lstnumnoe(:) = 0
        lstcoenoe(:) = 0.0
        nummaisur(:) = 0
!       Va permmettre de vérifier la bijectivité entre les noeuds de la surface et les POI1
!           Valeur négative pour ne pas être égale à un numéro de maille
        lstnummaipoi1(:) = -2
!
        nb_noeud_uniq = 0
        compte_maille = 0
        surfacetotale = 0.0
        do ii = 1, nbgrp
            call jelira(jexnom(magrma, grplmax(ii)), 'LONUTI', nb)
            call jeveuo(jexnom(magrma, grplmax(ii)), 'L', ldgm)
            do jj = ldgm, ldgm+nb-1
                imail = zi(jj)
!               Vérification que la maille n'est pas en double dans les groupes
                if (in_liste_entier(imail, nummaisur(1:compte_maille))) then
                    vmessk(2) = grplmax(ii)
                    vmessk(3) = int_to_char8(imail)
                    call utmess('F', 'AFFECARAELEM_12', nk=3, valk=vmessk, ni=1, vali=vmessi)
                end if
!               Calcul sur la maille
                compte_maille = compte_maille+1
                nummaisur(compte_maille) = imail
                call jelira(jexnum(connex, imail), 'LONMAX', nm)
                call jeveuo(jexnum(connex, imail), 'L', vi=noeuds)
                call calcul_cara_maille(coord, noeuds(1:nm), appui, maillesurf, maillecdg)
                surfacetotale = surfacetotale+maillesurf(1)
                valfongro = 1.0
                if (nfct .eq. 1) then
                    call fointe('F ', fongro, 3, ['X', 'Y', 'Z'], maillecdg, valfongro, iret)
                end if
                do kk = 1, nm
                    ok = in_liste_entier(noeuds(kk), lstnumnoe(1:nb_noeud_uniq), indx=ll)
                    if (.not. ok) then
                        nb_noeud_uniq = nb_noeud_uniq+1
                        lstnumnoe(nb_noeud_uniq) = noeuds(kk)
                        lstcoenoe(nb_noeud_uniq) = lamasse*maillesurf(2)*valfongro
                    else
                        lstcoenoe(ll) = lstcoenoe(ll)+lamasse*maillesurf(2)*valfongro
                    end if
                end do
            end do
        end do
!
!       Doit-on diviser par la surface totale
        if (repartition) then
            do ll = 1, nb_noeud_uniq
                lstcoenoe(ll) = lstcoenoe(ll)/surfacetotale
            end do
        end if
!
!       on éclate les GROUP_MA_POI1 pour vérifier que les noeuds sont dans le GROUP_MA
        call getvtx('MASS_REP', 'GROUP_MA_POI1', iocc=iocc, nbval=lmax, vect=grplmax, nbret=nbgrp)
        do ii = 1, nbgrp
            call jelira(jexnom(magrma, grplmax(ii)), 'LONUTI', nb)
            call jeveuo(jexnom(magrma, grplmax(ii)), 'L', ldgm)
            do jj = ldgm, ldgm+nb-1
                imail = zi(jj)
!               Si la maille n'est pas affectée
                if (zjdlm(imail) .eq. 0) then
                    vmessk(2) = grplmax(ii)
                    vmessk(3) = int_to_char8(imail)
                    call utmess('F', 'AFFECARAELEM_22', nk=3, valk=vmessk, ni=1, vali=vmessi)
                end if
!               Si la maille n'a pas 1 noeud
                call jelira(jexnum(connex, imail), 'LONMAX', nm)
                if (nm .ne. 1) then
                    vmessk(2) = grplmax(ii)
                    vmessk(3) = int_to_char8(imail)
                    call utmess('F', 'AFFECARAELEM_16', nk=3, valk=vmessk, ni=1, vali=vmessi)
                end if
                call jeveuo(jexnum(connex, imail), 'L', vi=noeuds)
                ok = in_liste_entier(noeuds(1), lstnumnoe(1:nb_noeud_uniq), indx=kk)
!               Si le noeud POI1 n'est pas dans la liste des noeuds de la surface ==> Pas bien
                if (.not. ok) then
                    vmessk(2) = grplmax(ii)
                    vmessk(3) = int_to_char8(imail)
                    call utmess('F', 'AFFECARAELEM_14', nk=3, valk=vmessk, ni=1, vali=vmessi)
                end if
!               Si on déjà mis une maille POI1 en face du noeud  ==> Pas bien
                if (lstnummaipoi1(kk) .ne. -2) then
                    vmessk(2) = grplmax(ii)
                    vmessk(3) = int_to_char8(imail)
                    vmessk(4) = int_to_char8(lstnummaipoi1(kk))
                    vmessk(5) = int_to_char8(noeuds(1))
                    call utmess('F', 'AFFECARAELEM_19', nk=5, valk=vmessk, ni=1, vali=vmessi)
                end if
                lstnummaipoi1(kk) = imail
            end do
        end do
!       La relation doit être bijective ==> on ne doit plus avoir lstnummaipoi1 = -2.
!           Tous les noeuds doivent avoir une maille POI1 en face
        do ll = 1, nb_noeud_uniq
            if (lstnummaipoi1(ll) .eq. -2) then
                vmessk(2) = int_to_char8(lstnumnoe(ll))
                call utmess('F', 'AFFECARAELEM_15', nk=2, valk=vmessk, ni=1, vali=vmessi)
            end if
        end do
!
        AS_DEALLOCATE(vi=nummaisur)
!
        if (ivr(3) .eq. 2) then
            write (ifm, 100) iocc, surfacetotale
        end if
        lamasse = 0.0
!       Par défaut on est dans le repère global, matrices symétriques
        irep = 1; isym = 1; eta = 0.0; zero(:) = 0.0
        discretm = 'M_T_D_N'; discretk = 'K_T_D_N'
        do ii = 1, nb_noeud_uniq
            iv = 1
            if (ivr(3) .eq. 2) then
                vmessk(1) = int_to_char8(lstnummaipoi1(ii))
                vmessk(2) = int_to_char8(lstnumnoe(ii))
                write (ifm, 110) vmessk(1), vmessk(2), lstcoenoe(ii)
            end if
            lamasse = lamasse+lstcoenoe(ii)
            call affdis(ndim, irep, eta, discretm, lstcoenoe(ii:), &
                        jdc, jdv, ivr, iv, ['K', 'M', 'A'], &
                        ncmp, ll, jdcinf, jdvinf, isym)
            call nocart(cart(ll), 3, ncmp, mode='NUM', nma=1, limanu=[lstnummaipoi1(ii)])
            call nocart(cartdi, 3, dimcar, mode='NUM', nma=1, limanu=[lstnummaipoi1(ii)])
!           On met 0 sur les raideurs
            iv = 1
            call affdis(ndim, irep, eta, discretk, zero, &
                        jdc, jdv, ivr, iv, ['K', 'M', 'A'], &
                        ncmp, ll, jdcinf, jdvinf, isym)
            call nocart(cart(ll), 3, ncmp, mode='NUM', nma=1, limanu=[lstnummaipoi1(ii)])
            call nocart(cartdi, 3, dimcar, mode='NUM', nma=1, limanu=[lstnummaipoi1(ii)])
        end do
        if (ivr(3) .eq. 2) then
            write (ifm, 120) lamasse
        end if
!
        AS_DEALLOCATE(vi=lstnumnoe)
        AS_DEALLOCATE(vr=lstcoenoe)
        AS_DEALLOCATE(vi=lstnummaipoi1)
    end do
!
!   Pour l'impression des valeurs affectées : Maille , Noeud , valeur
100 format(/, 'OCCURRENCE : ', i4, '. SURFACE TOTALE : ', 1pe12.5)
110 format('   MAILLE : ', A8, '. NOEUD : ', A8, '. MASSE : ', 1pe12.5)
120 format('==>POUR CETTE OCCURRENCE. MASSE TOTALE : ', 1pe12.5)
!
    call jedema()
!
999 continue
end subroutine
