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
subroutine cafono(load, loadLigrel, mesh, model, valeType)
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/getfac.h"
#include "asterc/r8dgrd.h"
#include "asterfort/affono.h"
#include "asterfort/alcart.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/char8_to_int.h"
#include "asterfort/char_crea_ligf.h"
#include "asterfort/char_nb_ligf.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisdg.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/nocart.h"
#include "asterfort/noligr.h"
#include "asterfort/reliem.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=8), intent(in) :: load, mesh, model
    character(len=19), intent(in) :: loadLigrel
    character(len=4), intent(in) :: valeType
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Treatment of load FORCE_NODALE
!
! --------------------------------------------------------------------------------------------------
!
! In  load             : load
! In  loadLigrel
! In  mesh             : mesh
! In  model            : model
! In  valeType         : affected value type (real, complex or function)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: keywordFact = 'FORCE_NODALE'
    integer(kind=8), parameter :: nmocl = 10
    integer(kind=8) :: nfono, n2dl, n3dl, n6dl, ncoq2d, nbcomp
    integer(kind=8) :: i, idgex, ii, in, ino, iret
    integer(kind=8) :: j, jj, jl, jnbno, jno
    integer(kind=8) :: jprnm, jval, jvalv, nangl, nbec, nbecf
    integer(kind=8) :: nbno, nbnoeu, nsurch, numel
    integer(kind=8) :: igrel, inema
    integer(kind=8) :: ntypel(nmocl), forimp(nmocl), ier
    real(kind=8) :: dgrd, valfor(nmocl)
    aster_logical :: verif, l_occu_void
    aster_logical :: lcolle
    character(len=8) :: nomn, typmcl(2), valfof(nmocl)
    character(len=16) :: motcle(nmocl), motcls(2)
    character(len=19) :: carte, modelLigrel
    character(len=24) :: liel, nomnoe, nomele, mesnoe
    integer(kind=8) :: nb_elem_late, nb_noel_maxi, jlgns, iexi
    integer(kind=8), pointer :: desgi(:) => null()
    character(len=8), pointer :: noms_noeuds(:) => null()
    character(len=8), pointer :: ncmp(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    call getfac(keywordFact, nfono)
    if (nfono .eq. 0) goto 999
!
    igrel = 0
    inema = 0
    motcls(1) = 'GROUP_NO'
    motcls(2) = 'NOEUD'
    typmcl(1) = 'GROUP_NO'
    typmcl(2) = 'NOEUD'
!
    verif = .true.
!
! - Count number of late elements
!
    call char_nb_ligf(mesh, keywordFact, 'Node', nb_elem_late, nb_noel_maxi)
!
! - Create <LIGREL> on late elements
!
    call char_crea_ligf(mesh, loadLigrel, nb_elem_late, nb_noel_maxi)
!
    call jenonu(jexnom('&CATA.TE.NOMTE', 'FORCE_NOD_2DDL'), n2dl)
    call jenonu(jexnom('&CATA.TE.NOMTE', 'FORCE_NOD_3DDL'), n3dl)
    call jenonu(jexnom('&CATA.TE.NOMTE', 'FORCE_NOD_6DDL'), n6dl)
    call jenonu(jexnom('&CATA.TE.NOMTE', 'FORCE_NOD_COQ2D'), ncoq2d)
    ntypel(1) = n2dl
    ntypel(2) = n2dl
    ntypel(3) = n3dl
    ntypel(4) = n6dl
    ntypel(5) = n6dl
    ntypel(6) = n6dl
!
! ---------------------------------------------------
!     RECUPERATION DES MOTS-CLES DDL POSSIBLES SOUS FORCE_NODALE
! ---------------------------------------------------
    motcle(1) = 'FX'
    motcle(2) = 'FY'
    motcle(3) = 'FZ'
    motcle(4) = 'MX'
    motcle(5) = 'MY'
    motcle(6) = 'MZ'
    motcle(7) = 'REP'
    motcle(8) = 'ALPHA'
    motcle(9) = 'BETA'
    motcle(10) = 'GAMMA'
    nbcomp = 10
!
! ---------------------------------------------------
! *** RECUPERATION DU DESCRIPTEUR GRANDEUR .PRNM
! *** DU MODELE
! ---------------------------------------------------
!
    call dismoi('NB_EC', 'FORC_R', 'GRANDEUR', repi=nbecf)
    if (nbecf .gt. 11) then
        call utmess('F', 'MODELISA2_65')
    else
        call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
        call jeveuo(modelLigrel//'.PRNM', 'L', jprnm)
    end if
!
    call dismoi('NB_EC', 'DEPL_R', 'GRANDEUR', repi=nbec)
    if (nbec .gt. 11) then
        call utmess('F', 'MODELISA_94')
    end if
!
    call jeveuo(loadLigrel//'.NBNO', 'E', jnbno)
    nomnoe = mesh//'.COORDO    .VALE'
    call jelira(nomnoe, 'LONMAX', nbnoeu)
    nbnoeu = nbnoeu/3
!
    mesnoe = '&&CAFONO.MES_NOEUDS'
    motcls(1) = 'GROUP_NO'
    motcls(2) = 'NOEUD'
    typmcl(1) = 'GROUP_NO'
    typmcl(2) = 'NOEUD'
!
! ---------------------------------------------------
!     ALLOCATION DE TABLEAUX DE TRAVAIL
! ---------------------------------------------------
!   OBJETS INTERMEDIAIRES PERMETTANT D'APPLIQUER LA REGLE DE SURCHARGE
!        -  VECTEUR (K8) CONTENANT LES NOMS DES NOEUDS
!        -  TABLEAU DES VALEURS DES DDLS DES FORCES IMPOSEES
!                         DIM NBNOEU * NBCOMP
!        -  VECTEUR (IS) CONTENANT LE DESCRIPTEUR GRANDEUR ASSOCIE AUX
!                         FORCES IMPOSEES PAR NOEUD
!
    AS_ALLOCATE(vk8=noms_noeuds, size=nbnoeu)
    if (valeType .eq. 'REEL') then
        call wkvect('&&CAFONO.VALDDLR', 'V V R', nbcomp*nbnoeu, jval)
    else
        call wkvect('&&CAFONO.VALDDLF', 'V V K8', nbcomp*nbnoeu, jval)
    end if
    AS_ALLOCATE(vi=desgi, size=nbnoeu)
!
    dgrd = r8dgrd()
    if (valeType .eq. 'FONC') then
        do i = 1, nbcomp*nbnoeu
            zk8(jval-1+i) = '&FOZERO'
        end do
    end if
    nsurch = 0
!
! --------------------------------------------------------------
!     BOUCLE SUR LES OCCURENCES DU MOT-CLE FACTEUR FORCE_NODALE
! --------------------------------------------------------------
!
    lcolle = .false.
    call jeexin(mesh//'.NOMNOE', ier)
    if (ier .ne. 0) then
        lcolle = .true.
    end if
    do i = 1, nfono
        do ii = 1, nbcomp
            forimp(ii) = 0
        end do
!
        if (valeType .eq. 'REEL') then
            do j = 1, 6
                call getvr8(keywordFact, motcle(j), iocc=i, scal=valfor(j), nbret=forimp(j))
            end do
!
            call getvr8(keywordFact, 'ANGL_NAUT', iocc=i, nbval=3, vect=valfor(8), &
                        nbret=nangl)
            if (nangl .ne. 0) then
!              --- REPERE UTILISATEUR ---
                valfor(7) = -1.d0
                forimp(7) = 1
                do ii = 1, min(3, abs(nangl))
                    valfor(7+ii) = valfor(7+ii)*dgrd
                    forimp(7+ii) = 1
                end do
            else
!              --- REPERE GLOBAL ---
                valfor(7) = 0.d0
            end if
!
        else if (valeType .eq. 'FONC') then
            do ii = 1, nbcomp
                valfof(ii) = '&FOZERO'
            end do
            do j = 1, 6
                call getvid(keywordFact, motcle(j), iocc=i, scal=valfof(j), nbret=forimp(j))
            end do
!
            call getvid(keywordFact, 'ANGL_NAUT', iocc=i, nbval=3, vect=valfof(8), &
                        nbret=nangl)
            if (nangl .ne. 0) then
!              --- REPERE UTILISATEUR ---
                valfof(7) = 'UTILISAT'
                forimp(7) = 1
                do ii = 1, min(3, abs(nangl))
                    forimp(7+ii) = 1
                end do
            else
!              --- REPERE GLOBAL ---
                valfof(7) = 'GLOBAL'
            end if
        end if
        if (nangl .lt. 0) then
            call utmess('A', 'MODELISA2_66')
        end if
!
!       ---------------------------
!       CAS DE GROUP_NO ET DE NOEUD
!       ---------------------------
!
        call reliem(' ', mesh, 'NO_NOEUD', keywordFact, i, &
                    2, motcls, typmcl, mesnoe, nbno)
        if (nbno .eq. 0) goto 110
        call jeveuo(mesnoe, 'L', jno)
!
        l_occu_void = .true.
        do jj = 1, nbno
            ino = char8_to_int(zk8(jno-1+jj), lcolle, mesh, "NOEUD")
            noms_noeuds(ino) = zk8(jno-1+jj)
            call affono(zr(jval), zk8(jval), desgi(ino), zi(jprnm-1+(ino-1)*nbec+1), nbcomp, &
                        valeType, zk8(jno-1+jj), ino, nsurch, forimp, &
                        valfor, valfof, motcle, verif, nbec)

            if (desgi(ino) .ne. 0) l_occu_void = .false.
        end do

        if (l_occu_void) then
            do ii = 1, 6
                if (forimp(ii) .ne. 0) then
                    call utmess('F', 'CHARGES2_46', sk=motcle(ii))
                end if
            end do
        end if
!
        call jedetr(mesnoe)
110     continue
    end do
!
!     -----------------------------------------------
!     AFFECTATION DU LIGREL ET STOCKAGE DANS LA CARTE
!              DIMENSIONS AUX VRAIES VALEURS
!     -----------------------------------------------
!
    liel = loadLigrel//'.LIEL'
    carte = load//'.CHME.FORNO'
!
    call jeexin(carte//'.DESC', iret)
!
    if (iret .eq. 0) then
        if (valeType .eq. 'REEL') then
            call alcart('G', carte, mesh, 'FORC_R')
        else if (valeType .eq. 'FONC') then
            call alcart('G', carte, mesh, 'FORC_F')
        else
            ASSERT(.false.)
        end if
    end if
!
    call jeveuo(carte//'.NCMP', 'E', vk8=ncmp)
    call jeveuo(carte//'.VALV', 'E', jvalv)
!
    ncmp(1) = 'FX'
    ncmp(2) = 'FY'
    ncmp(3) = 'FZ'
    ncmp(4) = 'MX'
    ncmp(5) = 'MY'
    ncmp(6) = 'MZ'
    ncmp(7) = 'REP'
    ncmp(8) = 'ALPHA'
    ncmp(9) = 'BETA'
    ncmp(10) = 'GAMMA'
!
    call jeveuo(loadLigrel//'.NBNO', 'E', jnbno)
    call jeexin(loadLigrel//'.LGNS', iexi)
    if (iexi .gt. 0) then
        call jeveuo(loadLigrel//'.LGNS', 'E', jlgns)
    else
        jlgns = 1
    end if
!
!     -----------------------------------------------
!     BOUCLE SUR TOUS LES NOEUDS DU MAILLAGE
!     -----------------------------------------------
!
    do ino = 1, nbnoeu
!
        if (desgi(ino) .ne. 0) then
!
            nomn = noms_noeuds(ino)
            in = char8_to_int(nomn, lcolle, mesh, "NOEUD")
            idgex = jprnm-1+(in-1)*nbec+1
!
            do i = 1, 6
                if (exisdg(zi(idgex), i)) then
                    numel = ntypel(i)
                end if
            end do
            if ((exisdg(zi(idgex), 6)) .and. (.not. (exisdg(zi(idgex), 4)))) then
                numel = ncoq2d
            end if
!
            igrel = igrel+1
            call jenuno(jexnum('&CATA.TE.NOMTE', numel), nomele)
            call noligr(loadLigrel, igrel, numel, in, &
                        1, inema, zi(jnbno), jlgns)
!
            call jeveuo(jexnum(liel, igrel), 'E', jl)
            if (valeType .eq. 'REEL') then
                do i = 1, nbcomp
                    zr(jvalv-1+i) = zr(jval-1+nbcomp*(ino-1)+i)
                end do
            else
                do i = 1, nbcomp
                    zk8(jvalv-1+i) = zk8(jval-1+nbcomp*(ino-1)+i)
                end do
            end if
!
!   ON CREE UNE CARTE POUR CHAQUE NOEUD AFFECTE ET ON NOTE TOUTES
!   LES COMPOSANTES (NBCOMP)
!
            call nocart(carte, -3, nbcomp, ligrel=liel, nma=1, &
                        limanu=[zi(jl)])
!
        end if
!
    end do
!
    AS_DEALLOCATE(vk8=noms_noeuds)
    AS_DEALLOCATE(vi=desgi)
    if (valeType .eq. 'REEL') then
        call jedetr('&&CAFONO.VALDDLR')
    else if (valeType .eq. 'FONC') then
        call jedetr('&&CAFONO.VALDDLF')
    end if
999 continue
    call jedema()
end subroutine
