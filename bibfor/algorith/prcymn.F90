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
subroutine prcymn(nomres, soumat, repmat)
    implicit none
!  P. RICHARD     DATE 11/03/91
!-----------------------------------------------------------------------
!  BUT : < PROJECTION DES MATRICES POUR CYCLIQUE MAC-NEAL >
!
!        PROJETER LES MATRICES MASSE ET RAIDEUR ET SORTIR LES SOUS
!        MATRICES POUR TRAITER LE CAS CYCLIQUE AVEC INTERFACES
!        DE MAC NEAL OU AUCUN ( = MAC NEAL SANS MODE ATTACHE)
!-----------------------------------------------------------------------
!
! NOMRES /I/ : NOM UTILISATEUR DU RESULTAT
! SOUMAT /I/ : NOM K24 DE LA FAMILLE DES SOUS-MATRICES
! REPMAT /I/ : NOM K24 DU REPERTOIRE DES NOMS DES SOUS-MATRICES
!
!
!
#include "jeveux.h"
#include "asterfort/amppr.h"
#include "asterfort/bmradi.h"
#include "asterfort/ctetax.h"
#include "asterfort/ctetgd.h"
#include "asterfort/dcapno.h"
#include "asterfort/dismoi.h"
#include "asterfort/flexib.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecreo.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/pmppr.h"
#include "asterfort/r8inir.h"
#include "asterfort/rsadpa.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/zerlag.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
!
!
!
    character(len=6) :: pgc
    character(len=8) :: nomres, basmod, intf, kbid, k8bid
    character(len=19) :: numddl, raid
    character(len=24) :: repmat, soumat, noeint, chamva
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ibid, iran, j, k
    integer(kind=8) :: ldk0aa, ldk0ai, ldk0aj, ldk0ii, ldk0ji, ldk0jj, ldkpaa
    integer(kind=8) :: ldkpai, ldkpaj, ldkpja, ldkpji, ldkpjj, ldm0ii, llcham
    integer(kind=8) :: llkge, llmge, llnoa, llnod, llnog
    integer(kind=8) :: ltcap, ltcdp, ltcgp, ltetax
    integer(kind=8) :: ltetgd, ltexa, ltexd, ltexg, ltflex, ltmat, ltvec
    integer(kind=8) :: nbdax, nbddr, nbmod, nbnoa, nbnod, nbnog, nbsec
    integer(kind=8) :: nbsma, nbv, neq, ntail, ntrian, numa, numd
    integer(kind=8) :: numg
    real(kind=8) :: xx
    integer(kind=8), pointer :: deeq(:) => null()
    integer(kind=8), pointer :: cycl_nuin(:) => null()
    character(len=24), pointer :: cycl_refe(:) => null()
    integer(kind=8), pointer :: cycl_desc(:) => null()
    integer(kind=8), pointer :: cycl_nbsc(:) => null()
    blas_int :: b_incx, b_incy, b_n
!-----------------------------------------------------------------------
    data pgc/'PRCYMN'/
!-----------------------------------------------------------------------
!
! --- RECUPERATION DES CONCEPTS AMONT
!
    call jemarq()
    call jeveuo(nomres//'.CYCL_REFE', 'L', vk24=cycl_refe)
    intf = cycl_refe(2) (1:8)
    basmod = cycl_refe(3) (1:8)
    call jelibe(nomres//'.CYCL_REFE')
!
!
    call dismoi('NUME_DDL', basmod, 'RESU_DYNA', repk=numddl)
!----ON AJOUT .NUME POUR OBTENIR LE NUME_EQUA
    numddl(15:19) = '.NUME'
!
    call dismoi('REF_RIGI_PREM', basmod, 'RESU_DYNA', repk=raid)
!
! --- RECUPERATION DES DIMENSIONS DU PROBLEME GENERALISE
!
    call jeveuo(nomres//'.CYCL_DESC', 'L', vi=cycl_desc)
    nbmod = cycl_desc(1)
    nbddr = cycl_desc(2)
    nbdax = cycl_desc(3)
    call jelibe(nomres//'.CYCL_DESC')
!
! --- ALLOCATION DU REPERTOIRE DES NOMS DES SOUS-MATRICES
!
!
!-- On a constaté que :
!--   Si on a des DDL d'axe, le code plante (erreur jeveux - numéro d'objet invalide)
!--   Si on a 2 ou 4 secteurs, le changement de variable introduit une singularité
!--
!-- Compte tenu de la faible utilisation de cette option, on s'arrête en erreur fatale
!-- dans ces deux cas, et on propose d'utiliser la méthode "CRAIG BAMPTON" à la place.
!
    if (nbdax .gt. 0) then
        call utmess('F', 'ALGORITH14_90')
    end if
!
!
! --- RECUPERATION DU NOMBRE DE SECTEURS
!
    call jeveuo(nomres//'.CYCL_NBSC', 'L', vi=cycl_nbsc)
    nbsec = cycl_nbsc(1)
    call jelibe(nomres//'.CYCL_NBSC')
!
    if (nbsec .eq. 2) then
        call utmess('F', 'ALGORITH14_92')
    end if
    if (nbsec .eq. 4) then
        call utmess('F', 'ALGORITH14_92')
    end if
!
    if (nbdax .gt. 0) then
        nbsma = 13
    else
        nbsma = 6
    end if
!
    call jecreo(repmat, 'V N K8')
    call jeecra(repmat, 'NOMMAX', nbsma)
!
! --- CREATION DE LA FAMILLE DES SOUS-MATRICES
!
    call jecrec(soumat, 'V V R', 'NU', 'DISPERSE', 'VARIABLE', &
                nbsma)
!
! --- ALLOCATION DES MATRICES (STOCKAGE PLEIN)
!
    ntail = nbmod*(nbmod+1)/2
!
    call jecroc(jexnom(repmat, 'K0II'))
    call jenonu(jexnom(repmat, 'K0II'), ibid)
    call jeecra(jexnum(soumat, ibid), 'LONMAX', ntail)
!
    call jecroc(jexnom(repmat, 'M0II'))
    call jenonu(jexnom(repmat, 'M0II'), ibid)
    call jeecra(jexnum(soumat, ibid), 'LONMAX', ntail)
!
    ntrian = nbddr*(nbddr+1)/2
    ntail = nbddr*nbddr
!
    call jecroc(jexnom(repmat, 'K0JI'))
    call jenonu(jexnom(repmat, 'K0JI'), ibid)
    call jeecra(jexnum(soumat, ibid), 'LONMAX', nbmod*nbddr)
!
    call jecroc(jexnom(repmat, 'K0JJ'))
    call jenonu(jexnom(repmat, 'K0JJ'), ibid)
    call jeecra(jexnum(soumat, ibid), 'LONMAX', ntrian)
!
    call jecroc(jexnom(repmat, 'KPLUSJI'))
    call jenonu(jexnom(repmat, 'KPLUSJI'), ibid)
    call jeecra(jexnum(soumat, ibid), 'LONMAX', nbmod*nbddr)
!
    call jecroc(jexnom(repmat, 'KPLUSJJ'))
    call jenonu(jexnom(repmat, 'KPLUSJJ'), ibid)
    call jeecra(jexnum(soumat, ibid), 'LONMAX', ntail)
!
    if (nbdax .gt. 0) then
!
        ntail = nbddr*nbdax
!
        call jecroc(jexnom(repmat, 'K0AI'))
        call jenonu(jexnom(repmat, 'K0AI'), ibid)
        call jeecra(jexnum(soumat, ibid), 'LONMAX', nbmod*nbdax)
!
        call jecroc(jexnom(repmat, 'K0AJ'))
        call jenonu(jexnom(repmat, 'K0AJ'), ibid)
        call jeecra(jexnum(soumat, ibid), 'LONMAX', ntail)
!
        call jecroc(jexnom(repmat, 'K0AA'))
        call jenonu(jexnom(repmat, 'K0AA'), ibid)
        call jeecra(jexnum(soumat, ibid), 'LONMAX', nbdax**2)
!
        call jecroc(jexnom(repmat, 'KPLUSAA'))
        call jenonu(jexnom(repmat, 'KPLUSAA'), ibid)
        call jeecra(jexnum(soumat, ibid), 'LONMAX', nbdax**2)
!
        call jecroc(jexnom(repmat, 'KPLUSAI'))
        call jenonu(jexnom(repmat, 'KPLUSAI'), ibid)
        call jeecra(jexnum(soumat, ibid), 'LONMAX', nbmod*nbdax)
!
        call jecroc(jexnom(repmat, 'KPLUSAJ'))
        call jenonu(jexnom(repmat, 'KPLUSAJ'), ibid)
        call jeecra(jexnum(soumat, ibid), 'LONMAX', ntail)
!
    end if
!
! --- TRAITEMENT DES PRODUITS MODES-MODES
!
    call jenonu(jexnom(repmat, 'K0II'), ibid)
    call jeveuo(jexnum(soumat, ibid), 'E', ldk0ii)
    call jenonu(jexnom(repmat, 'M0II'), ibid)
    call jeveuo(jexnum(soumat, ibid), 'E', ldm0ii)
!
    do i = 1, nbmod
        call rsadpa(basmod, 'L', 1, 'RIGI_GENE', i, &
                    0, sjv=llkge, styp=k8bid)
        zr(ldk0ii+i*(i-1)/2) = zr(llkge)
        call rsadpa(basmod, 'L', 1, 'MASS_GENE', i, &
                    0, sjv=llmge, styp=k8bid)
        zr(ldm0ii+i*(i-1)/2) = zr(llmge)
    end do
!
! --- ALLOCATION DES TABLEAUX DE TRAVAIL
!
    call wkvect('&&'//pgc//'.EXTR.DROI', 'V V I', nbddr, ltexd)
    call wkvect('&&'//pgc//'.EXTR.GAUC', 'V V I', nbddr, ltexg)
    if (nbdax .gt. 0) then
        call wkvect('&&'//pgc//'.EXTR.AXE', 'V V I', nbdax, ltexa)
    end if
!
! --- RECUPERATION NB EQUATIONS
!
    call dismoi('NB_EQUA', raid, 'MATR_ASSE', repi=neq)
    call jeveuo(numddl//'.DEEQ', 'L', vi=deeq)
!
! --- RECUPERATION DES NUMEROS D'INTERFACE DROITE ET GAUCHE
!
    call jeveuo(nomres//'.CYCL_NUIN', 'L', vi=cycl_nuin)
    numd = cycl_nuin(1)
    numg = cycl_nuin(2)
    numa = cycl_nuin(3)
!
! --- RECUPERATION DU NOMBRE DE NOEUDS DES INTERFACES
!
    noeint = intf//'.IDC_LINO'
!
    call jelira(jexnum(noeint, numd), 'LONMAX', nbnod)
    call jeveuo(jexnum(noeint, numd), 'L', llnod)
!
    call jelira(jexnum(noeint, numg), 'LONMAX', nbnog)
    call jeveuo(jexnum(noeint, numg), 'L', llnog)
!
    if (nbdax .gt. 0) then
        call jelira(jexnum(noeint, numa), 'LONMAX', nbnoa)
        call jeveuo(jexnum(noeint, numa), 'L', llnoa)
    end if
!
!
! --- RECUPERATION RANGS DDL INTERFACE DANS VECTEUR ASSEMBLE
!
! --- RECUPERATION  DEFORMEES ET DDL NOEUDS DROITE
!
    kbid = ' '
    call bmradi(basmod, kbid, '        ', numd, nbddr, &
                zi(ltexd), ibid)
!
! --- RECUPERATION DEFORMEES ET DDL NOEUDS GAUCHE
!
    kbid = ' '
    call bmradi(basmod, kbid, '        ', numg, nbddr, &
                zi(ltexg), ibid)
!
! --- RECUPERATION DEFORMEES EVENTUELS NOEUDS AXE
!
    if (nbdax .gt. 0) then
        kbid = ' '
        call bmradi(basmod, kbid, '        ', numa, nbdax, &
                    zi(ltexa), ibid)
    end if
!
! --- CALCUL DES TRANSPOSEES MATRICES TETA DE CHANGEMENT DE REPERE
!
!     ATTENTION LE PARTAGE DES ROUTINES DE CALCUL DE CES MATRICES
!     AVEC CRAIG-BAMPTON IMPOSE DE CALCULER LES TRANPOSEES DES
!     MATRICES TETA QUI APPARAISSENT DANS LA THEORIE DE MAC-NEAL
!
    call wkvect('&&'//pgc//'.TETGD', 'V V R', nbddr*nbddr, ltetgd)
    call ctetgd(basmod, numd, numg, nbsec, zr(ltetgd), &
                nbddr)
    if (nbdax .gt. 0) then
        call wkvect('&&'//pgc//'.TETAX', 'V V R', nbdax*nbdax, ltetax)
        call ctetax(basmod, numa, nbsec, zr(ltetax), nbdax)
    end if
!
    call wkvect('&&'//pgc//'.VECT', 'V V R', neq, ltvec)
!
! --- DETERMINATION DES TABLEAUX EXTRACTION SUR MODES PROPRES
!
! --- POUR K0JI
!
    call jenonu(jexnom(repmat, 'K0JI'), ibid)
    call jeveuo(jexnum(soumat, ibid), 'E', ldk0ji)
    call wkvect('&&'//pgc//'CDPHI', 'V V R', nbmod*nbddr, ltcdp)
!
    do i = 1, nbmod
        call dcapno(basmod, 'DEPL    ', i, chamva)
        call jeveuo(chamva, 'L', llcham)
        b_n = to_blas_int(neq)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(llcham), b_incx, zr(ltvec), b_incy)
        call zerlag(neq, deeq, vectr=zr(ltvec))
!
        do j = 1, nbddr
!
! ------- EXTRACTION DDL DROITE
!
            iran = zi(ltexd+j-1)
            xx = zr(ltvec+iran-1)
            call amppr(zr(ltcdp), nbddr, nbmod, [xx], 1, &
                       1, j, i)
        end do
        call jelibe(chamva)
    end do
!
    b_n = to_blas_int(nbmod*nbddr)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, -1.d0, zr(ltcdp), b_incx, zr(ldk0ji), &
               b_incy)
    call jedetr('&&'//pgc//'CDPHI')
!
! --- POUR KPLUSJI
!
    call jenonu(jexnom(repmat, 'KPLUSJI'), ibid)
    call jeveuo(jexnum(soumat, ibid), 'E', ldkpji)
    call wkvect('&&'//pgc//'CGPHI', 'V V R', nbmod*nbddr, ltcgp)
!
    do i = 1, nbmod
        call dcapno(basmod, 'DEPL    ', i, chamva)
        call jeveuo(chamva, 'L', llcham)
        b_n = to_blas_int(neq)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(llcham), b_incx, zr(ltvec), b_incy)
        call zerlag(neq, deeq, vectr=zr(ltvec))
!
        do j = 1, nbddr
!
! ------- EXTRACTION DDL GAUCHE
!
            iran = zi(ltexg+j-1)
            xx = zr(ltvec+iran-1)
            call amppr(zr(ltcgp), nbddr, nbmod, [xx], 1, &
                       1, j, i)
        end do
        call jelibe(chamva)
    end do
!
    call pmppr(zr(ltetgd), nbddr, nbddr, -1, zr(ltcgp), &
               nbddr, nbmod, 1, zr(ldkpji), nbddr, &
               nbmod)
!
    call jedetr('&&'//pgc//'CGPHI')
!
    if (nbdax .gt. 0) then
!
! ----- POUR K0AA ET KPLUSAA
!
        call wkvect('&&'//pgc//'CAPHI', 'V V R', nbmod*nbdax, ltcap)
        call jenonu(jexnom(repmat, 'K0AI'), ibid)
        call jeveuo(jexnum(soumat, ibid), 'E', ldk0ai)
        call jenonu(jexnom(repmat, 'KPLUSAI'), ibid)
        call jeveuo(jexnum(soumat, ibid), 'E', ldkpai)
!
        do i = 1, nbmod
            call dcapno(basmod, 'DEPL    ', i, chamva)
            call jeveuo(chamva, 'L', llcham)
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(llcham), b_incx, zr(ltvec), b_incy)
            call zerlag(neq, deeq, vectr=zr(ltvec))
!
            do j = 1, nbdax
!
! --------- EXTRACTION DDL AXE
!
                iran = zi(ltexa+j-1)
                xx = zr(ltvec+iran-1)
                call amppr(zr(ltcap), nbdax, nbmod, [xx], 1, &
                           1, j, i)
            end do
            call jelibe(chamva)
            call jedetr('&&'//pgc//'.VECT')
        end do
!
        call pmppr(zr(ltetax), nbdax, nbdax, -1, zr(ltcap), &
                   nbdax, nbmod, 1, zr(ldkpai), nbdax, &
                   nbmod)
        b_n = to_blas_int(nbmod*nbdax)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, -1.d0, zr(ltcap), b_incx, zr(ldk0ai), &
                   b_incy)
        call jedetr('&&'//pgc//'CAPHI')
!
    end if
!
! --- ALLOCATION TABLEAU DE TRAVAIL ET FLEXIBILITE COURANTE
!
    nbv = max(nbddr, nbmod)
    nbv = max(nbdax, nbv)
    ntail = nbv**2
    call wkvect('&&'//pgc//'.MAT.TRAV', 'V V R', ntail, ltmat)
    call wkvect('&&'//pgc//'.FLEX.RES', 'V V R', ntail, ltflex)
!
! --- CALCUL DES TERMES DE FLEXIBILITE
!
! --- POUR K0JJ
!
    call jenonu(jexnom(repmat, 'K0JJ'), ibid)
    call jeveuo(jexnum(soumat, ibid), 'E', ldk0jj)
!
    call flexib(basmod, nbmod, zr(ltflex), nbddr, nbddr, &
                numd, numd)
    b_n = to_blas_int(nbddr*nbddr)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, -1.d0, zr(ltflex), b_incx, zr(ltmat), &
               b_incy)
    k = 0
    do j = 1, nbddr
        do i = j, 1, -1
            zr(ldk0jj+k) = zr(ltmat-1+(j-1)*nbddr+i)
            k = k+1
        end do
    end do
!
    call flexib(basmod, nbmod, zr(ltflex), nbddr, nbddr, &
                numg, numg)
    call pmppr(zr(ltetgd), nbddr, nbddr, -1, zr(ltflex), &
               nbddr, nbddr, 1, zr(ltmat), nbddr, &
               nbddr)
    call pmppr(zr(ltmat), nbddr, nbddr, 1, zr(ltetgd), &
               nbddr, nbddr, 1, zr(ltflex), nbddr, &
               nbddr)
    call r8inir(nbddr*nbddr, 0.d0, zr(ltmat), 1)
    b_n = to_blas_int(nbddr*nbddr)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, -1.d0, zr(ltflex), b_incx, zr(ltmat), &
               b_incy)
    k = 0
    do j = 1, nbddr
        do i = j, 1, -1
            zr(ldk0jj+k) = zr(ldk0jj+k)+zr(ltmat-1+(j-1)*nbddr+i)
            k = k+1
        end do
    end do
!
! --- POUR KPLUSJJ
!
    call jenonu(jexnom(repmat, 'KPLUSJJ'), ibid)
    call jeveuo(jexnum(soumat, ibid), 'E', ldkpjj)
!
    call flexib(basmod, nbmod, zr(ltflex), nbddr, nbddr, &
                numg, numd)
    call pmppr(zr(ltetgd), nbddr, nbddr, -1, zr(ltflex), &
               nbddr, nbddr, 1, zr(ltmat), nbddr, &
               nbddr)
    call amppr(zr(ldkpjj), nbddr, nbddr, zr(ltmat), nbddr, &
               nbddr, 1, 1)
!
!
    if (nbdax .gt. 0) then
!
! --- POUR K0AJ ET KPLUSAJ
!
        call jenonu(jexnom(repmat, 'K0AJ'), ibid)
        call jeveuo(jexnum(soumat, ibid), 'E', ldk0aj)
        call jenonu(jexnom(repmat, 'KPLUSAJ'), ibid)
        call jeveuo(jexnum(soumat, ibid), 'E', ldkpaj)
!
        call flexib(basmod, nbmod, zr(ltflex), nbdax, nbddr, &
                    numa, numd)
        b_n = to_blas_int(nbdax*nbddr)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, -1.d0, zr(ltflex), b_incx, zr(ldk0aj), &
                   b_incy)
        call pmppr(zr(ltetax), nbdax, nbdax, -1, zr(ltflex), &
                   nbdax, nbddr, 1, zr(ltmat), nbdax, &
                   nbddr)
        call amppr(zr(ldkpaj), nbdax, nbddr, zr(ltmat), nbdax, &
                   nbddr, 1, 1)
!
        call flexib(basmod, nbmod, zr(ltflex), nbdax, nbddr, &
                    numa, numg)
        call pmppr(zr(ltflex), nbdax, nbddr, 1, zr(ltetgd), &
                   nbddr, nbddr, 1, zr(ltmat), nbdax, &
                   nbddr)
        call r8inir(nbdax*nbddr, 0.d0, zr(ltflex), 1)
        b_n = to_blas_int(nbdax*nbddr)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, -1.d0, zr(ltmat), b_incx, zr(ltflex), &
                   b_incy)
        call amppr(zr(ldk0aj), nbdax, nbddr, zr(ltflex), nbdax, &
                   nbddr, 1, 1)
!
!
! --- POUR KPLUSJA
!
        call jenonu(jexnom(repmat, 'KPLUSJA'), ibid)
        call jeveuo(jexnum(soumat, ibid), 'E', ldkpja)
!
        call flexib(basmod, nbmod, zr(ltflex), nbddr, nbdax, &
                    numg, numa)
        call pmppr(zr(ltetgd), nbddr, nbddr, -1, zr(ltflex), &
                   nbddr, nbdax, 1, zr(ldkpja), nbddr, &
                   nbdax)
!
!
! --- POUR K0AA ET KPLUSAA
!
        call jenonu(jexnom(repmat, 'K0AA'), ibid)
        call jeveuo(jexnum(soumat, ibid), 'E', ldk0aa)
        call jenonu(jexnom(repmat, 'KPLUSAA'), ibid)
        call jeveuo(jexnum(soumat, ibid), 'E', ldkpaa)
!
        call flexib(basmod, nbmod, zr(ltflex), nbdax, nbdax, &
                    numa, numa)
        b_n = to_blas_int(nbdax*nbdax)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, -1.d0, zr(ltflex), b_incx, zr(ldk0aa), &
                   b_incy)
        call pmppr(zr(ltetax), nbdax, nbdax, -1, zr(ltflex), &
                   nbdax, nbdax, 1, zr(ltmat), nbdax, &
                   nbdax)
        call amppr(zr(ldkpaa), nbdax, nbdax, zr(ltmat), nbdax, &
                   nbdax, 1, 1)
        call pmppr(zr(ltmat), nbdax, nbdax, 1, zr(ltetax), &
                   nbdax, nbdax, 1, zr(ltflex), nbdax, &
                   nbdax)
        call r8inir(nbdax*nbdax, 0.d0, zr(ltmat), 1)
        b_n = to_blas_int(nbdax*nbdax)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, -1.d0, zr(ltflex), b_incx, zr(ltmat), &
                   b_incy)
        call amppr(zr(ldk0aa), nbdax, nbdax, zr(ltmat), nbdax, &
                   nbdax, 1, 1)
!
!
    end if
!
    call jedetr('&&'//pgc//'.MAT.TRAV')
    call jedetr('&&'//pgc//'.TETGD')
    call jedetr('&&'//pgc//'.FLEX.RES')
    call jedetr('&&'//pgc//'.EXTR.DROI')
    call jedetr('&&'//pgc//'.EXTR.GAUC')
    call jedetr('&&'//pgc//'.VECT')
    if (nbdax .gt. 0) then
        call jedetr('&&'//pgc//'.TETAX')
        call jedetr('&&'//pgc//'.EXTR.AXE')
    end if
!
    call jedema()
end subroutine
