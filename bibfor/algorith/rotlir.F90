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
subroutine rotlir(nomres, sst1, intf1, lino1, codret, &
                  indin1, tramo1, ddla1, nbeq1, imast, &
                  numlia)
    implicit none
!    M. CORUS     DATE 02/02/10
!-----------------------------------------------------------------------
!  BUT:      < CALCUL DE LA TRACE DES MODES ORIENTES A L'INTERFACE >
!
!  CALCULER LES NOUVELLES MATRICE REDUITES DE LIAISON EN TENANT COMPTE
!  DE L'ORIENTATION DES SOUS-STRUCTURES.
!  ON DETERMINE LA MATRICE DE LIAISON, LES DIMENSIONS DE CES MATRICES
!  ET LE PRONO ASSOCIE
!
!  VERIFICATION DE LA COHERENCE DES INTERFACE EN VIS-A-VIS
!  GESTION DES LIAISONS INCOMPATIBLES
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! NOMRES   /I/: NOM UTILISATEUR DU RESULTAT
! SST1     /I/: NOM UTILISATEUR DE LA SOUS STRUCTURE
! INTF1    /I/: NOM UTILISATEUR DE L'INTERFACE
! LINO1 /I/: VECTEUR CONTENANT LA LISTE DES NOEUDS DE L'INTERFACE
!                 COURANTE
! CODRET   /I/: ENTIER PERMETTANT DE SAVOIR SI L'AUTRE "COTE" DE LA
!               MEME INTERFACE A DEJA ETE TRAITEE, POUR REORDONNER
!               LES NOEUDS
! NUMLIA   /I/: NUMERO DE LA LIAISON COURANTE
! INDIN1 /I/: VECTEUR CONTENANT LES INDICES ASSOCIES AUX DDL
!                 D'INTERFACE
! TRAMO1  /I/: MATRICE CONTENANT LA TRACE DES MODES ORIENTES
! DDLA1  /O/: NOMBRE DE DDL ACTIFS DE L'INTERFACE
! NBEQ1    /O/: NOMBRE DE MODES DANS LA BASE MODALE
! IMAST    /I/: ENTIER DETERMINANT SI ON CONSTRUIT :
!               -1 / -2 : LA MATRICE DE LIAISON STANDARD (C),
!                1 /  2 : LA MATRICE DE LIAISON PROJETEE (C.PHI)
!
!
!
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/bmnoin.h"
#include "asterfort/codent.h"
#include "asterfort/dismoi.h"
#include "asterfort/isdeco.h"
#include "asterfort/idensd.h"
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
#include "asterfort/rotati.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/int_to_char8.h"
#include "blas/ddot.h"
!
!
!
!   PARAMETER REPRESENTANT LE NOMBRE MAX DE COMPOSANTES DE LA GRANDEUR
!   SOUS-JACENTE TRAITEE
!
    character(len=8) :: nomres
    character(len=1) :: k1bid
    character(len=24) :: int1, indin1, lino1, maint1, restmo, tramo1, ordol, valk(2)
    character(len=19) :: kint, numeq1, numeq2
    character(len=8) :: sst1, intf1, lint1, bamo1, kbid, nmacr1, temp
    character(len=8) :: sst2, mailla, nomnoe
    character(len=4) :: nliai
    integer(kind=8) :: ibid, nbno1, llint1, nbeq1, lmod1, numlia, lonmod, i1, j1, k1, l1
    integer(kind=8) :: m1, n1, lindi1, lnoeu1, nbnoe, nbec, ipos1, ipos2, lmain1, nbcmpm
    integer(kind=8) :: leuler, lresmo, codret, lmacr1, nbddl1, imast, ddla1, lact1, iret
    integer(kind=8) :: lmarot, length, nbcmp, jnocmp, noer
    parameter(nbcmpm=300)
    integer(kind=8) :: deco(nbcmpm)
    real(kind=8) :: euler(3), rota(3, 3), norme, nortot
    logical :: nook
    blas_int :: b_incx, b_incy, b_n
!
!-----------C
!--       --C
!-- DEBUT --C
!--       --C
!-----------C
!
    call jemarq()
!
!-- RECUPERATION DU NOMBRE D'ENTIER CODES POUR LES DDL
    call dismoi('NB_EC', 'DEPL_R', 'GRANDEUR', repi=nbec)
!
!--------------------------------------------------------C
!--                                                    --C
!-- RECUPERATION DE LA TRACE DES MODES SUR L'INTERFACE --C
!--                                                    --C
!--------------------------------------------------------C
!
!-- NOM DE LA BASE MODALE ET NOMBRE DE MODES
    call mgutdm(nomres, sst1, ibid, 'NOM_BASE_MODALE', ibid, &
                bamo1)
    call dismoi('NB_MODES_TOT', bamo1, 'RESULTAT', repi=nbeq1)
!
!-- INTERFACE AMONT DE LA SOUS-STRUCTURE
    call mgutdm(nomres, sst1, ibid, 'NOM_LIST_INTERF', ibid, &
                lint1)
!
!-- NOMBRE DE NOEUDS DE L'INTERFACE
    int1 = lint1//'.IDC_LINO'
    call jenonu(jexnom(int1(1:13)//'NOMS', intf1), ibid)
    call jelira(jexnum(int1, ibid), 'LONMAX', nbno1)
!
!-- nombre max de composantes et noms des composantes
    call jelira(jexnom('&CATA.GD.NOMCMP', 'DEPL_R'), 'LONMAX', nbcmp)
    call jeveuo(jexnom('&CATA.GD.NOMCMP', 'DEPL_R'), 'L', jnocmp)
!
!-- LISTE DES NUMEROS DES NOEUDS DE L'INTERFACE
    call jenonu(jexnom(lint1//'.IDC_NOMS', intf1), ibid)
    call jeveuo(jexnum(lint1//'.IDC_LINO', ibid), 'L', llint1)
!
!-- NUMEROTATION DES NOEUDS DE L'INTERFACE DANS LES MAILLAGES INITAUX
    call wkvect(lino1, 'V V I', nbno1, lnoeu1)
    call bmnoin(bamo1, kbid, intf1, ibid, nbno1, &
                zi(lnoeu1), nbnoe)
!
!-- SI UNE AUTRE INTERFACE A DEJA ETE DEFINIE, ON REORDONNE LA LISTE
!-- COURANTE POUR QUE LES NOEUDS TOMBENT "EN FACE"
!
    if (codret .gt. 0) then
        temp = '&&OP0126'
        call codent(numlia, 'D', nliai)
        ordol = temp//'      .LINO.'//nliai
!
        call jeexin(ordol, iret)
!
        if (iret .ne. 0) then
            call jeveuo(ordol, 'L', l1)
            call wkvect('&&VECTEUR_NOEUDS_TEMP', 'V V I', nbno1, m1)
            call wkvect('&&VECTEUR_INDICES_TEMP', 'V V I', nbno1, n1)
            do i1 = 1, nbno1
                zi(m1+i1-1) = zi(lnoeu1+i1-1)
                do j1 = 1, nbno1
                    if (zi(l1+i1-1) .eq. zi(llint1+j1-1)) then
                        zi(n1+i1-1) = j1
                    end if
                end do
            end do
            do i1 = 1, nbno1
                zi(lnoeu1+i1-1) = zi(m1+zi(n1+i1-1)-1)
            end do
            call jedetr('&&VECTEUR_NOEUDS_TEMP')
            call jedetr('&&VECTEUR_INDICES_TEMP')
        end if
!
    end if
!
!-- RECUPERATION DES INDICES CORRESPONDANT AUX DDL D'INTERFACE
!-- DANS LA NUMEROTATION DES MAILLAGES INITIAUX
    call mgutdm(nomres, sst1, ibid, 'NOM_MACR_ELEM', ibid, &
                nmacr1)
!-- recuperation du nume_equa
    call dismoi('NUME_EQUA', nmacr1, 'NUME_DDL', repk=numeq1)
!
!-- recuperation du maillage
    call dismoi('NOM_MAILLA', nmacr1, 'NUME_DDL', repk=mailla)
!
!-- RECUPERATION DE LA NUMEROTATION DES EQUATIONS
    call jeveuo(jexnum(numeq1//'.PRNO', 1), 'L', lmacr1)
!
!-- REMPLISSAGE DES VECTEURS D'INDICES POUR REPERER LES DDL  D'INTEFACE
    nbddl1 = 6*nbno1
    call jeexin(indin1, iret)
    if (iret .eq. 0) then
        call wkvect(indin1, 'V V I', nbddl1, lindi1)
    else
        call jeveuo(indin1, 'L', lindi1)
    end if
!
!-- ON NE TRAITE QUE LES DDL DX DY DZ DRX DRY ET DRZ
!-- on renvoit des alarmes si d'autres ddls sont présents
    ddla1 = 0
    do i1 = 1, nbno1
        ipos1 = zi(lmacr1+(zi(lnoeu1+i1-1)-1)*(2+nbec))
        call isdeco(zi(lmacr1+(zi(lnoeu1+i1-1)-1)*(2+nbec)+2), deco, nbcmp)
        ipos2 = 0
        do k1 = 1, 6
            if (iret .eq. 0) then
                zi(lindi1+(i1-1)*6+k1-1) = (ipos1+ipos2)*deco(k1)
            end if
            ipos2 = ipos2+deco(k1)
            ddla1 = ddla1+deco(k1)
        end do
        do k1 = 7, nbcmp
            if (deco(k1) .eq. 1) then
                noer = zi(lnoeu1+i1-1)
                nomnoe = int_to_char8(noer)
                valk(1) = zk8(jnocmp-1+k1)
                valk(2) = nomnoe
                call utmess('A', 'ALGORITH12_35', nk=2, valk=valk)
            end if
        end do
    end do
!
!-- ALLOCATION DE LA PLACE POUR LES MATRICES TEMPORAIRES
    maint1 = '&&MATR_TEMP_MO_INT1'
    restmo = '&&VECT_MODE_NOEUD_6DDL'
    call wkvect(restmo, 'V V R', 6, lresmo)
!
!
!-- IL DOIT Y AVOIR MOYEN D'OPTIMISER LE PROCESS
!-- D'EXTRACTION / ROTATION / REMPLISSAGE
!-- MAIS LA, CA MARCHE...
!
!
!--
!-- CAS DE LA MATRICE IDENTITE SEULE, POUR CONSTRUIRE
!-- ENSUITE D'OBSERVATION C AVEC LIPSRB
!--
    if (imast .lt. 0) then
!
!-- DANS CE CAS, LA MATRICE EST UNE "IDENTITE", QUI EST ENSUITE
!-- MULTIPLIEE PAR UNE MATRICE DIAGONALE PAR BLOC, OU CHAQUE BLOC
!-- EST UNE MATRICE DE ROTATION. POUR GAGNER DE LA PLACE ET DU
!-- TEMPS, ON NE STOCKE QU'UN SEUL BLOC
!
        nbeq1 = ddla1
!
        call wkvect(maint1, 'V V R', nbddl1*ddla1, lmain1)
!
        if (ddla1 .eq. nbddl1) then
            do j1 = 1, ddla1
                zr(lmain1+(j1-1)*(ddla1+1)) = 1.0d0
            end do
        end if
        if (2*ddla1 .eq. nbddl1) then
            l1 = 0
            m1 = 0
            do j1 = 1, ddla1
                zr(lmain1+(j1-1)*nbddl1+l1+m1) = 1.0d0
                l1 = l1+1
                if (l1 .eq. 3) then
                    m1 = m1+6
                    l1 = 0
                end if
            end do
        end if
        if ((ddla1-nbddl1)*(2*ddla1-nbddl1) .ne. 0) then
            write (6, *) 'SEULS LES ELEMENTS PORTANT (DX,DY,DZ) OU ', &
                '(DX,DY,DZ,DRX,DRY,DRZ) SONT SUPPORTES'
        end if
!
!-- RECUPERATION DE LA MATRICE DE ROTATION DE LA SST ESCLAVE
        call jeveuo(jexnum(nomres//'      .MODG.LIDF', numlia), 'L', j1)
        if (sst1 .eq. zk8(j1)) then
            sst2 = zk8(j1+2)
        else
            sst2 = zk8(j1)
        end if
!
!-- CALCUL DE MATRICE DE ROTATION POUR LA SOUS STRUCTURE
!--  => IL FAUT CONSTRUIRE LA MATRICE INVERSE
        call jenonu(jexnom(nomres//'      .MODG.SSNO', sst2), ibid)
        call jeveuo(jexnum(nomres//'      .MODG.SSOR', ibid), 'L', leuler)
        do i1 = 1, 3
            euler(i1) = -zr(leuler+i1-1)
        end do
!
        call rotati(euler, rota)
!        CALL WKVECT('&&ROTLIR.MATR_ROTATION','V V R',
!     &               DDLA1*DDLA1,LMAROT)
        call wkvect('&&ROTLIR.MATR_ROTATION', 'V V R', 9, lmarot)
!
!-- REMPLISSAGE DE LA MATRICE DE ROTATION
!        DO 280 I1=1,DDLA1
        do i1 = 1, 3
            j1 = ((i1-1)/3)+1
            l1 = i1-(j1-1)*3
!
!          IBID=LMAROT+(I1-1)*DDLA1
            ibid = lmarot+(i1-1)*3
!
            zr(ibid+(j1-1)*3) = rota(1, l1)
            zr(ibid+(j1-1)*3+1) = rota(2, l1)
            zr(ibid+(j1-1)*3+2) = rota(3, l1)
!
        end do
!
!
!
!
    else if (imast .gt. 0) then
!--
!-- CAS DE LA TRACE DES MODES
!--
!
!-- CALCUL DE MATRICE DE ROTATION POUR LA SOUS STRUCTURE
        call jenonu(jexnom(nomres//'      .MODG.SSNO', sst1), ibid)
        call jeveuo(jexnum(nomres//'      .MODG.SSOR', ibid), 'L', leuler)
        do i1 = 1, 3
            euler(i1) = zr(leuler+i1-1)
        end do
        call rotati(euler, rota)
!
!-- ALLOCATION DE LA PLACE POUR LES MATRICES TEMPORAIRES
        call wkvect(maint1, 'V V R', nbddl1*nbeq1, lmain1)
!
!-- EXTRACTION ET ROTATION DE LA TRACE DES MODES SUR L'INTERFACE
!
        call jeveuo(jexnum(bamo1//'           .TACH', 1), 'L', lmod1)
!
! -- verification que les nume_equa des modes correspondent bien a celui du nume_ddl
        nook = .false.
        do i1 = 1, nbeq1
            kint = zk24(lmod1+i1-1) (1:19)
            call dismoi('NUME_EQUA', kint, 'CHAM_NO', repk=numeq2)
            if (.not. idensd('NUME_EQUA', numeq1, numeq2)) then
                call utmess('E', 'ALGORITH12_36', nk=1, valk=[bamo1], ni=1, &
                            vali=[i1])
                nook = .true.
            end if
        end do
        if (nook) then
            call utmess('F', 'ALGORITH12_37')
        end if
!
        do i1 = 1, nbeq1
            kint = zk24(lmod1+i1-1) (1:19)
            call jeveuo(kint//'.VALE', 'L', ibid)
!
            if (i1 .eq. 1) then
                call jelira(kint//'.VALE', 'LONMAX', lonmod)
            end if
!
            norme = 0.d0
            length = 0
            do j1 = 1, nbno1
!-- REMPLISSAGE TEMPORAIRE DE LA RESTRICTION DU MODE AUX 6 DDL
!-- DU NOEUD J1
                do k1 = 1, 6
                    if (zi(lindi1+(j1-1)*6+k1-1) .gt. 0) then
                        zr(lresmo+k1-1) = zr(ibid+zi(lindi1+(j1-1)*6+k1-1)-1)
                        length = length+1
                        norme = norme+zr(lresmo+k1-1)**2
                    else
                        zr(lresmo+k1-1) = 0.d0
                    end if
                end do
!
!-- ROTATION DU VECTEUR RESTRICTION
                do k1 = 1, 6
                    l1 = int(mod(k1-1, 3)+1)
                    m1 = int(int((k1-1)/3)*3)
                    zr(lmain1+(i1-1)*nbddl1+(j1-1)*6+k1-1) = rota(l1, 1)*zr(lresmo+m1)+rota(l1, 2&
                                                             &)*zr(lresmo+m1+1)+rota(l1, 3)*zr(l&
                                                             &resmo+m1+2)
                end do
            end do
!
!-- ON ANNULE BRUTALEMENT LES DEPLACEMENTS DES MODES
!--   DONT ON PENSE QUE L'INTERFACE EST FIXE
!-- SINON, ON GENERE DES PATHOLOGIE DANS LA RECHERCHE
!--   DES RELATIONS INDEPENDANTES
!
            b_n = to_blas_int(lonmod)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            nortot = sqrt(ddot(b_n, zr(ibid-1), b_incx, zr(ibid-1), b_incy))
!
            if (sqrt(norme)/length .lt. 100.d0*r8prem()*nortot) then
                do k1 = 1, nbddl1
                    zr(lmain1+(i1-1)*nbddl1+k1-1) = 0.d0
                end do
            end if
!
        end do
!
    else
!-- POUR EVITER DES CALCULS
        goto 999
    end if
!
!
!-- TRI DES MATRICES POUR ELIMINER LES LIGNES DES DDL NON REPRESENTES
!
    call wkvect(tramo1, 'V V R', ddla1*nbeq1, lact1)
    call codent(numlia, 'D0', k1bid)
!
    do j1 = 1, nbeq1
        ibid = 0
        do i1 = 1, nbddl1
            if (zi(lindi1+i1-1) .gt. 0) then
                zr(lact1+ddla1*(j1-1)+ibid) = zr(lmain1+nbddl1*(j1-1)+i1-1)
                ibid = ibid+1
            end if
        end do
    end do
!
!
!
    call jedetr(maint1)
!---------C
!--     --C
!-- FIN --C
!--     --C
!---------C
!
999 continue
    call jedetr(restmo)
    call jedema()
end subroutine
