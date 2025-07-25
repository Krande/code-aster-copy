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
subroutine mamodg(model, stolci, nomres, itxsto, itysto, &
                  itzsto, iprsto, iadirg, nbmo, max, &
                  may, maz, nbloc)
    implicit none
! ROUTINE NOUVEAU MODELE OPTIMISEE
! ROUTINE CALCULANT LA MASSE AJOUTEE SUR MODELE GENERALISE
! ARGUMENTS :
! IN : STOLCI : K19 : NOM CONCERNANT LES SD NUMEDDLGENE
! IN : NOMRES : K8 :NOM UTILISATEUR DU RESULTAT
! IN : MODEL : K2 : CHARACTER DISTINGUANT LE FLUIDE 2D ET 3D
! IN : MAX, MAY,MAZ : K19 : MATRICES AX, AY, AZ CALCULEES SUR
!                          L INTERFACE
! IN : ITXSTO,ITYSTO,ITZSTO,IPRSTO : ADR JEVEUX DES NOMS DES
!      CHAMPS DE DEPL_R STOCKEES PAR CMP ET DE LA PRESSION
!      CALCULEE SUR TOUS LES MODES
! IN : IADIRG : ADRESSE DU PREMIER ELEMENT D UN TABLEAU CONTENANT
!       LES RANGS GENERALISES DU COEFF DE MASSE AJOUTEE
! IN : NBMO : NOMBRE DE MODES TOTAL DANS LA BASE MODALE DES
!            SOUS-STRUCTURES - DEFORMEES STATIQUES + MODES
!             NORMAUX
!---------------------------------------------------------------------
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/getvtx.h"
#include "asterfort/infniv.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/mrmult.h"
#include "asterfort/mtdscr.h"
#include "blas/ddot.h"
!
    integer(kind=8) :: nbpres, imatx, imaty, itxsto, itysto, itzsto
    integer(kind=8) :: iprsto, imatz
    integer(kind=8) :: irang, jrang, i, j, iblo, ldblo, iadirg
    integer(kind=8) :: iblodi, nbloc, n1bloc, n2bloc, nbmo, nn
    integer(kind=8) :: ifm, niv, iret1, hc
    real(kind=8) :: mij, rx, ry, rz
    character(len=2) :: model
    character(len=8) :: repon
    character(len=8) :: nomres
    character(len=19) :: max, may, maz, stolci
    real(kind=8), pointer :: vectx(:) => null()
    real(kind=8), pointer :: vecty(:) => null()
    real(kind=8), pointer :: vectz(:) => null()
    integer(kind=8), pointer :: smdi(:) => null()
    integer(kind=8), pointer :: smde(:) => null()
    integer(kind=8), pointer :: indic(:) => null()
    integer(kind=4), pointer :: smhc(:) => null()
    real(kind=8), pointer :: pres(:) => null()
    real(kind=8), pointer :: tpx(:) => null()
    real(kind=8), pointer :: tpy(:) => null()
    real(kind=8), pointer :: tpz(:) => null()
    blas_int :: b_incx, b_incy, b_n
! ------------------------------------------------------------------
!----- ICI ON CALCULE LA MASSE AJOUTEE SUR UN MODELE GENERALISE ---
!
    call jemarq()
!
    call infniv(ifm, niv)
    call getvtx(' ', 'AVEC_MODE_STAT', scal=repon, nbret=nn)
    if (repon(1:3) .eq. 'NON') call jeveuo('&&DELAT.INDIC', 'L', vi=indic)
!
    call jeexin(stolci//'.SMHC', iret1)
    ASSERT(iret1 .gt. 0)
    call jeveuo(stolci//'.SMHC', 'L', vi4=smhc)
    call jeveuo(stolci//'.SMDI', 'L', vi=smdi)
    call jeveuo(stolci//'.SMDE', 'L', vi=smde)
!
    call jelira(zk24(iprsto) (1:19)//'.VALE', 'LONMAX', nbpres)
    AS_ALLOCATE(vr=vectx, size=nbpres)
    AS_ALLOCATE(vr=vecty, size=nbpres)
!
! --- RECUPERATION DES DESCRIPTEURS DE MATRICES ASSEMBLEES MAX ET MAY
! --- EVENTUELLEMENT MAZ
!
    call mtdscr(max)
    call jeveuo(max(1:19)//'.&INT', 'E', imatx)
    call mtdscr(may)
    call jeveuo(may(1:19)//'.&INT', 'E', imaty)
    if (model .eq. '3D') then
        call mtdscr(maz)
        call jeveuo(maz(1:19)//'.&INT', 'E', imatz)
        AS_ALLOCATE(vr=vectz, size=nbpres)
    end if
!
!     BOUCLE SUR LES BLOCS DE LA MATRICE ASSEMBLEE GENE
!
    do iblo = 1, nbloc
!
        call jecroc(jexnum(nomres//'           .UALF', iblo))
        call jeveuo(jexnum(nomres//'           .UALF', iblo), 'E', ldblo)
!-------------------------------------------------------------
!
!         BOUCLE SUR LES COLONNES DE LA MATRICE ASSEMBLEE
!
        n1bloc = 1
        n2bloc = smde(1)
!
!
        do i = n1bloc, n2bloc
            if (i .gt. nbmo) goto 10
            if (repon(1:3) .eq. 'NON') then
                if (indic(i) .ne. 1) goto 10
            end if
            call jeveuo(zk24(itxsto+i-1) (1:19)//'.VALE', 'L', vr=tpx)
            call jeveuo(zk24(itysto+i-1) (1:19)//'.VALE', 'L', vr=tpy)
            if (model .eq. '3D') then
                call jeveuo(zk24(itzsto+i-1) (1:19)//'.VALE', 'L', vr=tpz)
                call mrmult('ZERO', imatz, tpz, vectz, 1, &
                            .true._1)
            end if
!
!------MULTIPLICATIONS MATRICE MAX * CHAMNO MODX---------------------
!----------ET MATRICE MAY * CHAMNO MODY------------------------------
!
            call mrmult('ZERO', imatx, tpx, vectx, 1, &
                        .true._1)
            call mrmult('ZERO', imaty, tpy, vecty, 1, &
                        .true._1)
!
! RANG GENERALISE DU TERME DE MASSE CALCULEE : LIGNE
!
            irang = zi(iadirg+i-1)
            hc = smdi(i)
            if (i .gt. 1) hc = hc-smdi(i-1)
!
            do j = (i-hc+1), i
!
!----------------------------------------------------------------
! ICI ON CALCULE LA MASSE AJOUTEE SUR UN MODELE GENERALISE
!--------------------------------------------------------------
!-----------STOCKAGE DANS LA MATR_ASSE_GENE  ------
!
                if (repon(1:3) .eq. 'NON') then
                    if (indic(j) .ne. 1) goto 50
                end if
!
                call jeveuo(zk24(iprsto+j-1) (1:19)//'.VALE', 'L', vr=pres)
!
                b_n = to_blas_int(nbpres)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                rx = ddot(b_n, pres, b_incx, vectx, b_incy)
                b_n = to_blas_int(nbpres)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                ry = ddot(b_n, pres, b_incx, vecty, b_incy)
!
                if (model .eq. '3D') then
                    b_n = to_blas_int(nbpres)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    rz = ddot(b_n, pres, b_incx, vectz, b_incy)
                    mij = rx+ry+rz
                else
                    mij = rx+ry
                end if
50              continue
                if (repon(1:3) .eq. 'NON') then
                    if (indic(j) .ne. 1) mij = 0.d0
                end if
!
! RANG GENERALISE DU TERME DE MASSE: COLONNE
!
                jrang = zi(iadirg+j-1)
                iblodi = 1
!
                if (iblodi .ne. iblo) then
!
!                 CAS OU LE BLOC COURANT N EST PAS LE BON
!
                    call jelibe(jexnum(nomres//'           .UALF', iblo))
                    call jeveuo(jexnum(nomres//'           .UALF', iblodi), 'E', ldblo)
                    zr(ldblo+smdi(irang)+jrang-irang-1) = mij
                    if (niv .eq. 2) then
                        write (ifm, 350) irang, jrang, mij
                    end if
                    call jelibe(jexnum(nomres//'           .UALF', iblodi))
                    call jeveuo(jexnum(nomres//'           .UALF', iblo), 'E', ldblo)
!
                else
                    zr(ldblo+smdi(irang)+jrang-irang-1) = mij
                    if (niv .eq. 2) then
                        write (ifm, 350) irang, jrang, mij
                    end if
                end if
            end do
10          continue
        end do
    end do
!
350 format(18x, 'M', 2 i 4, 1x, '=', 1x, d 12.5)
!
!
!--MENAGE FINAL DES OBJETS DE TRAVAIL
!
    AS_DEALLOCATE(vr=vectz)
    AS_DEALLOCATE(vr=vectx)
    AS_DEALLOCATE(vr=vecty)
!
    call jedema()
end subroutine
