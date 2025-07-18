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
subroutine diag99(nomres)
    implicit none
#include "jeveux.h"
#include "asterfort/copisd.h"
#include "asterfort/copmod.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mrmult.h"
#include "asterfort/mtdscr.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rscrsd.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnoch.h"
#include "asterfort/vpgskp.h"
#include "asterfort/vtcrem.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "blas/ddot.h"
    character(len=8) :: nomres
!
!         DEFI_BASE_MODALE : DIAG_MASS
!
! CE MOT CLE PERMET DE DIAGONALISER LA MATRICE DE MASSE EN DEUX ETAPES
!
! 1- RETIRER AUX MODES STATIQUES LEUR CONTRIBUTION SUR LES
!    MODES DYNAMIQUES
!
! 2- ORTHOGONALISER LA FAMILLES DES MODES STATIQUES MODIFIES PAR
!    LE PROCEDE DE GRAAM-SCHMIDT
!
!----------------------------------------------------------------------
!
!
    integer(kind=8) :: iad, jiad, ier, idmode, lmasse, idstat
    integer(kind=8) :: jnsta, i, j, k, ieq, nbord
    integer(kind=8) :: nbmode, nbstat, neq, n1, iorne, iorol
    real(kind=8) :: alpha, r8scal
    complex(kind=8) :: cbid
    character(len=8) :: k8b, meca, stat
    character(len=14) :: nu
    character(len=24) :: masse, numddl, mailla
    character(len=19) :: chamol, chamne
    real(kind=8), pointer :: trav1(:) => null()
    real(kind=8), pointer :: trav2(:) => null()
    real(kind=8), pointer :: trav3(:) => null()
    integer(kind=8), pointer :: trav4(:) => null()
    real(kind=8), pointer :: vale(:) => null()
    integer(kind=8), pointer :: ordm(:) => null()
    integer(kind=8), pointer :: ords(:) => null()
    blas_int :: b_incx, b_incy, b_n
    cbid = dcmplx(0.d0, 0.d0)
!----------------------------------------------------------------------
    call jemarq()
!
!----------------------------------------------------------------------
! --- RECUPERATION DES MODES PROPRES
!-----------------------------------------------------------------------
!
    call getvid('DIAG_MASS', 'MODE_MECA', iocc=1, scal=meca, nbret=n1)
!
    call jelira(meca//'           .ORDR', 'LONUTI', nbmode)
    call jeveuo(meca//'           .ORDR', 'L', vi=ordm)
!
!
    call dismoi('REF_MASS_PREM', nomres, 'RESU_DYNA', repk=masse)
    call dismoi('NUME_DDL', nomres, 'RESU_DYNA', repk=numddl)
!
    nu = numddl(1:14)
!
    call dismoi('NOM_MAILLA', numddl, 'NUME_DDL', repk=mailla)
    call dismoi('NB_EQUA', masse, 'MATR_ASSE', repi=neq)
    call wkvect('&&DIAG99.MODE_MECA', 'V V R', nbmode*neq, idmode)
    call copmod(meca, bmodr=zr(idmode), numer=nu)
!
!-----------------------------------------------------------------------
! --- RECUPERATION DES MODES STATIQUES
!-----------------------------------------------------------------------
!
    call getvid('DIAG_MASS', 'MODE_STAT', iocc=1, scal=stat, nbret=n1)
!
    call jelira(stat//'           .ORDR', 'LONUTI', nbstat)
    call jeveuo(stat//'           .ORDR', 'L', vi=ords)
    call wkvect('&&DIAG99.MODE_STAT', 'V V R', nbstat*neq, idstat)
    call copmod(stat, bmodr=zr(idstat), numer=nu)
!
!-----------------------------------------------------------------------
! --- RECUPERATION DU DESCRIPTEUR DE LA MATRICE DE MASSE
!-----------------------------------------------------------------------
    call mtdscr(masse)
    call jeveuo(masse(1:19)//'.&INT', 'L', lmasse)
!
!-----------------------------------------------------------------------
! 1- RETIRER AUX MODES STATIQUES LEUR CONTRIBUTION AUX MODES PROPRES
! MODE STAT J =
! MODE STAT J - SOMME (T(MODE STAT J)*MASSE*MODE PROPRE I)*MODE PROPRE I
! OU T(VECTEUR) EST LA TRANSPOSEE DU VECTEUR
!-----------------------------------------------------------------------
    call wkvect('&&DIAG99.NEW_STAT', 'V V R', nbstat*neq, jnsta)
    AS_ALLOCATE(vr=trav1, size=neq)
    AS_ALLOCATE(vr=trav2, size=neq)
    AS_ALLOCATE(vr=trav3, size=nbstat)
    AS_ALLOCATE(vi=trav4, size=neq)
!
    do j = 1, nbstat
!
        trav1(:) = 0.d0
!
        do i = 1, nbmode
!
! --------- PRODUIT MASSE*MODE PROPRE I
            call mrmult('ZERO', lmasse, zr(idmode+(i-1)*neq), trav2, 1, &
                        .true._1)
!
! --------- (T(MODE STAT J)*MASSE*MODE PROPRE I)
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            r8scal = ddot(b_n, zr(idstat+(j-1)*neq), b_incx, trav2, b_incy)
!
! --------- PRODUIT (T(MODE STAT J)*MASSE*MODE PROPRE I)*MODE PROPRE I
! --------- PUIS
! --------- SOMME (T(MODE STAT J)*MASSE*MODE PROPRE I)*MODE PROPRE I
            do k = 1, neq
                trav1(1+(k-1)) = trav1(1+(k-1))+r8scal*zr(idmode+(i-1)*neq+(k-1))
            end do
        end do
!
        do k = 1, neq
            zr(jnsta+(j-1)*neq+(k-1)) = zr(idstat+(j-1)*neq+(k-1))-trav1(1+(k-1))
        end do
!
    end do
!
    do i = 1, neq
        trav4(i) = 1
    end do
    alpha = 0.717d0
!
    call vpgskp(neq, nbstat, zr(jnsta), alpha, lmasse, &
                2, trav1, trav4, trav3)
!
    nbord = nbmode+nbstat
    call rscrsd('G', nomres, 'MODE_MECA', nbord)
!
    iorne = 0
    do i = 1, nbmode
        iorol = ordm(i)
        iorne = iorne+1
!
        call rsexch('F', meca, 'DEPL', iorol, chamol, &
                    ier)
        call rsexch(' ', nomres, 'DEPL', iorne, chamne, &
                    ier)
        call copisd('CHAMP', 'G', chamol, chamne)
        call rsnoch(nomres, 'DEPL', iorne)
!
        call rsadpa(meca, 'L', 1, 'NUME_MODE', iorol, &
                    0, sjv=iad, styp=k8b, istop=0)
        call rsadpa(nomres, 'E', 1, 'NUME_MODE', iorne, &
                    0, sjv=jiad, styp=k8b)
        zi(jiad) = zi(iad)
!
        call rsadpa(meca, 'L', 1, 'FREQ', iorol, &
                    0, sjv=iad, styp=k8b, istop=0)
        call rsadpa(nomres, 'E', 1, 'FREQ', iorne, &
                    0, sjv=jiad, styp=k8b)
        zr(jiad) = zr(iad)
!
        call rsadpa(meca, 'L', 1, 'NORME', iorol, &
                    0, sjv=iad, styp=k8b, istop=0)
        call rsadpa(nomres, 'E', 1, 'NORME', iorne, &
                    0, sjv=jiad, styp=k8b)
        zk24(jiad) = zk24(iad)
!
        call rsadpa(meca, 'L', 1, 'OMEGA2', iorol, &
                    0, sjv=iad, styp=k8b, istop=0)
        call rsadpa(nomres, 'E', 1, 'OMEGA2', iorne, &
                    0, sjv=jiad, styp=k8b)
        zr(jiad) = zr(iad)
!
        call rsadpa(meca, 'L', 1, 'MASS_GENE', iorol, &
                    0, sjv=iad, styp=k8b, istop=0)
        call rsadpa(nomres, 'E', 1, 'MASS_GENE', iorne, &
                    0, sjv=jiad, styp=k8b)
        zr(jiad) = zr(iad)
!
        call rsadpa(meca, 'L', 1, 'RIGI_GENE', iorol, &
                    0, sjv=iad, styp=k8b, istop=0)
        call rsadpa(nomres, 'E', 1, 'RIGI_GENE', iorne, &
                    0, sjv=jiad, styp=k8b)
        zr(jiad) = zr(iad)
!
        call rsadpa(meca, 'L', 1, 'TYPE_MODE', iorol, &
                    0, sjv=iad, styp=k8b, istop=0)
        call rsadpa(nomres, 'E', 1, 'TYPE_MODE', iorne, &
                    0, sjv=jiad, styp=k8b)
        zk16(jiad) = zk16(iad)
!
        call rsadpa(meca, 'L', 1, 'MODELE', iorol, &
                    0, sjv=iad, styp=k8b, istop=0)
        call rsadpa(nomres, 'E', 1, 'MODELE', iorne, &
                    0, sjv=jiad, styp=k8b)
        zk8(jiad) = zk8(iad)
!
        call rsadpa(meca, 'L', 1, 'CHAMPMAT', iorol, &
                    0, sjv=iad, styp=k8b, istop=0)
        call rsadpa(nomres, 'E', 1, 'CHAMPMAT', iorne, &
                    0, sjv=jiad, styp=k8b)
        zk8(jiad) = zk8(iad)
!
        call rsadpa(meca, 'L', 1, 'CARAELEM', iorol, &
                    0, sjv=iad, styp=k8b, istop=0)
        call rsadpa(nomres, 'E', 1, 'CARAELEM', iorne, &
                    0, sjv=jiad, styp=k8b)
        zk8(jiad) = zk8(iad)
!
    end do
!
    do i = 1, nbstat
        iorol = ords(i)
        iorne = iorne+1
!
        call rsadpa(nomres, 'E', 1, 'NUME_MODE', iorne, &
                    0, sjv=jiad, styp=k8b)
        zi(jiad) = iorol+nbmode
!
        call rsadpa(nomres, 'E', 1, 'FREQ', iorne, &
                    0, sjv=jiad, styp=k8b)
        zr(jiad) = 0.d0
!
        call rsadpa(nomres, 'E', 1, 'OMEGA2', iorne, &
                    0, sjv=jiad, styp=k8b)
        zr(jiad) = 0.d0
!
        call rsadpa(nomres, 'E', 1, 'MASS_GENE', iorne, &
                    0, sjv=jiad, styp=k8b)
        zr(jiad) = 0.d0
!
        call rsadpa(nomres, 'E', 1, 'RIGI_GENE', iorne, &
                    0, sjv=jiad, styp=k8b)
        zr(jiad) = 0.d0
!
        call rsadpa(stat, 'L', 1, 'NOEUD_CMP', iorol, &
                    0, sjv=iad, styp=k8b, istop=0)
        call rsadpa(nomres, 'E', 1, 'NOEUD_CMP', iorne, &
                    0, sjv=jiad, styp=k8b)
        zk16(jiad) = zk16(iad)
!
        call rsadpa(stat, 'L', 1, 'TYPE_DEFO', iorol, &
                    0, sjv=iad, styp=k8b, istop=0)
        call rsadpa(nomres, 'E', 1, 'TYPE_DEFO', iorne, &
                    0, sjv=jiad, styp=k8b)
        zk16(jiad) = zk16(iad)
!
        call rsadpa(stat, 'L', 1, 'TYPE_MODE', iorol, &
                    0, sjv=iad, styp=k8b, istop=0)
        call rsadpa(nomres, 'E', 1, 'TYPE_MODE', iorne, &
                    0, sjv=jiad, styp=k8b)
        zk16(jiad) = zk16(iad)
!
        call rsadpa(stat, 'L', 1, 'MODELE', iorol, &
                    0, sjv=iad, styp=k8b, istop=0)
        call rsadpa(nomres, 'E', 1, 'MODELE', iorne, &
                    0, sjv=jiad, styp=k8b)
        zk8(jiad) = zk8(iad)
!
        call rsadpa(stat, 'L', 1, 'CHAMPMAT', iorol, &
                    0, sjv=iad, styp=k8b, istop=0)
        call rsadpa(nomres, 'E', 1, 'CHAMPMAT', iorne, &
                    0, sjv=jiad, styp=k8b)
        zk8(jiad) = zk8(iad)
!
        call rsadpa(stat, 'L', 1, 'CARAELEM', iorol, &
                    0, sjv=iad, styp=k8b, istop=0)
        call rsadpa(nomres, 'E', 1, 'CARAELEM', iorne, &
                    0, sjv=jiad, styp=k8b)
        zk8(jiad) = zk8(iad)
!
        call rsexch(' ', nomres, 'DEPL', iorne, chamol, &
                    ier)
        call vtcrem(chamol, masse, 'G', 'R')
        call jeveuo(chamol//'.VALE', 'E', vr=vale)
        do ieq = 1, neq
            vale(ieq) = zr(jnsta+(i-1)*neq+ieq-1)
        end do
        call rsnoch(nomres, 'DEPL', iorne)
    end do
!
    AS_DEALLOCATE(vr=trav1)
    AS_DEALLOCATE(vr=trav2)
    AS_DEALLOCATE(vr=trav3)
    AS_DEALLOCATE(vi=trav4)
    call jedetr('&&DIAG99.TRAV5')
    call jedetr('&&DIAG99.TRAV6')
    call jedetr('&&DIAG99.MODE_MECA')
    call jedetr('&&DIAG99.MODE_STAT')
!
    call jedema()
end subroutine
