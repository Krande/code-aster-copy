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

subroutine mdallr(resu1, resu2, basemo, nbmode, nbsauv, &
                  vecpr8, vecpc8, zcmplx)
!
!     ALLOCATION DES VECTEURS DE SORTIE (DONNEEES MODALES REELLES)
!     ------------------------------------------------------------------
! IN  : NOMRES : NOM DU CONCEPT RESULTAT
! IN  : NBMODE : NOMBRE DE MODES
! IN  : NBSAUV : NOMBRE DE PAS CALCULE (INITIAL COMPRIS)
! IN  : DATAx  : DONNEES MODALES AU COMPLET (x=I POUR ENTIER, x=K POUR
!                CHAR, x=R POUR REEL)
! ----------------------------------------------------------------------
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/detrsd.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/juveca.h"
#include "asterfort/nummo1.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnoch.h"
#include "asterfort/vpcrea.h"
#include "asterfort/vtcrem.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: nbmode, nbsauv, imode, ier, lvale, i, jrefa
    integer(kind=8) :: jdesc, igd, iarg, jrefe
    aster_logical :: lrefe, zcmplx
    character(len=8) :: resu1, resu2, matgen, k8b, basemo, typ
    character(len=14) :: nugene
    character(len=19) :: chamge
    real(kind=8) :: vecpr8(nbmode, *)
    complex(kind=8) :: vecpc8(nbmode, *)
!
    integer(kind=8) :: nbmax, ipar, ipar1, ipar2
    parameter(nbmax=50)
    character(len=24) :: kpar(nbmax)
!
    call jemarq()
!
    lrefe = .true.
    nugene = resu2//'.NUGEN'
    matgen = '&&MDALMA'
!
! CREATION DE LA NUMEROTATION GENERALISE SUPPORT
    call nummo1(nugene, basemo, nbmode, 'PLEIN')
!
! CREATION DE LA MATRICE GENERALISE SUPPORT
    call wkvect(matgen//'           .REFA', 'V V K24', 20, jrefa)
    zk24(jrefa-1+11) = 'MPI_COMPLET'
    zk24(jrefa-1+1) = basemo
    zk24(jrefa-1+2) = nugene
    zk24(jrefa-1+9) = 'MS'
    zk24(jrefa-1+10) = 'GENE'
!
! recuperation des parametres a garder dans le modele gene
    call getvtx(' ', 'NOM_PARA', nbval=nbmax, vect=kpar, nbret=ipar)
!
    do imode = 1, nbsauv
!        --- VECTEUR PROPRE ---
        call rsexch(' ', resu2, 'DEPL', imode, chamge, &
                    ier)
        if (ier .eq. 0) then
        else if (ier .eq. 100 .and. lrefe) then
            if (.not. zcmplx) then
                call vtcrem(chamge, matgen, 'G', 'R')
            else
                call vtcrem(chamge, matgen, 'G', 'C')
            end if
            ! GLUTE CAR ON A UTILISE VTCRE[ABM] POUR UN CHAM_GENE QUI A UN .REFE
            ! DE TAILLE 2 ET NON 4 COMME UN CHAM_NO ET PAS DE .DESC
            call wkvect(chamge//'.DESC', 'G V I', 2, jdesc)
            call dismoi("NUM_GD_SI", nugene, "NUME_DDL", repi=igd)
            zi(jdesc-1+1) = igd
            zi(jdesc-1+2) = 1
            call jeecra(chamge//'.DESC', 'DOCU', iarg, 'VGEN')
            call jeecra(chamge//'.REFE', 'DOCU', iarg, 'VGEN')
            call juveca(chamge//'.REFE', 2)

            call jeveuo(chamge//'.REFE', 'E', jrefe)
            zk24(jrefe) = basemo
            zk24(jrefe+1) = nugene
        else
            ASSERT(.false.)
        end if
        call jeecra(chamge//'.REFE', 'DOCU', cval='VGEN')
        call jeveuo(chamge//'.VALE', 'E', lvale)
        do ier = 1, nbmode
            if (.not. zcmplx) then
                zr(lvale+ier-1) = vecpr8(ier, imode)
            else
                zc(lvale+ier-1) = vecpc8(ier, imode)
            end if
        end do
        call rsnoch(resu2, 'DEPL', imode)
!
        do i = 1, ipar
            call rsadpa(resu1, 'L', 1, kpar(i), imode, &
                        1, sjv=ipar1, styp=typ, istop=0)
            call rsadpa(resu2, 'E', 1, kpar(i), imode, &
                        0, sjv=ipar2, styp=k8b)
            if (typ(1:1) .eq. 'I') then
                zi(ipar2) = zi(ipar1)
            else if (typ(1:1) .eq. 'R') then
                zr(ipar2) = zr(ipar1)
            else if (typ(1:2) .eq. 'K8') then
                zk8(ipar2) = zk8(ipar1)
            else if (typ(1:3) .eq. 'K16') then
                zk16(ipar2) = zk16(ipar1)
            else if (typ(1:3) .eq. 'K32') then
                zk32(ipar2) = zk32(ipar1)
            end if
        end do
    end do
!
    call vpcrea(0, resu2, ' ', ' ', ' ', &
                ' ', ier)
!
! --- MENAGE
    call detrsd('MATR_ASSE_GENE', matgen)
!
    call jedema()
!
end subroutine
