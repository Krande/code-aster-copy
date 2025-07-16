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
subroutine resu74(tran, nomres)
!
    use DynaGene_module
!
    implicit none
!
!     CETTE ROUTINE PERMET LA CONCATENATION DE DEUX CONCEPTS TRAN_GENE
!     CALCULES PAR DEUX COMMANDE DYNA_VIBRA//TRAN/GENE
!     TRAN TRONQUE ET NOMRES SONT COPIES DANS TRAN
! ----------------------------------------------------------------------
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jedupo.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/juveca.h"
#include "asterfort/wkvect.h"
#include "blas/dcopy.h"
    character(len=8) :: nomres, tran
!
! IN  : TRAN : PREMIER CONCEPT TRAN_GENE
! IN  : NOMRES : SECOND CONCEPT TRAN_GENE
!
!
!
!
    integer(kind=8) :: nbmode, nc, np, ni, nbsto1, nbinst
    integer(kind=8) :: nbnoli, iret
    real(kind=8) :: prec, tinit, prec2
    character(len=8) :: resu, crit
!     ------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, nbvint, jdesc
    real(kind=8), pointer :: inst1(:) => null()
    integer(kind=8), pointer :: ordr1(:) => null()
    real(kind=8), pointer :: v_bloc(:) => null()
    integer(kind=8), pointer :: v_blo2(:) => null()
    character(len=16) :: tran16, nomres16
    integer(kind=8) :: shift, last_bloc, i_bloc, n_bloc
    character(len=7) :: intk7
    type(DynaGene) :: dyna_gene
!
!
!
!-----------------------------------------------------------------------
    call jemarq()
    resu = '&&TMPRES'
!
    call jeveuo(tran//'           .DESC', 'E', jdesc)
    nbmode = zi(jdesc+1)
    ASSERT(nbmode .gt. 0)
    nbnoli = zi(jdesc+2)
    if (nbnoli .gt. 0) then
        nbvint = zi(jdesc+3)
    end if
!
!     --- RECUPERATION DE L'INSTANT DE REPRISE
!
    call getvtx('ETAT_INIT', 'CRITERE', iocc=1, scal=crit, nbret=nc)
    call getvr8('ETAT_INIT', 'PRECISION', iocc=1, scal=prec, nbret=np)
    if (nc .eq. 0) crit = 'RELATIF'
    if (np .eq. 0) prec = 1.d-6
!
!
!      --- RECHERCHE DU NUMERO D'ORDRE DE L'INSTANT DE REPRISE
!
!
    call getvr8('ETAT_INIT', 'INST_INIT', iocc=1, scal=tinit, nbret=ni)
!
    call dyna_gene%init(tran)
!
    if (ni .eq. 0) then
        call dyna_gene%get_values(dyna_gene%disc, dyna_gene%n_bloc, shift, nbsto1, vr=inst1)
        tinit = inst1(nbsto1)
    else
        call dyna_gene%get_values_by_disc(dyna_gene%disc, tinit, shift, nbsto1, vr=inst1)
    end if
!
    prec2 = prec
    if (crit(1:7) .eq. 'RELATIF') prec2 = prec*inst1(1)
    if (abs(tinit-inst1(1)) .le. prec2) then
        nbinst = 1+shift
! on prend la dernière valeur du bloc précédent
        call dyna_gene%get_values_by_index(dyna_gene%disc, nbinst, shift, vr=inst1)
        goto 202
    end if
    if (crit(1:7) .eq. 'RELATIF') prec2 = prec*inst1(nbsto1)
    if (abs(tinit-inst1(nbsto1)) .le. prec2) then
        nbinst = nbsto1+shift
        goto 202
    end if
    do i = 2, nbsto1-1
        if (crit(1:7) .eq. 'RELATIF') prec2 = prec*inst1(i)
        if (abs(tinit-inst1(i)) .le. prec2) then
            nbinst = i+shift
            goto 202
        end if
    end do
!
!   tran is thus to be truncated from i=1 up to i=nbinst
!   its total size is nbsto1
!
202 continue
!
!
!     --- Retrieval of DEPL, VITE, and ACCE fields ---
!     Note : fortran pointers are not used for nomres and resu fields since the
!            blas function dcopy is used to copy part of the corresponding variables
!            jeveux pointers referrals with zr(*+....) in dcopy are possible
!
!
! reduction taille du bloc de l'instant de reprise
!
    if (dyna_gene%n_bloc .eq. 0) then
        last_bloc = 1
        tran16 = tran//'        '
    else
        call dyna_gene%get_current_bloc(dyna_gene%disc, last_bloc)
        call codent(last_bloc, 'D0', intk7)
        tran16 = tran//'.'//intk7
    end if
!
    call juveca(tran16//'   .ORDR', (nbinst-shift))
    call juveca(tran16//'   .DISC', (nbinst-shift))
    call juveca(tran16//'   .PTEM', (nbinst-shift))
    call juveca(tran16//'   .DEPL', (nbinst-shift)*nbmode)
    call juveca(tran16//'   .VITE', (nbinst-shift)*nbmode)
    call juveca(tran16//'   .ACCE', (nbinst-shift)*nbmode)
    if (nbnoli .gt. 0) then
        call juveca(tran16//'.NL.VINT', (nbinst-shift)*nbvint)
    end if
!
    if (dyna_gene%n_bloc .eq. 0) then
! si 1 seul bloc, renomme '       ' -> '0000001'
        call dyna_gene%free()
        call codent(1, 'D0', intk7)
        tran16 = tran//'.'//intk7
        call jedupo(tran//'           .DEPL', 'G', tran16//'   .DEPL', .false._1)
        call jedetr(tran//'           .DEPL')
        call jedupo(tran//'           .VITE', 'G', tran16//'   .VITE', .false._1)
        call jedetr(tran//'           .VITE')
        call jedupo(tran//'           .ACCE', 'G', tran16//'   .ACCE', .false._1)
        call jedetr(tran//'           .ACCE')
        call jedupo(tran//'           .ORDR', 'G', tran16//'   .ORDR', .false._1)
        call jedetr(tran//'           .ORDR')
        call jedupo(tran//'           .DISC', 'G', tran16//'   .DISC', .false._1)
        call jedetr(tran//'           .DISC')
        call jedupo(tran//'           .PTEM', 'G', tran16//'   .PTEM', .false._1)
        call jedetr(tran//'           .PTEM')
!
        if (nbnoli .ne. 0) then
            call jedupo(tran//'        .NL.VINT', 'G', tran16//'.NL.VINT', .false._1)
            call jedetr(tran//'        .NL.VINT')
        end if
    else
! sinon suppression des blocs après l'instant de reprise
        call dyna_gene%get_current_bloc(dyna_gene%disc, last_bloc)
!
        do i_bloc = last_bloc+1, dyna_gene%n_bloc
            call codent(i_bloc, 'D0', intk7)
            tran16 = tran//'.'//intk7
            call jedetr(tran16//'   .ORDR')
            call jedetr(tran16//'   .DISC')
            call jedetr(tran16//'   .PTEM')
            call jedetr(tran16//'   .DEPL')
            call jedetr(tran16//'   .VITE')
            call jedetr(tran16//'   .ACCE')
            if (nbnoli .gt. 0) then
                call jedetr(tran16//'.NL.VINT')
            end if
        end do
        call dyna_gene%free()
    end if
!
!   déplacement les blocs de nomres à la suite de tran
!
    call jeexin(nomres//'           .BLOC', iret)
    if (iret .eq. 0) then
!       si 1 seul bloc, renomme '       ' -> '%07i'
        nomres16 = nomres//'        '
        n_bloc = 1
        call codent(last_bloc+n_bloc, 'D0', intk7)
        tran16 = tran//'.'//intk7
        call jedupo(nomres16//'   .DEPL', 'G', tran16//'   .DEPL', .false._1)
        call jedetr(nomres16//'   .DEPL')
        call jedupo(nomres16//'   .VITE', 'G', tran16//'   .VITE', .false._1)
        call jedetr(nomres16//'   .VITE')
        call jedupo(nomres16//'   .ACCE', 'G', tran16//'   .ACCE', .false._1)
        call jedetr(nomres16//'   .ACCE')
        call jedupo(nomres16//'   .ORDR', 'G', tran16//'   .ORDR', .false._1)
        call jedetr(nomres16//'   .ORDR')
        call jedupo(nomres16//'   .DISC', 'G', tran16//'   .DISC', .false._1)
        call jedetr(nomres16//'   .DISC')
        call jedupo(nomres16//'   .PTEM', 'G', tran16//'   .PTEM', .false._1)
        call jedetr(nomres16//'   .PTEM')
        if (nbnoli .ne. 0) then
            call jedupo(nomres16//'.NL.VINT', 'G', tran16//'.NL.VINT', .false._1)
            call jedetr(nomres16//'.NL.VINT')
        end if
    else
        call jelira(nomres//'           .BLOC', 'LONMAX', n_bloc)
        do i_bloc = 1, n_bloc
            call codent(i_bloc, 'D0', intk7)
            nomres16 = nomres//'.'//intk7
            call codent(last_bloc+i_bloc, 'D0', intk7)
            tran16 = tran//'.'//intk7
            call jedupo(nomres16//'   .DEPL', 'G', tran16//'   .DEPL', .false._1)
            call jedetr(nomres16//'   .DEPL')
            call jedupo(nomres16//'   .VITE', 'G', tran16//'   .VITE', .false._1)
            call jedetr(nomres16//'   .VITE')
            call jedupo(nomres16//'   .ACCE', 'G', tran16//'   .ACCE', .false._1)
            call jedetr(nomres16//'   .ACCE')
            call jedupo(nomres16//'   .ORDR', 'G', tran16//'   .ORDR', .false._1)
            call jedetr(nomres16//'   .ORDR')
            call jedupo(nomres16//'   .DISC', 'G', tran16//'   .DISC', .false._1)
            call jedetr(nomres16//'   .DISC')
            call jedupo(nomres16//'   .PTEM', 'G', tran16//'   .PTEM', .false._1)
            call jedetr(nomres16//'   .PTEM')
            if (nbnoli .ne. 0) then
                call jedupo(nomres16//'.NL.VINT', 'G', tran16//'.NL.VINT', .false._1)
                call jedetr(nomres16//'.NL.VINT')
            end if
        end do
    end if
!
!   mise à jour SD .BLOC
!
    call jeexin(tran//'           .BLOC', iret)
    if (iret .eq. 0) then
        call wkvect(tran//'           .BLOC', 'G V R', last_bloc+n_bloc, vr=v_bloc)
        call wkvect(tran//'           .BLO2', 'G V I', last_bloc+n_bloc, vi=v_blo2)
    else
        call juveca(tran//'           .BLOC', last_bloc+n_bloc)
        call jeveuo(tran//'           .BLOC', 'E', vr=v_bloc)
        call juveca(tran//'           .BLO2', last_bloc+n_bloc)
        call jeveuo(tran//'           .BLO2', 'E', vi=v_blo2)
    end if
!
!   redefinition taille de last_bloc
    call codent(last_bloc, 'D0', intk7)
    tran16 = tran//'.'//intk7
    call jeveuo(tran16//'   .DISC', 'L', vr=inst1)
    call jelira(tran16//'   .DISC', 'LONMAX', nbinst)
    call jeveuo(tran16//'   .ORDR', 'L', vi=ordr1)
    v_bloc(last_bloc) = inst1(nbinst)
    v_blo2(last_bloc) = ordr1(nbinst)
!   ajout blocs suivants
    if (n_bloc .eq. 1) then
        call codent(last_bloc+n_bloc, 'D0', intk7)
        tran16 = tran//'.'//intk7
        call jeveuo(tran16//'   .DISC', 'L', vr=inst1)
        call jelira(tran16//'   .DISC', 'LONMAX', nbinst)
        call jeveuo(tran16//'   .ORDR', 'L', vi=ordr1)
        v_bloc(last_bloc+n_bloc) = inst1(nbinst)
        v_blo2(last_bloc+n_bloc) = ordr1(nbinst)
    else
        call jeveuo(nomres//'           .BLOC', 'L', vr=inst1)
        call jeveuo(nomres//'           .BLO2', 'L', vi=ordr1)
        do i_bloc = 1, n_bloc
            v_bloc(last_bloc+i_bloc) = inst1(i_bloc)
            v_blo2(last_bloc+i_bloc) = ordr1(i_bloc)
        end do
    end if
!
!   menage nomres
!
    call jeexin(nomres//'           .BLOC', iret)
    if (iret .ne. 0) then
        call jedetr(nomres//'           .BLOC')
        call jedetr(nomres//'           .BLO2')
    end if
!
    if (nbnoli .ne. 0) then
        call jedetr(nomres//'        .NL.TYPE')
        call jedetr(nomres//'        .NL.INTI')
        call jedetr(nomres//'        .NL.VIND')
    end if
!
    call jedetr(nomres//'           .DESC')
!
!
!   --- Further cleanup
    call jeexin(nomres//'           .FDEP', iret)
    if (iret .ne. 0) then
        call jedetr(nomres//'           .FDEP')
        call jedetr(nomres//'           .FVIT')
        call jedetr(nomres//'           .FACC')
    end if
!
    call jedema()
!
end subroutine
