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
subroutine asmchc(matas)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/wkvect.h"
!
    character(len=*) :: matas
! OBJET :
!        TRAITEMENT DES CHARGES CINEMATIQUES DANS UNE MATRICE ASSEMBLEE
!        CALCUL DES OBJETS  .CCLL,  .CCVA,  .CCII
!-----------------------------------------------------------------------
! VAR  MATAS   K*19    : NOM DE LA MATR_ASSE
!-----------------------------------------------------------------------
    complex(kind=8) :: dcmplx
!----------------------------------------------------------------------
!     VARIABLES LOCALES
!----------------------------------------------------------------------
    character(len=1) :: base
    character(len=8) :: kbid
    character(len=14) :: nu
    character(len=19) :: mat, nomsto
    integer(kind=8) :: typmat, ielim, jelim, kdeb, kfin, nccva, ilig, jcol
    integer(kind=8) :: jsmhc, jvalm, jvalm2, jccva, jccll, nelim
    integer(kind=8) :: iret2, jnequ, ieq, k, deciel, nterm, neq, ier, imatd
    integer(kind=8) :: nblocm, jccii, decjel, iremp, keta
    aster_logical :: nonsym
    integer(kind=8), pointer :: elim(:) => null()
    integer(kind=8), pointer :: remplis(:) => null()
    integer(kind=8), pointer :: ccid(:) => null()
    character(len=24), pointer :: refa(:) => null()
    integer(kind=8), pointer :: smdi(:) => null()
    integer(kind=8), pointer :: nulg(:) => null()
!----------------------------------------------------------------------
    call jemarq()
    mat = matas
!     CALL CHEKSD('SD_MATR_ASSE',MAT,IRET)
    call jeexin(mat//'.CCVA', ier)
    ASSERT(ier .eq. 0)
    call jeexin(mat//'.CCID', ier)
    if (ier .eq. 0) goto 999
!
    call jelira(mat//'.REFA', 'CLAS', cval=base)
    call jeveuo(mat//'.REFA', 'E', vk24=refa)
    nu = refa(2) (1:14)
    call jeexin(nu//'.NUML.DELG', imatd)
    if (imatd .ne. 0) then
        call jeveuo(nu//'.NUML.NEQU', 'L', jnequ)
        call jeveuo(nu//'.NUML.NULG', 'L', vi=nulg)
    else
        call jeveuo(nu//'.NUME.NEQU', 'L', jnequ)
    end if
    neq = zi(jnequ)
!
!
!     -- ON DETRUIT LES OBJETS S'ILS EXISTENT DEJA :
    ASSERT(refa(3) .ne. 'ELIMF')
    call jedetr(mat//'.CCLL')
    call jedetr(mat//'.CCVA')
    call jedetr(mat//'.CCII')
!
!     -- CALCUL DE ELIM(*) ET NELIM :
!     -----------------------------------
!     ELIM    I(*)    : TABLEAU ENTIER DE DIM = NEQ DONNANT LES
!                       LES NUMEROS DES EQUATIONS A ELIMINER ET LEUR
!                       NUMERO D'ELIMINATION
!                       ZI(KKELI-1+IEQ) = / 0      -> PAS ELIMINE
!                                         / IELIM  -> ELIMINE
!     NELIM   I       : NOMBRE D'EQUATIONS DE LA MATRICE A ELIMINER
    AS_ALLOCATE(vi=elim, size=neq)
    call jeveuo(mat//'.CCID', 'L', vi=ccid)
    nelim = 0
    do ieq = 1, neq
        if (imatd .ne. 0) then
            keta = ccid(nulg(ieq))
        else
            keta = ccid(ieq)
        end if
        ASSERT(keta .eq. 1 .or. keta .eq. 0)
        if (keta .eq. 1) then
            nelim = nelim+1
            elim(ieq) = nelim
        else
            elim(ieq) = 0
        end if
    end do
!
!
!
    if (nelim .eq. 0) then
        refa(3) = 'ELIMF'
        goto 999
    end if
!     -----------------------------------------------------------------
!
    nomsto = nu//'.SMOS'
!
!
    call jeexin(nomsto//'.SMHC', iret2)
    ASSERT(iret2 .gt. 0)
    call jeveuo(nomsto//'.SMHC', 'L', jsmhc)
    call jeveuo(nomsto//'.SMDI', 'L', vi=smdi)
!
!
!
!     -- CALCUL DE .CCLL :
!     -----------------------------------------
    call wkvect(mat//'.CCLL', base//' V I ', 3*nelim, jccll)
!
    kfin = 0
    do jcol = 1, neq
        kdeb = kfin+1
        kfin = smdi(jcol)
        jelim = elim(jcol)
!
        if (jelim .ne. 0) then
            zi(jccll-1+3*(jelim-1)+1) = jcol
            do k = kdeb, kfin-1
                ilig = zi4(jsmhc-1+k)
                ielim = elim(ilig)
                if (ielim .eq. 0) zi(jccll-1+3*(jelim-1)+2) = zi(jccll-1+3*(jelim-1)+2)+1
            end do
!
        else
            do k = kdeb, kfin-1
                ilig = zi4(jsmhc-1+k)
                ielim = elim(ilig)
                if (ielim .ne. 0) zi(jccll-1+3*(ielim-1)+2) = zi(jccll-1+3*(ielim-1)+2)+1
            end do
        end if
    end do
!
!
!     -- CALCUL DE NCCVA ET .CCLL(3*(I-1)+3) :
!     -----------------------------------------
    deciel = 0
    do ielim = 1, nelim
        nterm = zi(jccll-1+3*(ielim-1)+2)
        zi(jccll-1+3*(ielim-1)+3) = deciel
        deciel = deciel+nterm
    end do
    nccva = max(deciel, 1)
!
!
!     -- RECUPERATION DE .VALM
!        CALCUL DE TYPMAT ET NONSYM :
!     ------------------------------------
    call jelira(jexnum(mat//'.VALM', 1), 'TYPE', cval=kbid)
    typmat = 1
    if (kbid(1:1) .eq. 'C') typmat = 2
    nonsym = .false.
    call jelira(mat//'.VALM', 'NMAXOC', nblocm)
    ASSERT(nblocm .eq. 1 .or. nblocm .eq. 2)
    if (nblocm .eq. 2) nonsym = .true.
    call jeveuo(jexnum(mat//'.VALM', 1), 'E', jvalm)
    if (nonsym) call jeveuo(jexnum(mat//'.VALM', 2), 'E', jvalm2)
!
!
!     -- ALLOCATION DE .CCVA ET .CCII :
!     ------------------------------------
    call wkvect(mat//'.CCVA', base//' V '//kbid(1:1), nccva, jccva)
    call wkvect(mat//'.CCII', base//' V I', nccva, jccii)
!
!
!     -- REMPLISSAGE DE .CCII ET .CCVA :
!     -----------------------------------------
!
    AS_ALLOCATE(vi=remplis, size=nelim)
    kfin = 0
    do jcol = 1, neq
        kdeb = kfin+1
        kfin = smdi(jcol)
        jelim = elim(jcol)
!
        if (jelim .ne. 0) then
            deciel = zi(jccll-1+3*(jelim-1)+3)
            do k = kdeb, kfin-1
                ilig = zi4(jsmhc-1+k)
                ielim = elim(ilig)
                if (ielim .eq. 0) then
                    remplis(jelim) = remplis(jelim)+1
                    iremp = remplis(jelim)
                    zi(jccii-1+deciel+iremp) = ilig
                    if (typmat .eq. 1) then
                        zr(jccva-1+deciel+iremp) = zr(jvalm-1+k)
                    else
                        zc(jccva-1+deciel+iremp) = zc(jvalm-1+k)
                    end if
                end if
            end do
!
        else
            do k = kdeb, kfin-1
                ilig = zi4(jsmhc-1+k)
                ielim = elim(ilig)
                decjel = zi(jccll-1+3*(ielim-1)+3)
                if (ielim .ne. 0) then
                    remplis(ielim) = remplis(ielim)+1
                    iremp = remplis(ielim)
                    zi(jccii-1+decjel+iremp) = jcol
                    if (typmat .eq. 1) then
                        if (nonsym) then
                            zr(jccva-1+decjel+iremp) = zr(jvalm2-1+k)
                        else
                            zr(jccva-1+decjel+iremp) = zr(jvalm-1+k)
                        end if
                    else
                        if (nonsym) then
                            zc(jccva-1+decjel+iremp) = zc(jvalm2-1+k)
                        else
                            zc(jccva-1+decjel+iremp) = zc(jvalm-1+k)
                        end if
                    end if
                end if
            end do
        end if
!
    end do
!
!
!---  "SIMPLIFICATION" DE .VALM : 1. SUR LA DIAGONALE ET 0. EN DEHORS
!---------------------------------------------------------------------
!     DANS LE CAS PARALLELE SUR N PROCS, IL EST A NOTER QU'UN N SE
!     TROUVERA SUR LA DIAGONALE ET QU'UN N L'EQUILIBRERA DANS LE
!     SECOND MEMBRE
    kfin = 0
    do jcol = 1, neq
        kdeb = kfin+1
        kfin = smdi(jcol)
        jelim = elim(jcol)
!
        if (jelim .ne. 0) then
            if (typmat .eq. 1) then
                zr(jvalm-1+kfin) = 1.d0
                if (nonsym) zr(jvalm2-1+kfin) = 1.d0
            else
                zc(jvalm-1+kfin) = dcmplx(1.d0, 0.d0)
                if (nonsym) zc(jvalm2-1+kfin) = dcmplx(1.d0, 0.d0)
            end if
        end if
!
        if (jelim .ne. 0) then
            do k = kdeb, kfin-1
                if (typmat .eq. 1) then
                    zr(jvalm-1+k) = 0.d0
                    if (nonsym) zr(jvalm2-1+k) = 0.d0
                else
                    zc(jvalm-1+k) = dcmplx(0.d0, 0.d0)
                    if (nonsym) zc(jvalm2-1+k) = dcmplx(0.d0, 0.d0)
                end if
            end do
!
        else
            do k = kdeb, kfin-1
                ilig = zi4(jsmhc-1+k)
                ielim = elim(ilig)
                if (ielim .ne. 0) then
                    if (typmat .eq. 1) then
                        zr(jvalm-1+k) = 0.d0
                        if (nonsym) zr(jvalm2-1+k) = 0.d0
                    else
                        zc(jvalm-1+k) = dcmplx(0.d0, 0.d0)
                        if (nonsym) zc(jvalm2-1+k) = dcmplx(0.d0, 0.d0)
                    end if
                end if
            end do
        end if
!
    end do
!
    refa(3) = 'ELIMF'
!
!
!
!
999 continue
    AS_DEALLOCATE(vi=remplis)
    AS_DEALLOCATE(vi=elim)
!     CALL CHEKSD('sd_matr_asse',MAT,IRET)
    call jedema()
end subroutine
