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

subroutine ualfcr(mataz, basz)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jevtbl.h"
#include "asterfort/jexnum.h"
#include "asterfort/smosli.h"
#include "asterfort/utmess.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
!
    character(len=*) :: mataz, basz
!     Creation de l'objet mataz.UALF pour contenir la factorisee ldlt
!     de la matrice mataz
!     rq : cette routine cree (si necessaire) le stockage morse de mataz
!     ------------------------------------------------------------------
! in  jxvar k19 mataz     : nom d'une s.d. matr_asse
! in        k1  basz      : base de creation pour .UALF
!                  si basz=' ' on prend la meme base que celle de .VALM
!
!-----------------------------------------------------------------------

    character(len=1) :: base, tyrc, basto
    character(len=4) :: kmpic
    character(len=14) :: nu
    character(len=19) :: stomor, stolci, matas
    integer(kind=8) ::  neq, nbloc, nblocm, nbloclc
    integer(kind=8) :: jsmhc, jschc, ico
    integer(kind=8) :: itbloc, jcollc, jcolm, ibloc, jualf, jvalm
    integer(kind=8) :: kterm, nbterm, iligm, iliglc, k1, n1, ii1, ii2
    integer(kind=8) :: ismdi, ismdim1, iscdi, kblocm, iret
    integer(kind=8) :: iliglc2, jcollc2, iblo2, iscdi2, nnz, iexi
    real(kind=8) :: rtbloc
    integer(kind=8), pointer :: smdi(:) => null()
    integer(kind=8), pointer :: smde(:) => null()
    integer(kind=8), pointer :: scde(:) => null()
    integer(kind=8), pointer :: scib(:) => null()
    integer(kind=8), pointer :: scdi(:) => null()
    integer(kind=8), pointer :: lc2m(:) => null()
    integer(kind=8), pointer :: m2lc(:) => null()
    real(kind=8), pointer :: digs(:) => null()
    real(kind=8), pointer :: dig2(:) => null()
    complex(kind=8), pointer :: digsc(:) => null()
    complex(kind=8), pointer :: digc2(:) => null()
    integer(kind=4), pointer :: copier(:) => null()
    integer(kind=8), pointer :: nbterm_lc(:, :) => null()
    integer(kind=8), pointer :: decal(:, :) => null()
    character(len=24), pointer :: refa(:) => null()
    logical :: sym
!     ------------------------------------------------------------------

    call jemarq()
    matas = mataz
    base = basz
    call dismoi('MPI_COMPLET', matas, 'MATR_ASSE', repk=kmpic)
    if (kmpic .ne. 'OUI') then
        call utmess('F', 'CALCULEL6_54')
    end if
    if (base .eq. ' ') call jelira(matas//'.VALM', 'CLAS', cval=base)

!   -- On detruit .UALF s'il existe deja :
!   --------------------------------------
    call jedetr(matas//'.UALF')

    call jeveuo(matas//'.REFA', 'L', vk24=refa)
    nu = refa(2) (1:14)
    stomor = nu//'.SMOS'
    stolci = nu//'.SLCS'

!   -- Si le stockage stolci n'est pas encore cree, on le fait :
!   -------------------------------------------------------------
    call jeexin(stolci//'.SCDE', iret)
    if (iret .eq. 0) then
        call jelira(stomor//'.SMDI', 'CLAS', cval=basto)
        rtbloc = jevtbl('TAILLE_BLOC')
        call smosli(stomor, stolci, basto, rtbloc)
    end if

    call jeveuo(stomor//'.SMDE', 'L', vi=smde)
    call jeveuo(stomor//'.SMDI', 'L', vi=smdi)
    call jeveuo(stomor//'.SMHC', 'L', jsmhc)
    nnz = smde(2)

    call jeveuo(stolci//'.SCDE', 'L', vi=scde)
    call jeveuo(stolci//'.SCDI', 'L', vi=scdi)
    call jeveuo(stolci//'.SCHC', 'L', jschc)
    call jeveuo(stolci//'.SCIB', 'L', vi=scib)
    neq = scde(1)
    itbloc = scde(2)
    nbloc = scde(3)

    call jeveuo(stolci//'.M2LC', 'L', vi=m2lc)
    call jeveuo(stolci//'.LC2M', 'L', vi=lc2m)

    call jelira(matas//'.VALM', 'NMAXOC', nblocm)
    ASSERT(nblocm .eq. 1 .or. nblocm .eq. 2)
    nbloclc = nblocm*nbloc
    sym = nblocm .eq. 1

!   -- reel ou complexe ?
    call jelira(matas//'.VALM', 'TYPE', cval=tyrc)
    ASSERT(tyrc .eq. 'R' .or. tyrc .eq. 'C')

!   1. Allocation de .UALF :
!   ------------------------
    call jecrec(matas//'.UALF', base//' V '//tyrc, 'NU', 'DISPERSE', 'CONSTANT', &
                nbloclc)
    call jeecra(matas//'.UALF', 'LONMAX', itbloc)
    do ibloc = 1, nbloclc
        call jecroc(jexnum(matas//'.UALF', ibloc))
    end do

!   2. Remplissage de .UALF :
!   -------------------------

!   -----------------------------------------------------------------------------------------
!   La permutation par .M2LC transfere des termes entre les deux moities de la matrice (I/S).
!   Il faut optimiser les jeveuo/jelibe sur les blocs de .UALF
!   -----------------------------------------------------------------------------------------

!   2.1 : On compte les termes qui vont s'assembler dans les differents blocs de .UALF
!         => l'objet nbterm_lc(kblocm,kbloclc)=n1 ;
!            n1 : nombre de termes de .VALM(kblocm) a recopier dans .UALF(kbloclc)
!   -----------------------------------------------------------------------------------
    allocate (nbterm_lc(nblocm, nbloclc))
    nbterm_lc(:, :) = 0
    do kblocm = 1, nblocm
        do jcollc = 1, neq
            ibloc = scib(jcollc)+nbloc*(kblocm-1)

            jcolm = lc2m(jcollc)
            ismdi = smdi(jcolm)
            if (jcolm .gt. 1) then
                ismdim1 = smdi(jcolm-1)
            else
                ismdim1 = 0
            end if
            nbterm = ismdi-ismdim1

            do kterm = 1, nbterm
                iligm = zi4(jsmhc-1+ismdim1+kterm)
                iliglc = m2lc(iligm)
!               -- si le terme reste du meme cote de la diagonale :
                if (iliglc .le. jcollc) then
                    nbterm_lc(kblocm, ibloc) = nbterm_lc(kblocm, ibloc)+1
                else
                    iliglc2 = jcollc
                    jcollc2 = iliglc
                    if (sym) then
                        iblo2 = scib(jcollc2)
                    else
                        if (kblocm .eq. 1) then
                            iblo2 = scib(jcollc2)+nbloc
                        else if (kblocm .eq. 2) then
                            iblo2 = scib(jcollc2)
                        else
                            ASSERT(.false.)
                        end if
                    end if
                    nbterm_lc(kblocm, iblo2) = nbterm_lc(kblocm, iblo2)+1
                end if
            end do
        end do
    end do

!   2.2  On note ou doivent etre recopies les termes de .VALM :
!   -----------------------------------------------------------

!   -- calcul de decal(kblocm,kbloclc) :
    allocate (decal(nblocm, nbloclc))
    decal(:, :) = 0
    ico = 0
    do kblocm = 1, nblocm
        do ibloc = 1, nbloclc
            decal(kblocm, ibloc) = ico
            ico = ico+nbterm_lc(kblocm, ibloc)
        end do
    end do

!   -- calcul du vecteur copier :
    AS_ALLOCATE(vi4=copier, size=2*nblocm*nnz)
    nbterm_lc(:, :) = 0
    do kblocm = 1, nblocm
        do jcollc = 1, neq
            iscdi = scdi(jcollc)

            jcolm = lc2m(jcollc)
            ismdi = smdi(jcolm)
            if (jcolm .gt. 1) then
                ismdim1 = smdi(jcolm-1)
            else
                ismdim1 = 0
            end if
            nbterm = ismdi-ismdim1

            do kterm = 1, nbterm
                iligm = zi4(jsmhc-1+ismdim1+kterm)
                iliglc = m2lc(iligm)
!               -- si le terme reste du meme cote de la diagonale :
                if (iliglc .le. jcollc) then
                    ibloc = scib(jcollc)+nbloc*(kblocm-1)
                    copier(2*decal(kblocm, ibloc)+2*nbterm_lc(kblocm, ibloc)+1) = ismdim1+kterm
                   copier(2*decal(kblocm, ibloc)+2*nbterm_lc(kblocm, ibloc)+2) = iscdi+iliglc-jcollc
                    nbterm_lc(kblocm, ibloc) = nbterm_lc(kblocm, ibloc)+1
                else
                    iliglc2 = jcollc
                    jcollc2 = iliglc
                    iscdi2 = scdi(jcollc2)
                    if (sym) then
                        iblo2 = scib(jcollc2)
                    else
                        if (kblocm .eq. 1) then
                            iblo2 = scib(jcollc2)+nbloc
                        else if (kblocm .eq. 2) then
                            iblo2 = scib(jcollc2)
                        else
                            ASSERT(.false.)
                        end if
                    end if
                    copier(2*decal(kblocm, iblo2)+2*nbterm_lc(kblocm, iblo2)+1) = ismdim1+kterm
                    copier(2*decal(kblocm, iblo2)+2*nbterm_lc(kblocm, iblo2)+2) = &
                        iscdi2+iliglc2-jcollc2
                    nbterm_lc(kblocm, iblo2) = nbterm_lc(kblocm, iblo2)+1
                end if
            end do
        end do
    end do

!   2.3 On recopie .VALM dans .UALF :
!   ----------------------------------
    do kblocm = 1, nblocm
        call jeveuo(jexnum(matas//'.VALM', kblocm), 'L', jvalm)
        do ibloc = 1, nbloclc
            call jeveuo(jexnum(matas//'.UALF', ibloc), 'E', jualf)
            n1 = nbterm_lc(kblocm, ibloc)
            do k1 = 1, n1
                ii1 = copier(2*decal(kblocm, ibloc)+2*(k1-1)+1)
                ii2 = copier(2*decal(kblocm, ibloc)+2*(k1-1)+2)
                if (tyrc .eq. 'R') then
                    zr(jualf-1+ii2) = zr(jvalm-1+ii1)
                else
                    zc(jualf-1+ii2) = zc(jvalm-1+ii1)
                end if
            end do
            call jelibe(jexnum(matas//'.UALF', ibloc))
        end do
        call jelibe(jexnum(matas//'.VALM', kblocm))
    end do

!   STOP
!   ASSERT(.false.)
    AS_DEALLOCATE(vi4=copier)
    deallocate (nbterm_lc)
    deallocate (decal)

!   3. L'objet .DIGS(1:neq) a deja ete rempli avec la numerotation de .VALM
!   Il faut le mettre dans la numerotation de .UALF :
!   -----------------------------------------------------------------
    call jeexin(matas//'.DIGS', iexi)
    if (iexi .gt. 0) then
        if (tyrc .eq. 'R') then
            call jeveuo(matas//'.DIGS', 'L', vr=digs)
            ASSERT(size(digs) .eq. 2*neq)
            AS_ALLOCATE(vr=dig2, size=neq)
            dig2(1:neq) = digs(1:neq)
            do k1 = 1, neq
                digs(m2lc(k1)) = dig2(k1)
            end do
            AS_DEALLOCATE(vr=dig2)
        else
            call jeveuo(matas//'.DIGS', 'L', vc=digsc)
            ASSERT(size(digsc) .eq. 2*neq)
            AS_ALLOCATE(vc=digc2, size=neq)
            digc2(1:neq) = digsc(1:neq)
            do k1 = 1, neq
                digsc(m2lc(k1)) = digc2(k1)
            end do
            AS_DEALLOCATE(vc=digc2)
        end if
    end if

    call jedema()
end subroutine
