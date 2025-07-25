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

subroutine mesomm(champ, long, vi, vr, vc, &
                  nbma, linuma)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/assert.h"
#include "asterfort/celver.h"
#include "asterfort/digdel.h"
#include "asterfort/dismoi.h"
#include "asterfort/jaexin.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelibe.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/nbelem.h"
#include "asterfort/nbgrel.h"
#include "asterfort/scalai.h"
#include "asterfort/utmess.h"
!
    character(len=*) :: champ
    integer(kind=8), intent(in) :: long
    integer(kind=8), intent(in), optional :: nbma
    integer(kind=8), intent(in), optional :: linuma(*)
    integer(kind=8), intent(out), optional :: vi(*)
    real(kind=8), intent(out), optional :: vr(*)
    complex(kind=8), intent(out), optional :: vc(*)
!     but :  faire la "somme" d'un cham_elem (ou d'un resuelem)
!                  ou d'une partie d'un cham_elem
!            (notion d'integrale du champ sur le modele)
!            la seule contrainte est que tous les type_element du ligrel
!            connaissent la grandeur avec la meme longueur cumulee :
!
!            l'exemple suivant sera traite par la routine, alors que
!            sa signification est probablement douteuse ...
!              tri3 :  e 2 iden 3 ....
!              seg2 :  e 3 iden 2 ....
!              poi1 :  e 6 iden 1 ....
!
! in  : champ  :  nom du champ a sommer
! in  : long   :  longueur des vecteurs vi vr ou vc
!
! in  : nbma   :  longueur de la liste linuma
! in  : linuma :  liste des numeros des mailles
! -- Si nbma et linuma ne sont pas fournis, on calcule sur TOUS les elements.
!
! out : vi     :  vecteur contenant la "somme" du champ si la grandeur
!                 est entiere.
! out : vr     :  vecteur contenant la "somme" du champ si la grandeur
!                 est reelle.
! out : vc     :  vecteur contenant la "somme" du champ si la grandeur
!                 est complexe.
!
! ----------------------------------------------------------------------
!
    integer(kind=8) :: longt, ncmpel, mode, j, igd
    real(kind=8) :: rzero
    character(len=4) :: typch, kmpic
    character(len=8) :: scal
    character(len=19) :: champ2, ligrel
    aster_logical :: first
    integer(kind=8) :: i, iavale, ibid, icoef, idecgr, iel, ier1, ier2
    integer(kind=8) :: im, inum, jceld, jligr, k, nbgr, nel, numel1, iexi, nbmail
    character(len=24), pointer :: celk(:) => null()
!
    call jemarq()
!
    champ2 = champ
    rzero = 0.0d0
    if (present(nbma)) then
        nbmail = nbma
    else
        nbmail = 0
    end if
!
!
!     1- ON CALCULE : TYPCH,LIGREL,IGD ET SCAL :
!     -----------------------------------------
!
    call jeexin(champ2//'.CELD', ier1)
    call jeexin(champ2//'.RESL', ier2)
    if (ier1+ier2 .eq. 0) then
        call utmess('F', 'CALCULEL3_73', sk=champ2)
    end if
!
!
    if (ier1 .gt. 0) then
        typch = 'CHML'
!       -- ON VERIFIE QUE LE CHAM_ELEM N'EST PAS TROP DYNAMIQUE :
        call celver(champ2, 'NBVARI_CST', 'STOP', ibid)
        call celver(champ2, 'NBSPT_1', 'STOP', ibid)
        call jeveuo(champ2//'.CELD', 'L', jceld)
    else
        typch = 'RESL'
        call jeveuo(champ2//'.DESC', 'L', jceld)
    end if
!
    call jeveuo(champ2//'.CELK', 'L', vk24=celk)
    ligrel = celk(1) (1:19)
!
    igd = zi(jceld-1+1)
    scal = scalai(igd)
!
    if (scal(1:1) .eq. 'R') then
        ASSERT(present(vr))
    else if (scal(1:1) .eq. 'C') then
        ASSERT(present(vc))
    else if (scal(1:1) .eq. 'I') then
        ASSERT(present(vi))
    else
        ASSERT(.false.)
    end if
!
!
!     2- ON VERIFIE LES LONGUEURS:
!     ----------------------------
    first = .true.
    nbgr = nbgrel(ligrel)
    do j = 1, nbgr
        mode = zi(jceld-1+zi(jceld-1+4+j)+2)
        if (mode .eq. 0) goto 10
        ncmpel = digdel(mode)
        icoef = max(1, zi(jceld-1+4))
        ncmpel = ncmpel*icoef
        if (first) then
            longt = ncmpel
        else
            if (longt .ne. ncmpel) then
                call utmess('F', 'CALCULEL3_54')
            end if
        end if
        first = .false.
10      continue
    end do
!
!     -- ON MET A ZERO LE VECTEUR "VSCAL":
!     ------------------------------------
    ASSERT(longt .le. long)
!
    do i = 1, long
        if (scal(1:1) .eq. 'I') then
            vi(i) = 0
        else if (scal(1:1) .eq. 'R') then
            vr(i) = rzero
        else if (scal(1:1) .eq. 'C') then
            vc(i) = dcmplx(rzero, rzero)
        else
            call utmess('F', 'CALCULEL3_74', sk=scal)
        end if
    end do
!
!        -- ON CUMULE :
!        --------------
    if (typch .eq. 'CHML') then
!        -- (CAS DES CHAM_ELEM):
        call jeveuo(champ2//'.CELV', 'L', iavale)
        if (nbmail .eq. 0) then
            do j = 1, nbgr
                mode = zi(jceld-1+zi(jceld-1+4+j)+2)
                if (mode .eq. 0) goto 50
                nel = nbelem(ligrel, j)
                idecgr = zi(jceld-1+zi(jceld-1+4+j)+8)
                do k = 1, nel
                    do i = 1, longt
                        if (scal(1:1) .eq. 'I') then
                            vi(i) = vi(i)+zi(iavale-1+idecgr+(k-1)*longt+i-1)
                        else if (scal(1:1) .eq. 'R') then
                            vr(i) = vr(i)+zr(iavale-1+idecgr+(k-1)*longt+i-1)
                        else if (scal(1:1) .eq. 'C') then
                            vc(i) = vc(i)+zc(iavale-1+idecgr+(k-1)*longt+i-1)
                        end if
                    end do
                end do
50              continue
            end do
        else
            call jeveuo(ligrel//'.LIEL', 'L', jligr)
            do im = 1, nbmail
                inum = 0
                do j = 1, nbgr
                    mode = zi(jceld-1+zi(jceld-1+4+j)+2)
                    if (mode .eq. 0) goto 79
                    nel = nbelem(ligrel, j)
                    idecgr = zi(jceld-1+zi(jceld-1+4+j)+8)
                    do k = 1, nel
                        iel = zi(jligr+inum+k-1)
                        if (iel .ne. linuma(im)) goto 70
                        do i = 1, longt
                            if (scal(1:1) .eq. 'I') then
                                vi(i) = vi(i)+zi(iavale-1+idecgr+(k-1)*longt+i-1)
                            else if (scal(1:1) .eq. 'R') then
                                vr(i) = vr(i)+zr(iavale-1+idecgr+(k-1)*longt+i-1)
                            else if (scal(1:1) .eq. 'C') then
                                vc(i) = vc(i)+zc(iavale-1+idecgr+(k-1)*longt+i-1)
                            end if
                        end do
                        goto 90
70                      continue
                    end do
79                  continue
                    inum = inum+nel+1
                end do
90              continue
            end do
        end if
!
    else if (typch .eq. 'RESL') then
!        -- (CAS DES RESUELEM):
        if (nbmail .le. 0) then
            numel1 = 0
            do j = 1, nbgr
                mode = zi(jceld-1+zi(jceld-1+4+j)+2)
                if (mode .eq. 0) goto 120
                call jaexin(jexnum(champ2//'.RESL', j), iexi)
                if (iexi .eq. 0) goto 120
                call jeveuo(jexnum(champ2//'.RESL', j), 'L', iavale)
                ncmpel = digdel(mode)
                nel = nbelem(ligrel, j)
                numel1 = numel1+nel
                do k = 1, nel
                    do i = 1, longt
                        if (scal(1:1) .eq. 'I') then
                            vi(i) = vi(i)+zi(iavale+(k-1)*ncmpel-1+i)
                        else if (scal(1:1) .eq. 'R') then
                            vr(i) = vr(i)+zr(iavale+(k-1)*ncmpel-1+i)
                        else if (scal(1:1) .eq. 'C') then
                            vc(i) = vc(i)+zc(iavale+(k-1)*ncmpel-1+i)
                        end if
                    end do
                end do
                call jelibe(jexnum(champ2//'.RESL', j))
120             continue
            end do
        else
            call jeveuo(ligrel//'.LIEL', 'L', jligr)
            do im = 1, nbmail
                inum = 0
                numel1 = 0
                do j = 1, nbgr
                    mode = zi(jceld-1+zi(jceld-1+4+j)+2)
                    if (mode .eq. 0) goto 149
                    call jaexin(jexnum(champ2//'.RESL', j), iexi)
                    if (iexi .eq. 0) goto 150
                    call jeveuo(jexnum(champ2//'.RESL', j), 'L', iavale)
                    ncmpel = digdel(mode)
                    nel = nbelem(ligrel, j)
                    numel1 = numel1+nel
                    do k = 1, nel
                        iel = zi(jligr+inum+k-1)
                        if (iel .ne. linuma(im)) goto 140
                        do i = 1, longt
                            if (scal(1:1) .eq. 'I') then
                                vi(i) = vi(i)+zi(iavale+(k-1)*ncmpel-1+i)
                            else if (scal(1:1) .eq. 'R') then
                                vr(i) = vr(i)+zr(iavale+(k-1)*ncmpel-1+i)
                            else if (scal(1:1) .eq. 'C') then
                                vc(i) = vc(i)+zc(iavale+(k-1)*ncmpel-1+i)
                            end if
                        end do
                        call jelibe(jexnum(champ2//'.RESL', j))
                        goto 160
140                     continue
                    end do
149                 continue
                    inum = inum+nel+1
150                 continue
                end do
160             continue
            end do
        end if
!
    end if
!
!
!     -- IL FAUT COMMUNIQUER LE RESULTAT ENTRE LES PROCS :
    call dismoi('MPI_COMPLET', champ, 'CHAMP', repk=kmpic)
    if (kmpic .eq. 'NON') then
        if (scal(1:1) .eq. 'I') then
            call asmpi_comm_vect('MPI_SUM', 'I', nbval=longt, vi=vi)
        else if (scal(1:1) .eq. 'R') then
            call asmpi_comm_vect('MPI_SUM', 'R', nbval=longt, vr=vr)
        else if (scal(1:1) .eq. 'C') then
            ASSERT(.false.)
        end if
    end if
!
    call jedema()
end subroutine
