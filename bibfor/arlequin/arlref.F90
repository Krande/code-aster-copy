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
subroutine arlref(elrefe, fami, nomte, ndim, nno, &
                  nnos, npg, jpoids, jcoopg, jvf, &
                  jdfde, jdfd2, jgano)
!
    use calcul_module, only: ca_iactif_, ca_jnolfp_, ca_jpnlfp_, ca_nblfpg_, ca_nbsav_
!
    implicit none
!
#include "MeshTypes_type.h"
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/assert.h"
#include "asterfort/elraca.h"
#include "asterfort/elref1.h"
#include "asterfort/indk32.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"

    character(len=*), intent(in), optional :: elrefe
    character(len=*), intent(in)    :: fami
    character(len=16), intent(in)   :: nomte
    integer(kind=8), intent(out), optional  :: ndim
    integer(kind=8), intent(out), optional  :: nno
    integer(kind=8), intent(out), optional  :: nnos
    integer(kind=8), intent(out), optional  :: npg
    integer(kind=8), intent(out), optional  :: jpoids
    integer(kind=8), intent(out), optional  :: jcoopg
    integer(kind=8), intent(out), optional  :: jvf
    integer(kind=8), intent(out), optional  :: jdfde
    integer(kind=8), intent(out), optional  :: jdfd2
    integer(kind=8), intent(out), optional  :: jgano
! but: recuperer des informations sur l'element de reference :
!      - dimension de l'espace, nombre de noeuds, de points de Gauss, ...
!      - poids des points de gauss  : jpoids
!      - coordonnees des points de gauss  : jcoopg
!      - valeurs des fonctions de forme : jvf
!      - valeurs des derivees 1eres des fonctions de forme : jdfde
!      - valeurs des derivees 2emes des fonctions de forme : jdfd2
!      - matrice de passage gauss -> noeuds : jgano
! ----------------------------------------------------------------------
!   in   elrefe : nom de l'elrefe (k8) (par defaut l'elrefe principal).
!        fami   : nom de la famille de points de gauss :
!                 'RIGI','MASS',...
!   out  ndim   : dimension de l'espace (=nb coordonnees)
!        nno    : nombre de noeuds du type_maille
!        nnos   : nombre de noeuds sommets du type_maille
!        npg    : nombre de points de gauss
!        jpoids : adresse dans zr du tableau poids(ipg)
!        jcoopg : adresse dans zr du tableau coopg(idim,ipg)
!        jvf    : adresse dans zr du tableau ff(ino,ipg)
!        jdfde  : adresse dans zr du tableau dff(idim,ino,ipg)
!        jdfd2  : adresse dans zr du tableau dff2(idim,jdim,ino,ipg)
!        jgano  : adresse dans zr de la matrice de passage
!                 gauss -> noeuds (dim= 2+nno*npg)
!                 remarque importante : les 2 1ers termes sont les
!                             dimensions de la matrice: nno et npg

!   -------------------------------------------------------------------

    character(len=8) :: elrf, famil, fapg(MT_NBFAMX)
    character(len=16) :: nofgpg
    character(len=32) :: noflpg
    integer(kind=8) :: nbfpg, nbpg(MT_NBFAMX), jvr, decal, ifam, lonfam
    integer(kind=8) :: nufpg, nufgpg, nuflpg, jdfd2l, jganol
    integer(kind=8) :: ndiml, nnosl, nnol, npgl, jpoidl, jcoopl, jvfl, jdfdel

!   -- pour faire des "save" et gagner du temps CPU :
    integer(kind=8) :: maxsav
    parameter(maxsav=5)
    integer(kind=8) :: addsav(5, 10), k1, k2, nusav
    character(len=32) :: nomsav(maxsav)
    save nomsav, addsav

! DEB ------------------------------------------------------------------

!   -- pour etre sur que elrefe est appele "sous" te0000
    ASSERT(ca_iactif_ .eq. 3)

    if (.not. present(elrefe)) then
        call elref1(elrf)
        ASSERT(elrf .ne. 'XXXXXXXX')
    else
        elrf = elrefe
    end if

    famil = fami
    noflpg = nomte//elrf//famil

!   -- pour gagner du temps, on regarde si la famille a ete sauvee:
!   ---------------------------------------------------------------
    nusav = indk32(nomsav, noflpg, 1, ca_nbsav_)
    if (nusav .gt. 0) then
        ndiml = addsav(nusav, 1)
        nnol = addsav(nusav, 2)
        nnosl = addsav(nusav, 3)
        npgl = addsav(nusav, 4)
        jpoidl = addsav(nusav, 5)
        jcoopl = addsav(nusav, 6)
        jvfl = addsav(nusav, 7)
        jdfdel = addsav(nusav, 8)
        jdfd2l = addsav(nusav, 9)
        jganol = addsav(nusav, 10)
        goto 40
    end if

!   -- calcul de nufpg :
!   --------------------
    nuflpg = indk32(zk32(ca_jpnlfp_), noflpg, 1, ca_nblfpg_)
    if (nuflpg .eq. 0) then
        call utmess('F', 'CALCUL_43', sk=noflpg)
    end if
    nufgpg = zi(ca_jnolfp_-1+nuflpg)
    if (nufgpg .eq. 0) then
        call utmess('F', 'CALCUL_7', sk=noflpg)
    end if
    call jenuno(jexnum('&CATA.TM.NOFPG', nufgpg), nofgpg)
    ASSERT(nofgpg(1:8) .eq. elrf)

! - Get list of integration schemes of geometric support
    call elraca(elrf, &
                nbfpg_=nbfpg, fapg_=fapg, nbpg_=nbpg, &
                ndim_=ndiml, nno_=nnol, nnos_=nnosl)

! - Get index for integration scheme
    nufpg = indik8(fapg, nofgpg(9:16), 1, nbfpg)
    ASSERT(nufpg .gt. 0)

    call jeveuo('&INEL.'//elrf//'.ELRA_R', 'L', jvr)

    decal = 0
    do ifam = 1, nufpg-1
        npgl = nbpg(ifam)
        lonfam = npgl
        lonfam = lonfam+npgl*ndiml
        lonfam = lonfam+npgl*nnol
        lonfam = lonfam+npgl*nnol*ndiml
        lonfam = lonfam+npgl*nnol*ndiml*ndiml
        lonfam = lonfam+2+npgl*nnol
        decal = decal+lonfam
    end do

    npgl = nbpg(nufpg)

    jpoidl = jvr+decal
    jcoopl = jpoidl+npgl
    jvfl = jcoopl+npgl*ndiml
    jdfdel = jvfl+npgl*nnol
    jdfd2l = jdfdel+npgl*nnol*ndiml
    jganol = jdfd2l+npgl*nnol*ndiml*ndiml

!   -- on sauvegarde les valeurs calculees :
!   ----------------------------------------
!   -- on decale tout le monde vers le bas:
    ca_nbsav_ = min(ca_nbsav_+1, maxsav)
    do k1 = ca_nbsav_-1, 1, -1
        nomsav(k1+1) = nomsav(k1)
        do k2 = 1, 10
            addsav(k1+1, k2) = addsav(k1, k2)
        end do
    end do

!   -- on recopie les nouvelles valeurs en position 1 :
    nomsav(1) = noflpg
    addsav(1, 1) = ndiml
    addsav(1, 2) = nnol
    addsav(1, 3) = nnosl
    addsav(1, 4) = npgl
    addsav(1, 5) = jpoidl
    addsav(1, 6) = jcoopl
    addsav(1, 7) = jvfl
    addsav(1, 8) = jdfdel
    addsav(1, 9) = jdfd2l
    addsav(1, 10) = jganol

40  continue

    if (present(ndim)) ndim = ndiml
    if (present(nnos)) nnos = nnosl
    if (present(nno)) nno = nnol
    if (present(npg)) npg = npgl
    if (present(jpoids)) jpoids = jpoidl
    if (present(jcoopg)) jcoopg = jcoopl
    if (present(jvf)) jvf = jvfl
    if (present(jdfde)) jdfde = jdfdel
    if (present(jdfd2)) jdfd2 = jdfd2l
    if (present(jgano)) jgano = jganol

end subroutine
