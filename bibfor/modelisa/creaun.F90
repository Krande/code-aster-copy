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
subroutine creaun(char, noma, nomo, nzocu, nnocu, &
                  lisnoe, poinoe, nbgdcu, coefcu, compcu, &
                  multcu, penacu)
!
!
    implicit none
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/exiscp.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/int_to_char8.h"
!
    character(len=8) :: char
    character(len=8) :: noma
    character(len=8) :: nomo
    integer(kind=8) :: nzocu, nnocu
    character(len=24) :: lisnoe
    character(len=24) :: poinoe
    character(len=24) :: nbgdcu
    character(len=24) :: coefcu
    character(len=24) :: compcu
    character(len=24) :: multcu
    character(len=24) :: penacu
!
! ----------------------------------------------------------------------
!
! ROUTINE LIAISON_UNILATERALE (CREATION SD)
!
! CONSTRUCTION FINALE DES VECTEURS ON OUBLIE LE CONCEPT DE ZONES
!
! ----------------------------------------------------------------------
!
!
! IN  CHAR   : NOM DU CONCEPT CHARGE
! IN  NOMA   : NOM DU MAILLAGE
! IN  NOMO   : NOM DU MODELE
! IN  NZOCU  : NOMBRE DE ZONES
! IN  NNOCU  : NOMBRE DE NOEUDS
! IN  POINOE : NOM DE L'OBJET CONTENANT LE VECTEUR D'INDIRECTION
!               DES NOEUDS
! IN  LISNOE : NOM DE L'OBJET CONTENANT LES NOEUDS
! IN  NBGDCU : NOM JEVEUX DE LA SD INFOS POINTEURS GRANDEURS DU MEMBRE
!              DE GAUCHE
! IN  COEFCU : NOM JEVEUX DE LA SD CONTENANT LES VALEURS DU MEMBRE
!              DE DROITE
! IN  COMPCU : NOM JEVEUX DE LA SD CONTENANT LES GRANDEURS DU MEMBRE
!              DE GAUCHE
! IN  MULTCU : NOM JEVEUX DE LA SD CONTENANT LES COEFFICIENTS DU MEMBRE
!              DE GAUCHE
! IN  PENACU : NOM JEVEUX DE LA SD CONTENANT LES COEFFICIENTS DE PENALITE
!
!
!
!
    integer(kind=8) :: nbgau
    character(len=24) :: deficu
    integer(kind=8) :: jmult, jnoe, jpoi, jnbgd
    integer(kind=8) :: jcoef, jncmp
    integer(kind=8) :: ino, icmp, izone
    character(len=24) :: noeucu
    character(len=24) :: valk(2)
    integer(kind=8) :: jnoeu
    integer(kind=8) :: numnd, exist(1), nbsup
    integer(kind=8) :: nbno, nbcmp
    integer(kind=8) :: jdebcp, jdebnd
    character(len=8) :: cmp, k8bla, nomno
    integer(kind=8) :: cptd, ncmpg, cptnd
    character(len=24) :: cmpgcu, ndimcu, coegcu, coedcu, poincu
    integer(kind=8) :: jcmpg, jdim, jcoefg, jcoefd, jpoin, jpena
    integer(kind=8) :: ifm, niv
    character(len=8), pointer :: cmpg(:) => null()
    character(len=8), pointer :: coefd(:) => null()
    character(len=8), pointer :: coefg(:) => null()
    integer(kind=8), pointer :: indir(:) => null()
    real(kind=8), pointer :: cpena(:) => null()
!
! ----------------------------------------------------------------------
!
    call jemarq()
    call infniv(ifm, niv)
!
! --- INITIALISATIONS
!
    deficu = char(1:8)//'.UNILATE'
    k8bla = ' '
    call jeveuo(multcu, 'L', jmult)
    call jeveuo(poinoe, 'L', jpoi)
    call jeveuo(lisnoe, 'L', jnoe)
    call jeveuo(nbgdcu, 'L', jnbgd)
    call jeveuo(compcu, 'L', jncmp)
    call jeveuo(coefcu, 'L', jcoef)
    call jeexin(penacu, jpena)
    if (jpena .ne. 0) then
        call jeveuo(penacu, 'L', jpena)
    end if
!
!
! --- CALCUL DU NOMBRE TOTAL DE GRANDEURS A GAUCHE
!
    nbgau = 0
    do izone = 1, nzocu
        nbno = zi(jpoi+izone)-zi(jpoi+izone-1)
        nbcmp = zi(jnbgd+izone)-zi(jnbgd+izone-1)
        nbgau = nbgau+nbno*nbcmp
    end do
!
! --- CREATION DES VECTEURS DEFINITIFS
!
    noeucu = deficu(1:16)//'.LISNOE'
    call wkvect(noeucu, 'G V I', nnocu, jnoeu)
!
! --- CREATION DES VECTEURS TEMPORAIRES
!
    AS_ALLOCATE(vi=indir, size=nnocu+1)
    AS_ALLOCATE(vk8=cmpg, size=nbgau)
    AS_ALLOCATE(vk8=coefg, size=nbgau)
    AS_ALLOCATE(vk8=coefd, size=nnocu)
    if (jpena .ne. 0) then
        AS_ALLOCATE(vr=cpena, size=nnocu)
    end if
    indir(1) = 1
!
! ---
!
    cptnd = 1
    cptd = 1
    ncmpg = 1
!
    do izone = 1, nzocu
!
        nbno = zi(jpoi+izone)-zi(jpoi+izone-1)
        jdebnd = zi(jpoi+izone-1)
        nbcmp = zi(jnbgd+izone)-zi(jnbgd+izone-1)
        jdebcp = zi(jnbgd+izone-1)
!
        do ino = 1, nbno
!
            numnd = zi(jnoe-1+jdebnd+ino-1)
            nbsup = 0
!
            do icmp = 1, nbcmp
!
                cmp = zk8(jncmp-1+jdebcp+icmp-1)
!
                call exiscp(cmp, k8bla, nomo, 1, 'NUM', &
                            k8bla, [numnd], exist)
!
                if (exist(1) .eq. 1) then
                    if (niv .ge. 2) then
                        nomno = int_to_char8(numnd)
                        valk(1) = nomno
                        valk(2) = cmp
                        call utmess('I', 'UNILATER_58', nk=2, valk=valk)
                    end if
                    cmpg(ncmpg) = cmp
                    coefg(ncmpg) = zk8(jmult-1+jdebcp+icmp-1)
                    ncmpg = ncmpg+1
                else
                    nbsup = nbsup+1
                    nomno = int_to_char8(numnd)
                    valk(1) = nomno
                    valk(2) = cmp
                    call utmess('I', 'UNILATER_75', nk=2, valk=valk)
                end if
!
            end do
!
            zi(jnoeu-1+cptnd) = numnd
            coefd(cptd) = zk8(jcoef+izone-1)
            indir(cptnd+1) = indir(cptnd)+nbcmp-nbsup
            if (jpena .ne. 0) then
                cpena(cptd) = zr(jpena+izone-1)
            end if
!
            cptd = cptd+1
            cptnd = cptnd+1
!
        end do
    end do
!
    cptd = cptd-1
    cptnd = cptnd-1
    ncmpg = ncmpg-1
!
    ASSERT(cptd .eq. nnocu)
    ASSERT(cptnd .eq. nnocu)
!
! --- QUELQUES INFOS DIMENSIONS
!
    ndimcu = deficu(1:16)//'.NDIMCU'
    call jeveuo(ndimcu, 'E', jdim)
    zi(jdim) = nnocu
    zi(jdim+1) = ncmpg
!
! --- LISTE DES POINTEURS DES NOEUDS
!
    poincu = deficu(1:16)//'.POINOE'
    call wkvect(poincu, 'G V I', nnocu+1, jpoin)
    do ino = 1, nnocu+1
        zi(jpoin-1+ino) = indir(ino)
    end do
!
! --- LISTE DES NOMS DE COMPOSANTES A GAUCHE
!
    cmpgcu = deficu(1:16)//'.CMPGCU'
    call wkvect(cmpgcu, 'G V K8', ncmpg, jcmpg)
    do icmp = 1, ncmpg
        zk8(jcmpg-1+icmp) = cmpg(icmp)
    end do
!
! --- LISTE DES COEFFICIENTS A DROITE ET A GAUCHE
!
    coegcu = deficu(1:16)//'.COEFG'
    coedcu = deficu(1:16)//'.COEFD'
    call wkvect(coegcu, 'G V K8', ncmpg, jcoefg)
    call wkvect(coedcu, 'G V K8', nnocu, jcoefd)
!
    do icmp = 1, ncmpg
        zk8(jcoefg-1+icmp) = coefg(icmp)
    end do
!
    do icmp = 1, cptd
        zk8(jcoefd-1+icmp) = coefd(icmp)
    end do
!
! --- LISTE DES COEFFICIENTS DE PENALITE
!
    if (jpena .ne. 0) then
        call wkvect(deficu(1:16)//'.COEFPE', 'G V R', nnocu, jpena)
        zr(jpena:jpena+nnocu-1) = cpena(1:nnocu)
    end if
!
! --- NETTOYAGE
!
    AS_DEALLOCATE(vi=indir)
    AS_DEALLOCATE(vk8=cmpg)
    AS_DEALLOCATE(vk8=coefg)
    AS_DEALLOCATE(vk8=coefd)
    AS_DEALLOCATE(vk8=coefd)
    if (jpena .ne. 0) then
        AS_DEALLOCATE(vr=cpena)
    end if
!
! ======================================================================
    call jedema()
!
end subroutine
