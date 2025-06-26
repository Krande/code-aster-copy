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
subroutine ldlt_renum(numeDof1Z, numeDof2Z_, perm_, permJvBase_)
!
    implicit none
!
#include "asterc/getres.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/infbav.h"
#include "asterfort/infmue.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/nueffe.h"
#include "asterfort/promor.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
!
    character(len=*), intent(in) :: numeDof1Z
    character(len=*), intent(in), optional :: numeDof2Z_
    character(len=24), intent(in), optional :: perm_
    character(len=1), intent(in), optional :: permJvBase_

! --------------------------------------------------------------------------------------------------
! Fonctions :
! ------------
!  Cette routine sert a renumeroter une matrice pour que la factorisation LDLT
!  soit plus efficace. La méthode de renumerotation que l'on essaye d'utiliser est 'RCMK'
!
!  1) Si numeDof2Z_ est absent :
!     Ajouter les objets .M2LC et .LC2M dans la sd_nume_ddl numeDof1Z
!
!     Ces 2 objets etablissent la correspondance entre les numerotations du stockage Morse (.VALM)
!     et du stockage Ligne de Ciel (.UALF). Ces deux vecteurs sont reciproques l'un de l'autre :
!       * ieqlc=.M2LC(ieqm)
!       * ieqm =.LC2M(ieqlc)
!     Ces 2 objets sont crees sur la meme base (G/V) que la base de numeDof1

!  2) Si numeDof2Z_ est present :
!     * Calculer un nouveau nume_ddl numeDof2Z_ correspondant au nume_ddl numeDof1Z mais avec
!       la numerotation 'RCMK'.
!     * Retourner egalement l'objet perm_ qui etablit la permutation entre numeDof1Z et numeDof2Z_.
!     numeDof2Z_ est cree sur la base volatile ('V')
!     perm_ est cree sur la base 'permJvBase_'
!
! --------------------------------------------------------------------------------------------------
!
! (o) In/Jxvar  numeDof1Z  : sd_nume_ddl
! (f) In/Jxout  numeDof2Z_  : sd_nume_ddl
! (f) In/Jxout  perm_  : K24 : OJB V I
! (f) In        permJvBase_  : K1 : base pour la creation de perm_
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24), parameter :: renumRCMK = "RCMK"
    character(len=14) :: numeDof1, numeDof2, numeDofTrav
    character(len=8) :: nogd
    character(len=1) :: jvBase
  integer(kind=8) :: nbLigr, k, neq2, ncmp1, ncmp2, nbno1, nbno2, iad1, iad2, icmp, ieq1, ieq2, nbec
    integer(kind=8) :: nec1, nec2, n1, n2, ino, neq, iexi, ifm, niv
    character(len=24), pointer :: listLigr(:) => null()
    character(len=24), pointer :: refn(:) => null()
    character(len=19) :: ligrName
    character(len=8) :: modeLoc, result, model
    character(len=16) :: command, resultType
    logical :: non_renum, limpr
    integer(kind=8), pointer :: prno1(:) => null()
    integer(kind=8), pointer :: prno2(:) => null()
    integer(kind=8), pointer :: nueq1(:) => null()
    integer(kind=8), pointer :: nueq2(:) => null()
    integer(kind=8), pointer :: m2lc(:) => null()
    integer(kind=8), pointer :: lc2m(:) => null()
    character(len=16), save :: commandSave = ' ', resultSave = ' '
    integer(kind=8), save :: neqSav = 0
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

    numeDof1 = numeDof1Z
    call dismoi('NB_EQUA', numeDof1, 'NUME_DDL', repi=neq)

! - Manage output in mess file
    call getres(result, resultType, command)
    limpr = ASTER_TRUE
    call infniv(ifm, niv)
    if (command .eq. commandSave .and. result .eq. resultSave) then
        if (neq .eq. neqSav .and. niv .eq. 1) then
            limpr = .false.
        end if
    end if
    neqSav = neq
    resultSave = result
    commandSave = command
    if (.not. limpr) then
        call infmue()
    end if

    if (.not. present(numeDof2Z_)) then
        call jelira(numeDof1//'.SMOS.SMDE', 'CLAS', cval=jvBase)
        call jedetr(numeDof1//'.SLCS.M2LC')
        call jedetr(numeDof1//'.SLCS.LC2M')
        call wkvect(numeDof1//'.SLCS.M2LC', jvBase//' V I', neq, vi=m2lc)
        call wkvect(numeDof1//'.SLCS.LC2M', jvBase//' V I', neq, vi=lc2m)
    else
        numeDof2 = numeDof2Z_
        ASSERT(present(perm_))
        ASSERT(present(permJvBase_))
        call jedetr(perm_)
        call wkvect(perm_, permJvBase_//' V I', neq, vi=m2lc)
    end if

!   -- On ne renumerote pas (avec RCMK) si :
!       * C'est un NUME_DDL_GENE (pas de noeuds)
!       * On est dans la commande MACR_ELEM_STAT (RCMK a deja ete fait)
!       * C'est un le NUME_DDL d'une matrice "reduite" avec ELIM_LAGR,
!         (car on a supprime certains ddls)
!   -------------------------------------------------------------------
    non_renum = (command .eq. 'MACR_ELEM_STAT')
    call jenuno(jexnum(numeDof1//'.NUME.LILI', 1), ligrName)
    non_renum = non_renum .or. (ligrName .eq. '&SOUSSTR')
    call jeveuo(numeDof1//'.NUME.REFN', 'L', vk24=refn)
    non_renum = non_renum .or. (refn(4) .eq. 'ELIM_LAGR')
    if (non_renum) then
        if (present(numeDof2Z_)) then
            call copisd('NUME_DDL', 'V', numeDof1, numeDof2)
            do k = 1, neq
                m2lc(k) = k
            end do
        else
            do k = 1, neq
                m2lc(k) = k
                lc2m(k) = k
            end do
        end if
        goto 999
    end if

!   -- Si on est dans la commande CALC_ERREUR, il faut utiliser l'argument modelocz de nueffe :
!   -------------------------------------------------------------------------------------------
    if (command .eq. 'CALC_ERREUR') then
        modeLoc = 'DDL_NOZ1'
    else
        modeLoc = ' '
    end if

!   -- Renumerotation RCMK :
!   ---------------------------------------------------------------
    if (present(numeDof2Z_)) then
        numeDofTrav = numeDof2
    else
        numeDofTrav = '&&ldlt_RENUM.N'
    end if
    call copisd('NUME_EQUA', 'V', numeDof1//'.NUME', numeDofTrav//'.NUME')

!   -- il faut ignorer .LILI(1) -> '&MAILLA' :
    call jelira(numeDof1//'.NUME.LILI', 'NOMUTI', nbLigr)
    ASSERT(nbLigr .ge. 2)
    nbLigr = nbLigr-1
    AS_ALLOCATE(vk24=listLigr, size=nbLigr)
    do k = 1, nbLigr
        call jenuno(jexnum(numeDof1//'.NUME.LILI', 1+k), ligrName)
        listLigr(k) = ligrName
    end do
    call jedetr(numeDofTrav//'     .ADNE')
    call jedetr(numeDofTrav//'     .ADLI')

! - Affreuse glute (tant qu'on stockera le modèle dans NUME_EQUA/REFN)
    model = refn(4) (1:8)
    call nueffe(nbLigr, listLigr, 'VV', numeDofTrav, renumRCMK, model, &
                modeLocZ_=modeLoc)
    AS_DEALLOCATE(vk24=listLigr)

! -- Pour etablir la correspondance entre les deux numerotations,
!      il faut comparer les objets .PRNO et .NUEQ
    call dismoi('NB_EQUA', numeDofTrav, 'NUME_DDL', repi=neq2)
    ASSERT(neq .eq. neq2)

    call dismoi('NOM_GD', numeDof1, 'NUME_DDL', repk=nogd)
    call dismoi('NB_EC', nogd, 'GRANDEUR', repi=nec1)
    call dismoi('NOM_GD', numeDofTrav, 'NUME_DDL', repk=nogd)
    call dismoi('NB_EC', nogd, 'GRANDEUR', repi=nec2)
    ASSERT(nec1 .eq. nec2)
    nbec = nec1

    call jelira(numeDof1//'.NUME.PRNO', 'NMAXOC', n1)
    call jelira(numeDofTrav//'.NUME.PRNO', 'NMAXOC', n2)
    ASSERT(n1 .eq. n2)

    call jeveuo(numeDof1//'.NUME.NUEQ', 'L', vi=nueq1)
    call jeveuo(numeDofTrav//'.NUME.NUEQ', 'L', vi=nueq2)

    do k = 1, n1
        call jelira(jexnum(numeDof1//'.NUME.PRNO', k), 'LONMAX', iexi)
        if (iexi .eq. 0) cycle
        call jeveuo(jexnum(numeDof1//'.NUME.PRNO', k), 'L', vi=prno1)
        call jeveuo(jexnum(numeDofTrav//'.NUME.PRNO', k), 'L', vi=prno2)
        nbno1 = size(prno1)/(nbec+2)
        nbno2 = size(prno2)/(nbec+2)
        ASSERT(nbno1 .eq. nbno2)
        do ino = 1, nbno1
            ncmp1 = prno1((ino-1)*(2+nbec)+2)
            ncmp2 = prno2((ino-1)*(2+nbec)+2)
            ASSERT(ncmp1 .eq. ncmp2)
            iad1 = prno1((ino-1)*(2+nbec)+1)
            iad2 = prno2((ino-1)*(2+nbec)+1)
            do icmp = 1, ncmp1
                ieq1 = nueq1(iad1-1+icmp)
                ieq2 = nueq2(iad2-1+icmp)
                ASSERT(ieq1 .ge. 1 .and. ieq1 .le. neq)
                ASSERT(ieq2 .ge. 1 .and. ieq2 .le. neq)
                m2lc(ieq1) = ieq2
            end do
        end do
    end do

    do k = 1, neq
        ASSERT(m2lc(k) .ge. 1 .and. m2lc(k) .le. neq)
    end do

    if (.not. present(numeDof2Z_)) then
        do k = 1, neq
            lc2m(m2lc(k)) = k
        end do
        call detrsd('NUME_EQUA', numeDofTrav//'.NUME')
    else
        call promor(numeDofTrav, 'V')
    end if

999 continue

    if (.not. limpr) then
        call infbav()
    end if
    call jedema()
end subroutine
