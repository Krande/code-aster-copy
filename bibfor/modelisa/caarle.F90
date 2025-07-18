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
subroutine caarle(numdlz, iocc, lisrez, chargz)
!
    implicit none
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/arlcou.h"
#include "asterfort/arllec.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: iocc
    character(len=*) :: numdlz, chargz, lisrez
!
! ----------------------------------------------------------------------
!
!              CREATION D'UN CHARGEMENT DE TYPE ARLEQUIN
!
! ----------------------------------------------------------------------
!
! I/O  CHARGE : SD CHARGE
!
! SD EN ENTREE :
! ==============
!
! .CHME.MODEL.NOMO       : NOM DU MODELE ASSOCIE A LA CHARGE
! .TYPE                  : TYPE DE LA CHARGE
!
! SD EN SORTIE ENRICHIE PAR :
! ===========================
! .CHME.LIGRE     : LIGREL DE CHARGE
! .CHME.CIMPO     : CARTE COEFFICIENTS IMPOSES
! .CHME.CMULT     : CARTE COEFFICIENTS MULTIPLICATEURS
!
!
    character(len=8) :: charge
    character(len=14) :: numddl
    character(len=19) :: lisrel
    character(len=24) :: typmai
    character(len=10) :: noma, nomb, nom1, nom2
    character(len=8) :: nomo, mail, model(3), cine(3)
    character(len=8) :: k8bid
    integer(kind=8) :: dime, nocc, iop
    integer(kind=8) :: nbtyp
    integer(kind=8) :: zocc, ibid, i
    integer(kind=8) :: jtypm, jlgrf
    character(len=16) :: motfac, option
!
    data typmai/'&&CAARLE.NOMTM'/
!
! ----------------------------------------------------------------------
    call jemarq()
!
    numddl = numdlz
    charge = chargz
    lisrel = lisrez
    zocc = iocc
    motfac = 'LIAISON_ELEM'
!
! --- IMPOSE-T-ON UNE CHARGE ARLEQUIN ?
!
    call getfac(motfac, nocc)
    call getvtx(motfac, 'OPTION', iocc=zocc, scal=option, nbret=iop)
    if (option .ne. '3D_POU_ARLEQUIN') then
        call utmess('F', 'MODELISA6_39', sk=option)
    end if
!
    call getfac(motfac, nocc)
    if (nocc .eq. 0) goto 999
!
! --- LECTURE NOMS DU MODELE ET DU MAILLAGE
!
    call getvid(' ', 'MODELE', iocc=0, nbval=1, scal=nomo, &
                nbret=ibid)
    call jeveuo(nomo(1:8)//'.MODELE    .LGRF', 'L', jlgrf)
    mail = zk8(jlgrf)
!
! --- STRUCTURES DE DONNEES
!
    noma = charge(1:8)//'.A'
    nomb = charge(1:8)//'.B'
    nom1 = charge(1:8)//'.1'
    nom2 = charge(1:8)//'.2'
!
! --- CREATION D'UN VECTEUR CONTENANT LE NOM DES TYPES DE MAILLES
!
    call jelira('&CATA.TM.NOMTM', 'NOMMAX', nbtyp, k8bid)
    call wkvect(typmai, 'V V K8', nbtyp, jtypm)
    do i = 1, nbtyp
        call jenuno(jexnum('&CATA.TM.NOMTM', i), zk8(jtypm-1+i))
    end do
!
! --- LECTURE ET VERIFICATION DES MAILLES DES MODELES
!
    call arllec(motfac, zocc, nomo, noma, nomb, &
                model, cine, dime)
!
! --- CALCUL DES EQUATIONS DE COUPLAGE
!
    call arlcou(mail, zocc, nomo, typmai, noma, &
                nomb, cine, dime, lisrel, charge)
!
! --- DESALLOCATION GROUPES 1 ET 2
!
    call jedetr(nom1)
    call jedetr(nom2)
    call jedetr(noma//'.GROUPEMA')
    call jedetr(nomb//'.GROUPEMA')
!
! --- DESALLOCATION AUTRES OBJETS ARLEQUIN
!
    call jedetr(typmai)
!
999 continue
!
    call jedema()
!
end subroutine
