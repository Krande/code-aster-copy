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

subroutine algocp(ds_measure, resoco, numedd, matass)
!
! person_in_charge: mickael.abbas at edf.fr
!
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/cfcpem.h"
#include "asterfort/cfcpes.h"
#include "asterfort/cfcpma.h"
#include "asterfort/cfdisd.h"
#include "asterfort/cfecrd.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmrvai.h"
    type(NL_DS_Measure), intent(inout) :: ds_measure
    character(len=24) :: resoco
    character(len=19) :: matass
    character(len=14) :: numedd
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (METHODES DISCRETES - RESOLUTION)
!
! ALGO. POUR CONTACT    : PENALISATION
! ALGO. POUR FROTTEMENT : SANS
!
! ----------------------------------------------------------------------
!
!
! RESOLUTION DE : C.DU + K ATA.DU  = F- K ATA(E-U)
!                 A(U+DU)      <= E (= POUR LES LIAISONS ACTIVES)
!
! AVEC E = JEU COURANT (CORRESPONDANT A U/I/N)
!
!      C = ( K  BT ) MATRICE DE RIGIDITE INCLUANT LES LAGRANGE
!          ( B  0  )
!
!      U = ( DEPL )
!          ( LAM  )
!
!      F = ( DL  ) DANS LA PHASE DE PREDICTION
!          ( DUD )
!
!      F = ( L - QT.SIG - BT.LAM  ) AU COURS D'UNE ITERATION DE NEWTON
!          (           0          )
!
! IO  ds_measure       : datastructure for measure and statistics management
! IN  RESOCO : SD DE TRAITEMENT NUMERIQUE DU CONTACT
! IN  NUMEDD : NUME_DDL
! IN  MATASS : NOM DE LA MATRICE DU PREMIER MEMBRE ASSEMBLEE
!
!
!
!
    integer(kind=8) :: ifm, niv
    character(len=24) :: afmu
    integer(kind=8) :: jafmu
    integer(kind=8) :: nbliai, neq, nbliac
    integer(kind=8) :: iter
    integer(kind=8) :: lmat
    character(len=19) :: matrcf
!
! ----------------------------------------------------------------------
!
    call jemarq()
    call infdbg('CONTACT', ifm, niv)
!
! --- AFFICHAGE
!
    if (niv .ge. 2) then
        write (ifm, *) '<CONTACT><CALC> ALGO_CONTACT   : PENALISATION'
        write (ifm, *) '<CONTACT><CALC> ALGO_FROTTEMENT: SANS'
    end if
!
! --- LECTURE DES STRUCTURES DE DONNEES DE CONTACT
!
    afmu = resoco(1:14)//'.AFMU'
    matrcf = resoco(1:14)//'.MATR'
    call jeveuo(afmu, 'E', jafmu)
!
! --- RECUPERATION DU DESCRIPTEUR DE LA MATRICE GLOBALE
!
    call jeveuo(matass//'.&INT', 'E', lmat)
!
! --- INITIALISATION DES VARIABLES
!
    nbliai = cfdisd(resoco, 'NBLIAI')
    neq = cfdisd(resoco, 'NEQ')
    nbliac = cfdisd(resoco, 'NBLIAC')
    iter = 1
!
! --- CREATION DU SECOND MEMBRE AFMU = -E_N*AT*JEU
!
    call cfcpes(resoco, jafmu)
!
    if (nbliac .eq. 0) then
        goto 999
    end if
!
! --- CALCUL DE LA MATRICE DE CONTACT PENALISEE ELEMENTAIRE [E_N*AT]
!
    call cfcpem(resoco, nbliai)
!
! --- CALCUL DE LA MATRICE DE CONTACT PENALISEE GLOBALE [E_N*AT*A]
!
    call cfcpma(resoco, neq, nbliai, numedd, matrcf)
!
999 continue
!
! --- ETAT DES VARIABLES DE CONTROLE DU CONTACT
!
    call cfecrd(resoco, 'NBLIAC', nbliac)
!
! --- SAUVEGARDE DES INFOS DE DIAGNOSTIC
!
    call nmrvai(ds_measure, 'Cont_Algo ', input_count=iter)
    call nmrvai(ds_measure, 'Cont_NCont', input_count=nbliac)
!
    call jedema()
!
end subroutine
