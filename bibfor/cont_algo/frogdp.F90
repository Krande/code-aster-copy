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

subroutine frogdp(ds_measure, resoco, numedd, matass, resigr)
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
#include "asterfort/cffpfo.h"
#include "asterfort/cffpm1.h"
#include "asterfort/cffpm2.h"
#include "asterfort/cffrot.h"
#include "asterfort/cfmata.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mrmult.h"
#include "asterfort/mtdscr.h"
#include "asterfort/nmrvai.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    real(kind=8) :: resigr
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
! ALGO. POUR FROTTEMENT : PENALISATION
!
! ----------------------------------------------------------------------
!
!
! RESO. DE : C.DU + KC ACT.AC.DU + KG AGT.AG.DU = F - KG AGT.AG (E-U)
!            AC. (U+DU)      <= E  (= POUR LES LIAISONS ACTIVES)
!
! AVEC E = JEU COURANT (CORRESPONDANT A U/I/N)
!
!     AC = MATRICE DE CONTACT
!
!    AG  = MATRICE DE FROTTEMENT POUR LES NOEUDS GLISSANTS
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
! IN  RESIGR : RESI_GLOB_RELA
! ON UTILISE UNIQUEMENT LE VECTEUR AFMU CAR LES DONNEES DE ATMU SONT
! NECESSAIRE POUR LE CALCUL DE LA MATRICE TANGENTE QUI SE FAIT
! A L'AIDE DU VECTEUR AFMU
!
!
!
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: nmult
    integer(kind=8) :: ieq
    integer(kind=8) :: neq, nbliac, nbliai, ndim
    integer(kind=8) :: lmat
    integer(kind=8) :: lmaf1, iter
    integer(kind=8) :: nesmax
    character(len=14) :: numef1, numef2, nutemp
    character(len=19) :: maf1, maf2, matemp, mact
    character(len=14) :: numecf
    character(len=19) :: matrcf, fro1, fro2
    character(len=19) :: atmu, afmu
    integer(kind=8) :: jatmu, jafmu
    character(len=19) :: depl0
    real(kind=8), pointer :: vale(:) => null()
!
! ----------------------------------------------------------------------
!
    call jemarq()
    call infdbg('CONTACT', ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<CONTACT><CALC> ALGO_CONTACT   : PENALISATION'
        write (ifm, *) '<CONTACT><CALC> ALGO_FROTTEMENT: PENALISATION'
    end if
!
! --- LECTURE DES STRUCTURES DE DONNEES DE CONTACT
!
    atmu = resoco(1:14)//'.ATMU'
    afmu = resoco(1:14)//'.AFMU'
    call jeveuo(atmu, 'E', jatmu)
    call jeveuo(afmu, 'E', jafmu)
!
! --- MATRICES DE FROTTEMENT
!
    maf1 = '&&FROGDP.MAF1'
    numef1 = '&&FROGDP.NUF1'
    fro1 = resoco(1:14)//'.FRO1'
    numef2 = '&&FROGDP.NUF2'
    maf2 = '&&FROGDP.MAF2'
    fro2 = resoco(1:14)//'.FRO2'
!
! --- ACCES AUX CHAMPS DE TRAVAIL
! --- DEPL0: INCREMENT DE DEPLACEMENT CUMULE DEPUIS DEBUT
! ---        DU PAS DE TEMPS SANS CORRECTION DU CONTACT
!
    depl0 = resoco(1:14)//'.DEP0'
    call jeveuo(depl0(1:19)//'.VALE', 'L', vr=vale)
!
! --- RECUPERATION DU DESCRIPTEUR DE LA MATRICE GLOBALE
!
    call jeveuo(matass//'.&INT', 'E', lmat)
!
! --- INITIALISATIONS DES VARIABLES
!
    nbliai = cfdisd(resoco, 'NBLIAI')
    neq = cfdisd(resoco, 'NEQ')
    ndim = cfdisd(resoco, 'NDIM')
    nesmax = cfdisd(resoco, 'NESMAX')
    nbliac = cfdisd(resoco, 'NBLIAC')
    iter = 1
!
! --- CREATION DU SECOND MEMBRE ATMU = -E_N.[Ac]T.{JEU}
!
    call cfcpes(resoco, jatmu)
!
    if (nbliac .eq. 0) then
        goto 999
    end if
!
! --- CALCUL DES COEFFICIENTS DE LAGRANGE MU POUR LE FROTTEMENT
!
    call cffpfo(resoco, nbliai, nbliac, ndim)
!
! --- CALCUL DE LA MATRICE DE CONTACT PENALISEE "ELEMENTAIRE" [E_N*AcT]
!
    call cfcpem(resoco, nbliai)
!
! --- CALCUL DE LA MATRICE DE CONTACT PENALISEE "GLOBALE" [E_N*AcT*Ac]
!
    mact = '&&FROGDP.MACT'
    call cfcpma(resoco, neq, nbliai, numedd, mact)
!
! --- CREATION DE LA MATRICE FRO1 = E_T*AaT
!
    call cffpm1(resoco, nbliai, ndim, nesmax)
!
! --- CREATION DE LA MATRICE MAF1 = E_T*AaT*Aa
! --- CETTE MATRICE NE SERT QU'AU CALCUL DU SECOND MEMBRE
!
    nmult = ndim-1
    call cfmata(resoco, neq, nbliai, nmult, numedd, &
                fro1, numef1, maf1)
!
! --- RECUPERATION DU SECOND MEMBRE E_T*AaT*Aa * DELTA -> AFMU
!
    call mtdscr(maf1)
    call jeveuo(maf1//'.&INT', 'L', lmaf1)
    call mrmult('ZERO', lmaf1, vale, zr(jafmu), 1, &
                .true._1)
!
! --- CREATION DE FRO2 = E_T*AT
!
    call cffpm2(resoco, resigr, nbliai, nbliac, ndim)
!
! --- CREATION DE LA SECONDE PARTIE DE LA MATRICE DE FROTTEMENT MAF2
!
    nmult = 1
    call cfmata(resoco, neq, nbliai, nmult, numedd, &
                fro2, numef2, maf2)
!
! --- CALCUL DE LA MATRICE TANGENTE MAFROT = MACT+MAF1+MAF2
!
    matemp = '&&FROGDP.MATP'
    nutemp = '&&FROGDP.NUTP'
    matrcf = resoco(1:14)//'.MATR'
    numecf = '&&FROGDP.NUFR'
    call cffrot(maf1, '-', maf2, matemp, nutemp)
    call cffrot(mact, '+', matemp, matrcf, numecf)
!
! --- CALCUL DES FORCES DE CONTACT (AT.MU) ET FROTTEMENT (AF.MU)
!
    do ieq = 1, neq
        zr(jafmu-1+ieq) = zr(jafmu-1+ieq)+zr(jatmu-1+ieq)
        zr(jatmu-1+ieq) = 0.d0
    end do
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
! ======================================================================
!
end subroutine
