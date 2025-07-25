! --------------------------------------------------------------------
! Copyright (C) 2005 IFP - MARTIN GUITTON         WWW.CODE-ASTER.ORG
! Copyright (C) 2007 - 2025 - EDF R&D - www.code-aster.org
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
subroutine algogl(ds_measure, sdcont_defi, sdcont_solv, solveu, matass, &
                  noma, ctccvg)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cfacat.h"
#include "asterfort/cfaduc.h"
#include "asterfort/cfatmu.h"
#include "asterfort/cfdisd.h"
#include "asterfort/cfdisi.h"
#include "asterfort/cfecrd.h"
#include "asterfort/cffact.h"
#include "asterfort/cfimp2.h"
#include "asterfort/cfmajc.h"
#include "asterfort/cfpeti.h"
#include "asterfort/cfreso.h"
#include "asterfort/cftabl.h"
#include "asterfort/elpiv1.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mtdsc3.h"
#include "asterfort/nmrvai.h"
#include "blas/daxpy.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    type(NL_DS_Measure), intent(inout) :: ds_measure
    character(len=8) :: noma
    character(len=24) :: sdcont_defi, sdcont_solv
    character(len=19) :: solveu, matass
    integer(kind=8) :: ctccvg
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (METHODES DISCRETES - RESOLUTION)
!
! ALGO. POUR CONTACT    : CONTRAINTES TOUTES ACTIVES CAR GLISSIERE
! ALGO. POUR FROTTEMENT : SANS
! ALGO GLISSIERE
!
! ----------------------------------------------------------------------
!
!
! RESOLUTION DE : C.DU + AT.MU  = F
!                 A(U+DU)      <= E (POUR LES LIAISONS ACTIVES)
!
! AVEC E = JEU COURANT (CORRESPONDANT A U/I/N)
!
!      A = MATRICE DE CONTACT
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
! IN  DEFICO : SD DE DEFINITION DU CONTACT
! IN  RESOCO : SD DE TRAITEMENT NUMERIQUE DU CONTACT
! IN  SOLVEU : SD SOLVEUR
! IN  MATASS : NOM DE LA MATRICE DU PREMIER MEMBRE ASSEMBLEE
! IN  NOMA   : NOM DU MAILLAGE
! OUT CTCCVG : CODE RETOUR CONTACT DISCRET
!                -1 : PAS DE CALCUL DU CONTACT DISCRET
!                 0 : CAS DU FONCTIONNEMENT NORMAL
!                 1 : NOMBRE MAXI D'ITERATIONS
!                 2 : MATRICE SINGULIERE
!
!
!
!
    integer(kind=8) :: ifm, niv
    aster_logical :: lechec
    integer(kind=8) :: ieq, iter
    integer(kind=8) :: llliai, llliac
    integer(kind=8) :: indic, indfac, ajliai, spliai, spavan
    integer(kind=8) :: neq, nbliac, nbliai, ndim, nesmax
    real(kind=8) :: rho, xjvmax
    character(len=1) :: typeaj
    character(len=19) :: macont
    integer(kind=8) :: ldscon, lmat
    character(len=19) :: ddeplc, ddepl0, ddelt
    integer(kind=8) :: itemax, isto, itemul
    real(kind=8), pointer :: vddelt(:) => null()
    real(kind=8), pointer :: ddep0(:) => null()
    real(kind=8), pointer :: ddepc(:) => null()
    blas_int :: b_incx, b_incy, b_n
!
! ----------------------------------------------------------------------
!
    call jemarq()
    call infdbg('CONTACT', ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<CONTACT><CALC> ALGO_CONTACT   : CONT. ACTIVES'
    end if
!
! --- ACCES AUX CHAMPS DE TRAVAIL
! --- DDEPL0: INCREMENT DE SOLUTION SANS CORRECTION DU CONTACT
! --- DDEPLC: INCREMENT DE SOLUTION APRES CORRECTION DU CONTACT
! --- DDELT : INCREMENT DE SOLUTION ITERATION DE CONTACT
!
    ddepl0 = sdcont_solv(1:14)//'.DEL0'
    ddeplc = sdcont_solv(1:14)//'.DELC'
    ddelt = sdcont_solv(1:14)//'.DDEL'
    call jeveuo(ddepl0(1:19)//'.VALE', 'L', vr=ddep0)
    call jeveuo(ddeplc(1:19)//'.VALE', 'E', vr=ddepc)
    call jeveuo(ddelt(1:19)//'.VALE', 'E', vr=vddelt)
!
! --- PREPARATION DE LA MATRICE DE CONTACT
!
    macont = sdcont_solv(1:14)//'.MATC'
    call mtdsc3(macont)
    call jeecra(macont(1:19)//'.REFA', 'DOCU', cval='ASSE')
    call jeveuo(macont(1:19)//'.&INT', 'E', ldscon)
!
! --- RECUPERATION DU DESCRIPTEUR DE LA MATRICE GLOBALE
!
    call jeveuo(matass//'.&INT', 'L', lmat)
!
! --- INITIALISATION DES VARIABLES
!
    nbliai = cfdisd(sdcont_solv, 'NBLIAI')
    neq = cfdisd(sdcont_solv, 'NEQ')
    ndim = cfdisd(sdcont_solv, 'NDIM')
    nesmax = cfdisd(sdcont_solv, 'NESMAX')
    nbliac = cfdisd(sdcont_solv, 'NBLIAC')
    itemul = 2
    itemax = itemul*nbliai
    typeaj = 'A'
    xjvmax = 0.d0
!
! --- GESTION DE LA FACTORISATION
!
    isto = cfdisi(sdcont_defi, 'STOP_SINGULIER')
    ajliai = 0
    spliai = 0
    if (nbliac .gt. 0) then
        indic = 1
    else
        indic = 0
    end if
    indfac = 1
!
    iter = 1
!
    if (niv .ge. 2) then
        write (ifm, 210) itemax
    end if
!
! ======================================================================
!                    REPRISE DE LA BOUCLE PRINCIPALE
! ======================================================================
!
40  continue
!
! --- MISE A JOUR DE LA SOLUTION ITERATION DE CONTACT
!
    do ieq = 1, neq
        vddelt(ieq) = ddep0(ieq)-ddepc(ieq)
    end do
!
! --- RESOLUTION MATRICIELLE POUR DES LIAISONS ACTIVES
!
    if (nbliac .ne. 0) then
!
        spavan = spliai
!
! ----- CALCUL DE [-A.C-1.AT] COLONNE PAR COLONNE (A PARTIR DE INDFAC)
!
        call cfacat(indic, nbliac, ajliai, spliai, indfac, &
                    sdcont_defi, sdcont_solv, solveu, lmat, xjvmax)
!
! ----- DETECTION DES PIVOTS NULS
!
        call elpiv1(xjvmax, indic, nbliac, ajliai, spliai, &
                    spavan, noma, sdcont_defi, sdcont_solv)
!
! ----- ON A SUPPRIME UNE LIAISON
!
        if (indic .eq. -1) then
            goto 150
        end if
!
! ----- FACTORISATION LDLT DE [-A.C-1.AT]
!
        call cffact(ldscon, isto, nbliac, indfac, lechec)
!
! ----- LA MATRICE DE CONTACT EST-ELLE SINGULIERE ?
!
        if (lechec) then
            ctccvg = 2
            goto 999
        end if
!
! ----- CALCUL DU SECOND MEMBRE : {JEU(DEPTOT) - A.DDEPL0} -> {MU}
!
        call cfaduc(sdcont_solv, nbliac)
!
! ----- RESOLUTION : [-A.C-1.AT].{MU} = {JEU(DEPTOT) - A.DDEPL0}
!
        call cfreso(sdcont_solv, ldscon, nbliac)
!
! ----- MISE A JOUR DU VECTEUR SOLUTION ITERATION DE CONTACT
!
        call cfmajc(sdcont_solv, neq, nbliac)
    end if
!
! --- VERIFICATION SI L'ENSEMBLE DES LIAISONS SUPPOSEES EST TROP PETIT
!
    call cfpeti(sdcont_solv, neq, nbliai, nbliac, rho, &
                llliai, llliac)
!
! --- ACTUALISATION DE DDEPLC = DDEPLC + RHO .DDELT
!
    b_n = to_blas_int(neq)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, rho, vddelt, b_incx, ddepc, &
               b_incy)
!
! --- MODIFICATIONS DES LIAISONS
!
    if (rho .lt. 1.d0) then
!
! ----- SI AU MOINS UNE LIAISON SUPPOSEE NON ACTIVE EST VIOLEE
! ----- ON AJOUTE A L'ENSEMBLE DES LIAISONS ACTIVES LA PLUS VIOLEE
!
        ASSERT(llliai .gt. 0)
        call cftabl(indic, nbliac, ajliai, spliai, sdcont_solv, &
                    typeaj, llliac, llliai)
        call cfimp2(sdcont_defi, sdcont_solv, noma, llliai, 'ACT')
    else
!
! ----- ON A CONVERGE
!
        goto 160
    end if
!
! --- ON PASSE A L'ITERATION DE CONTRAINTES ACTIVES SUIVANTE
!
150 continue
!
    iter = iter+1
!
! --- A-T-ON DEPASSE LE NOMBRE D'ITERATIONS DE CONTACT AUTORISE ?
!
    if (iter .gt. itemax+1) then
        ctccvg = 1
        goto 999
    end if
!
    goto 40
!
! ======================================================================
!                            ON A CONVERGE
! ======================================================================
!
160 continue
!
! --- CALCUL DES FORCES DE CONTACT (AT.MU)
!
    call cfatmu(neq, nbliac, sdcont_solv)
!
! --- CODE RETOUR
!
    ctccvg = 0
!
! --- ETAT DES VARIABLES DE CONTROLE DU CONTACT
!
    call cfecrd(sdcont_solv, 'NBLIAC', nbliac)
!
    if (niv .ge. 2) then
        write (ifm, 220) iter
    end if
!
999 continue
!
! --- SAUVEGARDE DES INFOS DE DIAGNOSTIC
!
    call nmrvai(ds_measure, 'Cont_Algo ', input_count=iter)
    call nmrvai(ds_measure, 'Cont_NCont', input_count=nbliac)
!
    call jedema()
!
210 format(' <CONTACT><CALC> DEBUT DES ITERATIONS (MAX: ', i6, ')')
220 format(' <CONTACT><CALC> FIN DES ITERATIONS (NBR: ', i6, ')')
!
end subroutine
