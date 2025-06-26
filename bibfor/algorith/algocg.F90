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
subroutine algocg(ds_measure, defico, resoco, solveu, matass, &
                  ctccvg)
!
! person_in_charge: mickael.abbas at edf.fr
!
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/cfdisd.h"
#include "asterfort/cfdisi.h"
#include "asterfort/cfdisr.h"
#include "asterfort/cfecrd.h"
#include "asterfort/cfgcac.h"
#include "asterfort/cfgccj.h"
#include "asterfort/cfgcin.h"
#include "asterfort/cfgcpc.h"
#include "asterfort/cfgcpr.h"
#include "asterfort/cfgcrl.h"
#include "asterfort/cfgcsg.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jedupo.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmrvai.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
    type(NL_DS_Measure), intent(inout) :: ds_measure
    character(len=24) :: defico, resoco
    character(len=19) :: matass, solveu
    integer(kind=8) :: ctccvg
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (METHODES DISCRETES - RESOLUTION)
!
! ALGO. POUR CONTACT    : GRADIENT CONJUGUE PROJETE
! ALGO. POUR FROTTEMENT : SANS
!
! ----------------------------------------------------------------------
!
!
!
! RESOLUTION DE : C.DU + AT.MU  = F
!                 A(U+DU)      <= E (= POUR LES LIAISONS ACTIVES)
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
    aster_logical :: conjug
    integer(kind=8) :: iliai, iter, premax
    integer(kind=8) :: neq, nbliac, nbliai
    integer(kind=8) :: gcpmax
    character(len=16) :: precon, search, pceffe
    character(len=19) :: sgradm, sgradp, sgrprm, sgrprp, mum
    integer(kind=8) :: jsgram, jsgrap, jsgprm, jsgprp, jmum
    character(len=19) :: mu
    integer(kind=8) :: jmu
    character(len=19) :: ddeplc, ddelt
    real(kind=8) :: tole, coefrs
    real(kind=8) :: ninf, ninfpc, alpha, epsi
    real(kind=8), pointer :: vddelt(:) => null()
    real(kind=8), pointer :: ddepc(:) => null()
    blas_int :: b_incx, b_incy, b_n
!
! ----------------------------------------------------------------------
!
    call jemarq()
    call infdbg('CONTACT', ifm, niv)
!
! --- AFFICHAGE
!
    if (niv .ge. 2) then
        write (ifm, *) '<CONTACT><CALC> ALGO_CONTACT   : GRADIENT CONJUGUE PROJETE'
        write (ifm, *) '<CONTACT><CALC> ALGO_FROTTEMENT: SANS'
    end if
!
! --- INITIALISATION DES VARIABLES
!
    nbliai = cfdisd(resoco, 'NBLIAI')
    neq = cfdisd(resoco, 'NEQ')
    ctccvg = 0
!
! --- LECTURE DES STRUCTURES DE DONNEES DE CONTACT
!
    mu = resoco(1:14)//'.MU'
    sgradm = resoco(1:14)//'.SGDM'
    sgradp = resoco(1:14)//'.SGDP'
    sgrprm = resoco(1:14)//'.SGPM'
    sgrprp = resoco(1:14)//'.SGPP'
    mum = resoco(1:14)//'.MUM'
    call jeveuo(mu, 'E', jmu)
    call jeveuo(sgradm, 'E', jsgram)
    call jeveuo(sgradp, 'E', jsgrap)
    call jeveuo(sgrprm, 'E', jsgprm)
    call jeveuo(sgrprp, 'E', jsgprp)
!
! --- ACCES AUX CHAMPS DE TRAVAIL
! --- DDEPLC: INCREMENT DE SOLUTION APRES CORRECTION DU CONTACT
! --- DDELT : INCREMENT DE SOLUTION ITERATION DE CONTACT
!
    ddeplc = resoco(1:14)//'.DELC'
    ddelt = resoco(1:14)//'.DDEL'
!
! --- INITIALISATION DES VECTEURS DE TRAVAIL
!
    call jedupo(mu, 'V', mum, .false._1)
    call jeveuo(mum, 'E', jmum)
!
! ======================================================================
!                             INITIALISATIONS
! ======================================================================
!
    iter = 1
    conjug = .false.
    tole = r8prem()
    ninfpc = 0.d0
!
! --- RECUPERATION DU CRITERE DE CONVERGENCE
!
    epsi = cfdisr(defico, 'RESI_ABSO')
    coefrs = cfdisr(defico, 'COEF_RESI')
    gcpmax = 10*nbliai
    premax = cfdisi(defico, 'ITER_PRE_MAXI')
    if (cfdisi(defico, 'ITER_GCP_MAXI') .ne. 0) then
        gcpmax = max(gcpmax, cfdisi(defico, 'ITER_GCP_MAXI'))
    end if
    if (cfdisi(defico, 'PRE_COND') .eq. 1) then
        precon = 'DIRICHLET'
    else
        precon = 'SANS'
    end if
    if (cfdisi(defico, 'RECH_LINEAIRE') .eq. 1) then
        search = 'NON_ADMISSIBLE'
    else
        search = 'ADMISSIBLE'
    end if
!
    if (niv .ge. 2) then
        write (ifm, 9010) gcpmax
    end if
!
! --- INITIALISATION A PARTIR DU CHAMP DE MULTIPLICATEURS INITIAL
!
    call cfgcin(resoco, matass, solveu, neq, nbliai)
!
! ======================================================================
!                    REPRISE DE LA BOUCLE PRINCIPALE
! ======================================================================
!
40  continue
!
    if (niv .eq. 2) then
        write (ifm, *) '<CONTACT><CALC> --------------------------------'
        write (ifm, *) '<CONTACT><CALC> ITERATION DE GCP = ', iter
    end if
!
! --- CALCUL DU SOUS-GRADIENT
!
    call cfgcsg(resoco, neq, nbliai, tole, ninf)
!
! --- A-T-ON CONVERGE ?
!
    if (niv .eq. 2) then
        write (ifm, 9060) ninf, epsi
    end if
!
    if (ninf .lt. epsi) then
        goto 160
    end if
!
! --- PRECONDITIONNEMENT UNIQUEMENT AU VOISINAGE DE LA SOLUTION
! --- LE VOISINAGE EST DEFINI PAR COEF_RESI
!
    if (iter .eq. 1) then
        ninfpc = coefrs*ninf
    end if
    pceffe = precon
    if (ninfpc .gt. 0.d0) then
        if (ninf .gt. ninfpc) then
            pceffe = 'SANS'
        end if
    end if
!
! --- PRECONDITIONNEMENT
!
    call cfgcpc(resoco, matass, solveu, neq, nbliai, &
                pceffe, tole, premax, epsi)
!
! --- CONJUGAISON
!
    call cfgccj(resoco, nbliai, conjug)
!
! --- RECHERCHE LINEAIRE: PAS D'AVANCEMENT
!
    call cfgcrl(resoco, neq, nbliai, matass, solveu, &
                alpha)
!
! --- PROJECTION DU PAS D'AVANCEMENT
!
    call cfgcpr(resoco, matass, solveu, neq, nbliai, &
                search, alpha)
!
! --- ACTUALISATION DE {DDEPLC} = {DDEPLC} - ALPHA . {DDELT}
!
    call jeveuo(ddelt(1:19)//'.VALE', 'L', vr=vddelt)
    call jeveuo(ddeplc(1:19)//'.VALE', 'E', vr=ddepc)
    b_n = to_blas_int(neq)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, -alpha, vddelt, b_incx, ddepc, &
               b_incy)
!
! --- ON VERIFIE SI L'ETAT DE CONTACT A CHANGE (ON NE CONJUGUE PAS)
!
    conjug = .true.
    do iliai = 1, nbliai
        if (((zr(jmum-1+iliai) .le. tole) .and. (zr(jmu-1+iliai) .gt. tole)) .or. &
            ((zr(jmum-1+iliai) .gt. tole) .and. (zr(jmu-1+iliai) .le. tole))) then
            conjug = .false.
            if (niv .eq. 2) then
                write (ifm, *) '<CONTACT><CALC>'//&
     &        ' CHANGEMENT DE L''ETAT DE CONTACT'
            end if
            goto 100
        end if
    end do
100 continue
!
! --- MISE A JOUR DES GRADIENTS ET DES DIRECTIONS DE RECHERCHE
!
    b_n = to_blas_int(nbliai)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, zr(jsgrap), b_incx, zr(jsgram), b_incy)
    b_n = to_blas_int(nbliai)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, zr(jsgprp), b_incx, zr(jsgprm), b_incy)
    b_n = to_blas_int(nbliai)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, zr(jmu), b_incx, zr(jmum), b_incy)
!
! --- ON PASSE A L'ITERATION SUIVANTE
!
    iter = iter+1
!
! --- A-T-ON DEPASSE LE NOMBRE D'ITERATIONS DE CONTACT AUTORISE ?
!
    if (iter .ge. gcpmax) then
        ctccvg = 1
        goto 160
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
! --- ACTIVATION DES LIAISONS ET CALCUL DE LA FORCE DE CONTACT
!
    call cfgcac(resoco, tole, neq, nbliai, nbliac)
!
! --- ETAT DES VARIABLES DE CONTROLE DU CONTACT
!
    call cfecrd(resoco, 'NBLIAC', nbliac)
!
    if (niv .ge. 2) then
        write (ifm, 9020) iter
    end if
!
! --- SAUVEGARDE DES INFOS DE DIAGNOSTIC
!
    call nmrvai(ds_measure, 'Cont_Algo ', input_count=iter)
    call nmrvai(ds_measure, 'Cont_NCont', input_count=nbliac)
!
    call jedema()
!
9010 format(' <CONTACT><CALC> DEBUT DES ITERATIONS (MAX: ', i8, ')')
9020 format(' <CONTACT><CALC> FIN DES ITERATIONS (NBR: ', i8, ')')
9060 format(' <CONTACT><CALC> NORME INFINIE DU RESIDU : ', 1pe12.5, ' (CRITERE: ', 1pe12.5, ')')
end subroutine
