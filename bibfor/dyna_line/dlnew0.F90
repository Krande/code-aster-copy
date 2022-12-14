! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
! aslint: disable=W1504
!
subroutine dlnew0(result, force0, force1, iinteg, neq,&
                  istoc, iarchi, nbexci, nondp, nmodam,&
                  lamort, limped, lmodst, imat, masse,&
                  rigid, amort, nchar, nveca, liad,&
                  lifo, modele, mate, mateco, carele,&
                  charge, infoch, fomult, numedd, depla,&
                  vitea, accea, dep0, vit0, acc0,&
                  fexte, famor, fliai, depl1, vite1,&
                  acce1, psdel, fammo, fimpe, fonde,&
                  vien, vite, vita1, mltap, a0,&
                  a2, a3, a4, a5, a6,&
                  a7, a8, c0, c1, c2,&
                  c3, c4, c5, nodepl, novite,&
                  noacce, matres, maprec, solveu, criter,&
                  chondp, vitini, vitent, valmod, basmod,&
                  veanec, vaanec, vaonde, veonde, dt,&
                  theta, tempm, temps, iforc2, tabwk1,&
                  tabwk2, archiv, nbtyar, typear, numrep,&
                  ds_energy, kineLoad)
!
!     CALCUL MECANIQUE TRANSITOIRE PAR INTEGRATION DIRECTE
!     AVEC METHODES IMPLICITES :                  - THETA-WILSON
!                                                 - NEWMARK
! ----------------------------------------------------------------------
!  IN  : IINTEG    : ENTIER INDIQUANT LA METHODE D'INTEGRATION
!  IN  : NEQ       : NOMBRE D'EQUATIONS
!  IN  : ISTOC     : PILOTAGE DU STOCKAGE DES RESULTATS
!  IN  : IARCHI    : PILOTAGE DE L'ARCHIVAGE DES RESULTATS
!  IN  : NBEXCI    : NOMBRE D'EXCITATIONS
!  IN  : NONDP     : NOMBRE D'ONDES PLANES
!  IN  : NMODAM    : NOMBRE D'AMORTISSEMENTS MODAUX
!  IN  : LAMORT    : LOGIQUE INDIQUANT SI IL Y A AMORTISSEMENT
!  IN  : LIMPED    : LOGIQUE INDIQUANT SI
!  IN  : LMODST    : LOGIQUE INDIQUANT SI MODE STATIQUE
!  IN  : IMAT      : TABLEAU D'ADRESSES POUR LES MATRICES
!  IN  : MASSE     : MATRICE DE MASSE
!  IN  : NCHAR     : NOMBRE D'OCCURENCES DU MOT CLE CHARGE
!  IN  : NVECA     : NOMBRE D'OCCURENCES DU MOT CLE VECT_ASSE
!  IN  : LIAD      : LISTE DES ADRESSES DES VECTEURS CHARGEMENT (NVECT)
!  IN  : LIFO      : LISTE DES NOMS DES FONCTIONS EVOLUTION (NVECT)
!  IN  : MODELE    : NOM DU MODELE
!  IN  : MATE      : NOM DU CHAMP DE MATERIAU
!  IN  : CARELE    : CARACTERISTIQUES DES POUTRES ET COQUES
!  IN  : CHARGE    : LISTE DES CHARGES
!  IN  : INFOCH    : INFO SUR LES CHARGES
!  IN  : FOMULT    : LISTE DES FONC_MULT ASSOCIES A DES CHARGES
!  IN  : NUMEDD    : NUME_DDL DE LA MATR_ASSE RIGID
!  IN  : SOLVEU    : NOM DU SOLVEUR
!  IN  : CRITER    :
!  IN  : CHONDP    : NOMS DES ONDES PLANES
!  VAR : DEP0      : TABLEAU DES DEPLACEMENTS A L'INSTANT N
!  VAR : VIT0      : TABLEAU DES VITESSES A L'INSTANT N
!  VAR : ACC0      : TABLEAU DES ACCELERATIONS A L'INSTANT N
!  IN  : DT        : PAS DE TEMPS
!  IN  : THETA     : PARAMETRE DU SCHEMA TEMPOREL
!  IN  : TEMPM     : INSTANT PRECEDENT
!  IN  : TEMPS     : INSTANT COURANT
! IN  NUMREP : NUMERO DE REUSE POUR LA TABLE PARA_CALC
!
! CORPS DU PROGRAMME
!
!
    use NonLin_Datastructure_type
!
    implicit none
!
! DECLARATION PARAMETRES D'APPELS
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/assert.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dlarch.h"
#include "asterfort/dlfext.h"
#include "asterfort/enerca.h"
#include "asterfort/fimped.h"
#include "asterfort/fmodam.h"
#include "asterfort/fointe.h"
#include "asterfort/fondpl.h"
#include "asterfort/forcdy.h"
#include "asterfort/fteta.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/newacc.h"
#include "asterfort/newdep.h"
#include "asterfort/newvit.h"
#include "asterfort/nmarpc.h"
#include "asterfort/resoud.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsorac.h"
#include "asterfort/vtcopy.h"
#include "asterfort/vtcreb.h"
#include "asterfort/wkvect.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/dscal.h"
    integer :: nbexci, nondp, nmodam, iinteg, neq
    integer :: istoc, iarchi, imat(3), nchar, nveca, liad(*)
    integer :: iforc2, archiv, nbtyar, mltap(nbexci)
!
    real(kind=8) :: depla(neq), vitea(neq), accea(neq), dep0(*), vit0(*)
    real(kind=8) :: acc0(*), fexte(*), famor(*), fliai(*), depl1(neq)
    real(kind=8) :: vite1(neq), acce1(neq), psdel(neq), fammo(neq), fimpe(neq)
    real(kind=8) :: fonde(neq), vien(neq), vite(neq), vita1(neq), a0, a2, a3
    real(kind=8) :: a4, a5, a6, a7, a8, c0, c1, c2, c3, c4, c5, tabwk1(neq)
    real(kind=8) :: tabwk2(neq), dt, theta, tempm, temps
!
    aster_logical :: lamort, limped, lmodst
    type(NL_DS_Energy), intent(inout) :: ds_energy
!
    character(len=8) :: nodepl(nbexci), novite(nbexci), noacce(nbexci)
    character(len=8) :: masse, rigid, amort
    character(len=8) :: result
    character(len=19) :: force0, force1
    character(len=8) :: chondp(nondp)
    character(len=8) :: matres
    character(len=16) :: typear(nbtyar)
    character(len=19) :: solveu
    character(len=19) :: maprec
    character(len=24) :: criter, kineLoad
    character(len=24) :: modele, mate, mateco, carele, charge, infoch, fomult, numedd
    character(len=24) :: vitini, vitent, valmod, basmod
    character(len=24) :: lifo(*)
    character(len=24) :: veanec, vaanec, vaonde, veonde
!
!
! DECLARATION VARIABLES LOCALES
!
    character(len=6) :: nompro
    parameter ( nompro = 'DLNEW0' )
!
    integer :: iforc1
    integer :: nbexre, item2(1), iret, lvale, ibid, i
    integer :: ltps0, ltps1, nbinst, ifnobi, ifcibi, alarm
    integer :: iexci, ieq, iresu
!
    real(kind=8) :: coefd, coefv, coefa, prec, eps0, alpha
    integer :: numrep
    character(len=8) :: k8bid
    character(len=16) :: typa(6)
    character(len=19) :: chsol, cham19, chamno, chamn2, k19bid
    character(len=19) :: masse1, amort1, rigid1
    character(len=24) :: veccor, vecond
    complex(kind=8) :: cbid
    character(len=8), pointer :: listresu(:) => null()
    real(kind=8), pointer :: coef_rre(:) => null()
    real(kind=8), pointer :: forc0(:) => null()
    real(kind=8), pointer :: nlval1(:) => null()
    real(kind=8), pointer :: nlval2(:) => null()
!     -----------------------------------------------------------------
!
!
! --- NOM DES STRUCTURES DE TRAVAIL
!
    chsol = '&&'//nompro//'.SOLUTION '
    veccor = '&&VECCOR'
    vecond = '&&VECOND'
    chamno = '&&'//nompro//'.CHAMNO'
    call jeexin(chamno(1:19)//'.REFE', iret)
    if (iret .eq. 0) then
        call vtcreb(chamno, 'V', 'R', nume_ddlz = numedd, nb_equa_outz = neq)
    endif
    chamn2 = '&&'//nompro//'.CHAMN2'
    call jeexin(chamn2(1:19)//'.REFE', iret)
    if (iret .eq. 0) then
        call vtcreb(chamn2, 'V', 'R', nume_ddlz = numedd, nb_equa_outz = neq)
    endif
!====
! 2. DEPLACEMENT, VITESSE ET ACCELERATIONS A
!====
    depla(1:neq) = 0.d0
    vitea(1:neq) = 0.d0
    accea(1:neq) = 0.d0
!
    if (lmodst) then
!
        do iexci = 1, nbexci
!
            if (mltap(iexci) .eq. 1) then
                call fointe('F ', nodepl(iexci), 1, ['INST'], [temps],&
                            coefd, ieq)
                call fointe('F ', novite(iexci), 1, ['INST'], [temps],&
                            coefv, ieq)
                call fointe('F ', noacce(iexci), 1, ['INST'], [temps],&
                            coefa, ieq)
                do ieq = 1, neq
                    depla(ieq) = depla(ieq) + psdel(ieq)*coefd
                    vitea(ieq) = vitea(ieq) + psdel(ieq)*coefv
                    accea(ieq) = accea(ieq) + psdel(ieq)*coefa
                end do
            endif
!
        end do
!
    endif
!
!====
! 3.
!====
    do ieq = 1, neq
        vite(ieq) = vit0(ieq)
    end do
    if (lmodst) then
        do 32 , ieq = 1,neq
        vien(ieq) = vitea(ieq)
 32     continue
    endif
    if (limped) then
        call fimped(modele, mateco, numedd, neq, vitini,&
                    vitent, veccor, veanec, vaanec, tempm,&
                    fimpe)
    endif
    if (nondp .ne. 0) then
        call fondpl(modele, mate, mateco, numedd, neq,&
                    chondp, nondp, vecond, veonde, vaonde,&
                    temps, fonde)
    endif
!
    if (nmodam .ne. 0) then
        if (lmodst) then
            do 33 , ieq = 1,neq
            vita1(ieq) = vit0(ieq) + vitea(ieq)
 33         continue
            call fmodam(neq, vita1, valmod, basmod, fammo)
        else
            call fmodam(neq, vit0, valmod, basmod, fammo)
        endif
    endif
!
!====
! 4. CALCUL DU SECOND MEMBRE F*
!====
    call jeveuo(force0(1:19)//'.VALE', 'E', vr=forc0)
    call jeveuo(force1(1:19)//'.VALE', 'E', iforc1)
!
    call dlfext(nveca, nchar, temps, neq, liad,&
                lifo, charge, infoch, fomult, modele,&
                mate, mateco, carele, numedd, zr(iforc1))
!
    if (nondp .ne. 0) then
        do 43 , ieq = 1,neq
        zr(iforc1+ieq-1) = zr(iforc1+ieq-1) - fonde(ieq)
 43     continue
    endif
    if (ds_energy%l_comp) then
        do ieq = 1, neq
            fexte(ieq)=fexte(ieq+neq)
            fexte(ieq+neq)=zr(iforc1+ieq-1)
        end do
    endif
!
    if (limped) then
        do 41 , ieq = 1,neq
        zr(iforc1+ieq-1) = zr(iforc1+ieq-1) - fimpe(ieq)
 41     continue
        if (ds_energy%l_comp) then
            do ieq = 1, neq
                fliai(ieq)=fliai(ieq+neq)
                fliai(ieq+neq)=fimpe(ieq)
            end do
        endif
    endif
!
    if (nmodam .ne. 0) then
        do 42 , ieq = 1,neq
        zr(iforc1+ieq-1) = zr(iforc1+ieq-1) - fammo(ieq)
 42     continue
        if (ds_energy%l_comp) then
            do ieq = 1, neq
                famor(ieq)=famor(ieq+neq)
                famor(ieq+neq)=fammo(ieq)
            end do
        endif
    endif
!
!
!   Chargement venant d'un RESU a TEMPS
!
    call getfac('EXCIT_RESU', nbexre)
    if (nbexre .ne. 0) then
        call jeveuo('&&COMDLT.COEF_RRE', 'L', vr=coef_rre)
        call jeveuo('&&COMDLT.LISTRESU', 'L', vk8=listresu)
        prec=1.d-9
        eps0 =1.d-12
        do iresu = 1, nbexre
            if (abs(temps) .gt. eps0) then
                call rsorac(listresu(iresu), 'INST', 0, temps, k8bid,&
                            cbid, prec, 'RELATIF', item2, 1,&
                            ibid)
            else
                call rsorac(listresu(iresu), 'INST', 0, temps, k8bid,&
                            cbid, eps0, 'ABSOLU', item2, 1,&
                            ibid)
            endif
            if (ibid .gt. 0) then
                call rsexch('F', listresu(iresu), 'DEPL', item2(1), cham19,&
                            iret)
                call vtcopy(cham19, chamno, 'F', iret)
                call jeveuo(chamno//'.VALE', 'L', lvale)
!
            else
                call wkvect('&&DLNEW0.XTRAC', 'V V R8', neq, lvale)
                call jelira(listresu(iresu)//'           .ORDR', 'LONUTI', nbinst)
!
!        --- INTERPOLATION LINEAIRE ---
                do i = 1, nbinst-1
!
                    call rsadpa(listresu(iresu), 'L', 1, 'INST', i,&
                                0, sjv=ltps0, styp=k8bid)
                    call rsadpa(listresu(iresu), 'L', 1, 'INST', i+1,&
                                0, sjv=ltps1, styp=k8bid)
                    if (i .eq. 1 .and. temps .lt. zr(ltps0)) then
                        call rsexch('F', listresu(iresu), 'DEPL', i, cham19,&
                                    iret)
                        call vtcopy(cham19, chamno, 'F', iret)
                        call jeveuo(chamno//'.VALE', 'L', lvale)
                        goto 213
                    endif
                    if (temps .ge. zr(ltps0) .and. temps .lt. zr(ltps1)) then
                        alpha = (temps - zr(ltps0)) / (zr(ltps1) - zr( ltps0))
                        call rsexch('F', listresu(iresu), 'DEPL', i, cham19,&
                                    iret)
                        call vtcopy(cham19, chamno, 'F', iret)
                        call jeveuo(chamno//'.VALE', 'L', vr=nlval1)
                        call rsexch('F', listresu(iresu), 'DEPL', i+1, cham19,&
                                    iret)
                        call vtcopy(cham19, chamn2, 'F', iret)
                        call jeveuo(chamn2//'.VALE', 'L', vr=nlval2)
                        call dcopy(neq, nlval1, 1, zr(lvale), 1)
                        call dscal(neq, (1.d0-alpha), zr(lvale), 1)
                        call daxpy(neq, alpha, nlval2, 1, zr(lvale),&
                                   1)
                        goto 213
                    endif
                    if (i .eq. nbinst-1 .and. temps .ge. zr(ltps1)) then
                        call rsexch('F', listresu(iresu), 'DEPL', i+1, cham19,&
                                    iret)
                        call vtcopy(cham19, chamno, 'F', iret)
                        call jeveuo(chamno//'.VALE', 'L', lvale)
                        goto 213
                    endif
                end do
213             continue
            endif
            do ieq = 1, neq
                zr(iforc2+ieq-1) = zr(lvale+ieq-1)*coef_rre(iresu)
                zr(iforc1+ieq-1) = zr(iforc1+ieq-1) + zr(lvale+ieq-1)* coef_rre(iresu)
            end do
            if (ibid .gt. 0) then
                call jelibe(cham19//'.VALE')
            else
                call jelibe(cham19//'.VALE')
                call jedetr('&&DLNEW0.XTRAC')
            endif
        end do
        if (ds_energy%l_comp) then
            do ieq = 1, neq
                fexte(ieq+neq)=fexte(ieq+neq)+ zr(iforc2+ieq-1)
            end do
        endif
    endif
!
    if (iinteg .eq. 2) then
        call dcopy(neq, zr(iforc1), 1, zr(iforc2), 1)
        call fteta(theta, neq, forc0, zr(iforc1))
    endif
!
!====
! 5. FORCE DYNAMIQUE F*
!====
    call forcdy(imat(2), imat(3), lamort, neq, c0,&
                c1, c2, c3, c4, c5,&
                dep0, vit0, acc0, tabwk1, tabwk2,&
                zr(iforc1))
!
!====
! 6.  RESOLUTION DU PROBLEME K*  . U*  =  P*
!           --- RESOLUTION AVEC FORCE1 COMME SECOND MEMBRE ---
!====
    call resoud(matres, maprec, solveu, kineLoad, 0,&
                force1, chsol, 'V', [0.d0], [cbid],&
                criter, .true._1, 0, iret)
    call copisd('CHAMP_GD', 'V', chsol(1:19), force1(1:19))
    call jeveuo(force1(1:19)//'.VALE', 'E', iforc1)
    call detrsd('CHAMP_GD', chsol)
    call dcopy(neq, zr(iforc1), 1, depl1, 1)
!
!====
! 7. CALCUL DES DEPLACEMENTS,VITESSES ET ACCELERATIONS
!====
    if (iinteg .eq. 2) then
!
        call newacc(neq, a4, a5, a6, dep0,&
                    vit0, acc0, depl1, acce1)
        call newvit(neq, a7, a7, vit0, acc0,&
                    vite1, acce1)
        call newdep(neq, a8, dt, dep0, vit0,&
                    acc0, depl1, acce1)
!
    else if (iinteg.eq.1) then
!
        call newacc(neq, a0, -a2, -a3, dep0,&
                    vit0, acc0, depl1, acce1)
        call newvit(neq, a6, a7, vit0, acc0,&
                    vite1, acce1)
!
    endif
!
!====
! 8. CALCUL DES ENERGIES
!====
!
    if (ds_energy%l_comp) then
        masse1=masse//'           '
        amort1=amort//'           '
        rigid1=rigid//'           '
        ASSERT(kineLoad .eq. ' ')
        call wkvect('FNODABID', 'V V R', 2*neq, ifnobi)
        call wkvect('FCINEBID', 'V V R', 2*neq, ifcibi)
        call enerca(k19bid, dep0, vit0, depl1, vite1,&
                    masse1, amort1, rigid1, fexte, famor,&
                    fliai, zr(ifnobi), zr(ifcibi), lamort, .true._1,&
                    .false._1, ds_energy, '&&DLNEWI')
        call jedetr('FNODABID')
        call jedetr('FCINEBID')
    endif
!
!====
! 9. TRANSFERT DES NOUVELLES VALEURS DANS LES ANCIENNES
!====
!
    call dcopy(neq, depl1, 1, dep0, 1)
    call dcopy(neq, vite1, 1, vit0, 1)
    call dcopy(neq, acce1, 1, acc0, 1)
!
    if (iinteg .eq. 2) then
        call dcopy(neq, zr(iforc2), 1, forc0, 1)
    else
        call dcopy(neq, zr(iforc1), 1, forc0, 1)
    endif
!
!
!====
! 11. ARCHIVAGE EVENTUEL DANS L'OBJET SOLUTION
!====
    if (archiv .eq. 1) then
!
        istoc = 0
        alarm = 1
!
        if (lmodst) then
!
            typa(1) = 'DEPL_ABSOLU'
            typa(2) = 'VITE_ABSOLU'
            typa(3) = 'ACCE_ABSOLU'
            typa(4) = '    '
            typa(5) = '    '
            typa(6) = '    '
            do 101 , ieq = 1,neq
            depla(ieq) = depla(ieq) + dep0(ieq)
            vitea(ieq) = vitea(ieq) + vit0(ieq)
            accea(ieq) = accea(ieq) + acc0(ieq)
101         continue
            call dlarch(result, neq, istoc, iarchi, ' ',&
                        alarm, temps, nbtyar, typa, masse,&
                        depla, vitea, accea, fexte( neq+1), famor(neq+1),&
                        fliai(neq+1))
!
        endif
!
        call dlarch(result, neq, istoc, iarchi, ' ',&
                    alarm, temps, nbtyar, typear, masse,&
                    dep0, vit0, acc0, fexte(neq+1), famor(neq+1),&
                    fliai(neq+1))
!
    endif
!
    call nmarpc(ds_energy, numrep, temps)
!
end subroutine
