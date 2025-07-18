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
subroutine modint(ssami, raiint, nddlin, nbmod, shift, &
                  matmod, masse, raide, neq, coint, &
                  noddli, nnoint, vefreq, switch)
    implicit none
!    M. CORUS     DATE 05/02/10
!-----------------------------------------------------------------------
!
!  BUT:      < CALCUL DES MODES D'INTERFACE >
!
!  ON EXTRAIT, A PARTIR DE LA SOUS MATRICE ASSOCIEE A L'INTERFACE, LA
!  CONNECTIVITE DU TREILLIS DE POUTRE SOUS JACENT. ON DETERMINE LE
!  NOMBRE DE PARTIES INDEPENDANTES DE L'INTERFACE, ET ON CALCULE, POUR
!  CHAQUE PARTIE, LES PREMIERS MODES PROPRES, EN PRENANT SOIN DE BIEN
!  CAPTER LES MODES DE CORPS RIGIDE. ON CONSTRUIT ENSUITE, SUR LA BASE
!  DE CES MODES, UN SOUS ESPACE POUR PROJETER LES MATRICES DU PROBLEME
!  COMPLET, ET ON CALCULE, SUR CE SOUS ESPACE, LES MODES D'INTERFACE.
!
!-----------------------------------------------------------------------
!  IN  : SSAMI   : MATRICE DE MASSE DU MODELE D'INTERFACE
!  IN  : RAIINT  : MATRICE DE RAIDEUR DU MODELE D'INTERFACE
!  IN  : NDDLIN  : NOMBRE D'EQUATIONS DU NUME_DDL D'INTERFACE
!  IN  : NBMOD   : NOMBRE DE MODES D'INTERFACE DEMANDE
!  IN  : SHIFT   : VALEUR DE FREQUENCE POUR LE DECALAGE DE LA RAIDEUR
!  OUT : MATMOD  : MATRICE DES MODES D'INTERFACE
!  IN  : MASSE   : MATRICE DE MASSE DU MODELE COMPLET
!  IN  : RAIDE   : MATRICE DE RAIDEUR DU MODELE COMPLET
!  IN  : NEQ     : NOMBRE D'EQUATIONS DU NUME_DDL COMPLET
!  IN  : COINT   : DEFINITION DE LA CONNECTIVITE DE L'INTERFACE
!  IN  : NODDLI  : DEFINITION DES DDL PORTES PAR LES NOEUDS D'INTERFACE
!  IN  : NNOINT  : NOMBRE DE NOEUD A L'INTERFACE
!  OUT : VEFREQ  : NOM DU VECTEUR CONTENANT LES FREQUENCES PROPRES
!  IN  : SWITCH  : 1 SI ON DOIT RELEVER LES MODES SUR TOUTE LA SST
!                  0 SI ON NE CONSERVE QUE LES MODES SUR L'INTERFACE
!
!-----------------------------------------------------------------------
!
!
!
!
!     ------------------------------------------------------------------
!
!-- VARIABLES EN ENTREES / SORTIE
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/isnnem.h"
#include "asterc/matfpe.h"
#include "asterc/r8vide.h"
#include "asterc/r8pi.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvr8.h"
#include "asterfort/intdis.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetc.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/infniv.h"
#include "asterfort/nmop45.h"
#include "asterfort/omega2.h"
#include "asterfort/mrmult.h"
#include "asterfort/mtcmbl.h"
#include "asterfort/mtdefs.h"
#include "asterfort/mtdscr.h"
#include "asterfort/preres.h"
#include "asterfort/resoud.h"
#include "asterfort/utmess.h"
#include "asterfort/vpcres.h"
#include "asterfort/vpleci.h"
#include "asterfort/wkvect.h"
#include "asterfort/rsadpa.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "blas/ddot.h"
#include "blas/dggev.h"
    integer(kind=8) :: nddlin, nbmod, nnoint, neq, switch, jfreq
    real(kind=8) :: shift
    character(len=6) :: k6bid
    character(len=8) :: modes
    character(len=19) :: masse, raide, ssami, raiint
    character(len=24) :: coint, noddli, matmod, vefreq
!
!
!-- VARIABLES DE LA ROUTINE
    aster_logical :: l_hpp
    integer(kind=8) :: lmatmo, i1, j1, k1, m1, lmakry, nsekry, nsekry2, nsekry3, nbborn
    integer(kind=8) :: lmatk, lmatm, lmapro, nbrss, lkpro, lmatrm, lmatrk, lwork
    integer(kind=8) :: limped, lmatma, iret, nbvect, ibid, no, nbsst, lindin, coeff, lvp
    integer(kind=8) :: ifm, niv, mode_symetrique, maxitr, nbvec2
    integer(kind=4) :: info
    real(kind=8) :: temp, rbid, norm, lambda, comlin(2), swork(1), max, omecor
    real(kind=8) :: bande(2), freq1, freq2, alpha, tolsor, precsh, fcorig, precdc
    real(kind=8) :: mval, normx, valx, avalx, r8bid, vrbid(1), pi
    complex(kind=8) :: cbid
    character(len=1) :: listyp(2), k1bid
    character(len=4) :: mod45
    character(len=8) :: k8bid, method, sdstab, arret
    character(len=16) :: typres, k16bid, optiof, sturm, modrig, stoper, typcal
    character(len=19) :: lismat(2), imped, solveu, nume91, nume, prno, eigsol, k19bid
    character(len=19) :: raide2, masse2
    character(len=24) :: valk, k24bid
!
    integer(kind=8), pointer :: v_ind_lag(:) => null()
    integer(kind=8), pointer :: delg(:) => null()
    integer(kind=8), pointer :: ddl_actif_int(:) => null()
    integer(kind=8), pointer :: v_ind_f_pro(:) => null()
    real(kind=8), pointer :: matr_mod_red(:) => null()
    real(kind=8), pointer :: matr_work_dggev(:) => null()
    real(kind=8), pointer :: vect_alphai(:) => null()
    real(kind=8), pointer :: vect_alphar(:) => null()
    real(kind=8), pointer :: vect_beta(:) => null()
    real(kind=8), pointer :: v_f_pro(:) => null()
    real(kind=8), pointer :: vale(:) => null()
    blas_int :: b_lda, b_ldb, b_ldvl, b_ldvr, b_lwork, b_n
    blas_int :: b_incx, b_incy
!
!-- DEBUT --C
!
    cbid = dcmplx(0.d0, 0.d0)
    pi = r8pi()
    call jemarq()
    call infniv(ifm, niv)
!
!------------------------------------------------------------C
!--                                                        --C
!-- CONSTRUCTION DES MATRICES D'IMPEDANCE DYNAMIQUE K+MU*M --C
!--            ET DE MASSE DU MODELE D'INTERFACE           --C
!--                                                        --C
!------------------------------------------------------------C
!
!-- ON DESACTIVE LE TEST FPE
    call matfpe(-1)
!
    call mtdscr(ssami)
    call jeveuo(ssami(1:19)//'.&INT', 'L', lmatma)
!
    imped = '&&MOIN93.RAID_SHIFT'
    call mtdefs(imped, raiint, 'V', ' ')
    lismat(1) = raiint
    lismat(2) = ssami
    if (switch .eq. 1) then
        call getvr8('MODE_INTERF', 'SHIFT', iocc=1, scal=rbid, nbret=ibid)
        shift = -(rbid*2.d0*pi)**2
    end if
!
    comlin(1) = 1.d0
    comlin(2) = shift
    listyp(1) = 'R'
    listyp(2) = 'R'
!
    call mtcmbl(2, listyp, comlin, lismat, imped, &
                ' ', nume91, 'ELIM1')
    call mtdscr(imped)
    call jeveuo(imped(1:19)//'.&INT', 'E', limped)
    call dismoi('SOLVEUR', ssami, 'MATR_ASSE', repk=solveu)
    call dismoi('NOM_NUME_DDL', ssami, 'MATR_ASSE', repk=nume91)
    call preres(solveu, 'V', iret, '&&OP0091.MATPRE', imped, &
                ibid, 1)
!
    if (iret .eq. 2) then
        valk = imped
        call utmess('F', 'ALGELINE4_37', sk=valk)
    end if
!
!-------------------------------------------------------------------C
!--                                                               --C
!-- RECUPERATION DU NOMBRE DE STRUCTURES DISJOINTES A L'INTERFACE --C
!--      AINSI QUE LES CONNECTIVITES PAR DECOMPOSITION QR         --C
!--                                                               --C
!-------------------------------------------------------------------C
!
    call intdis(coint, nnoint, noddli, '&&MODINT.INTERFACES_SST ', nbsst)
!
    call jeveuo('&&MODINT.INTERFACES_SST ', 'L', lindin)
!
!------------------------------------------------------C
!--                                                  --C
!-- CALCUL DES MODES DU MODELE D'INTERFACE (ARNOLDI) --C
!--                                                  --C
!------------------------------------------------------C
!
!-- ESTIMATION DU NOMBRE DE MODES A CALCULER PAR SOUS STRUCURE
!
    norm = dble(nbmod/nbsst)
    temp = 6.d0/nbsst
!
    coeff = 3
    if (norm .gt. 7) then
        nbvect = coeff*(int(norm)+2*(int(temp)+1))
    else
        nbvect = coeff*(6+int(norm)+2)
    end if
    nsekry = int(nbvect*nbsst/coeff)
    nsekry2 = nsekry+5
    coeff = int(nbvect/coeff)
    write (ifm, *) '------------------------------------------------',&
     &'------------------------'
    write (ifm, *) ' VOUS AVEZ DEMANDE', nbmod, ' MODES'
    write (ifm, *) ' LA TAILE DU SOUS ESPACE RETENU EST', nsekry
    write (ifm, *) '------------------------------------------------',&
     &'------------------------'
!--
!-- Appel a nmop45 pour le calcul des modes du modele d'interface
!--
!
!
! --- CREATION DE LA SD EIGENSOLVER PARAMETRANT LE CALCUL MODAL
! --- UN GEP SYM REEL RESOLU VIA SORENSEN
!
! BANDE MODALE EN DUR
    bande(1) = 0.d0
    bande(2) = 1.d0
    modes = '&&MODEST'
    eigsol = '&&MODINT.EIGSOL'
!
    k1bid = 'R'
    k8bid = ''
    k16bid = ''
    k19bid = ''
    k24bid = ''
    ibid = isnnem()
    r8bid = r8vide()
    typres = 'DYNAMIQUE'
    method = 'SORENSEN'
! MATR_A
    raide2 = imped
! MATR_B
    masse2 = ssami
! OPTION MODALE EN DUR
    optiof = 'PLUS_PETITE'
    nbborn = 1
! DIM_SOUS_ESPACE EN DUR
    nbvect = 0
! COEF_SOUS_ESPACE EN DUR
    nbvec2 = 2
! NMAX_ITER_SHIFT EN DUR
    nbrss = 5
! PARA_ORTHO_SOREN EN DUR
    alpha = 0.717d0
! NMAX_ITER_SOREN EN DUR
    maxitr = 200
! PREC_SOREN EN DUR
    tolsor = 0.d0
! CALC_FREQ/FLAMB/PREC_SHIFT EN DUR
    precsh = 5.d-2
! SEUIL_FREQ/CRIT EN DUR
    fcorig = 1.d-2
    omecor = omega2(fcorig)
! VERI_MODE/PREC_SHIFT EN DUR
    precdc = 5.d-2
! STOP_BANDE_VIDE EN DUR
    arret = 'NON'
! STURM EN DUR
    sturm = 'NON'
! OPTION MODE RIGIDE EN DUR
    modrig = 'SANS'
! OPTION STOP_ERREUR EN DUR
    stoper = 'NON'
! TYPE DE CALCUL: 'CALIBRATION' OU 'TOUT'.
    typcal = 'TOUT'
    call vpcres(eigsol, typres, raide2, masse2, k19bid, &
                optiof, method, modrig, arret, k19bid, &
                stoper, sturm, typcal, k1bid, k16bid, &
                nsekry2, nbvect, nbvec2, nbrss, nbborn, &
                ibid, ibid, ibid, ibid, maxitr, &
                bande, precsh, omecor, precdc, r8bid, &
                r8bid, r8bid, r8bid, r8bid, tolsor, &
                alpha)
!
! --- CALCUL MODAL PROPREMENT DIT
!
    l_hpp = ASTER_TRUE
    mod45 = 'VIBR'
    sdstab = '&&DUMMY2'
    call nmop45(eigsol, l_hpp, mod45, modes, sdstab)
    call vpleci(eigsol, 'I', 1, k24bid, r8bid, &
                nsekry2)
    call detrsd('EIGENSOLVER', eigsol)
!
!   -- on examine les modes calcules pour savoir ou tronquer sans couper
!      un sous-espace propre en 2 :
    nsekry3 = 0
    do j1 = nsekry, nsekry2-1
        call rsadpa(modes, 'L', 1, 'FREQ', j1, &
                    0, sjv=jfreq)
        freq1 = zr(jfreq)
        call rsadpa(modes, 'L', 1, 'FREQ', j1+1, &
                    0, sjv=jfreq)
        freq2 = zr(jfreq)
        if (abs(freq1-freq2)/(abs(freq1)+abs(freq2)) .gt. 1.e-6) then
! on peut "couper" a j1 :
            nsekry3 = j1
            exit
        end if
    end do
    ASSERT(nsekry3 .ne. 0)
    if (nsekry3 .gt. nsekry) nsekry = nsekry3
!
!
    call wkvect('&&MODINT.SE_KRYLOV', 'V V R', neq*nsekry, lmakry)
    call jeveuo('&&MOIN93.V_IND_LAG', 'L', vi=v_ind_lag)
    call jeveuo('&&MOIN93.DDL_ACTIF_INT', 'L', vi=ddl_actif_int)
!restaurant trélazé
    do j1 = 1, nsekry
! construction du nom du mode "a la main"
        call codent(j1-1, 'D0', k6bid)
        call jeveuo(modes//'.001.'//k6bid//'.VALE', 'L', vr=vale)
        b_n = to_blas_int(6*nnoint)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        norm = ddot(b_n, vale, b_incx, vale, b_incy)
        norm = sqrt(norm)
! normalisation dans L2
        do i1 = 1, nddlin
            m1 = ddl_actif_int(i1)
            zr(lmakry+(j1-1)*neq+v_ind_lag(1+(i1-1)*2)-1) = vale(m1)/norm
            zr(lmakry+(j1-1)*neq+v_ind_lag(1+(i1-1)*2+1)-1) = vale(m1)/norm
        end do
    end do
!
    call jedetc('G', modes, 1)
    call jedetc('G', '&&DUMMY2', 1)
!
    no = max(nsekry, nbvect)
    AS_ALLOCATE(vr=v_f_pro, size=no)
    AS_ALLOCATE(vi=v_ind_f_pro, size=no)
!
!-- RELEVE STATIQUE DU SOUS ESPACE DE KRYLOV SUR LE MODELE COMPLET
!
    call dismoi('SOLVEUR', raide, 'MATR_ASSE', repk=solveu)
    call resoud(raide, '&&MOIN93.MATPRE', solveu, ' ', nsekry, &
                ' ', ' ', ' ', zr(lmakry), [cbid], &
                ' ', .true._1, 0, iret)
!
!---------------------------------------------C
!--                                         --C
!-- PROJECTION DE LA MASSE ET DE LA RAIDEUR --C
!--                                         --C
!---------------------------------------------C
!
    call wkvect('&&MODINT.M_PROJ_TEMP', 'V V R', neq*nsekry, lmatrm)
    call wkvect('&&MODINT.K_PROJ_TEMP', 'V V R', neq*nsekry, lmatrk)
    call wkvect('&&MODINT.M_PROJ', 'V V R', nsekry**2, lmapro)
    call wkvect('&&MODINT.K_PROJ', 'V V R', nsekry**2, lkpro)
!
!-- MISE A 0 DES DDL DE LAGRANGE
    call dismoi('NOM_NUME_DDL', raide, 'MATR_ASSE', repk=nume)
    call dismoi('NUME_EQUA', nume, 'NUME_DDL', repk=prno)
    call jeveuo(prno//'.DELG', 'L', vi=delg)
!
    do i1 = 1, neq
        if (delg(i1) .lt. 0) then
            do j1 = 1, nsekry
                zr(lmakry+(j1-1)*neq+i1-1) = 0.d0
            end do
        end if
    end do
!
    call jeveuo(masse(1:19)//'.&INT', 'L', lmatm)
    call mrmult('ZERO', lmatm, zr(lmakry), zr(lmatrm), nsekry, &
                .true._1)
    call jeveuo(raide(1:19)//'.&INT', 'L', lmatk)
    call mrmult('ZERO', lmatk, zr(lmakry), zr(lmatrk), nsekry, &
                .true._1)
!
    do j1 = 1, nsekry
        b_n = to_blas_int(neq)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        zr(lmapro+(j1-1)*nsekry+j1-1) = ddot( &
                                        b_n, zr(lmakry+(j1-1)*neq), b_incx, &
                                        zr(lmatrm+(j1-1)*neq), b_incy &
                                        )
        b_n = to_blas_int(neq)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        zr(lkpro+(j1-1)*nsekry+j1-1) = ddot( &
                                       b_n, zr(lmakry+(j1-1)*neq), b_incx, zr(lmatrk+(j1-1)*neq), &
                                       b_incy &
                                       )
        do i1 = 1, j1-1
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            zr(lmapro+(j1-1)*nsekry+i1-1) = ddot( &
                                            b_n, zr(lmakry+(i1-1)*neq), b_incx, &
                                            zr(lmatrm+(j1-1)*neq), b_incy &
                                            )
            zr(lmapro+(i1-1)*nsekry+j1-1) = zr(lmapro+(j1-1)*nsekry+i1-1)
!
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            zr(lkpro+(j1-1)*nsekry+i1-1) = ddot( &
                                           b_n, zr(lmakry+(i1-1)*neq), b_incx, &
                                           zr(lmatrk+(j1-1)*neq), b_incy &
                                           )
            zr(lkpro+(i1-1)*nsekry+j1-1) = zr(lkpro+(j1-1)*nsekry+i1-1)
        end do
    end do
!
!-------------------------------------------------C
!--                                             --C
!-- RESOLUTION DU PB AUX VALEURS PROPRES REDUIT --C
!--                                             --C
!-------------------------------------------------C
!
    AS_ALLOCATE(vr=vect_alphar, size=nsekry)
    AS_ALLOCATE(vr=vect_alphai, size=nsekry)
    AS_ALLOCATE(vr=vect_beta, size=nsekry)
    AS_ALLOCATE(vr=matr_mod_red, size=nsekry**2)
!
    b_ldvr = to_blas_int(nsekry)
    b_ldvl = to_blas_int(1)
    b_ldb = to_blas_int(nsekry)
    b_lda = to_blas_int(nsekry)
    b_n = to_blas_int(nsekry)
    b_lwork = to_blas_int(-1)
    call dggev('N', 'V', b_n, zr(lkpro), b_lda, &
               zr(lmapro), b_ldb, vect_alphar, vect_alphai, vect_beta, &
               vrbid, b_ldvl, matr_mod_red, b_ldvr, swork, &
               b_lwork, info)
    lwork = int(swork(1))
    AS_ALLOCATE(vr=matr_work_dggev, size=lwork)
    b_ldvr = to_blas_int(nsekry)
    b_ldvl = to_blas_int(1)
    b_ldb = to_blas_int(nsekry)
    b_lda = to_blas_int(nsekry)
    b_n = to_blas_int(nsekry)
    b_lwork = to_blas_int(lwork)
    call dggev('N', 'V', b_n, zr(lkpro), b_lda, &
               zr(lmapro), b_ldb, vect_alphar, vect_alphai, vect_beta, &
               vrbid, b_ldvl, matr_mod_red, b_ldvr, matr_work_dggev, &
               b_lwork, info)
!-- ON REACTIVE LE TEST FPE
    call matfpe(1)
!
!-- CLASSEMENT DES FREQUENCES PROPRES
    temp = 1.d+16
    do i1 = 1, nsekry
        if (abs(vect_beta(i1)) .gt. 0) then
            lambda = vect_alphar(i1)/vect_beta(i1)
            v_f_pro(i1) = (sqrt(abs(lambda)))/2/pi
        else
            v_f_pro(i1) = temp
        end if
    end do
!
!
    call jeexin(vefreq, iret)
    if (iret .eq. 0) then
        call wkvect(vefreq, 'V V R', nbmod, lvp)
    else
! appel depuis modexp => il faut redimensionner l'objet
        call jedetr(vefreq)
        call wkvect(vefreq, 'V V R', nbmod, lvp)
    end if
!
    do i1 = 1, nbmod
        temp = 1.d+16
        do j1 = 1, nsekry
            if (v_f_pro(j1) .lt. temp) then
                temp = v_f_pro(j1)
                v_ind_f_pro(i1) = j1
            end if
        end do
        v_f_pro(1+v_ind_f_pro(i1)-1) = 1.d+16
        lambda = vect_alphar(1+v_ind_f_pro(i1)-1)/vect_beta(1+v_ind_f_pro(i1)-1)-0*shift
        zr(lvp+i1-1) = (sqrt(abs(lambda)))/2/pi
    end do
!
!-------------------------------------------------------C
!--                                                   --C
!-- RESTITUTION DES MODES D'INTERFACE SUR LE MAILLAGE --C
!--                                                   --C
!-------------------------------------------------------C
!
    call wkvect(matmod, 'V V R', neq*nbmod, lmatmo)
    do j1 = 1, nbmod
        do k1 = 1, nsekry
            temp = matr_mod_red(1+(v_ind_f_pro(j1)-1)*nsekry+k1-1)
            do i1 = 1, neq
                zr(lmatmo+(j1-1)*neq+i1-1) = zr( &
                                             lmatmo+(j1-1)*neq+i1-1)+temp*zr(lmakry+(k1-1)*neq+i1&
                                             &-1 &
                                             )
            end do
        end do
    end do
!
!-- normalisation des modes et mise a zero des ddls de Lagrange :
    do j1 = 1, nbmod
        normx = -1.d0
        do i1 = 1, neq
            if (delg(i1) .lt. 0) then
                zr(lmatmo+(j1-1)*neq+i1-1) = 0.d0
            else
                valx = zr(lmatmo+(j1-1)*neq+i1-1)
                avalx = abs(valx)
                if (avalx .gt. normx) then
                    if (avalx/normx .le. 1.0000001) then
                        mode_symetrique = 1
                    else
                        mode_symetrique = 0
                    end if
                    normx = avalx
                    mval = valx
                end if
            end if
        end do
        ASSERT(normx .gt. 0.d0)
!
!       -- Si le mode est symetrique, on choisit celui qui a son premier "max" >0:
        if (mode_symetrique .eq. 1) then
            do i1 = 1, neq
                valx = zr(lmatmo+(j1-1)*neq+i1-1)
                avalx = abs(valx)
                if (avalx .gt. 0.9999999d0*normx) then
                    if (mval*valx .lt. 0.d0) mval = -mval
                    exit
                end if
            end do
        end if
    end do
!
!---------------------------------------C
!--                                   --C
!-- DESTRUCTION DES OBJETS DE TRAVAIL --C
!--                                   --C
!---------------------------------------C
!
    call jedetr('&&MODINT.M_PROJ_TEMP')
    call jedetr('&&MODINT.K_PROJ_TEMP')
    call jedetr('&&MODINT.M_PROJ')
    call jedetr('&&MODINT.K_PROJ')
    call jedetr('&&MODINT.INTERFACES_SST')
!
    AS_DEALLOCATE(vr=vect_alphar)
    AS_DEALLOCATE(vr=vect_alphai)
    AS_DEALLOCATE(vr=vect_beta)
    AS_DEALLOCATE(vr=matr_mod_red)
    AS_DEALLOCATE(vr=matr_work_dggev)
!
!
    call detrsd('MATR_ASSE', imped)
!
    call jedetr('&&MODINT.SE_KRYLOV')
!
    AS_DEALLOCATE(vr=v_f_pro)
    AS_DEALLOCATE(vi=v_ind_f_pro)
!
!-- FIN --C
!
    call jedema()
end subroutine
