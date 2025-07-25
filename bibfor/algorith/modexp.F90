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
subroutine modexp(modgen, sst1, indin1, lino1, nbmod, &
                  numlia, tramod, modet, solveu)
!-------------------------------------------------------------C
!--                                                         --C
!--  ROUTINE QUI REALISE L'EXPANSION DES MODES D'INTERFACE  --C
!--    DE L'INTERFACE ESCLAVE VERS L'INTERFACE MAITRESSE    --C
!--                                                         --C
!-------------------------------------------------------------C
!--   VARIABLES E/S  :
!--   MODGEN   /IN/  : NOM DU MODELE GENERALISE
!--   SST1     /IN/  : NOM DE LA SOUS STRUCTURE
!--   INDIN1   /IN/  : NOM DU VECTEUR CONTENANT LES DDL D'INTERFACE
!--   LINO1    /IN/  : NOM DU VECTEUR CONTENANT LES NOEUDS D'INTERFACE
!--   NBMOD    /IN/  : NOMBRE DE MODE A ETENDRE
!--   NUMLIA   /IN/  : NUMERO DE LA LIAISON ASSOCIEE A L'INTERFACE
!--   TRAMOD   /IN/  : TRACE DES MODES SUR L'INTERFACE ESCLAVE
!--   MODET    /OUT/ : MODES ETENDUS SUR L'INTERFACE MAITRE
!
    implicit none
!
!
!
#include "jeveux.h"
#include "asterc/matfpe.h"
#include "asterc/r8pi.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/conint.h"
#include "asterfort/ddllag.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedetc.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/infniv.h"
#include "asterfort/modint.h"
#include "asterfort/preres.h"
#include "asterfort/resoud.h"
#include "asterfort/tri.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "blas/dgelss.h"
!
!
!
!
    integer(kind=8) :: nbmod, lbid, i1, j1, k1, isst1, ibid, nbeq1, nl, nc
    integer(kind=8) :: lmast, numlia, nbno, lnres, lmodet, sizeco, connec, lconnc, nbec
    integer(kind=8) :: lprno, ipos1, lcphi, nbddl, lnoint, lindin, llino, lindno, lipos
    integer(kind=8) :: ik, lddld, linlag, lintrf, linddl, nddlin, nbvect, ltramo, lmatmo
    integer(kind=8) :: lclin, lwork, jwork, lphiex, lcpet, ifm, niv
    integer(kind=4) :: info, rank
    real(kind=8) :: shift, swork(1), pi
    complex(kind=8) :: cbid
    character(len=4) :: k4bid
    character(len=8) :: modgen, sst1
    character(len=14) :: nume, nume91
    character(len=19) :: raide, masse, solveu, prno, ssami, raiint
    character(len=24) :: coint, noddli, matmod, vefreq, indin1, lino1, tramod
    character(len=24) :: modet
    integer(kind=8) :: iret
    integer(kind=8), pointer :: vect_clefs(:) => null()
    integer(kind=8), pointer :: vect_num(:) => null()
    blas_int :: b_lda, b_ldb, b_lwork, b_m, b_n, b_nrhs
    cbid = dcmplx(0.d0, 0.d0)
    pi = r8pi()
    call infniv(ifm, niv)
!
!---------------------------------------------------C
!--                                               --C
!-- CALCUL DES MODES D'INTERFACE POUR L'EXPANSION --C
!--                                               --C
!---------------------------------------------------C
!
!-- RECUPERATION DES MATRICES DE MASSE ET RAIDEUR
    call jenonu(jexnom(modgen//'      .MODG.SSNO', sst1), isst1)
    call jeveuo(jexnum(modgen//'      .MODG.SSME', isst1), 'L', ibid)
    call jeveuo(zk8(ibid)//'.MAEL_RAID_REFE', 'L', lbid)
    raide = zk24(lbid+1) (1:19)
    call jeveuo(zk8(ibid)//'.MAEL_MASS_REFE', 'L', lbid)
    masse = zk24(lbid+1) (1:19)
!
!-- RECUPERATION DU NUME_DDL
    call dismoi('NOM_NUME_DDL', masse(1:8), 'MATR_ASSE', repk=nume)
    call dismoi('NUME_EQUA', nume, 'NUME_DDL', repk=prno)
    call jeveuo(jexnum(prno//'.PRNO', 1), 'L', lprno)
!
    call dismoi('NB_EQUA', raide, 'MATR_ASSE', repi=nbeq1)
    call dismoi('NB_EC', 'DEPL_R', 'GRANDEUR', repi=nbec)
!
!-- REMPLISSAGE DES VECTEURS D'INDICES ASSOCIES AUX DDL D'INTERFACE
!--      POUR LA CREATION DES MATRICES
    call jelira(indin1, 'LONMAX', nbddl)
    call jeveuo(indin1, 'L', lindin)
    noddli = '&&MOIN93.NOEUDS_DDL_INT'
    call jelira(lino1, 'LONMAX', nbno)
    call jeveuo(lino1, 'L', llino)
!
    call wkvect(noddli, 'V V I', 9*nbno, lnoint)
    call wkvect('&&MOIN93.V_IND_LAG', 'V V I', 2*nbddl, linlag)
    call wkvect('&&MOIN93.DDL_ACTIF_INT', 'V V I', nbddl, lintrf)
    call wkvect('&&MOIN93.V_IND_DDL_INT', 'V V I', nbddl, linddl)
!
!-- LA LISTE DES NOEUDS D'INTERFACE DOIT ETRE ORDONNEE
    AS_ALLOCATE(vi=vect_clefs, size=nbno)
    AS_ALLOCATE(vi=vect_num, size=nbno)
!
    do k1 = 1, nbno
        vect_clefs(k1) = k1
        vect_num(k1) = zi(llino+k1-1)
    end do
    call tri(vect_num, vect_clefs, 1, nbno)
    do k1 = 1, nbno
        ik = vect_clefs(k1)
        zi(lnoint+(k1-1)) = zi(llino+ik-1)
        zi(lnoint+(k1-1)+nbno) = zi(lindin+(ik-1)*6)
        zi(lnoint+(k1-1)+2*nbno) = 6
        do j1 = 1, 6
            zi(lnoint+(k1-1)+((2+j1)*nbno)) = j1
        end do
    end do
    AS_DEALLOCATE(vi=vect_clefs)
    AS_DEALLOCATE(vi=vect_num)
!
    call wkvect('&&MOIN93.IS_DDL_INTERF', 'V V I', nbeq1, lddld)
    k1 = 1
    do i1 = 1, nbddl
        if (zi(lindin+i1-1) .gt. 0) zi(lddld+zi(lindin+i1-1)-1) = 1
        ipos1 = zi(lindin+i1-1)
        if (ipos1 .gt. 0) then
            zi(linddl+k1-1) = ipos1
            call ddllag(nume, ipos1, nbeq1, zi(linlag+(k1-1)*2), zi(linlag+(k1-1)*2+1))
            zi(lintrf+k1-1) = i1
            k1 = k1+1
        end if
    end do
    nddlin = k1-1
!
!-- CONSTRUCTION DES MATRICES DE MASSE ET DE RAIDEUR DU PROBLEME
!--    D'INTERFACE
    sizeco = 36*nbno
    coint = '&&MODEXP.CONNEC_INTERF'
    nume91 = '&&NUME91'
    raiint = '&&RAID91'
    ssami = '&&MASS91'
    call wkvect(coint, 'V V I', sizeco, lconnc)
    call wkvect('&&MOIN93.IND_NOEUD', 'V V I', zi(lnoint+nbno-1), lindno)
    call wkvect('&&MOIN93.IPOS_DDL_INTERF', 'V V I', nbno, lipos)
    call conint(nume, raide, coint, connec, noddli, &
                nbno, nume91, raiint, ssami)
!
!-- CALCUL DES MODES DU MODELE D'INTERFACE
    call getvr8(' ', 'SHIFT', scal=shift, nbret=ibid)
    shift = -((shift*2.d0*pi)**2)
    matmod = '&&MODEXP.MATRICE_MODES'
    vefreq = '&&MODEXP.VECTEUR_FREQ'
!
    call codent(numlia, 'D0', k4bid)
!-- MOUVEMENTS DE L'INTERFACE ESCLAVE A ETENDRE
    call jeveuo('&&OP0091.MAS'//k4bid, 'L', lmast)
    call jelira('&&OP0091.MAS'//k4bid, 'LONMAX', ibid)
    nl = int(ibid/nbmod)
!
!-- MATRICE D'OBSERVATION
    call jelira(tramod, 'LONMAX', ibid)
    nc = int(ibid/nl)
    call jeveuo(tramod, 'L', ltramo)
!
!   Par defaut, on double la taille du sous espace de recherche
    nbvect = 2*nbmod
!
!-- Factorisation de la matrice de raideur
    call preres(solveu, 'V', ibid, '&&OP0091.MATPRE', raide, &
                ibid, 1)
!
!-- ON BOUCLE POUR AVOIR UNE EXPANSION CORRECTE. TANT QUE C'EST PAS BON,
!-- ON ENRICHIT LA BASE DES MODES D'INTERFACE
500 continue
!  IL FAUDRA TRAVAILLER CE POINT POUR RENDRE L'EXPANSION PLUS ROBUSTE
    call modint(ssami, raiint, nddlin, nbvect, shift, &
                matmod, masse, raide, nbeq1, coint, &
                noddli, nbno, vefreq, 0)
    call jeveuo(matmod, 'L', lmatmo)
!
!
!---------------------------C
!--                       --C
!-- EXPANSION DES DONNEES --C
!--                       --C
!---------------------------C
!
!-- TAILLES DES MATRICES
!--
!--  C      : NL x NC
!--  MATMOD : NC x NBVECT
!--  CPHI   : NL x NBVECT
!
!
!-- PRODUIT C*PHI_mast
!
    call matfpe(-1)
    call wkvect('&&MODEXP.CPHI', 'V V R', nl*nbvect, lcphi)
!
    do j1 = 1, nbvect
        do k1 = 1, nc
            do i1 = 1, nl
                zr(lcphi+(j1-1)*nl+i1-1) = zr( &
                                           lcphi+(j1-1)*nl+i1-1)+zr(ltramo+(k1-1)*nl+(i1-1))*zr(&
                                           &lmatmo+(j1-1)*nc+(k1-1) &
                                           )
            end do
        end do
    end do
!      CALL DGEMM('N','N',NL,NBVECT,NC,1.,ZR(LTRAMO),
!     &            NL,ZR(LMATMO),NC,0.,ZR(LCPHI),NL)
!
!-- EXPANSION DE DONNEES
    ibid = max(nl, nbvect)
    call wkvect('&&MODEXP.COMB_LIN', 'V V R', ibid*nbmod, lclin)
!
!-- RECOPIE A LA MAIN POUR ETRE COMPATIBLE AVEC LES TAILLES
    do j1 = 1, nbmod
        do i1 = 1, nl
            zr(lclin+(j1-1)*max(nl, nbvect)+(i1-1)) = zr(lmast+(j1-1)*nl+(i1-1))
        end do
    end do
!
    call wkvect('&&MODEXP.VEC_VAL_SING', 'V V R', min(nl, nbvect), ibid)
!
    b_ldb = to_blas_int(max(nl, nbvect))
    b_lda = to_blas_int(nl)
    b_m = to_blas_int(nl)
    b_n = to_blas_int(nbvect)
    b_nrhs = to_blas_int(nbmod)
    b_lwork = to_blas_int(-1)
    call dgelss(b_m, b_n, b_nrhs, zr(lcphi), b_lda, &
                zr(lclin), b_ldb, zr(ibid), -1.0d0, rank, &
                swork, b_lwork, info)
!
    lwork = int(swork(1))
    call wkvect('&&MODEXP.MAT_AXB', 'V V R', lwork, jwork)
!
    b_ldb = to_blas_int(max(nl, nbvect))
    b_lda = to_blas_int(nl)
    b_m = to_blas_int(nl)
    b_n = to_blas_int(nbvect)
    b_nrhs = to_blas_int(nbmod)
    b_lwork = to_blas_int(lwork)
    call dgelss(b_m, b_n, b_nrhs, zr(lcphi), b_lda, &
                zr(lclin), b_ldb, zr(ibid), 1.0d-12, rank, &
                zr(jwork), b_lwork, info)
!
!-- MOUVEMENTS DE L'INTERFACE
    call wkvect('&&MODEXP.PHI_EXP', 'V V R', nc*nbmod, lphiex)
!
    do j1 = 1, nbmod
        do k1 = 1, nbvect
            do i1 = 1, nc
                zr(lphiex+(j1-1)*nc+i1-1) = zr( &
                                            lphiex+(j1-1)*nc+i1-1)+zr(lmatmo+(k1-1)*nc+(i1-1))*z&
                                            &r(lclin+(j1-1)*nbvect+(k1-1) &
                                            )
            end do
        end do
    end do
!
!-- PROJECTION POUR VERIFIER L'EXPANSION
    call wkvect('&&MODEXP.C_PHI_EXP', 'V V R', nl*nbmod, lcpet)
!
    do j1 = 1, nbmod
        do k1 = 1, nc
            do i1 = 1, nl
                zr(lcpet+(j1-1)*nl+i1-1) = zr( &
                                           lcpet+(j1-1)*nl+i1-1)+zr(ltramo+(k1-1)*nl+(i1-1))*zr(&
                                           &lphiex+(j1-1)*nc+(k1-1) &
                                           )
            end do
        end do
    end do
!
!-- VERIFICATION DE LA QUALITE DE L'EXPANSION
!
    call wkvect('&&MODEXP.NORM_RESIDU', 'V V R', nbmod, lnres)
!
!-- NORME DE LA DIFFERENCE
    swork(1) = 0
    do j1 = 1, nbmod
        do i1 = 1, nl
            zr(lnres+j1-1) = zr(lnres+j1-1)+(zr(lmast+(j1-1)*nl+(i1-1))-zr(lcpet+(j1-1)*nl+(i1-1&
                             &)))**2
        end do
        zr(lnres+j1-1) = sqrt(zr(lnres+j1-1))/nl
        swork(1) = max(swork(1), zr(lnres+j1-1))
    end do
!
    call matfpe(1)
!
    if (swork(1) .gt. 1.d-1) then
!
        if (nbvect .gt. nddlin/2) then
!-- ON N'ARRIVE PAS A AVOIR UNE EXPANSION CORRECTE
            write (ifm, *) '*----------------------------*'
            write (ifm, *) '*'
            write (ifm, *) 'Residu de l''expansion :', swork(1)
            write (ifm, *) '  Taille max. du sous espace atteinte'
            write (ifm, *) '  Probleme mal conditionne'
            write (ifm, *) '*'
            write (ifm, *) '*----------------------------*'
            ASSERT(.false.)
        end if
!
        write (ifm, *) 'Residu de l''expansion :', swork(1)
        nbvect = nbvect+nbmod
!
        call jedetr('&&OP0091.MAT_SM1XUT')
        call jedetr('&&OP0091.MAT_VXSM1XUT')
        call jedetr('&&OP0091.MAT_S')
        call jedetr('&&OP0091.MAT_U')
        call jedetr('&&OP0091.MAT_V')
        call jedetr('&&MODINT.INTERFACES_SST')
        call jedetr('&&MODEXP.CPHI')
        call jedetr('&&MODEXP.COMB_LIN')
        call jedetr('&&MODEXP.VEC_VAL_SING')
        call jedetr('&&MODEXP.MAT_AXB')
        call jedetr('&&MODEXP.PHI_EXP')
        call jedetr('&&MODEXP.MATRICE_MODES')
        call jedetr('&&MODEXP.C_PHI_EXP')
        call jedetr('&&MODEXP.NORM_RESIDU')
        goto 500
    else
        write (ifm, *) 'Residu de l''expansion :', swork(1)
    end if
!
!-------------------------C
!--                     --C
!-- RELEVEMENT STATIQUE --C
!--                     --C
!-------------------------C
!
    call wkvect(modet, 'V V R', nbeq1*nbmod, lmodet)
!
    if (nc .ne. nddlin) then
        ASSERT(.false.)
    end if
!
!-- EXPANSION STATIQUE
!
    do j1 = 1, nbmod
        do i1 = 1, nddlin
            zr(lmodet+(j1-1)*nbeq1+zi(linlag+(i1-1)*2)-1) = zr(lphiex+(j1-1)*nc+(i1-1))
            zr(lmodet+(j1-1)*nbeq1+zi(linlag+(i1-1)*2+1)-1) = zr(lphiex+(j1-1)*nc+(i1-1))
        end do
    end do
!
    call resoud(raide, '&&MOIN93.MATPRE', solveu, ' ', nbmod, &
                ' ', ' ', ' ', zr(lmodet), [cbid], &
                ' ', .true._1, 0, iret)
!
!------------C
!-- MENAGE --C
!------------C
    call jedetr('&&MODEXP.CPHI')
    call jedetr('&&MODEXP.COMB_LIN')
    call jedetr('&&MODEXP.PHI_EXP')
    call jedetr('&&MODEXP.VEC_VAL_SING')
    call jedetr('&&MODEXP.MAT_AXB')
    call jedetr('&&MODEXP.CONNEC_INTERF')
    call jedetr('&&MOIN93.NOEUDS_DDL_INT')
    call jedetr('&&MOIN93.V_IND_LAG')
    call jedetr('&&MOIN93.DDL_ACTIF_INT')
    call jedetr('&&MOIN93.V_IND_DDL_INT')
    call jedetr('&&MOIN93.IS_DDL_INTERF')
    call jedetr('&&MOIN93.IND_NOEUD')
    call jedetr('&&MOIN93.IPOS_DDL_INTERF')
    call jedetr('&&MODL91      .MODG.SSNO')
    call jedetr('&&MODL91      .MODG.SSME')
    call jedetr('&&MODINT.INTERFACES_SST')
    call jedetr('&&MODEXP.C_PHI_EXP')
    call jedetr('&&MODEXP.NORM_RESIDU')
!
    call detrsd('MATR_ASSE', '&&RAID91')
    call detrsd('MATR_ASSE', '&&MASS91')
    call jedetc('V', '&&NUME91', 1)
    call jedetr('&&MODEXP.MATRICE_MODES')
    call jedetr('&&MODEXP.VECTEUR_FREQ')
!
end subroutine
