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
subroutine traint(resgen, modgen, numlia, sst1, sst2, &
                  intf1, intf2, nbmod, nl, nc)
    implicit none
!
!-------------------------------------------------------------C
!--             ROUTINE XXXXX1    M. CORUS - AOUT 2011      --C
!--             CALCUL DES TRAVAUX AUX INTERFACES           --C
!--                                                         --C
!-------------------------------------------------------------C
!--   VARIABLES E/S  :
!--   RESGEN   /IN/  : NOM DE LA SD RESULTAT CONTENANT LES MODES REDUITS
!--   MODGEN   /IN/  : NOM DU MODELE GENERALISE
!--   NUMLIA   /IN/  : NUMERO DE LA LIAISON ASSOCIEE A L'INTERFACE
!--   SST1     /IN/  : NOM DE LA SOUS STRUCTURE 1
!--   SST2     /IN/  : NOM DE LA SOUS STRUCTURE 2
!--   INTF1    /IN/  : NOM DE L'INTERFACE DE LA SST 1
!--   INTF1    /IN/  : NOM DE L'INTERFACE DE LA SST 2
!--   NBMOD    /IN/  : NOMBRE DE MODE DU MODELE REDUIT
!--   NL       /IN/  : NB DE LIGNES DE LA MATRICE D'OBSERVATION
!--   NC       /IN/  : NB DE COLONNES DE LA MATRICE D'OBSERVATION
!
!
!
#include "jeveux.h"
#include "asterfort/codent.h"
#include "asterfort/ddllag.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvis.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/lceqvn.h"
#include "asterfort/mrmult.h"
#include "asterfort/mtdscr.h"
#include "asterfort/regeec.h"
#include "asterfort/wkvect.h"
#include "asterfort/zerlag.h"
#include "blas/ddot.h"
!
!
!
    character(len=4) :: k4bid
    character(len=8) :: modgen, resgen, sst1, sst2, intf1, intf2, rest1
    character(len=8) :: rest2, mraid
    character(len=19) :: nume91
    character(len=24) :: indin1, indin2
    integer(kind=8) :: i1, ibid, nbsst, j1, k1, l1, imast1, nbeq, lmast, imast2
    integer(kind=8) :: nbeq1, nbeq2, numlia, lefi1, lefi2, nbmod, lnusst, isst1
    integer(kind=8) :: isst2, llint1, lobs2, lobs1, llint2, nbddl1, nbddl2, tach1
    integer(kind=8) :: nc, tach2, lomeg, lmod1, lmod2, lbid, lag1, lag2
    integer(kind=8) :: lesc, leff1, leff2, ideeq, nbddl, llint
    integer(kind=8) :: lcopy1, lcopy2, unit, ldepma, ldepsl, nl
    real(kind=8) :: travm, travk, trvint
    integer(kind=8), pointer :: num_sst_maitre(:) => null()
    real(kind=8), pointer :: trav_interf(:) => null()
    integer(kind=8), pointer :: num_sst_esclave(:) => null()
    character(len=8), pointer :: nom_sst(:) => null()
    integer(kind=8), pointer :: lipr(:) => null()
    integer(kind=8), pointer :: matrice_mass(:) => null()
    real(kind=8), pointer :: tr_mod_mast_pro(:) => null()
    integer(kind=8), pointer :: matrice_raid(:) => null()
    blas_int :: b_incx, b_incy, b_n
!
    call getvis(' ', 'UNITE', scal=unit, nbret=ibid)
    i1 = numlia
!
    call jelira(modgen//'      .MODG.SSNO', 'NOMMAX', nbsst)
!
    call jeveuo('&&OP0091.NOM_SST', 'E', vk8=nom_sst)
    call jeveuo('&&OP0091.NUME_SST', 'E', lnusst)
    call jeveuo('&&OP0091.MATRICE_MASS', 'E', vi=matrice_mass)
    call jeveuo('&&OP0091.MATRICE_RAID', 'E', vi=matrice_raid)
    call jeveuo('&&OP0091.TRAV_INTERF', 'E', vr=trav_interf)
    call jeveuo('&&OP0091.NUM_SST_ESCLAVE', 'E', vi=num_sst_esclave)
    call jeveuo('&&OP0091.NUM_SST_MAITRE', 'E', vi=num_sst_maitre)
    call jeveuo('&&OP0091.PULSA_PROPRES', 'L', lomeg)
!
    call jeveuo(modgen//'      .MODG.LIPR', 'L', vi=lipr)
!
!-- RECHERCHE DE LA SOUS STRUCTURE ESCLAVE DE LA LIAISON
!-- EN CHERCHANT LE FACTEUR MULTIPLICATEUR DANS LA MATRICE DE
!-- LIAISON LAGRANGE / LAGRANGE
    call jeveuo(jexnum(modgen//'      .MODG.LIMA', lipr(1+(i1-1)*9+8)), 'L', ibid)
!
!-- L'INTERFACE MAITRE EST CELLE DONT LE NUMERO EST NEGATIF
    if (int(zr(ibid)) .eq. -2) then
        imast1 = 1
        imast2 = -2
    else
        imast2 = 2
        imast1 = -1
    end if
!
!-- SOUS STRUCTURE 1
    indin1 = '&&VEC_DDL_INTF_'//intf1
    nbeq1 = lipr(1+(i1-1)*9+1)
!-- SOUS STRUCTURE 2
    indin2 = '&&VEC_DDL_INTF_'//intf2
    nbeq2 = lipr(1+(i1-1)*9+4)
!
!-- RESTITUTION DES MODES SUR LES DIFFERENTES SOUS STRUCTURES
!
    if (i1 .eq. 1) then
        nom_sst(1) = sst1
        nom_sst(2) = sst2
        zi(lnusst) = 2
        zi(lnusst+1) = 1
        zi(lnusst+2) = 2
        rest1 = '&&910001'
        rest2 = '&&910002'
        call regeec(rest1, resgen, sst1)
        call jedetr('&&MODE_ETENDU_REST_ELIM')
        call regeec(rest2, resgen, sst2)
        call jedetr('&&MODE_ETENDU_REST_ELIM')
    else
        ibid = zi(lnusst)
        isst1 = 0
        isst2 = 0
        do j1 = 1, ibid
            if (sst1 .eq. nom_sst(j1)) isst1 = j1
            if (sst2 .eq. nom_sst(j1)) isst2 = j1
        end do
        if (isst1 .eq. 0) then
            zi(lnusst) = zi(lnusst)+1
            zi(lnusst+zi(lnusst)) = zi(lnusst)
            nom_sst(1+zi(lnusst)-1) = sst1
            call codent(zi(lnusst), 'D0', k4bid)
            rest1 = '&&91'//k4bid
            call regeec(rest1, resgen, sst1)
            call jedetr('&&MODE_ETENDU_REST_ELIM')
        else
            call codent(zi(lnusst+isst1), 'D0', k4bid)
            rest1 = '&&91'//k4bid
        end if
        if (isst2 .eq. 0) then
            zi(lnusst) = zi(lnusst)+1
            zi(lnusst+zi(lnusst)) = zi(lnusst)
            nom_sst(1+zi(lnusst)-1) = sst2
            call codent(zi(lnusst), 'D0', k4bid)
            rest2 = '&&91'//k4bid
            call regeec(rest2, resgen, sst2)
            call jedetr('&&MODE_ETENDU_REST_ELIM')
        else
            call codent(zi(lnusst+isst2), 'D0', k4bid)
            rest2 = '&&91'//k4bid
        end if
    end if
!
!-- RECUPERATION DE LA SST ESCLAVE EN FONCTION DU NUMERO DE LA LIAISON
    do j1 = 1, zi(lnusst)
        if (imast1 .eq. 1) then
            if (sst1 .eq. nom_sst(j1)) then
                num_sst_esclave(i1) = zi(lnusst+j1)
            end if
            if (sst2 .eq. nom_sst(j1)) then
                num_sst_maitre(i1) = zi(lnusst+j1)
            end if
        end if
        if (imast2 .eq. 2) then
            if (sst2 .eq. nom_sst(j1)) then
                num_sst_esclave(i1) = zi(lnusst+j1)
            end if
            if (sst1 .eq. nom_sst(j1)) then
                num_sst_maitre(i1) = zi(lnusst+j1)
            end if
        end if
    end do
!
!-- OBSERVATION DE L'INTERFACE MAITRE
!
    call jeveuo('&&LIPSRB.TR_MOD_MAST_PRO', 'L', vr=tr_mod_mast_pro)
!
!-- RECUPERATION D'INFOS (TAILLES, DDL INTERF, MODES RESTITUES, ETC.)
    call jelira(indin1, 'LONMAX', nbddl1)
    call jelira(indin2, 'LONMAX', nbddl2)
    call jeveuo(indin1, 'L', llint1)
    call jeveuo(indin2, 'L', llint2)
    call jeveuo(jexnum(rest1//'           .TACH', 1), 'L', tach1)
    call jeveuo(jexnum(rest2//'           .TACH', 1), 'L', tach2)
    call jelira(zk24(tach1) (1:19)//'.VALE', 'LONMAX', nbeq1)
    call jelira(zk24(tach2) (1:19)//'.VALE', 'LONMAX', nbeq2)
!
!-- ALLOCATION DES OBJETS TEMPORAIRE POUR LE CALCUL DES TRAVAUX
    call codent(i1, 'D0', k4bid)
    call wkvect('&&OP0091.OBS_01', 'V V R', nbddl1, lobs1)
    call wkvect('&&OP0091.OBS_02', 'V V R', nbddl2, lobs2)
    call wkvect('&&OP0091.SLA'//k4bid, 'V V R', nl*nbmod, lesc)
    call wkvect('&&OP0091.MAS'//k4bid, 'V V R', nl*nbmod, lmast)
!
!
!  POUR LE STOCKAGE DES DEPLACEMENTS IMPOSES A L'INTERFACE
    if (imast1 .eq. 1) then
        call wkvect('&&OP0091.DEPL_IMPO_'//k4bid, 'V V R', nbmod*nbeq1, ldepsl)
        call wkvect('&&OP0091.MAST_IMPO_'//k4bid, 'V V R', nbmod*nbeq2, ldepma)
        llint = llint1
        nbddl = nbddl1
        nbeq = nbeq1
        call jenonu(jexnom(modgen//'      .MODG.SSNO', sst1), isst1)
!
    else
        call wkvect('&&OP0091.DEPL_IMPO_'//k4bid, 'V V R', nbmod*nbeq2, ldepsl)
        call wkvect('&&OP0091.MAST_IMPO_'//k4bid, 'V V R', nbmod*nbeq1, ldepma)
!
        llint = llint2
        nbddl = nbddl2
        nbeq = nbeq2
        call jenonu(jexnom(modgen//'      .MODG.SSNO', sst2), isst1)
!
    end if
!
    call wkvect('&&OP0091.LAGRANGES_1', 'V V I', nbddl, lag1)
    call wkvect('&&OP0091.LAGRANGES_2', 'V V I', nbddl, lag2)
!
    call jeveuo(jexnum(modgen//'      .MODG.SSME', isst1), 'L', ibid)
    call jeveuo(zk8(ibid)//'.MAEL_RAID_REFE', 'L', lbid)
    mraid = zk24(lbid+1) (1:8)
    call dismoi('NOM_NUME_DDL', mraid, 'MATR_ASSE', repk=nume91)
    do k1 = 1, nbddl
        if (zi(llint+k1-1) .gt. 0) then
            call ddllag(nume91, zi(llint+k1-1), nbeq, zi(lag1+k1-1), zi(lag2+k1-1))
        end if
    end do
!
!  STOCKAGES DES EFFORTS
    call wkvect('&&OP0091.MODE_SST1', 'V V R', nbeq1, lmod1)
    call wkvect('&&OP0091.MODE_SST2', 'V V R', nbeq2, lmod2)
    call wkvect('&&OP0091.MODE_SST1_EFF', 'V V R', nbeq1, leff1)
    call wkvect('&&OP0091.MODE_SST2_EFF', 'V V R', nbeq2, leff2)
!
    call wkvect('&&OP0091.MODE_SST1_EFFWI', 'V V R', nbeq1, lefi1)
    call wkvect('&&OP0091.MODE_SST2_EFFWI', 'V V R', nbeq2, lefi2)
!
    call wkvect('&&OP0091.MODE_SST1_COPY', 'V V R', nbeq1, lcopy1)
    call wkvect('&&OP0091.MODE_SST2_COPY', 'V V R', nbeq2, lcopy2)
!--
!-- RECHERCHE DES MATRICES DE MASSE ET DE RAIDEUR
!--
    call jenonu(jexnom(modgen//'      .MODG.SSNO', sst1), isst1)
    call jenonu(jexnom(modgen//'      .MODG.SSNO', sst2), isst2)
!-- SST1
    call jeveuo(jexnum(modgen//'      .MODG.SSME', isst1), 'L', ibid)
    call jeveuo(zk8(ibid)//'.MAEL_RAID_REFE', 'L', lbid)
    call mtdscr(zk24(lbid+1) (1:19))
    call jeveuo(zk24(lbid+1) (1:19)//'.&INT', 'L', matrice_raid(isst1))
    call jeveuo(zk8(ibid)//'.MAEL_MASS_REFE', 'L', lbid)
    call mtdscr(zk24(lbid+1) (1:19))
    call jeveuo(zk24(lbid+1) (1:19)//'.&INT', 'L', matrice_mass(isst1))
!-- SST2
    call jeveuo(jexnum(modgen//'      .MODG.SSME', isst2), 'L', ibid)
    call jeveuo(zk8(ibid)//'.MAEL_RAID_REFE', 'L', lbid)
    call mtdscr(zk24(lbid+1) (1:19))
    call jeveuo(zk24(lbid+1) (1:19)//'.&INT', 'L', matrice_raid(isst2))
    call jeveuo(zk8(ibid)//'.MAEL_MASS_REFE', 'L', lbid)
    call mtdscr(zk24(lbid+1) (1:19))
    call jeveuo(zk24(lbid+1) (1:19)//'.&INT', 'L', matrice_mass(isst2))
!
!
!--
!-- BOUCLE SUR LES MODES PROPRES
!--
    do j1 = 1, nbmod
!--
!-- SST1
!--
!-- RECUPERATION DU MODE RESTITUE
        call jeveuo(zk24(tach1+j1-1) (1:19)//'.VALE', 'L', ibid)
!
!-- RECOPIE DANS UN VECTEUR DE TRAVAIL
        call lceqvn(nbeq1, zr(ibid), zr(lcopy1))
!
!-- RECHERCHE DU MACRO ELEMENT ASSOCIE A LA SST
        call jeveuo(jexnum(modgen//'      .MODG.SSME', isst1), 'L', ibid)
        call jeveuo(zk8(ibid)//'      .NUME.DEEQ', 'L', ideeq)
!
!-- ANULATION DES LAGRANGES
        call zerlag(nbeq1, zi(ideeq), vectr=zr(lcopy1))
!
!-- EXTRACTION DES COMPOSANTES ASSOCIEES A L'INTERFACE
        ibid = 0
        do k1 = 1, nbddl1
            if (zi(llint1+k1-1) .gt. 0) then
                zr(lmod1+zi(llint1+k1-1)-1) = zr(lcopy1+zi(llint1+k1-1)-1)
                zr(lobs1+ibid) = zr(lcopy1+zi(llint1+k1-1)-1)
                ibid = ibid+1
            end if
        end do
!
!-- CALCUL DU TRAVAIL
        call mrmult('ZERO', matrice_mass(isst1), zr(lcopy1), zr(leff1), 1, &
                    .true._1)
        call lceqvn(nbeq1, zr(leff1), zr(lefi1))
        b_n = to_blas_int(nbeq1)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        travm = ddot(b_n, zr(lmod1), b_incx, zr(leff1), b_incy)
        call mrmult('ZERO', matrice_raid(isst1), zr(lcopy1), zr(leff1), 1, &
                    .true._1)
        b_n = to_blas_int(nbeq1)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        travk = ddot(b_n, zr(lmod1), b_incx, zr(leff1), b_incy)
        trvint = travk-(zr(lomeg+j1-1)**2)*travm
!--
!-- SST2
!--
!-- RECUPERATION DU MODE RESTITUE
        call jeveuo(zk24(tach2+j1-1) (1:19)//'.VALE', 'L', ibid)
!
!-- RECOPIE DANS UN VECTEUR DE TRAVAIL
        call lceqvn(nbeq2, zr(ibid), zr(lcopy2))
!
!-- RECHERCHE DU MACRO ELEMENT ASSOCIE A LA SST
        call jeveuo(jexnum(modgen//'      .MODG.SSME', isst2), 'L', ibid)
        call jeveuo(zk8(ibid)//'      .NUME.DEEQ', 'L', ideeq)
!
!-- ANULATION DES LAGRANGES
        call zerlag(nbeq2, zi(ideeq), vectr=zr(lcopy2))
!
!-- EXTRACTION DES COMPOSANTES ASSOCIEES A L'INTERFACE
        ibid = 0
        do k1 = 1, nbddl2
            if (zi(llint2+k1-1) .gt. 0) then
                zr(lmod2+zi(llint2+k1-1)-1) = zr(lcopy2+zi(llint2+k1-1)-1)
                zr(lobs2+ibid) = zr(lcopy2+zi(llint2+k1-1)-1)
                ibid = ibid+1
            end if
        end do
!
!--
!-- CALCUL DE L'EFFORT RESIDUEL
!--
        call mrmult('ZERO', matrice_mass(isst2), zr(lcopy2), zr(leff2), 1, &
                    .true._1)
        call lceqvn(nbeq2, zr(leff2), zr(lefi2))
        b_n = to_blas_int(nbeq2)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        travm = ddot(b_n, zr(lmod2), b_incx, zr(leff2), b_incy)
!
        call mrmult('ZERO', matrice_raid(isst2), zr(lcopy2), zr(leff2), 1, &
                    .true._1)
        b_n = to_blas_int(nbeq2)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        travk = ddot(b_n, zr(lmod2), b_incx, zr(leff2), b_incy)
!
        trvint = trvint+travk-(zr(lomeg+j1-1)**2)*travm
        if (zr(lomeg+j1-1) .gt. 1) trvint = trvint/zr(lomeg+j1-1)
!
        trav_interf(1+nbmod*(i1-1)+j1-1) = trvint
        write (unit, *) 'MODE ', j1, ' - TRAVAIL INTERFACE =', trvint
!
!--
!-- PROJECTION DES MOUVEMENTS DE L'INTERFACE MAITRE
!--
        if (imast1 .eq. -1) then
            lbid = lobs1
            ibid = lobs2
        else
            lbid = lobs2
            ibid = lobs1
        end if
!
        do k1 = 1, nl
!-- STOCKAGE DE LA PROJECTION DES MVTS MAITRES SUR L'INT. ESCLAVE
            do l1 = 1, nc
                zr(lesc+(k1-1)+(j1-1)*nl) = zr( &
                                            lesc+(k1-1)+(j1-1)*nl)+tr_mod_mast_pro(1+(k1-1)+(l1-&
                                            &1)*nl)*zr(lbid+l1-1 &
                                            )
            end do
!-- STOCKAGE DES MVTS DE L'INT. ESCLAVE POUR EXPANSION
            zr(lmast+(k1-1)+(j1-1)*nl) = zr(ibid+k1-1)
        end do
!
!-- DEPLACEMENT IMPOSE DE L'INTERFACE ESCLAVE
!
        l1 = 0
        do k1 = 1, nbddl
            if (zi(llint+k1-1) .gt. 0) then
                zr(ldepsl+zi(lag1+k1-1)-1+nbeq*(j1-1)) = zr(lesc+l1+(j1-1)*nl)
                zr(ldepsl+zi(lag2+k1-1)-1+nbeq*(j1-1)) = zr(lesc+l1+(j1-1)*nl)
                l1 = l1+1
            end if
        end do
!
!--
!-- CALCUL DU TRAVAIL ASSOCIE AU DEPLACEMENT DIFERENTIEL DE L'INTERFACE
!--
!-- PAS AU POINT DU TOUT. FAUT TROUVER UNE NORME AD HOC
!-- JE COMMENTE TANT QUE J'AI PAS TROUVE
!
        if (1 .gt. 2) then
            travk = 0.d0
            l1 = 0
            if (imast1 .eq. -1) then
!
                do k1 = 1, nbddl
                    if (zi(llint+k1-1) .gt. 0) then
                        travk = travk+( &
                                zr(lobs2+l1)-zr(lesc+l1+(j1-1)*nl))*(zr(leff2+zi(llint+k1-1)-1)-&
                                &zr(lomeg+j1-1)*zr(lefi2+zi(llint+k1-1)-1) &
                                )
                        l1 = l1+1
                    end if
                end do
!
            else
                do k1 = 1, nbddl
                    if (zi(llint+k1-1) .gt. 0) then
                        travk = travk+( &
                                zr(lobs1+l1)-zr(lesc+l1+(j1-1)*nl))*(zr(leff1+zi(llint+k1-1)-1)-&
                                &zr(lomeg+j1-1)*zr(lefi1+zi(llint+k1-1)-1) &
                                )
                        l1 = l1+1
                    end if
                end do
!
            end if
            if (zr(lomeg+j1-1) .gt. 1) travk = travk/zr(lomeg+j1-1)
!
            write (unit, *) '                - DECOLLEMENT =', travk
            write (unit, *) ' '
!
        end if
!
!
    end do
!
end subroutine
