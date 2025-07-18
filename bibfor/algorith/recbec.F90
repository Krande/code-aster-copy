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
subroutine recbec(nomres, typesd, basmod, modcyc, numsec)
    implicit none
!-----------------------------------------------------------------------
!
!  BUT:  < RESTITUTION CRIAG-BAMPTON ECLATEE >
!
!   RESTITUER LES RESULTATS ISSUS D'UN CALCUL CYCLIQUE
!          AVEC DES INTERFACES DE TYPE CRAIG-BAMPTON
!     => RESULTAT COMPOSE DE TYPE MODE_MECA ALLOUE PAR LA
!   ROUTINE
!-----------------------------------------------------------------------
!
! NOMRES  /I/: NOM UT DU CONCEPT RESULTAT A REMPLIR
! TYPESD  /I/: NOM K16 DU TYPE DE LA STRUCTURE DE DONNEE
! BASMOD  /I/: NOM UT DE LA BASE MODALE EN AMONT DU CALCUL CYCLIQUE
! MODCYC  /I/: NOM UT DU RESULTAT ISSU DU CALCUL CYCLIQUE
! NUMSEC  /I/: NUMERO DU SECTEUR  SUR LEQUEL RESTITUER
!
!
!
!
#include "jeveux.h"
#include "asterc/r8depi.h"
#include "asterfort/bmnodi.h"
#include "asterfort/ctetgd.h"
#include "asterfort/dismoi.h"
#include "asterfort/genecy.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mtdscr.h"
#include "asterfort/mtexis.h"
#include "asterfort/ordr8.h"
#include "asterfort/recbbn.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rscrsd.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnoch.h"
#include "asterfort/utmess.h"
#include "asterfort/vtcrem.h"
#include "asterfort/wkvect.h"
#include "blas/daxpy.h"
!
    character(len=8) :: nomres, basmod, modcyc, kbid, k8b
    character(len=16) :: depl, typesd, typsup(1)
    character(len=19) :: chamva, numddl, matrix, mass
    character(len=24) :: tetgd, nomvec
    character(len=24) :: valk(2)
    complex(kind=8) :: dephc
    real(kind=8) :: para(2), depi, fact, genek, beta
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iad, ibid, icomp, iddi, idi, idia
    integer(kind=8) :: idiam, idicou, ier, ii, inum, iorc, iormo
    integer(kind=8) :: j, jj, ldfre, ldkge, ldmge, ldom2, ldomo
    integer(kind=8) :: ldotm, ldtyd, llcham, llmoc
    integer(kind=8) :: lmass, ltetgd, ltora, ltord
    integer(kind=8) :: ltorf, ltorg, ltorto, ltveco, ltvere, ltvezt, mdiapa
    integer(kind=8) :: nbdax, nbddg, nbddr, nbdia, nbmoc, nbmod, nbmor
    integer(kind=8) :: nborc, nbsec, nbtmp, neq, numa, numd, numg
    integer(kind=8) :: numsec
    real(kind=8) :: aaa, bbb, betsec
    integer(kind=8), pointer :: cycl_diam(:) => null()
    integer(kind=8), pointer :: cycl_nuin(:) => null()
    integer(kind=8), pointer :: cycl_desc(:) => null()
    integer(kind=8), pointer :: cycl_nbsc(:) => null()
    real(kind=8), pointer :: cycl_freq(:) => null()
    blas_int :: b_incx, b_incy, b_n
!-----------------------------------------------------------------------
    data depl/'DEPL            '/
    data typsup/'MODE_MECA       '/
!-----------------------------------------------------------------------
!
    call jemarq()
!
    depi = r8depi()
    ltora = 1
!
!----------------VERIFICATION DU TYPE DE STRUCTURE RESULTAT-------------
!
    if (typesd .ne. typsup(1)) then
        valk(1) = typesd
        valk(2) = typsup(1)
        call utmess('F', 'ALGORITH14_4', nk=2, valk=valk)
    end if
!
!--------------------------RECUPERATION DU .DESC------------------------
!
    call jeveuo(modcyc//'.CYCL_DESC', 'L', vi=cycl_desc)
    nbmod = cycl_desc(1)
    nbddr = cycl_desc(2)
    nbdax = cycl_desc(3)
!
!-------------------RECUPERATION DU NOMBRE DE SECTEUR-------------------
!
    call jeveuo(modcyc//'.CYCL_NBSC', 'L', vi=cycl_nbsc)
    nbsec = cycl_nbsc(1)
    mdiapa = int(nbsec/2)*(1-nbsec+(2*int(nbsec/2)))
!
!------------------RECUPERATION DES NOMBRES DE DIAMETRES MODAUX---------
!
    call jeveuo(modcyc//'.CYCL_DIAM', 'L', vi=cycl_diam)
    call jelira(modcyc//'.CYCL_DIAM', 'LONMAX', nbdia)
    nbdia = nbdia/2
!
!-----------------RECUPERATION DU NOMBRE DE DDL PHYSIQUES---------------
!
    call dismoi('NUME_DDL', basmod, 'RESU_DYNA', repk=numddl)
    call dismoi('REF_RIGI_PREM', basmod, 'RESU_DYNA', repk=matrix)
    call dismoi('NB_EQUA', numddl, 'NUME_DDL', repi=neq)
!
!-------------RECUPERATION DES FREQUENCES-------------------------------
!
    call jeveuo(modcyc//'.CYCL_FREQ', 'L', vr=cycl_freq)
!
!----------------RECUPERATION MATRICE DE MASSE--------------------------
!
    call dismoi('REF_MASS_PREM', basmod, 'RESU_DYNA', repk=mass, arret='C')
    call mtexis(mass, ier)
    if (ier .eq. 0) then
        valk(1) = mass(1:8)
        call utmess('F', 'ALGORITH12_39', sk=valk(1))
    end if
    call mtdscr(mass)
    call jeveuo(mass(1:19)//'.&INT', 'E', lmass)
!
!------------------ALLOCATION DES VECTEURS DE TRAVAIL-------------------
!
    call wkvect('&&RECBEC.VEC.TRAVC', 'V V C', neq, ltvezt)
    call wkvect('&&RECBEC.VEC.COMP', 'V V C', neq, ltveco)
    call wkvect('&&RECBEC.VEC.REEL', 'V V R', neq, ltvere)
!
!-----------------RECUPERATION DES NUMERO D'INTERFACE-------------------
!
    call jeveuo(modcyc//'.CYCL_NUIN', 'L', vi=cycl_nuin)
    numd = cycl_nuin(1)
    numg = cycl_nuin(2)
    numa = cycl_nuin(3)
!
!---------------RECUPERATION DU NUMERO ORDRE DES DEFORMEES--------------
!
    nomvec = '&&RECBEC.ORD.DEF.DR'
    call wkvect(nomvec, 'V V I', nbddr, ltord)
    kbid = ' '
    call bmnodi(basmod, kbid, '       ', numd, nbddr, &
                zi(ltord), ibid)
    nomvec = '&&RECBEC.ORD.DEF.GA'
    call wkvect(nomvec, 'V V I', nbddr, ltorg)
    kbid = ' '
    call bmnodi(basmod, kbid, '       ', numg, nbddr, &
                zi(ltorg), ibid)
!
    if (nbdax .gt. 0) then
        nomvec = '&&RECBEC.ORD.DEF.AX'
        call wkvect(nomvec, 'V V I', nbdax, ltora)
        kbid = ' '
        call bmnodi(basmod, kbid, '       ', numa, nbdax, &
                    zi(ltora), ibid)
    end if
!
!--------------------CLASSEMENT DES MODES PROPRES-----------------------
!               COMPTAGE DU NOMBRE DE MODES PHYSIQUES
!
    nbmoc = 0
    nbmor = 0
    do iddi = 1, nbdia
        nbtmp = cycl_diam(1+nbdia+iddi-1)
        nbmoc = nbmoc+nbtmp
        idia = cycl_diam(iddi)
        if (idia .eq. 0 .or. idia .eq. mdiapa) then
            nbmor = nbmor+nbtmp
        else
            nbmor = nbmor+2*nbtmp
        end if
    end do
    call wkvect('&&RECBEC.ORDRE.FREQ', 'V V I', nbmoc, ltorf)
    call wkvect('&&RECBEC.ORDRE.TMPO', 'V V I', nbmoc, ltorto)
    call ordr8(cycl_freq, nbmoc, zi(ltorto))
!
!
!-----------------ALLOCATION STRUCTURE DE DONNEES-----------------------
!
    call rscrsd('G', nomres, typesd, nbmor)
!
!-------DETERMINATION DES FUTUR NUMERO ORDRES DE MODES REELS------------
!
    nborc = 0
    do ii = 1, nbmoc
        iormo = zi(ltorto+ii-1)
        icomp = 0
        idicou = 0
        do jj = 1, nbdia
            icomp = icomp+cycl_diam(1+nbdia+jj-1)
            if (icomp .ge. iormo .and. idicou .eq. 0) idicou = jj
        end do
        nborc = nborc+1
        zi(ltorf+iormo-1) = nborc
        idiam = cycl_diam(idicou)
        if (idiam .ne. 0 .and. idiam .ne. mdiapa) nborc = nborc+1
    end do
    call jedetr('&&RECBEC.ORDRE.TMPO')
!
!---------------------RECUPERATION DES MODES COMPLEXES------------------
!
    call jeveuo(modcyc//'.CYCL_CMODE', 'L', llmoc)
!
!--------------CALCUL DU TETA DE CHANGEMENT DE BASE GAUCHE DROITE-------
!
    tetgd = '&&RECBEC.TETGD'
    call wkvect(tetgd, 'V V R', nbddr*nbddr, ltetgd)
    call ctetgd(basmod, numd, numg, nbsec, zr(ltetgd), &
                nbddr)
!
!--------------------------RESTITUTION----------------------------------
!
    nbddg = nbmod+nbddr+nbdax
    icomp = 0
    inum = 0
!
! --- BOUCLE SUR LES DIAMETRES NODAUX
!
    do idi = 1, nbdia
!
! ----- CALCUL DU DEPHASAGE DU SECTEUR DEMANDE
!
        idiam = cycl_diam(idi)
        beta = (depi/nbsec)*idiam
        betsec = (numsec-1)*beta
        aaa = cos(betsec)
        bbb = sin(betsec)
        dephc = dcmplx(aaa, bbb)
!
! ----- BOUCLE SUR LES MODE DU DIAMETRE COURANT
!
        do i = 1, cycl_diam(1+nbdia+idi-1)
!
            icomp = icomp+1
            inum = inum+1
            iorc = zi(ltorf+icomp-1)
!
! ------- DETERMINATION DU NUMERO DE DIAMETRE MODAL
!
            iad = llmoc+((icomp-1)*nbddg)
!
! ------- CALCUL MODE COMPLEXE SECTEUR DE BASE
!
            call recbbn(basmod, nbmod, nbddr, nbdax, tetgd, &
                        zi(ltord), zi(ltorg), zi(ltora), zc(iad), zc(ltveco), &
                        neq, beta)
!
! ------- CALCUL MASSE GENERALISEE
!
            call genecy(zc(ltveco), zc(ltveco), neq, lmass, para, &
                        nbsec, beta, beta, zc(ltvezt))
!
            do j = 1, neq
                zc(ltveco+j-1) = zc(ltveco+j-1)*dephc
                zr(ltvere+j-1) = dble(zc(ltveco+j-1))
            end do
!
! ------- RESTITUTION DU MODE PROPRE REEL (PARTIE RELLE)
!
            call rsexch(' ', nomres, depl, inum, chamva, &
                        ier)
            call vtcrem(chamva, matrix, 'G', 'R')
            call jeveuo(chamva//'.VALE', 'E', llcham)
!
! ------- COMMUN POUR MODE_MECA ET BASE_MODALE
!
            call rsadpa(nomres, 'E', 1, 'FREQ', inum, &
                        0, sjv=ldfre, styp=k8b)
            call rsadpa(nomres, 'E', 1, 'RIGI_GENE', inum, &
                        0, sjv=ldkge, styp=k8b)
            call rsadpa(nomres, 'E', 1, 'MASS_GENE', inum, &
                        0, sjv=ldmge, styp=k8b)
            call rsadpa(nomres, 'E', 1, 'OMEGA2', inum, &
                        0, sjv=ldom2, styp=k8b)
            call rsadpa(nomres, 'E', 1, 'NUME_MODE', inum, &
                        0, sjv=ldomo, styp=k8b)
            call rsadpa(nomres, 'E', 1, 'TYPE_MODE', inum, &
                        0, sjv=ldotm, styp=k8b)
!
            fact = 1.d0/(para(1)**0.5d0)
            genek = (cycl_freq(icomp)*depi)**2
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call daxpy(b_n, fact, zr(ltvere), b_incx, zr(llcham), &
                       b_incy)
            zr(ldfre) = cycl_freq(icomp)
            zr(ldkge) = genek
            zr(ldmge) = 1.d0
            zr(ldom2) = genek
            zi(ldomo) = iorc
            zk16(ldotm) = 'MODE_DYN'
!
! ------- SPECIFIQUE A BASE_MODALE
!
            call rsadpa(nomres, 'E', 1, 'TYPE_DEFO', inum, &
                        0, sjv=ldtyd, styp=k8b)
            zk16(ldtyd) = 'PROPRE          '
!
            call rsnoch(nomres, depl, inum)
!
! ------- EVENTUELLE RESTITUTION DE LA PARTIE IMAGINAIRE
!
            if (idiam .ne. 0 .and. idiam .ne. mdiapa) then
!
                do j = 1, neq
                    zr(ltvere+j-1) = dimag(zc(ltveco+j-1))
                end do
                iorc = iorc+1
                inum = inum+1
!
                call rsexch(' ', nomres, depl, inum, chamva, &
                            ier)
                call vtcrem(chamva, matrix, 'G', 'R')
                call jeveuo(chamva//'.VALE', 'E', llcham)
!
                call rsadpa(nomres, 'E', 1, 'FREQ', inum, &
                            0, sjv=ldfre, styp=k8b)
                call rsadpa(nomres, 'E', 1, 'RIGI_GENE', inum, &
                            0, sjv=ldkge, styp=k8b)
                call rsadpa(nomres, 'E', 1, 'MASS_GENE', inum, &
                            0, sjv=ldmge, styp=k8b)
                call rsadpa(nomres, 'E', 1, 'OMEGA2', inum, &
                            0, sjv=ldom2, styp=k8b)
                call rsadpa(nomres, 'E', 1, 'NUME_MODE', inum, &
                            0, sjv=ldomo, styp=k8b)
                call rsadpa(nomres, 'E', 1, 'TYPE_MODE', inum, &
                            0, sjv=ldotm, styp=k8b)
!
                fact = 1.d0/(para(2)**0.5d0)
                genek = (cycl_freq(icomp)*depi)**2
                b_n = to_blas_int(neq)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call daxpy(b_n, fact, zr(ltvere), b_incx, zr(llcham), &
                           b_incy)
                zr(ldfre) = cycl_freq(icomp)
                zr(ldkge) = genek
                zr(ldmge) = 1.d0
                zr(ldom2) = genek
                zi(ldomo) = iorc
                zk16(ldotm) = 'MODE_DYN'
!
                call rsadpa(nomres, 'E', 1, 'TYPE_DEFO', inum, &
                            0, sjv=ldtyd, styp=k8b)
                zk16(ldtyd) = 'PROPRE          '
!
                call rsnoch(nomres, depl, inum)
!
            end if
!
        end do
!
    end do
!
    call jedetr('&&RECBEC.VEC.TRAVC')
    call jedetr('&&RECBEC.VEC.COMP')
    call jedetr('&&RECBEC.VEC.REEL')
    call jedetr('&&RECBEC.ORD.DEF.DR')
    call jedetr('&&RECBEC.ORD.DEF.GA')
    call jedetr('&&RECBEC.ORDRE.FREQ')
    call jedetr('&&RECBEC.TETGD')
    if (nbdax .gt. 0) call jedetr('&&RECBEC.ORD.DEF.AX')
!
    call jedema()
end subroutine
