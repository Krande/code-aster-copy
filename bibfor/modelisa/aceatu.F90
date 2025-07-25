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

subroutine aceatu(noma, nomo, nbepo, ntyele, ivr, nbocc)
    use cara_elem_parameter_module
    implicit none
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/aceat2.h"
#include "asterfort/aceat3.h"
#include "asterfort/acemmt.h"
#include "asterfort/getvem.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
!
    integer(kind=8) :: nbepo, ntyele(*), nbocc(*), ivr(*), ifm
    character(len=8) :: noma, nomo
!     AFFE_CARA_ELEM
!     AFFECTATION DES CARACTERISTIQUES GEOMETRIQUES POUR LES TUYAUX
!     PAR CREATION D'UNE CARTE : NOMU//'.CAORTU' CONTENANT :
!        PGL1, PGL2, PGL3 : MATRICES DE PASSAGE REPERE LOCAL / GLOBAL
!        ICOUDE           : =0 SI ELEMENT DROIT =1 SI COUDE
!        DN1N2            : DISTANCE ENTRE LES DEUX SOMMETS
!        RCOURB           : RAYON DE COURBURE DU COUDE
!        ANGCOU           : ANGLE DU COUDE
!        ANGZZK           : ANGLE OMEGA ENTRE LA NORMALE AU PLAN DU
!                           COUDE ET LE VECTEUR ZK (GENERATRICE)
! ----------------------------------------------------------------------
! IN  : NOMA   : NOM DU MAILLAGE
! IN  : NOMO   : NOM DU MODELE
! IN  : NBEPO  : NOMBRE DE TYPES D'ELEMENTS
! IN  : NTYELE : NUMEROS DES TYPES ELEMENTS
! IN  : IVR    : (3) = INFO 2
! IN  : IFM    : FICHIER MESSAGES
! IN  : NBOCC  : NBOCC(4) NB OCCURENCES ORIENTATION
! ----------------------------------------------------------------------
!
    integer(kind=8) :: iext1, iext2, ima, inn, ioc, jcozk, jdco, jdgn, jdno, jdme
    integer(kind=8) ::   jma, jnozk
    integer(kind=8) :: nbext2, nbpart, nbtuy, ncar, ni1, ni2, nj, nj1, nj2, nng
    integer(kind=8) ::  nutyel, nval, ixma, j
    integer(kind=8) ::  nno, nbtuy4, nbext1, jzkpar, ibid
    integer(kind=8) :: ier, nbmail
    real(kind=8) :: val(3), epsi
    character(len=8) :: nomu, car, crit
    character(len=16) :: concep, cmd, nunoel
    character(len=24) :: mlggno, mlgcoo, mlgcnx, modmai, nomlu
    integer(kind=8), pointer :: eltuy(:) => null()
    integer(kind=8), pointer :: lismapart(:) => null()
    integer(kind=8), pointer :: lisnopart(:) => null()
    integer(kind=8), pointer :: mmt(:) => null()
    integer(kind=8), pointer :: nbmapart(:) => null()
    integer(kind=8), pointer :: noext1(:) => null()
    integer(kind=8), pointer :: noext2(:) => null()
    integer(kind=8), pointer :: notuy(:) => null()
    integer(kind=8), pointer :: sens(:) => null()
!     ------------------------------------------------------------------
!
    call jemarq()
    ier = 0
    call getres(nomu, concep, cmd)
!
!   Reconstruction des noms jeveux du concept maillage associe
    mlgcnx = noma//'.CONNEX'
    mlggno = noma//'.GROUPENO'
    mlgcoo = noma//'.COORDO    .VALE'
    call jelira(noma//'.TYPMAIL', 'LONMAX', nbmail)
    call jeveuo(mlgcoo, 'L', jdco)
!
    modmai = nomo//'.MAILLE'
    call jeexin(modmai, ixma)
    if (ixma .ne. 0) call jeveuo(modmai, 'L', jdme)
!
!   Comptage des MET3SEG3
    nbtuy = 0
    do ima = 1, nbmail
        nutyel = zi(jdme+ima-1)
        do j = 1, nbepo
            if (nutyel .eq. ntyele(j)) then
                call jenuno(jexnum('&CATA.TE.NOMTE', nutyel), nunoel)
                if ((nunoel .eq. 'MET3SEG3') .or. (nunoel .eq. 'MET6SEG3') .or. &
                    (nunoel .eq. 'MET3SEG4')) then
                    nbtuy = nbtuy+1
                end if
            end if
        end do
    end do
!
!   Stockage des éléments MET3SEG3 et des noeuds
    AS_ALLOCATE(vi=notuy, size=nbtuy*4)
    AS_ALLOCATE(vi=eltuy, size=nbtuy)
    nno = 0
    nbtuy = 0
    do ima = 1, nbmail
        nutyel = zi(jdme+ima-1)
        do j = 1, nbepo
            if (nutyel .eq. ntyele(j)) then
                call jenuno(jexnum('&CATA.TE.NOMTE', nutyel), nunoel)
                if ((nunoel .eq. 'MET3SEG3') .or. (nunoel .eq. 'MET6SEG3')) then
                    nno = 3
                    nbtuy = nbtuy+1
                    eltuy(nbtuy) = ima
                    call jeveuo(jexnum(mlgcnx, ima), 'L', jdno)
                    notuy(3*nbtuy-2) = zi(jdno)
                    notuy(3*nbtuy-1) = zi(jdno+1)
                    notuy(3*nbtuy) = zi(jdno+2)
                end if
            end if
        end do
    end do
!
    nbtuy4 = 0
    do ima = 1, nbmail
        nutyel = zi(jdme+ima-1)
        do j = 1, nbepo
            if (nutyel .eq. ntyele(j)) then
                call jenuno(jexnum('&CATA.TE.NOMTE', nutyel), nunoel)
                if (nunoel .eq. 'MET3SEG4') then
                    nno = 4
                    nbtuy4 = nbtuy4+1
                    eltuy(nbtuy4) = ima
                    call jeveuo(jexnum(mlgcnx, ima), 'L', jdno)
                    notuy(4*nbtuy4-3) = zi(jdno)
                    notuy(4*nbtuy4-2) = zi(jdno+1)
                    notuy(4*nbtuy4-1) = zi(jdno+2)
                    notuy(4*nbtuy4) = zi(jdno+3)
                end if
            end if
        end do
    end do
!
    if (nbtuy4 .ne. 0) then
        if (nbtuy .ne. 0) then
            call utmess('F', 'MODELISA_27')
        else
            nbtuy = nbtuy4
        end if
    end if
!
!   Comptage des parties connexes. HYPOTHESE : LES MAILLES SONT TOUTES ORIENTÉES DANS LE MÊME SENS
    nbext1 = 0
    nbext2 = 0
    do ima = 1, nbtuy
        iext1 = 0
        iext2 = 0
        ni1 = notuy(nno*(ima-1)+1)
        ni2 = notuy(nno*(ima-1)+2)
        do jma = 1, nbtuy
            if (jma .ne. ima) then
                nj1 = notuy(nno*(jma-1)+1)
                nj2 = notuy(nno*(jma-1)+2)
                if (ni1 .eq. nj2) then
                    iext1 = 1
                end if
                if (ni2 .eq. nj1) then
                    iext2 = 1
                end if
            end if
        end do
        if (iext1 .eq. 0) then
            nbext1 = nbext1+1
        end if
        if (iext2 .eq. 0) then
            nbext2 = nbext2+1
        end if
    end do
    if (nbext1 .ne. nbext2) then
        call utmess('F', 'MODELISA10_4')
    end if
    nbpart = nbext1
    ifm = ivr(4)
    if (ivr(3) .eq. 2) then
        write (ifm, *) 'NOMBRE DE PARTIES CONNEXES DE TUYAU : ', nbpart
    end if
!
!   Vérification et stockage des parties connexes.
    AS_ALLOCATE(vi=sens, size=nbpart)
    AS_ALLOCATE(vi=nbmapart, size=nbpart)
    AS_ALLOCATE(vi=noext1, size=nbpart)
    AS_ALLOCATE(vi=noext2, size=nbpart)
    AS_ALLOCATE(vi=lismapart, size=nbpart*nbtuy)
    AS_ALLOCATE(vi=lisnopart, size=nbpart*nbtuy*nno)
    call wkvect('&&ACEATU.ZKPART', 'V V I', nbpart*nbtuy*nno, jzkpar)
    call aceat2(nbtuy, eltuy, notuy, nbpart, noext1, noext2, nbmapart, lismapart, lisnopart, nno)
!
!   Lecture de MODI_METRIQUE
    AS_ALLOCATE(vi=mmt, size=nbmail)
    call acemmt(noma, mmt)
!
!   Lecture du mot-clef GENE_TUYAU
    inn = 0
!   Valeurs par défaut cohérentes avec le catalogue
    epsi = 1.d-4
    crit = 'RELATIF'
!   Pour ne pas passer des variables non-initialisées en argument
    jnozk = 1
    jcozk = 1
    if (nbocc(ACE_ORIENTATION) .ne. 0) then
        call wkvect('&&ACEATU.LISNOZK', 'V V I', nbocc(ACE_ORIENTATION), jnozk)
        call wkvect('&&ACEATU.LISCOZK', 'V V R', 3*nbocc(ACE_ORIENTATION), jcozk)
        do ioc = 1, nbocc(ACE_ORIENTATION)
!           Un seul noeud permis
            call getvem(noma, 'GROUP_NO', 'ORIENTATION', 'GROUP_NO', ioc, 1, nomlu, nj)
            call getvtx('ORIENTATION', 'CARA', iocc=ioc, scal=car, nbret=ncar)
            call getvr8('ORIENTATION', 'VALE', iocc=ioc, nbval=3, vect=val, nbret=nval)
            call getvr8('ORIENTATION', 'PRECISION', iocc=ioc, scal=epsi, nbret=ibid)
            if (ibid .eq. 0) then
                epsi = 1.d-4
                crit = 'RELATIF'
            end if
            call getvtx('ORIENTATION', 'CRITERE', iocc=ioc, scal=crit, nbret=ibid)
            if (car .eq. 'GENE_TUY') then
                if (nj .gt. 0) then
                    if (nj .eq. 1) then
                        call jeveuo(jexnom(mlggno, nomlu), 'L', jdgn)
                        call jelira(jexnom(mlggno, nomlu), 'LONUTI', nng)
                        if (nng .eq. 1) then
                            inn = inn+1
                            zi(jnozk-1+inn) = zi(jdgn)
                            zr(jcozk-1+3*inn-2) = val(1)
                            zr(jcozk-1+3*inn-1) = val(2)
                            zr(jcozk-1+3*inn) = val(3)
                        else
                            ier = 1
                            goto 999
                        end if
                    else
                        ier = 1
                        goto 999
                    end if
                end if
            end if
        end do
    end if
    call aceat3(noma, nomu, nbtuy, nbpart, nbmapart, &
                lismapart, lisnopart, ivr, inn, &
                zi(jnozk), zr(jcozk), sens, zr(jdco), epsi, &
                crit, nno, mmt)
!
999 continue
    if (ier .ne. 0) then
        call utmess('F', 'MODELISA_28')
    end if
!
!   Ménage
    AS_DEALLOCATE(vi=notuy)
    AS_DEALLOCATE(vi=eltuy)
    AS_DEALLOCATE(vi=sens)
    AS_DEALLOCATE(vi=nbmapart)
    AS_DEALLOCATE(vi=noext1)
    AS_DEALLOCATE(vi=noext2)
    AS_DEALLOCATE(vi=lismapart)
    AS_DEALLOCATE(vi=lisnopart)
    call jedetr('&&ACEATU.ZKPART')
    AS_DEALLOCATE(vi=mmt)
    call jedetr('&&ACEATU.LISNOZK')
    call jedetr('&&ACEATU.LISCOZK')
!
    call jedema()
end subroutine
