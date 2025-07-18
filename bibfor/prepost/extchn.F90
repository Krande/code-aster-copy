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

subroutine extchn(nchmno, nnoeud, numnd, ncmp, nbn, &
                  nbc, indic, nsschn, mcf, iocc)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8dgrd.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/exchnn.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecreo.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/matrot.h"
#include "asterfort/numek8.h"
#include "asterfort/rvchn1.h"
#include "asterfort/rvchn2.h"
#include "asterfort/rvrecu.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/char8_to_int.h"
!
    integer(kind=8) :: nbn, nbc, numnd(*), iocc
    character(len=6) :: indic
    character(len=8) :: nnoeud(*), ncmp(*)
    character(len=19) :: nchmno, nsschn
    character(len=*) :: mcf
!
!     OPERATION REALISEE
!     ------------------
!
!       CONSTRUCTION DE LA SD ASSOCIEE A L' EXTRACTION SUR UN CHAM_NO
!
!     ARGUMENTS
!     ---------
!
!       NCHMNO (IN) : NOM DE LA SD DE TYPE CHAM_NO SUR LAQUELLE
!                     ON EXTRAIT
!
!       NNOEUD (IN) : TABLEAU DES NOMS DE NOEUDS SUR LESQUELS
!                     ON EXTRAIT
!
!       NUMND  (IN) : TABLEAU DES NUMEROS DE NOEUDS SUR LESQUELS
!                     ON EXTRAIT
!
!       NCMP   (IN) : TABLEAU DES NOMS DE COMPOSANTES SUR LESQUELS
!                     ON EXTRAIT
!
!       NBN    (IN) : NBR DE NOEUDS DE NNOEUD
!
!       NBC    (IN) : NBR DE COMPOSANTES DE NCMP
!
!       INDIC  (IN) : INDIQUE SI LES NOEUDS MIS EN JEU SONT DONNES
!                     PAR NOMS ('NOMME', ALORS NNOEUDS EST UTILISE)
!                     OU PAR NUMEROS ('NUMERO', ALORS NUMND EST UTILISE)
!
!       NSSCHN (IN) : NOM DE LA SD DE TYPE SOUS_CHAM_NO CONSTRUITE
!
!                      .VALE     : OJB V R8 --> SGT DE VALEURS
!                                  DOCU = 'CHNO'
!
!                      .PADR     : OJB V I  --> POINTEUR DES NOEUDS
!                                               SUR .VALE
!
!                      .PCMP     : OJB V I  --> TABLE DES CMP ACTIVES
!                                               POUR L' EXTRACTION
!
!                      .NOMA     : OJB E K8 --> NOM DU MAILLAGE
!
!                      .NUGD     : OJB E I  --> NUMERO DE LA GRANDEUR
!
!                      .ERRE     : XD  V I  --> VERIFICATION DE LA
!                                               DEFINITION DES CMP SUR
!                                               LES NOEUDS
!
!     ORGANISATION DE LA SD SOUS_CHAM_NO
!     ----------------------------------
!
!       SI PADR(IN) = 0 ALORS
!       --              -----
!
!           LE NOEUD NUMERO IN N' EST PAS CONCERNE PAR L' EXTRACTION
!
!       SINON
!       -----
!
!           PADR(IN) EST L' ADRESSE DU SOUS-SGT DE VALEURS
!           DANS VALE, ASSOCIEE AU NOEUD NUMERO IN
!
!       FSI
!       ---
!       SI PCMP(ICMP) = 0 ALORS
!       --                -----
!
!           LA CMP NUMERO ICMP N' EST PAS CONCERNE PAR L' EXTRACTION
!
!       SINON
!       -----
!
!           PCMP(ICMP) EST L' ADRESSE DE LA VALEUR DE LA CMP
!           NUMERO ICMP (DANS LE CATALOGUE DES GRANDEURS) DANS
!           TOUS LES SOUS-SGT ASSOCIES AUX NOEUDS
!
!       FSI
!       ---
!
!       INDLOC EST LE NUMERO LOCALE D' UN NOEUD DE L' EXTRACTRATION
!       IL LUI CORRESPOND UN VECTEUR V DANS LA XD '.ERRE' :
!
!          V(JLOCCMP) = 0 <=> LA CMP NUMERO JLOCCMP POUR L' EXTRACTION
!                             EST DEFINIE SUR CE NOEUD
!
!
!***********************************************************************
!
!   FONCTIONS EXTERNES
!   ------------------
!
!   NOMS ET ADRESSES DES OJB ASSOCIES AUX CHAM_NO
!   ---------------------------------------------
!
    integer(kind=8) :: avalch, anueq
    character(len=3) :: type
    character(len=24) :: nvalch, nnueq, nomvec
!
!   NOMS ET ADREESES DES OJB ASSOCIES AUX SOUS_CHAM_NO
!   --------------------------------------------------
!
    integer(kind=8) :: apadr, apcmp, apval, anugd, aperr, anoma
    character(len=24) :: npadr, npcmp, npval, nnugd, nperr, nnoma
!
!   NOMS ET ADRESSES DES OJB ASSOCIES AUX NUMEEQUA
!   ----------------------------------------------
!
    integer(kind=8) :: aprno
    character(len=24) :: nprno
    character(len=19) :: nprof
!
!   VARIABLES ASSOCIEES A LA GRANDEUR
!   ---------------------------------
!
    integer(kind=8) :: gd, acmpgd, adesgd, nbtcmp, nbec
!
!   VARIABLES ASSOCIEES AU MAILLAGE
!   -------------------------------
!
    integer(kind=8) :: nbtnd
    character(len=8) :: nmaila
!
!   VARIABLES COMPLEMENTAIRES
!   -------------------------
!
    integer(kind=8) :: anumcp, i, ind, n1, ibid
    real(kind=8) :: angl(3), pgl(3, 3), orig(3), axez(3)
    real(kind=8) :: zero, xnormz, epsi
    aster_logical :: utili
    character(len=8) :: repere
    character(len=24) :: nomjv
! ----------------------------------------------------------------------
!
    call jemarq()
!
    nvalch = nchmno//'.VALE'
    ibid = 0
    zero = 0.0d0
    epsi = 1.0d-6
!
    call jelira(nvalch, 'TYPE', cval=type)
    ASSERT((type(1:1) .eq. 'R') .or. (type(1:1) .eq. 'C'))
    if (type(1:1) .eq. 'R') then
        call jeveuo(nvalch, 'L', avalch)
    else if (type(1:1) .eq. 'C') then
        nomvec = 'EXTCHN.VECTEUR'
        call rvrecu(mcf, iocc, nchmno, nomvec)
        call jeveuo(nomvec, 'L', avalch)
    end if
!
!   EST-CE UN REPERE "UTILISATEUR"
!   ------------------------------
!
    utili = .false.
    if (mcf(1:6) .eq. 'ACTION') then
        call getvtx(mcf, 'REPERE', iocc=iocc, scal=repere, nbret=n1)
        if (repere .eq. 'UTILISAT') then
            utili = .true.
            nomjv = '&&EXTCHN.NEW_CHAMP'
            call getvr8(mcf, 'ANGL_NAUT', iocc=iocc, nbval=3, vect=angl, &
                        nbret=n1)
            angl(1) = angl(1)*r8dgrd()
            angl(2) = angl(2)*r8dgrd()
            angl(3) = angl(3)*r8dgrd()
            call matrot(angl, pgl)
            call rvchn1(nchmno, nomjv, nbn, numnd, pgl)
            call jeveuo(nomjv, 'L', avalch)
        else if (repere .eq. 'CYLINDRI') then
            utili = .true.
            nomjv = '&&EXTCHN.NEW_CHAMP'
            call getvr8(mcf, 'ORIGINE', iocc=iocc, nbval=3, vect=orig, &
                        nbret=n1)
            call getvr8(mcf, 'AXE_Z', iocc=iocc, nbval=3, vect=axez, &
                        nbret=n1)
            xnormz = zero
            do i = 1, 3
                xnormz = xnormz+axez(i)*axez(i)
            end do
            if (xnormz .lt. epsi) then
                call utmess('F', 'PREPOST_38')
            end if
            xnormz = 1.0d0/sqrt(xnormz)
            do i = 1, 3
                axez(i) = axez(i)*xnormz
            end do
            call rvchn2(nchmno, nomjv, nbn, numnd, orig, &
                        axez)
            call jeveuo(nomjv, 'L', avalch)
        end if
    end if
!
!   RECUPERATION DE L' INFORMATION
!   ------------------------------
!
    nnoma = nsschn//'.NOMA'
    nnugd = nsschn//'.NUGD'
!
    call dismoi('NUM_GD', nchmno, 'CHAM_NO', repi=gd)
    call dismoi('NOM_MAILLA', nchmno, 'CHAM_NO', repk=nmaila)
!
    call wkvect(nnoma, 'V V K8', 1, anoma)
!
    zk8(anoma) = nmaila
!
    call wkvect(nnugd, 'V V I', 1, anugd)
!
    zi(anugd) = gd
!
!   RECUPERATION DE LA DESCRIPTION DE LA GRANDEUR DANS LE CATALOGUE
!   ---------------------------------------------------------------
!
    call jeveuo(jexnum('&CATA.GD.NOMCMP', gd), 'L', acmpgd)
    call jelira(jexnum('&CATA.GD.NOMCMP', gd), 'LONMAX', nbtcmp)
    call jeveuo(jexnum('&CATA.GD.DESCRIGD', gd), 'L', adesgd)
!
    nbec = zi(adesgd+3-1)
!
!   RECUPERATION DU NOMBRE TOTAL DE NOEUD DANS LE MAILLAGE
!   ------------------------------------------------------
!
    call dismoi('NB_NO_MAILLA', nmaila, 'MAILLAGE', repi=nbtnd)
!
!   TRADUCTION DES NOMS DE NOEUDS ET DE CMP EN NUMEROS
!   --------------------------------------------------
!
    if (indic .eq. 'NOMME') then
!
        do i = 1, nbn, 1
!
            numnd(i) = char8_to_int(nnoeud(i))
!
        end do
!
    end if
!
    call jecreo('&&EXTRCHNNUMCP', 'V V I')
    call jeecra('&&EXTRCHNNUMCP', 'LONMAX', nbc)
    call jeveuo('&&EXTRCHNNUMCP', 'E', anumcp)
!
    do i = 1, nbc, 1
!
        zi(anumcp+i-1) = 0
!
    end do
!
    call numek8(zk8(acmpgd), ncmp, nbtcmp, nbc, zi(anumcp))
!
!   CREATIONS DES OJB DE LA SD SOUS_CHAM_NO
!   ---------------------------------------
!
    npadr = nsschn//'.PADR'
    npcmp = nsschn//'.PCMP'
    npval = nsschn//'.VALE'
    nperr = nsschn//'.ERRE'
!
    call jecreo(npadr, 'V V I')
    call jeecra(npadr, 'LONMAX', nbtnd)
    call jeveuo(npadr, 'E', apadr)
!
    call jecreo(npcmp, 'V V I')
    call jeecra(npcmp, 'LONMAX', nbtcmp)
    call jeveuo(npcmp, 'E', apcmp)
!
    call jecreo(npval, 'V V R')
    call jeecra(npval, 'LONMAX', nbc*nbn)
    call jeecra(npval, 'DOCU', cval='CHNO')
    call jeveuo(npval, 'E', apval)
!
    call jecrec(nperr, 'V V I', 'NU', 'DISPERSE', 'VARIABLE', &
                nbn)
!
!   REMPLISSAGE DU POINTEUR D' ADRESSES
!   -----------------------------------
!
    do i = 1, nbtnd, 1
!
        zi(apadr+i-1) = 0
!
    end do
!
!
    do i = 1, nbn, 1
!
        zi(apadr+numnd(i)-1) = nbc*(i-1)+1
!
    end do
!
!   REMPLISSAGE DU POINTEUR DES POSITIONS DES CMP
!   ---------------------------------------------
!
    do i = 1, nbtcmp, 1
!
        zi(apcmp+i-1) = 0
!
    end do
!
    do i = 1, nbc, 1
!
!CC          ZI(APCMP+ZI(ANUMCP+I-1)-1) = I
        zi(apcmp+i-1) = zi(anumcp+i-1)
!
    end do
!
!   REMPISSAGE DU SGT DE VALEURS
!   ----------------------------
!
!
!          /* CAS D' UN CHAM_NO A REPRESENTATION PAR NOEUD */
!          /* CONTENUE DANS LE PRNO (OC 1) DU NUME_EQUA    */
!
!        RECUPERATION DU NOM DU NUME_EQUA
!        --------------------------------
!
    call dismoi('NUME_EQUA', nchmno, 'CHAM_NO', repk=nprof)
!
!        ACCES AU DESCRIPTEUR DU CHAM_NO
!        -------------------------------
!
    nprno = nprof//'.PRNO'
    nnueq = nprof//'.NUEQ'
    call jeveuo(nnueq, 'L', anueq)
!
    call jeveuo(jexnum(nprno, 1), 'L', aprno)
!
!        RECUPERATION DES VALEURS
!        ------------------------
!
!
    do i = 1, nbn
!
        call jecroc(jexnum(nperr, i))
        call jeecra(jexnum(nperr, i), 'LONMAX', nbc)
        call jeveuo(jexnum(nperr, i), 'E', aperr)
!
        ind = numnd(i)
!
        call exchnn(zi(aprno+(ind-1)*(2+nbec)+1-1), 0, zi(anumcp), nbc, zr(avalch), &
                    zi(anueq), .true._1, zr(apval+nbc*(i-1)+1-1), zi(aperr))
!
    end do
!
!
!   DESTRUCTION DES OJB TEMPORAIRES
!   -------------------------------
!
    if (utili) call jedetr(nomjv)
    call jedetr('&&EXTRCHNNUMCP')
    if (type(1:1) .eq. 'C') call jedetr(nomvec)
!
    call jedema()
end subroutine
