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

subroutine extche(nchme2, nmaile, nummai, ncmp, nbm, &
                  nbc, indic, nssche, mcf, iocc, &
                  nbnac, nnoeud)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8dgrd.h"
#include "asterfort/assert.h"
#include "asterfort/celcel.h"
#include "asterfort/celver.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exchem.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jecrec.h"
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
#include "asterfort/rvche1.h"
#include "asterfort/rvche2.h"
#include "asterfort/rvrecu.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/char8_to_int.h"
!
    integer(kind=8) :: nbm, nbc, nummai(*), iocc, nbnac, nnoeud(*)
    character(len=6) :: indic
    character(len=8) :: nmaile(*), ncmp(*)
    character(len=19) :: nchmel, nssche, nchme2
    character(len=*) :: mcf
!
!   OPERATION REALISEE
!   ------------------
!
!     CONSTRUCTION DE LA SD ASSOCIEE A L' EXTRACTION SUR UN CHAM_ELEM
!
!   ARGUMENTS
!   ---------
!
!     NCHME2 (IN) : NOM DE LA SD DE TYPE CHAM_ELEM SUR LAQUELLE
!                   ON EXTRAIT
!
!     NMAILE (IN) : TABLEAU DES NOMS DE MAILLE SUR LESQUELS
!                   ON EXTRAIT
!     NUNMAI (IN) : TABLEAU DES NUMEROS DE MAILLE SUR LESQUELS
!                   ON EXTRAIT
!     NCMP   (IN) : TABLEAU DES NOMS DE COMPOSANTES SUR LESQUELS
!                   ON EXTRAIT
!     NBM    (IN) : NBR DE MAILLE DE NMAILE
!     NBC    (IN) : NBR DE COMPOSANTES DE NCMP
!     INDIC  (IN) : INDICATEUR DU MODE DE PASSAGE DES MAILLES
!                   'NOMME'  <=> PAR NOMS,    ALORS NMAILE EST UTILISE
!                   'NUMERO' <=> PAR NUMEROS, ALORS NUMMAI EST UTILISE
!
!     NSSCHE (IN) : NOM DE LA SD DE TYPE SOUS_CHAM_NO CONSTRUITE
!
!                    .CELV     : OJB V R8 --> SGT DE VALEURS
!                                DOCU = 'CHLM'
!
!                    .PADR     : OJB V I  --> POINTEUR DES NOEUDS
!                                             SUR .CELV
!
!                    .PCMP     : OJB V I  --> TABLE DES CMP ACTIVES
!                                               POUR L' EXTRACTION
!
!                    .PNBN     : OJB V I  --> TABLE DES NBR DE NOEUDS
!                                             PAR MAILLES
!                    .PNCO     : OJB V I  --> TABLE DES NBR DE COUCHES
!                                             PAR MAILLES
!
!                    .PNSP     : OJB V I  --> TABLE DES NBR DE
!                                             SOUS POINTS PAR MAILLES
!
!                    .NOMA     : OJB E K8 --> NOM DU MAILLAGE
!
!                    .NUGD     : OJB E I  --> NUMERO DE LA GRANDEUR
!
!                    .ERRE     : XD  V I  --> VERIFICATION DE LA DEFINI-
!                                             TION DES CMP SUR LES
!                                             NOEUDS DES MAILLES
!
!     DESCRIPTION DE LA SD SOUS_CHAM_ELEM
!     -----------------------------------
!
!       SI PADR(IM) = 0 ALORS
!       --              -----
!
!           LA MAILLE NUMERO IM N' EST PAS CONCERNE PAR L' EXTRACTION
!
!       SINON
!       -----
!
!           PADR(IM) EST L' ADRESSE DU SOUS-SGT DE VALEURS
!           DANS VALE, ASSOCIEE A LA MAILLE NUMERO IM
!
!       FSI
!       ---
!
!       SI PNBN(IM) = 0 ALORS
!       --              -----
!
!           LA MAILLE NUMERO IM N' EST PAS CONCERNE PAR L' EXTRACTION
!
!       SINON
!       -----
!
!           PNBN(IM) EST LE NOMBRE DE NOEUDS DE LA MAILLE NUMERO IM
!
!       FSI
!       ---
!
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
!       MEME CONVENTION (0 OU <>0) POUR PNCO ET PNSP
!
!       IMALOC EST LE NUMERO LOCALE D' UNE MAILLE DE L' EXTRACTION
!       IL LUI CORRESPOND UN VECTEUR V DANS LA XD '.ERRE' :
!
!          V(JLOCCMP) = 0 <=> LA CMP NUMERO JLOCCMP POUR L' EXTRACTION
!                             EST DEFINIE SUR TOUS LES NOEUDS DE LA
!                             MAILLE
!
!*********************************************************************
!
!   FONCTIONS EXTERNES
!   ------------------
!
!
!   -------------------------
!
!
!   NOMS ET ADRESSES DES OJB ASSOCIES AU CHAM_ELEM
!   ---------------------------------------------
!
    character(len=3) :: type
    character(len=24) :: ndesc, nvale, ncelk, nomvec
    integer(kind=8) :: jceld, avale, acelk
!
!   NOMS ET ADRESSES DES OJB ASSOCIES AU LIGREL SOUS-JACENT
!   -------------------------------------------------------
!
    character(len=24) :: nrepe, nnoma
    character(len=19) :: nligrl
    character(len=8) :: nmaila, nomgd
    integer(kind=8) :: arepe, anoma
!
!   NOMS ET ADRESSES DES OJB ASSOCIES AU SOUS CHAM_ELEM
!   --------------------------------------------------
!
    character(len=24) :: npnbn, npadr, npcmp, nvalcp, nnugd, nperr, npnco, npnsp
    integer(kind=8) :: apnbn, apadr, apcmp, avalcp, anugd, aperr, apnco, apnsp
!
!   ENTIER REPERANT UNE MAILLE DANS UN LIGREL
!   -----------------------------------------
!
    integer(kind=8) :: grel, posm
!
!   ADRESSES DES SEGMENTS DE VALEURS DANS LE '.CELV'
!   ------------------------------------------------
!
    integer(kind=8) :: agrel, asgtm
!
!   ADRESSES LIEES AUX MODE LOCAUX
!   ------------------------------
!
    integer(kind=8) :: amodlo, mod
!
!   ADRESSE DE NUMERO DE CMP CONCERNEES PAR L' EXTRACTION
!   -----------------------------------------------------
!
    integer(kind=8) :: anumcp
!
!   DIVERS
!   ------
!
    integer(kind=8) :: i, m, nbscal, numm, nbtmai, gd, acmpgd, nbtcmp, adesgd
    integer(kind=8) :: nbval, nbn, nbco, nbsp, n1, kk
    real(kind=8) :: angl(3), pgl(3, 3), orig(3), axez(3)
    real(kind=8) :: zero, xnormz, epsi
    aster_logical :: utili
    character(len=8) :: repere
    character(len=24) :: nomjv, nomaux
!
!================= FIN DES DECLARATIONS ============================
!
!
!   CONSTRUCTIONS DES NOMS ET RECUPERATIONS DES SEGMENTS DE VALEURS
!   ---------------------------------------------------------------
!
    call jemarq()
    nchmel = nchme2
!
!
!     -- ON VERIFIE QUE LE CHAM_ELEM N'EST PAS TROP DYNAMIQUE :
    call celcel('NBVARI_CST', nchmel, 'V', '&&EXTCHE.CHAMEL1')
    nchmel = '&&EXTCHE.CHAMEL1'
    call celver(nchmel, 'NBSPT_1', 'COOL', kk)
    if (kk .eq. 1) then
        call dismoi('NOM_GD', nchmel, 'CHAMP', repk=nomgd)
        call utmess('I', 'PREPOST_36', sk=nomgd)
        call celcel('PAS_DE_SP', nchmel, 'V', '&&EXTCHE.CHAMEL2')
        nchmel = '&&EXTCHE.CHAMEL2'
    end if
!
!
    ndesc = nchmel//'.CELD'
    nvale = nchmel//'.CELV'
    ncelk = nchmel//'.CELK'
    zero = 0.0d0
    epsi = 1.0d-6
    call dismoi('NOM_MAILLA', nchmel, 'CHAMP', repk=nmaila)
!
!   TRADUCTION DES NOMS DE MAILLES EN NUMERO (SI NECESSAIRE)
!   --------------------------------------------------------
!
    if (indic .eq. 'NOMME') then
!
        do m = 1, nbm, 1
!
            nummai(m) = char8_to_int(nmaile(m))
!
        end do
!
    end if
!
    call jelira(nvale, 'TYPE', cval=type)
    ASSERT((type(1:1) .eq. 'R') .or. (type(1:1) .eq. 'C'))
    if (type(1:1) .eq. 'R') then
        call jeveuo(nvale, 'L', avale)
    else
        nomvec = 'EXTCHE.VECTEUR'
        call rvrecu(mcf, iocc, nchmel, nomvec)
        call jeveuo(nomvec, 'L', avale)
    end if
    call jeveuo(ndesc, 'L', jceld)
    call jeveuo(ncelk, 'L', acelk)
!
!   EST-CE UN REPERE "UTILISATEUR"
!   ------------------------------
!
    utili = .false.
    if (mcf(1:6) .eq. 'ACTION') then
        call getvtx(mcf, 'REPERE', iocc=iocc, scal=repere, nbret=n1)
        if (repere .eq. 'UTILISAT') then
            utili = .true.
            nomjv = '&&EXTCHE.NEW_CHAMP'
            call getvr8(mcf, 'ANGL_NAUT', iocc=iocc, nbval=3, vect=angl, &
                        nbret=n1)
            angl(1) = angl(1)*r8dgrd()
            angl(2) = angl(2)*r8dgrd()
            angl(3) = angl(3)*r8dgrd()
            call matrot(angl, pgl)
            call rvche1(nchmel, nomjv, nbm, nummai, pgl)
            call jeveuo(nomjv, 'L', avale)
        else if (repere .eq. 'CYLINDRI') then
            utili = .true.
            nomjv = '&&EXTCHE.NEW_CHAMP'
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
            call rvche2(nchmel, nomjv, nbm, nummai, orig, &
                        axez, nbnac, nnoeud)
            call jeveuo(nomjv, 'L', avale)
        end if
    end if
!
    nnugd = nssche//'.NUGD'
    npnbn = nssche//'.PNBN'
    npnco = nssche//'.PNCO'
    npnsp = nssche//'.PNSP'
    npadr = nssche//'.PADR'
    npcmp = nssche//'.PCMP'
    nvalcp = nssche//'.VALE'
    nperr = nssche//'.ERRE'
!
    nomaux = zk24(acelk+1-1)
    nligrl = nomaux(1:19)
!
    nrepe = nligrl//'.REPE'
    nnoma = nligrl//'.LGRF'
!
    call jeveuo(nnoma, 'L', anoma)
    call jeveuo(nrepe, 'L', arepe)
!
    nmaila = zk8(anoma)
!
    call wkvect(nssche//'.NOMA', 'V V K8', 1, anoma)
!
    zk8(anoma) = nmaila
!
!   RECUPERATION DES NBR TOTAUX DE MAILLE ET DE GRELS
!   -------------------------------------------------
!
    call jelira(nmaila//'.CONNEX', 'NMAXOC', nbtmai)
!
!
!   CONSERVATION DU NUMERO DE LA GRANDEUR
!   -------------------------------------
!
    call wkvect(nnugd, 'V V I', 1, anugd)
!
    gd = zi(jceld+1-1)
!
    zi(anugd) = gd
!
!   TRADUCTION DES NOMS DE CMP EN NUMEROS
!   -------------------------------------
!
    call jeveuo(jexnum('&CATA.GD.NOMCMP', gd), 'L', acmpgd)
    call jelira(jexnum('&CATA.GD.NOMCMP', gd), 'LONMAX', nbtcmp)
    call jeveuo(jexnum('&CATA.GD.DESCRIGD', gd), 'L', adesgd)
!
    call wkvect('&&EXTRCHE.NUMCP', 'V V I', nbc, anumcp)
!
    call numek8(zk8(acmpgd), ncmp, nbtcmp, nbc, zi(anumcp))
!
!   CONSTRUCTION DU POINTEUR SUR LES POSITIONS DE CMP
!   -------------------------------------------------
!
    call wkvect(npcmp, 'V V I', nbtcmp, apcmp)
    do i = 1, nbc, 1
        zi(apcmp+zi(anumcp+i-1)-1) = i
    end do
!
!   CONSTRUCTION DES POINTEURS SUR : NBR DE NDS, DE COUCHES ET DE SSPT
!   ------------------------------------------------------------------
!
    call wkvect(npnbn, 'V V I', nbtmai, apnbn)
    call wkvect(npnco, 'V V I', nbtmai, apnco)
    call wkvect(npnsp, 'V V I', nbtmai, apnsp)
!
    nbval = 0
!
    do m = 1, nbm, 1
!
        numm = nummai(m)
        grel = zi(arepe+2*(numm-1)+1-1)
        agrel = zi(jceld-1+zi(jceld-1+4+grel)+8)
        mod = zi(jceld-1+zi(jceld-1+4+grel)+2)
!
!
        call jelira(jexnum(nmaila//'.CONNEX', numm), 'LONMAX', nbn)
!
        if (mod .ne. 0) then
!
            call jeveuo(jexnum('&CATA.TE.MODELOC', mod), 'L', amodlo)
!
            nbco = zi(amodlo+4-1)
            nbco = max(nbco, nbco-10000)/nbn
            nbco = max(1, nbco)
            nbsp = max(1, zi(jceld-1+4))
!
        else
!
            nbn = 0
            nbco = 0
            nbsp = 0
!
        end if
!
        zi(apnbn+numm-1) = nbn
        zi(apnco+numm-1) = nbco
        zi(apnsp+numm-1) = nbsp
        nbval = nbval+nbn*nbc*nbco*nbsp
!
    end do
!
!   CONSTRUCTION DU POINTEUR SUR LES ADR DES SGT EXTRAITS PAR MAILLES
!   -----------------------------------------------------------------
!
    call wkvect(npadr, 'V V I', nbtmai, apadr)
!
    numm = nummai(1)
!
    zi(apadr+numm-1) = 1
!
    do m = 2, nbm, 1
!
        numm = nummai(m-1)
!
        zi(apadr+nummai(m)-1) = zi(apadr+numm-1)+nbc*zi(apnbn+numm-1)*zi(apnsp+numm-1)*zi(apn&
                                &co+numm-1)
!
    end do
!
!   CONSTRUCTION DE LA TABLE DES VALEURS DES CMP EXTRAITES
!   ------------------------------------------------------
!
    call wkvect(nvalcp, 'V V R', nbval, avalcp)
    call jeecra(nvalcp, 'DOCU', cval='CHLM')
!
    call jecrec(nperr, 'V V I', 'NU', 'DISPERSE', 'VARIABLE', &
                nbm)
!
    do m = 1, nbm, 1
!
        call jecroc(jexnum(nperr, m))
        call jeecra(jexnum(nperr, m), 'LONMAX', nbc)
        call jeveuo(jexnum(nperr, m), 'E', aperr)
!
        numm = nummai(m)
!
        grel = zi(arepe+2*(numm-1)+1-1)
        posm = zi(arepe+2*(numm-1)+2-1)
!
        agrel = zi(jceld-1+zi(jceld-1+4+grel)+8)
        mod = zi(jceld-1+zi(jceld-1+4+grel)+2)
!
        if ((mod .ne. 0) .and. (grel .ne. 0)) then
!
            call jeveuo(jexnum('&CATA.TE.MODELOC', mod), 'L', amodlo)
!
            nbscal = zi(amodlo+3-1)
            nbsp = zi(apnsp+numm-1)
            asgtm = agrel+(posm-1)*nbscal*nbsp-1
!
            call exchem(zi(amodlo), zi(anumcp), nbc, nbsp, zr(avale+asgtm), &
                        zr(avalcp+zi(apadr+numm-1)-1), zi(aperr))
!
        end if
!
    end do
!
    if (utili) call jedetr(nomjv)
    call jedetr('&&EXTRCHE.NUMCP')
    if (type(1:1) .eq. 'C') call jedetr(nomvec)
!
    call detrsd('CHAM_ELEM', '&&EXTCHE.CHAMEL1')
    call detrsd('CHAM_ELEM', '&&EXTCHE.CHAMEL2')
    call jedema()
end subroutine
