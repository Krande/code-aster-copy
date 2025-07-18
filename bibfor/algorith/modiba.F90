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

subroutine modiba(nomres, basemo, basefl, numvit, newres, &
                  itypfl, imasse, nuor, nbnuor, numo, &
                  nbmfl)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8depi.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/freqom.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/pteddl.h"
#include "asterfort/refdcp.h"
#include "asterfort/rscrsd.h"
#include "asterfort/rsexch.h"
#include "asterfort/titre.h"
#include "asterfort/vpnor2.h"
#include "asterfort/vpnorm.h"
#include "asterfort/vprecu.h"
#include "asterfort/vpstor.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
    integer(kind=8) :: numvit, itypfl, imasse
    integer(kind=8) :: nbnuor, nuor(*), nbmfl, numo(*)
    character(len=8) :: nomres, basemo
    character(len=19) :: basefl
    aster_logical :: newres
!     MODIFICATION D'UNE BASE MODALE DE TYPE MODE_MECA APRES UN CALCUL
!     DE COUPLAGE FLUIDE-STRUCTURE
!     APPELANT : OP0149, MODI_BASE_MODALE
! ----------------------------------------------------------------------
! IN  : NOMRES : NOM DU CONCEPT MODE_MECA DE SORTIE
! IN  : BASEMO : NOM DU CONCEPT MODE_MECA D'ENTREE
! IN  : BASEFL : NOM DU CONCEPT MELASFLU
! IN  : NUMVIT : NUMERO DE LA VITESSE DU FLUIDE POUR LAQUELLE ON RETIENT
!                LES NOUVELLES CARACTERISTIQUES MODALES
! IN  : NEWRES : INDICATEUR BOOLEEN
!       NEWRES = .TRUE.  CREATION D'UN NOUVEAU CONCEPT EN SORTIE
!                        => NOMRES <> BASEMO
!       NEWRES = .FALSE. MODIFICATION DU CONCEPT D'ENTREE
!                        => NOMRES =  BASEMO
! IN  : ITYPFL : INDICE CARACTERISTIQUE DE LA CONFIGURATION ETUDIEE
!       ITYPFL = 1  FAISCEAU_TRANS   ITYPFL = 2  GRAPPE
!       ITYPFL = 3  FAISCEAU_AXIAL   ITYPFL = 4  COQUE_COAX
! IN  : IMASSE : INDICE CARACTERISTIQUE LORSQUE ITYPFL = 4
! IN  : NUOR   : LISTE DES NUMEROS D'ORDRE DES MODES RETENUS POUR
!                RECONSTRUIRE LA BASE MODALE
! IN  : NBNUOR : NOMBRE DE MODES RETENUS POUR RECONSTRUIRE LA BASE
! IN  : NUMO   : LISTE DES NUMEROS D'ORDRE DES MODES PERTURBES PAR LE
!                COUPLAGE FLUIDE-STRUCTURE
! IN  : NBMFL  : NOMBRE DE MODES PERTURBES PAR LE COUPLAGE
!
!
    integer(kind=8) :: iddl(6), neq, nbmode, j, i
    integer(kind=8) :: lmod, iret, ivit, numod, imas
    integer(kind=8) :: ifac, ifre, ieq, k, icm, iprec
    integer(kind=8) :: lmat(2), lddl, lvali, lvalr, lvalk
    integer(kind=8) :: npari, nparr, npark
    integer(kind=8) :: nbpari, nbparr, nbpark, nbpara
    parameter(nbpari=1, nbparr=15, nbpark=1, nbpara=17)
    real(kind=8) :: frequ, amort, omeg2, masg, rigg
    real(kind=8) :: factx, facty, factz, depi, xmastr(3)
    character(len=1) :: typmod
    character(len=19) :: numeq
    character(len=16) :: norm
    character(len=19) :: nomcha
    character(len=24) :: chamfl, kvec, nopara(nbpara)
    character(len=24) :: kvali, kvalr, kvalk
    aster_logical :: lnorm
    real(kind=8), pointer :: coef_mode(:) => null()
    real(kind=8), pointer :: fact(:) => null()
    integer(kind=8), pointer :: deeq(:) => null()
    real(kind=8), pointer :: vmasg(:) => null()
    real(kind=8), pointer :: freq(:) => null()
    real(kind=8), pointer :: vale(:) => null()
!
    data iddl/1, 2, 3, 4, 5, 6/
    data nopara/&
     &  'NUME_MODE', 'NORME',&
     &  'FREQ',&
     &  'OMEGA2', 'AMOR_REDUIT',&
     &  'MASS_GENE', 'RIGI_GENE', 'AMOR_GENE',&
     &  'MASS_EFFE_DX', 'MASS_EFFE_DY', 'MASS_EFFE_DZ',&
     &  'FACT_PARTICI_DX', 'FACT_PARTICI_DY', 'FACT_PARTICI_DZ',&
     &  'MASS_EFFE_UN_DX', 'MASS_EFFE_UN_DY', 'MASS_EFFE_UN_DZ'/
!
!     ------------------------------------------------------------------
!
    call jemarq()
    depi = r8depi()
    lmat(1) = 0
    lmat(2) = 0
    lnorm = .false.
!   --- xmastr = [0,0,0] pour forcer les masses unitaires nulles sur les 3 directions
    xmastr = [0.d0, 0.d0, 0.d0]
!
!     --- CREATION DU CONCEPT MODE_MECA DE SORTIE LE CAS ECHEANT ---
!
    if (newres) then
        call rscrsd('G', nomres, 'MODE_MECA', nbnuor)
        call refdcp(basemo, nomres)
    end if
!
!     --- PARAMETRES SOUS ECOULEMENT ---
!
    call jeveuo(basefl//'.FREQ', 'L', vr=freq)
    call jeveuo(basefl//'.MASG', 'L', vr=vmasg)
    call jeveuo(basefl//'.FACT', 'L', vr=fact)
!
    typmod = 'R'
    kvec = '&&MODIBA.VECT'
    kvali = '&&MODIBA.PARA_I'
    kvalr = '&&MODIBA.PARA_R'
    kvalk = '&&MODIBA.PARA_K'
    call vprecu(basemo, 'DEPL', nbnuor, nuor, kvec, &
                nbpara, nopara(1), kvali, kvalr, kvalk, &
                neq, nbmode, typmod, npari, nparr, &
                npark)
    ASSERT(npari .eq. nbpari)
    ASSERT(nparr .eq. nbparr)
    ASSERT(npark .eq. nbpark)
    call jeveuo(kvec, 'E', lmod)
    call jeveuo(kvali, 'E', lvali)
    call jeveuo(kvalr, 'E', lvalr)
    call jeveuo(kvalk, 'E', lvalk)
!
!     --- ON RECUPERE UN NUME_DDL ---
!
    call rsexch('F', basemo, 'DEPL', nuor(1), nomcha, &
                iret)
    call dismoi('NUME_EQUA', nomcha, 'CHAM_NO', repk=numeq)
    call jeveuo(numeq//'.DEEQ', 'L', vi=deeq)
!
!     --- CAS DU COUPLAGE ---
!
    ivit = 1
    if (itypfl .eq. 3) ivit = numvit
    chamfl(1:13) = basefl(1:8)//'.C01.'
!
    do j = 1, nbmfl
        do i = 1, nbnuor
            numod = nuor(i)
!
            if (numo(j) .eq. numod) then
                imas = nbmfl*(ivit-1)+j
                ifac = nbmfl*(ivit-1)+3*(j-1)
                ifre = 2*nbmfl*(numvit-1)+2*(j-1)
                frequ = freq(ifre+1)
                amort = freq(1+ifre+1)
                omeg2 = (depi*frequ)**2
                masg = vmasg(imas)
                rigg = omeg2*masg
                factx = fact(ifac+1)
                facty = fact(ifac+2)
                factz = fact(ifac+3)
                if (amort .le. 0.d0) amort = 1.d-06
!
!           --- FREQUENCE ---
                zr(lvalr+i-1) = freqom(omeg2)
!           --- OMEGA2 ---
                zr(lvalr+nbnuor+i-1) = omeg2
!           --- AMOR_REDUIT ---
                zr(lvalr+nbnuor*2+i-1) = amort
!           --- MASS_GENE , RIGI_GENE ---
                zr(lvalr+nbnuor*3+i-1) = masg
                zr(lvalr+nbnuor*4+i-1) = rigg
!           --- MASS_EFFE_D... ---
                zr(lvalr+nbnuor*6+i-1) = factx*factx/masg
                zr(lvalr+nbnuor*7+i-1) = facty*facty/masg
                zr(lvalr+nbnuor*8+i-1) = factz*factz/masg
!           --- FACT_PARTICI_D... ---
                zr(lvalr+nbnuor*9+i-1) = factx/masg
                zr(lvalr+nbnuor*10+i-1) = facty/masg
                zr(lvalr+nbnuor*11+i-1) = factz/masg
!
                if (itypfl .eq. 3 .or. (itypfl .eq. 4 .and. imasse .ne. 0)) then
                    lnorm = .true.
                    zk24(lvalk+i-1) = 'SANS_CMP: LAGR'
                    write (chamfl(14:19), '(2I3.3)') numod, numvit
                    call jeveuo(chamfl(1:19)//'.VALE', 'L', vr=vale)
                    icm = 0
                    do ieq = 1, neq
                        do k = 1, 6
                            if (deeq(1+(2*ieq)-1) .eq. iddl(k)) then
                                icm = icm+1
                                zr(lmod+neq*(i-1)+ieq-1) = vale(icm)
                                goto 30
                            end if
                        end do
30                      continue
                    end do
                end if
            end if
        end do
    end do
!
!     --- ON NORMALISE 'SANS_CMP: LAGR'
!
    iprec = 0
    if (lnorm) then
        norm = 'AVEC_CMP'
        call wkvect('&&MODIBA.POSITION.DDL', 'V V I', neq, lddl)
        call pteddl('CHAM_NO', nomcha, 1, 'LAGR    ', neq, &
                    list_equa=zi(lddl))
        do ieq = 0, neq-1
            zi(lddl+ieq) = 1-zi(lddl+ieq)
        end do
        AS_ALLOCATE(vr=coef_mode, size=nbmode)
!        --- ON NORMALISE LES DEFORMEES
        call vpnorm(norm, 'OUI', lmat(1), neq, nbmode, &
                    zi(lddl), zr(lmod), zr(lvalr), xmastr, 0, &
                    0, coef_mode)
!        --- ON STOCKE LES DEFORMEES
        call vpstor(-1, typmod, nomres, nbnuor, neq, &
                    zr(lmod), zc(1), nbnuor, nbpari, nbparr, &
                    nbpark, nopara, '    ', zi(lvali), zr(lvalr), &
                    zk24(lvalk), iprec)
!        --- ON NORMALISE LES AUTRES CHAMPS
        call vpnor2(nomres, nbmode, nuor, coef_mode)
        AS_DEALLOCATE(vr=coef_mode)
    else
!        --- ON STOCKE LES DEFORMEES
        call vpstor(-1, typmod, nomres, nbnuor, neq, &
                    zr(lmod), zc(1), nbnuor, nbpari, nbparr, &
                    nbpark, nopara, '    ', zi(lvali), zr(lvalr), &
                    zk24(lvalk), iprec)
    end if
!
!     --- TITRE ASSOCIE AU CONCEPT ---
!
    call titre()
!
! --- MENAGE
!
    call jedetr('&&MODIBA.VECT')
    call jedetr('&&MODIBA.PARA_I')
    call jedetr('&&MODIBA.PARA_R')
    call jedetr('&&MODIBA.PARA_K')
    call jedetr('&&MODIBA.POSITION.DDL')
!
    call jedema()
end subroutine
