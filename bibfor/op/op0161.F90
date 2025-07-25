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
subroutine op0161()
    implicit none
!     COMBINAISON FOURIER
!     ------------------------------------------------------------------
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterc/r8dgrd.h"
#include "asterfort/assert.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/refode.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsagsd.h"
#include "asterfort/rscrsd.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnoch.h"
#include "asterfort/rsorac.h"
#include "asterfort/wkvect.h"
    complex(kind=8) :: cbid
    character(len=8) :: k8b, resu, resuin, modele, mate, carele
    character(len=16) :: concep, nomcmd, nsymb
    character(len=24) :: nomch
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: ian, ibid, ich, ior, iordr, iret, jangl
    integer(kind=8) :: jcara, jcham, jcoe, jjan, jmat, jmod, jnch
    integer(kind=8) :: jnha, jnmo, jordr, jpara, jtch, jtmo, k, tord(1)
    integer(kind=8) :: n1, n2, n3, nbangl, nbcham, nbordr, nbtrou
!
    real(kind=8) :: angle, rbid
!-----------------------------------------------------------------------
    call jemarq()
    call infmaj()
!
    call getres(resu, concep, nomcmd)
!
    call getvid(' ', 'RESULTAT', scal=resuin, nbret=n1)
    call rsorac(resuin, 'LONUTI', ibid, rbid, k8b, &
                cbid, rbid, k8b, tord, 1, &
                nbtrou)
    nbordr = tord(1)
    call wkvect('&&OP0161.NUME_ORDRE', 'V V I', nbordr, jordr)
    call rsorac(resuin, 'TOUT_ORDRE', 0, rbid, k8b, &
                cbid, rbid, k8b, zi(jordr), nbordr, &
                ibid)
!
    call getvtx(' ', 'NOM_CHAM', nbval=0, nbret=n2)
    nbcham = -n2
    call wkvect('&&OP0161.CHAMP', 'V V K16', nbcham, jcham)
    call getvtx(' ', 'NOM_CHAM', nbval=nbcham, vect=zk16(jcham), nbret=n2)
!
    call getvr8(' ', 'ANGLE', nbval=0, nbret=n3)
    nbangl = -n3
    call wkvect('&&OP0161.ANGLE', 'V V R', nbangl, jangl)
    call getvr8(' ', 'ANGLE', nbval=nbangl, vect=zr(jangl), nbret=n3)
!
    call rscrsd('G', resu, concep, nbangl)
!
    call wkvect('&&OP0161.NOM_CHAMP', 'V V K24', nbordr, jnch)
    call wkvect('&&OP0161.TYP_CHAMP', 'V V K8', nbordr, jtch)
    call wkvect('&&OP0161.NUM_HARMO', 'V V I', nbordr, jnha)
    call wkvect('&&OP0161.COEFFICIE', 'V V R', nbordr, jcoe)
!
    do ich = 1, nbcham
        nsymb = zk16(jcham+ich-1)
        k = 0
        do ior = 0, nbordr-1
            iordr = zi(jordr+ior)
            call rsexch(' ', resuin, nsymb, iordr, nomch, &
                        iret)
            if (iret .eq. 0) then
                k = k+1
                zk24(jnch+k-1) = nomch
                call rsadpa(resuin, 'L', 1, 'TYPE_MODE', iordr, &
                            0, sjv=jtmo, styp=k8b)
                zk8(jtch+k-1) = zk8(jtmo)
                call rsadpa(resuin, 'L', 1, 'NUME_MODE', iordr, &
                            0, sjv=jnmo, styp=k8b)
                zi(jnha+k-1) = zi(jnmo)
                zr(jcoe+k-1) = 1.d0
                call rsadpa(resuin, 'L', 1, 'MODELE', iordr, &
                            0, sjv=jmod, styp=k8b)
                modele = zk8(jmod)
                call rsadpa(resuin, 'L', 1, 'CHAMPMAT', iordr, &
                            0, sjv=jmat, styp=k8b)
                mate = zk8(jmat)
                call rsadpa(resuin, 'L', 1, 'CARAELEM', iordr, &
                            0, sjv=jcara, styp=k8b)
                carele = zk8(jcara)
            end if
        end do
!
        do ian = 1, nbangl
!
!     STOCKAGE DU NOM DU MODELE
!     -------------------------
            call rsadpa(resu, 'E', 1, 'MODELE', ian, &
                        0, sjv=jpara, styp=k8b)
            zk8(jpara) = modele
!
!     STOCKAGE DU NOM DU CHAMP MATERIAU
!     ---------------------------------
            call rsadpa(resu, 'E', 1, 'CHAMPMAT', ian, &
                        0, sjv=jpara, styp=k8b)
            zk8(jpara) = mate
!
!     STOCKAGE DU NOM DE LA CARACTERISTIQUE ELEMENTAIRE
!     -------------------------------------------------
            call rsadpa(resu, 'E', 1, 'CARAELEM', ian, &
                        0, sjv=jpara, styp=k8b)
            zk8(jpara) = carele
!
        end do
!
        if (k .ne. 0) then
            do ian = 1, nbangl
                call rsexch(' ', resu, nsymb, ian, nomch, &
                            iret)
                if (iret .eq. 110) then
                    call rsagsd(resu, 0)
                    call rsexch(' ', resu, nsymb, ian, nomch, &
                                iret)
                else if (iret .eq. 100) then
                else
                    ASSERT(.false.)
                end if
                angle = zr(jangl+ian-1)*r8dgrd()
                call refode(k, angle, zk24(jnch), zi(jnha), zk8(jtch), &
                            zr(jcoe), 'G', nomch)
                call rsnoch(resu, nsymb, ian)
                call rsadpa(resu, 'E', 1, 'ANGLE', ian, &
                            0, sjv=jjan, styp=k8b)
                zr(jjan) = zr(jangl+ian-1)
            end do
        end if
    end do
!
    call jedema()
end subroutine
