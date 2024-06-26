# Med file write issue


```
medfwrap.dll!MFIFVOP(const char * name, const unsigned int bidon, const int * access, const __int64 * const major, const __int64 * const minor, const __int64 * const release, const __int64 * const len) Line 89
	at D:\a\build\libmed_1719319141648\work\src\cfi\filecf.c(89)
medfwrap.dll!MFIFVOP.t93p.t3v.t94p.t95p.t96p.t97p.t98p() Line 0
	at D:\a\build\libmed_1719319141648\work\src\fi\medfile.f(0)
medfwrap.dll!MFIVOP() Line 31
	at D:\a\build\libmed_1719319141648\work\src\fi\medfile.f(31)
bibfor.dll!AS_MFIVOP() Line 44
	at C:\Work\code\code-aster-src\bibfor\echange\as_mfivop.F90(44)
bibfor.dll!AS_MED_MODULE::AS_MED_OPEN() Line 78
	at C:\Work\code\code-aster-src\bibfor\echange\as_med_module.F90(78)
bibfor.dll!IRMHDF() Line 213
	at C:\Work\code\code-aster-src\bibfor\prepost\irmhdf.F90(213)
bibfor.dll!IRMAIL() Line 251
	at C:\Work\code\code-aster-src\bibfor\prepost\irmail.F90(251)
bibfor.dll!IRCAME() Line 214
	at C:\Work\code\code-aster-src\bibfor\prepost\ircame.F90(214)
bibfor.dll!IRCNME() Line 122
	at C:\Work\code\code-aster-src\bibfor\prepost\ircnme.F90(122)
bibfor.dll!IRCHME() Line 261
	at C:\Work\code\code-aster-src\bibfor\prepost\irchme.F90(261)
bibfor.dll!IREMED() Line 291
	at C:\Work\code\code-aster-src\bibfor\prepost\iremed.F90(291)
bibfor.dll!IRMFAC() Line 353
	at C:\Work\code\code-aster-src\bibfor\prepost\irmfac.F90(353)
bibfor.dll!OP0039() Line 213
	at C:\Work\code\code-aster-src\bibfor\op\op0039.F90(213)
bibfor.dll!EX0000() Line 247
	at C:\Work\code\code-aster-src\bibfor\supervis\ex0000.F90(247)

```