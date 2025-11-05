# Checking pointers

`PARACHECKNAME` is a pointer that is not associated with a target. 
It is used in the `IRCHOR` routine in the `bibfor/prepost/irchor.F90` file.


```
forrtl: severe (408): fort: (7): Attempt to use pointer PARACHECKNAME when it is not associated with a target

Image              PC                Routine            Line        Source             
bibfor.dll         00007FF8D235B18A  IRCHOR                    314  irchor.F90
bibfor.dll         00007FF8D23DDA56  IRMFAC                    291  irmfac.F90
bibfor.dll         00007FF8D1F76B92  OP0039                    212  op0039.F90
bibfor.dll         00007FF8D277C48E  EX0000                    247  ex0000.F90
bibfor.dll         00007FF8D277CB4B  EXECOP                     71  execop.F90
bibfor.dll         00007FF8D277CE9A  EXPASS                     41  expass.F90
bibcxx.dll         00007FF8DA3CB40A  Unknown               Unknown  Unknown
bibcxx.dll         00007FF8DA3D0082  Unknown               Unknown  Unknown
bibcxx.dll         00007FF8DA3CFD13  Unknown               Unknown  Unknown
bibcxx.dll         00007FF8DA3CFB3A  Unknown               Unknown  Unknown
bibcxx.dll         00007FF8DA3CF9DA  Unknown               Unknown  Unknown
bibcxx.dll         00007FF8D9D1ED25  Unknown               Unknown  Unknown
python311.dll      00007FF8DC13E459  Unknown               Unknown  Unknown
python311.dll      00007FF8DC0F6B4A  Unknown               Unknown  Unknown
python311.dll      00007FF8DC201B42  Unknown               Unknown  Unknown
python311.dll      00007FF8DC1FCF29  Unknown               Unknown  Unknown
python311.dll      00007FF8DC2000EE  Unknown               Unknown  Unknown
python311.dll      00007FF8DC0F6FFD  Unknown               Unknown  Unknown
python311.dll      00007FF8DC0F9281  Unknown               Unknown  Unknown
python311.dll      00007FF8DC0F97AD  Unknown               Unknown  Unknown
python311.dll      00007FF8DC0F6C59  Unknown               Unknown  Unknown
python311.dll      00007FF8DC0F6E38  Unknown               Unknown  Unknown
python311.dll      00007FF8DC201CED  Unknown               Unknown  Unknown
python311.dll      00007FF8DC1FDF6D  Unknown               Unknown  Unknown
python311.dll      00007FF8DC1F81F1  Unknown               Unknown  Unknown
python311.dll      00007FF8DC1F2D79  Unknown               Unknown  Unknown
python311.dll      00007FF8DC1F072F  Unknown               Unknown  Unknown
python311.dll      00007FF8DC13E048  Unknown               Unknown  Unknown
python311.dll      00007FF8DC0F6729  Unknown               Unknown  Unknown
python311.dll      00007FF8DC201B42  Unknown               Unknown  Unknown
python311.dll      00007FF8DC1FCF29  Unknown               Unknown  Unknown
python311.dll      00007FF8DC1F81F1  Unknown               Unknown  Unknown
python311.dll      00007FF8DC27A6FE  Unknown               Unknown  Unknown
python311.dll      00007FF8DC27A7C8  Unknown               Unknown  Unknown
python311.dll      00007FF8DC27A398  Unknown               Unknown  Unknown
python311.dll      00007FF8DC2772E5  Unknown               Unknown  Unknown
python311.dll      00007FF8DC073ADA  Unknown               Unknown  Unknown
python311.dll      00007FF8DC074552  Unknown               Unknown  Unknown
python311.dll      00007FF8DC07493A  Unknown               Unknown  Unknown
python.exe         00007FF64FB31490  Unknown               Unknown  Unknown
KERNEL32.DLL       00007FF9E2D0257D  Unknown               Unknown  Unknown
ntdll.dll          00007FF9E416AF28  Unknown               Unknown  Unknown

Process finished with exit code 408
```



