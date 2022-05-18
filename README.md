# Decomp

[![CC BY-NC-SA 4.0][cc-shield]][cc]
[![stability-beta](https://img.shields.io/badge/stability-beta-33bbff.svg)](https://github.com/mkenney/software-guides/blob/master/STABILITY-BADGES.md#beta)


Decomp project for university course [Znanstveno računanje 1](https://www.pmf.unizg.hr/math/predmet/znarac1), academic year 2019/2020, University of Zagreb, Faculty of Science, Department of Mathematics.


## Kratki opis 
Program za numeričko rješavanje Poissonove jednadžbe uz Dirichletove rubne uvjete, na nepravilnoj domeni u obliku slova 'L'. Druge derivacije aproksimiramo centralnim diferencijama, a na domeni 'L' uvodimo ekvidistantnu mrežu i dekompozicjiu na manje pravokutne poddomene koje se ne preklapaju. Dobiveni sustav linearnih jednadžbi u matrici organiziramo na nain da prvo poredamo čvorove unutar poddomene, a zatim čvorove koji se nalaze na rubovima. Takvim postupkom dolazimo do *algoritma blok-Gaussovih eliminacija*, a detaljnije o algoritmu i rezultatima može se naći u [tekstu](https://github.com/sopetra/decomp/blob/main/Iterativne%20metode%20za%20sustave%20-%20dekompozicija%20domene.pdf).

## Pokretanje programa
Potrebne biblioteke:
1. `LAPACK` (Linear Algebra Package). Dokumentacija je dostupna na [linku](http://www.netlib.org/lapack/).
2. `BLAS` (Basic Linear ALgebra Subprograms). Dokumentacija je dostupna na [linku](http://www.netlib.org/blas/).
3. `f2c.c` prebacuje Fortran kod u C. 
Testiranje algoritma se postiže pokretanjem `domena.c`.


## Licence
  
 [Decomp](https://github.com/sopetra/decomp) © 2021 by [Petra Sočo](https://github.com/sopetra) is licensed under [Attribution-NonCommercial-ShareAlike 4.0 International][cc].

[![CC 4.0][cc-image]][cc]


[cc]: https://creativecommons.org/licenses/by-nc-sa/4.0/?ref=chooser-v1
[cc-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-shield]: https://img.shields.io/badge/License-CC%20BY--SA%204.0-lightgrey.svg


License can be found under [License](LICENSE).
