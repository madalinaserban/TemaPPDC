# TemaPPDC Serban Elena Madalina 10LF203

Pentru tema de laborator am realizat 4 sortari : Shell Sort,Bitonic sort,Bucket sort si Odd Even sort implementate MPI in limbajul c.
Am generat un sir de 1000000 de elemnte de tip double pentru a face sortarile.

##Shell Sort
MPI_Shell_sort este o implementare a algoritmului Shell Sort folosind MPI pentru a realiza sortarea paralela pe mai multe procesoare si primeste ca parametrii : 
a: tabloul de sortat
n: dimensiunea tabloului
root:  numarul de elemente din tablou
comm: comunicatorul MPI

Algoritmul functioneaza prin impartirea vectorului initial in sub-secvente mai mici si sortarea fiecareia. Pe masura ce sub-secventele sunt sortate, distanta intre elementele comparate se reduce treptat, pana cand intregul vector este sortat.

Functia imparte array-ul de sortat in bucati egale ai le distribuie catre toate procesele folosind MPI_Scatter. Apoi, procesele sorteaza bucatile lor local folosind algoritmul Shell Sort si se aduna folosind MPI_Allreduce. La final, datele sortate sunt adunate inapoi in array-ul initial folosind MPI_Gather.
Se returneaza MPI_SUCCESS daca totul a decurs fara erori si folosim MPI_Is_Sorted pentru a verifica daca vectorul final este sortat.


##Bucket Sort

MPI_Sort_bucket este o implementare paralela a algoritmului Bucket Sort si primeste ca parametrii : 
n: numarul de elemente din tablou.
a: tabloul de sortat
max: valoarea maxima posibila din tablou.
root: radacina procesului comunicarii MPI
comm: comunicatorul MPI folosit pentru a coordona operatiunile intre procese.

Alocam spatiul necesar pentru bucket, se face broadcast-ul array-ului a si se scaneaza array-ul pentru a colecta elementele din bucket-ul corespunzator rank-ului procesului curent. Apoi, se sorteaza bucket-ul folosind algoritmul Merge Sort. 
La final, se aduna elementele din fiecare bucket si se stocheaza in array-ul a.
Se returneaza MPI_SUCCESS daca totul a decurs fara erori si folosim MPI_Is_Sorted pentru a verifica daca vectorul final este sortat.

##Odd Even sort
MPI_Sort_oddEven este o implementare paralela a algoritmului de sortare odd-even folosind MPI si primeste ca parametrii :
n:  numarul de elemente din tablou
a: tabloul de sortat
root: procesul radacina al comunicarii MPI
comm - comunicatorul MPI folosit pentru comunicarea intre procese

Fiecare proces primeste un array prin MPI_Scatter si il sorteaza folosind mergeSort. Apoi, se executa un numar de iteratii odd-even pentru a asigura ordinea crescatoare a elementelor. In fiecare iteratie, procesele comunica intre ele prin MPI_Exchange. Dupa finalizarea tuturor iteratiilor, elementele sunt adunate inapoi in array-ul initial folosind MPI_Gather.
Se returneaza MPI_SUCCESS daca totul a decurs fara erori si folosim MPI_Is_Sorted pentru a verifica daca vectorul final este sortat.

##Bitonic Sort

MPI_Bitonic_sort primeste ca parametrii:

n: numarul de elemente din array-ul de sortat
a: array-ul de sortat
root: procesul radacina al comunicarii MPI
comm: comunicatorul MPI folosit pentru a coordona operatiunile intre procese

Array-ul de sortat este impartit in bucati egale, care sunt trimise fiecarui proces folosind functia MPI_Scatter. Fiecare proces sorteaza bucata primita folosind functia bitonicSort.

Functia bitonicSort imparte bucata in doua jumatati si apoi sorteaza fiecare jumatate, cele doua sunt unificate prin apelarea functiei bitonicMerge.

BitonicMerge compara elementele din bucata data si le schimba pozitiile daca este necesar. Apoi, bucata este impartita in doua jumatati, fiecare dintre ele fiind sortata. Aceste doua jumatati sunt apoi unificate prin apelarea functiei bitonicMerge.
Dupa finalizarea tuturor iteratiilor, elementele sunt adunate inapoi in array-ul initial folosind MPI_Gather.
Se returneaza MPI_SUCCESS daca totul a decurs fara erori si folosim MPI_Is_Sorted pentru a verifica daca vectorul final este sortat.


Excel :https://docs.google.com/spreadsheets/d/1z9kRgRTT7xQVowTdB47o5Lofl511KNR4DxyNOs_2Tp8/edit?usp=sharing


 
