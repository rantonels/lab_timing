# Esperienza "Timing"

Questo è il programma di analisi per l'esperienza "Timing".

##Guida all'analisi

###Prerequisiti

Installare git.

```bash
sudo apt-get install git
```

Scaricare e compilare root.


###Come scaricare

per scaricare la prima volta, da terminale:

```bash
cd ~
git clone https://github.com/rantonels/lab_timing.git
```

nella vostra home dovreste trovarvi una cartella "timing" oppure "lab_timing" con il progetto dentro. Quando volete aggiornare il progetto all'ultima versione, **NON** ripetete la procedura sopra. Piuttosto, fate ```cd``` nella cartella del progetto, ed eseguite

```bash
git pull
```

###Compilazione

I programmi in C++ per l'analisi vanno ricompilati. (Non si possono distribuire direttamente i binari, per motivi ovvi). Entrate nella cartella del progetto ed eseguite:

```bash
make
```

###Analisi Dati

**NB: lo script di analisi non è ancora pronto**


##Struttura progetto

Nella cartella del progetto trovate le seguenti cartelle:

###bin/

I binari c++ compilati. Questa cartella è vuota e viene popolata al momento della compilazione.

*OH NO! DOVE È FINITO IL MIO a.out?*

è qui dentro, zio.

###build/

Oggetti temporanei della compilazione.

###src/

sorgenti ed header c++, script python, etc.

###data/

Cartella dedicata esclusivamente ai **dati grezzi**. Tutti e soli i **dati grezzi** vanno inseriti qui. Questo vuol dire che eventuali dati convertiti, istogrammi, prodotti dell'analisi, grafici, latex, etc. vanno messi da qualche altra parte.

###tmp/

Cartella per i file temporanei dell'analisi. Ad esempio, file prodotti dal programma numero 3 che servono al programma numero 5. Qui possono anche starci eventuali file di debug o test o log prodotti dai programmi. Tenete conto che un file messo in questa cartella non è al sicuro, ed è a rischio di essere cancellato.

###out/

Risultati definitivi dell'analisi dati. Questa cartella può ad esempio contenere gli include definitivi per il .tex.



##Guida ai singoli programmi

###comptonfit

Esegue il fit iterativo del profilo Compton su di un file istogramma.

```bash
bin/comptonfit XXX
```

carica i dati da XXX, esegui il fit e salva in tmp/XXX.cfit

```bash
bin/comptonfit test
```

esegui un test con dati generati (che salva in tmp/randomcompton) e stampa a terminale i risultati

###datafile

datafile.h è un header per caricare file di dati testuali in array c++. Questo header è preferibile perché è in grado di rimuovere righe vuote e righe di commenti (prefissate con il cancelletto). È usato da comptonfit.

###root2csv

Apre un file root a quattro canali, genera gli istogrammi e li esporta in file di dati testuali.

Per esaminare gli istogrammi di un file:

```bash
bin/root2csv XXX
```

carica il file di root XXX e stampa a terminale gli isto dei quattro canali (non salva niente)

```bash
bin/root2csv XXX YYY
```

carica XXX, genera gli isto e li salva nei file YYY_h0, YYY_h1, eccetera. Non stampa niente a terminale.

**Nota Bene:** Come specificato sopra, è bene che XXX sia nella cartella data/ e YYY nella cartella tmp/. I dati grezzi e qualsiasi rielaborazione di essi devono rimanere separati. È importante poter cancellare tutti i prodotti dell'analisi senza paura di intaccare i dati di laboratorio. Inoltre, in caso di corruzione/modifica dei dati, è molto facile ripristinarli dalla cartella Dropbox.
