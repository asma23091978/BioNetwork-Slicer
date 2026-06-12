# Pulizia fatta sui file

## File Prolog con le reazioni

Il file che hai fornito conteneva all'inizio:

```prolog
/* To generate the .dot file */
%:-["<PERCORSO_LOCALE_PRIVATO>/params_fixed_context_30-01_nul.pl"].
/* For the slicing */
:- ["params.pl"].
```

Per la versione pubblica ho rimosso le righe di commento e il percorso privato locale.

L'inizio pulito del file ora è:

```prolog
/*
Biological reaction model for hgd0_04-02.
...
*/

:- ["params.pl"].

myreactions([
```

## Doppioni rimossi

Non ho inserito copie duplicate dello stesso file Prolog. È presente una sola versione pulita:

```text
prolog/hgd0_04_02_reactions.pl
```

## File non inclusi perché mancanti

Non ho creato artificialmente questi file, perché devono essere quelli veri del progetto:

```text
BioReSolvePositive.pl
output/hgd0_04-02.dot
```
