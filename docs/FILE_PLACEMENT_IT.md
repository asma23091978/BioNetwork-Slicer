# Dove mettere i file

Questa è la struttura pulita per GitHub.

## File già inseriti

```text
slicing_script.py
requirements.txt
prolog/hgd0_04_02_reactions.pl
scripts/repo_check.py
README.md
.gitignore
```

Il file `prolog/hgd0_04_02_reactions.pl` è il file con `myreactions([...])`. È stato pulito: il percorso locale privato è stato rimosso.

## File ancora da aggiungere per eseguire il progetto

### 1. Motore Prolog principale

Deve stare qui:

```text
BioReSolvePositive.pl
```

cioè nella cartella principale, accanto a `slicing_script.py`.

Questo file deve contenere il predicato `main_do(...)`, perché il Python esegue query come:

```prolog
main_do(ordslice, EKs, ListReactNumbR, ListEpos, ListCS).
main_do(negslice, EKs, ListReactNumbR, ListEpos, ListEneg, ListCS).
```

### 2. File grafo DOT

Deve stare qui:

```text
output/hgd0_04-02.dot
```

Il nome deve essere esattamente `hgd0_04-02.dot`, perché nel Python c'è:

```python
dot_name = "hgd0_04-02"
G = nx.DiGraph(nx.nx_pydot.read_dot("output/" + dot_name + ".dot"))
```

## File da non caricare su GitHub

```text
params.pl
Results/
tmp-slice*.txt
tmp-legenda*.txt
tmp-slicingrun*.txt
table_results*.csv
```

Questi file vengono generati automaticamente durante l'esecuzione.
