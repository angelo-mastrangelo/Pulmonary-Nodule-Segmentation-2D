# Segmentazione 2D di Noduli Polmonari in TC: Approccio Rule-Based e Analisi dei Limiti

Questo repository contiene il progetto di **Image Processing** sviluppato per il Corso di Laurea Magistrale in **Artificial Intelligence & Data Science** presso il Politecnico di Bari. Il lavoro affronta la sfida della segmentazione automatica dei noduli polmonari in immagini di Tomografia Computerizzata (TC) utilizzando un framework interamente deterministico.

## 🎯 Obiettivi del Progetto
L'obiettivo principale non è competere con le moderne architetture di deep learning, ma indagare il potenziale e i limiti strutturali di un sistema basato su regole.
* **Trasparenza Algoritmica**: Sviluppo di una pipeline senza ricorso a reti neurali o modelli black-box.
* **Interpretabilità Clinica**: Ogni classificazione è spiegabile matematicamente tramite feature geometriche e radiologiche.
* **Analisi Critica**: Misurazione scientifica del "punto di saturazione" dell'approccio 2D causato dalle ambiguità anatomiche.

## 🛠️ Metodologia e Pipeline
L'algoritmo segue una sequenza logica di moduli per ridurre progressivamente lo spazio delle soluzioni:

1. **Pre-processing**: Conversione DICOM in Unità Hounsfield (HU), windowing polmonare e normalizzazione.
2. **Isolamento Parenchima**: Estrazione dell'aria polmonare, analisi delle componenti connesse e uso dell'inviluppo convesso per recuperare i noduli juxtapleurici adesi alla pleura.
3. **Estrazione Candidati**: Generazione massiva di potenziali lesioni (solidi, GGO, calcificazioni) per massimizzare la sensibilità.
4. **Feature Extraction & Scoring**: Calcolo di vettori di proprietà per ogni candidato:
    * **Area**: Dimensione fisica della lesione.
    * **Solidity & Eccentricità**: Fondamentali per discriminare forme sferiche da strutture tubolari o vasi.
    * **Halo Score**: Misura il contrasto locale tra il candidato e il parenchima circostante.
5. **Filtro Pseudo-3D (2.5D)**: Analisi della coerenza spaziale sulle slice adiacenti per penalizzare i vasi sanguigni che, in sezione 2D, appaiono identici ai noduli.

## 📊 Risultati Sperimentali
Il sistema è stato testato sul dataset **LIDC-IDRI**, utilizzando il criterio di **Best-Match Dice** rispetto alle annotazioni di 4 radiologi esperti.

* **Sui casi rilevati**: Dice Score di **0.7859**, Precisione **0.7695** e Sensibilità **0.8596**.
* **Performance Globali**: Dice Score medio di **0.6175**, con una mediana di **0.7962**.

### Caso di Successo vs Limite 2D
L'analisi evidenzia che la causa primaria di fallimento è l'ambiguità intrinseca del dominio 2D, dove strutture vascolari mimano perfettamente la geometria tumorale.

![Confronto Successo vs Limite](03_success_vs_limit_clean.png)

## 🚀 Sviluppi Futuri
* **Transizione Volumetrica (3D)**: Analisi nativa di volumi voxel per eliminare l'ambiguità del taglio trasversale.
* **Integrazione Radiomica**: Estrazione di feature di texture avanzate per distinguere patologie a parità di forma.
* **Architettura Ibrida**: Uso del Deep Learning per la fase di detection, mantenendo l'approccio deterministico per l'estrazione dei contorni.
