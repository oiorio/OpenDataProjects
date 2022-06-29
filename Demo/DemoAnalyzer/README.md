Nell'attuale cartella ci sono: 
 A) Tre file cruciali:
 Demo/DemoAnalyzer/BuildFile.xml
 Demo/DemoAnalyzer/src/AnalyzeTT.cc
 Demo/DemoAnalyzer/analyzett_cfg.py
 B) numerosi file ancillari chiamati mclist*.py

Cosa fare:
0) Fai il setup come per lanciare demoanalyzer_cfg.py
1) Copia i file nella stessa cartella del tuo DemoAnalyzer
2) Usa il comando "scram b " per compilare
3) Fai "cmsenv" di nuovo

Se non ti ha dato errori nelle fasi 1,2,3, puoi testare il comando:

edmPluginDump | grep AnalyzerTT

che dovrebbe darti come output:

AnalyzerTT
AnalyzerTT

A questo punto puoi provare a girare

cmsRun analyzett_cfg.py

che dovrebbe produrti un file di 550 eventi con il contenuto di muoni, elettroni, jet etc
