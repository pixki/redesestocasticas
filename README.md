# redesestocasticas
Scripts y programas desarrollados para la materia de redes estocásticas

S-ALOHA.py   
==========

Genera matrices de transición de acuerdo al modelo de ALOHA ranurado de población finita sin colas. La salida puede ser "pipeada" a otro script/proceso o guardada como archivo CSV (Excel, Calc, ...)


dtmc-csim.py
============

Calcula la distribución estacionaria de una cadena de Markov que modele el protocolo  S-ALOHA, dada su matriz de transición. Este es el caso particular de resolverlo mediante la simulación de la cadena.


dtmc-dir.py
===========

Calcula la distribución estacionaria de una cadena de Markov que modele el protocolo  S-ALOHA, dada su matriz de transición. Emplea el método directo para calcular dicha distribución estacionaria.


dtmc-gs.py
==========

Calcula la distribución estacionaria de una cadena de Markov que modele el protocolo  S-ALOHA, dada su matriz de transición. Utiliza el método de Gauss-Seidel.


dtmc-throughput.py
==================

Combina todos los programas anteriores, pero para graficar el throughput unicamente. No muestra soluciones particulares.


graph_cov.py
============

Grafica el Coeficiente de Variación de la distribución Hiperexponencial como función de los parámetros \lambda_1 y \lambda_2 en un espacio tridimensional.

phasedist.py
============

Concentra las tres distribuciones tipo fase vistas en clase (Exponencial, Erlang e Hiperexponencial), grafica y muestrea valores experimentales con éstas distribuciones.


search-cov-he.py
================

Utilidad para encontrar valores de \lambda_1, \lambda_2, p que satisfagan cierto valor de CoV. Realiza una búsqueda de fuerza bruta en un conjunto reducido de valores. Si no se encuentran parámetros que satisfagan el CoV deseado, entonces se debe modificar el código para ampliar el espacio de búsqueda.


ctmc.py
=======
Para solucionar cadenas de Markov en tiempo continuo.