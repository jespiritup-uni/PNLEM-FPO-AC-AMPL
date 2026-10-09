# Código y Formulaciones PNLEM para FPO-AC en AMPL/KNITRO
Este repositorio contiene el código fuente y datos para reproducir las simulaciones 
de la tesis "EVALUACION COMPARATIVA DEL DESEMPEÑO COMPUTACIONAL DE CUATRO FORMULACIONES
 DE PROGRAMACIÓN NO LINEAL ENTERA MIXTA PARA FLUJO DE POTENCIA ÓPTIMO".

## Requisitos
- AMPL (versión 1.9.29  o superior)
- Solver KNITRO (versión 16.0.0 )

## Cómo ejecutar
1. Abrir la consola de AMPL.
2. Subir y elegir la carpeta a ejecutar (1.1 Modelo para FPO Coordenada Polar; 1.2 Modelo para FPO Coordenada Retangular; 1.3 Modelo para FPO Flujo en Ramos; 1.4 Modelo para FPO Injeccion de Corrientes)
3. En el script que esta dentro de la carpeta según el modelo elegido de nombre (FPO_CP;FPO_CR;FPO_FR;FPO_IC) y con extension ".dat" elegir que sistema de prueba IEEE desea ejecutar por defecto ejecutara el sistema de 14 barras "include sistema14nos_nesta.dat;" si desea probar otro sistema solo descomentar borrando el "#" y comentar los otros sistemas con "#".
4. Ejecutar el script de nombre (FPO_CP;FPO_CR;FPO_FR;FPO_IC) y con extension ".run" con el siguiente comando: "include FPO_CP.run", si es para coordenadas polares; "include FPO_CR.run", si es para coordenadas rectangualares, "include FPO_FR", si es para flujo en ramos,"include FPO_IC", si es para inyeccion de corrientes.
5. Los resultados saldran en la ventada de salida del AMPL.