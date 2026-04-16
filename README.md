# 🏗️ 3D Finite Element Analysis (FEA) Solver for I-Beam Structures

Este proyecto presenta un solver de **Elementos Finitos (MEF)** desarrollado en MATLAB para el análisis estructural de pórticos tridimensionales con perfiles metálicos de **Doble T**. 

El software es capaz de calcular desplazamientos, reacciones y esfuerzos internos en estructuras espaciales sometidas a cargas puntuales, desplazamientos impuestos y cargas distribuidas variables.

## 🚀 Características Principales
* **Análisis 3D Completo:** Uso de matrices de rotación $12 \times 12$ para orientar vigas en cualquier dirección del espacio.
* **Cargas Complejas:** Implementación de cargas trapezoidales variables discretizadas mediante funciones de forma.
* **Perfil Doble T:** Cálculo automático de propiedades mecánicas ($A, I_y, I_z, J$) para perfiles I-Beam.
* **Visualización:** Generación automática de la estructura deformada vs. original y mapa de ocupación de la matriz (Skyline).

## 📊 Resultados Visuales
![Estructura Deformada](images/deformada_resultado.png)
*Comparativa entre la geometría original y la deformada elástica bajo carga.*


## 🛠️ Instrucciones de Uso
1. Abre MATLAB.
2. Asegúrate de tener el archivo `main_structural_solver.m` en tu directorio de trabajo.
3. Ejecuta el script. El solver procesará automáticamente el ensamblaje de la matriz de rigidez global ($42 \times 42$) y mostrará los resultados en la consola y en figuras 3D.

## 👨‍💻 Autor
Desarrollado por **[Tu Nombre]** - Estudiante de Ingeniería Aeroespacial.
Especializado en cálculo estructural y métodos numéricos.
