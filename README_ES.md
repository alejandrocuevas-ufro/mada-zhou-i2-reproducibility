# Archivos de reproducibilidad — I² bivariado de Zhou–Dendukuri en `mada`

Este repositorio acompaña la auditoría metodológica del I² bivariado de Zhou–Dendukuri informado por `summary.reitsma()` en `mada` 0.5.12.

## Qué permite reproducir

1. El I² que actualmente entrega `mada`.
2. El mismo cálculo sustituyendo únicamente las varianzas de los efectos combinados por los componentes de `Psi`.
3. El efecto de corregir únicamente los denominadores intraestudio.
4. La reconstrucción completa de Zhou–Dendukuri.
5. Los cinco reanálisis publicados incluidos en el manuscrito y el ejemplo de TCCS isquémico con I² próximo a cero.
6. La simulación Monte Carlo utilizada en la Figura 1.

## Ejecución

El análisis está fijado a **`mada` 0.5.12**. Si no está instalada esa versión:

```r
install.packages("remotes")
remotes::install_version("mada", version = "0.5.12", repos = "https://cloud.r-project.org")
```

Desde la carpeta raíz del repositorio:

```bash
Rscript run_all.R
```

Los resultados quedarán en `results/`.

El primer script registra además `sessionInfo()` y guarda el código exacto de `summary.reitsma()` correspondiente a la versión instalada. Esto permite documentar de forma verificable cuál fue la implementación auditada.

## Comprobación de los artículos publicados

Antes de calcular el I² corregido, `03_reanalyse_published.R` reproduce el valor informado en cada artículo. El script comprueba automáticamente que la diferencia respecto del valor publicado —informado con un decimal— sea ≤0,15 puntos porcentuales.

## Datos

`data/published_reanalyses_2x2.csv` contiene únicamente datos agregados TP, FP, TN y FN recuperados de publicaciones. El origen exacto de cada conjunto se encuentra en `data/data_provenance.csv`. No se incluyen datos individuales de pacientes.

## Simulación

La semilla es `20260918`. Se realizan 500 réplicas para K = 6, 10, 20, 40, 80 y 160, manteniendo constantes los parámetros poblacionales y los tamaños de los grupos con y sin enfermedad.

Los archivos de `expected/` son referencias de desarrollo. Antes del envío del manuscrito deben sustituirse, cuando corresponda, por los resultados obtenidos al ejecutar los scripts directamente con R y `mada` 0.5.12.

## Antes de depositar el repositorio

La versión que se archive con DOI debería contener:

- todos los archivos de `data/` y `scripts/`;
- los resultados finales producidos por `Rscript run_all.R`;
- `results/sessionInfo.txt`;
- la figura final utilizada en el manuscrito; y
- una versión etiquetada del repositorio que coincida con la versión sometida del artículo.
