# A-Locality-Aware-Bruck-Allgather-replica

# Resultados de Comparación**

## Parámetros de Ejecución
- **Número de Procesos (`P`)**: Representa la cantidad total de procesos en el sistema.
- **Número de Procesos por Nodo (`P_local`)**: Indica cuántos procesos se ejecutan en cada nodo.
- **Número de Regiones (`Regiones`)**: Define cuántas particiones locales (regiones) existen en el sistema.

## Comparación de Rendimiento

### Caso 1: 16 Procesos (P=16, P_local=8, Regiones=2)
- **Bruck promedio**: 0.00901163 s
- **Locality-Aware Bruck promedio**: 0.0042912 s
- **Aceleración**: 2.10003x

### Caso 2: 8 Procesos (P=8, P_local=4, Regiones=2)
- **Bruck promedio**: 0.00697801 s
- **Locality-Aware Bruck promedio**: 0.00351423 s
- **Aceleración**: 1.98565x

### Caso 3: 4 Procesos (P=4, P_local=2, Regiones=2)
- **Bruck promedio**: 0.00868393 s
- **Locality-Aware Bruck promedio**: 0.00395475 s
- **Aceleración**: 2.19582x

### Caso 4: 2 Procesos (P=2, P_local=1, Regiones=2)
- **Bruck promedio**: 0.00343597 s
- **Locality-Aware Bruck promedio**: 0.00220968 s
- **Aceleración**: 1.55496x

## Resultados de Verificación

Los algoritmos **Bruck** y **Locality-Aware Bruck** fueron verificados contra la referencia estándar **MPI_Allgather**. Los resultados mostraron que ambos métodos dieron los mismos resultados en cada ejecución.

### Verificación
- **Bruck**: OK
- **Locality-Aware Bruck**: OK

