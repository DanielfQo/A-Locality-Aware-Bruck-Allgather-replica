#include <mpi.h>
#include <vector>
#include <cstring>
#include <iostream>
#include <string>
#include <algorithm>
#include <functional>

const int N = 2;         // elementos por proceso
const int WARMUP = 20;   // iteraciones de calentamiento
const int ITERS  = 1000;   // iteraciones cronometradas
const bool CHECK = true; // validar resultados

void bruck_allgather(MPI_Comm comm, int id, int p, int* data, const int* init_data, int n) {
    std::memcpy(data, init_data, n * sizeof(int));

    int dist = 1, step = 0;
    while (dist < p) {
        int sendto = (id - dist + p) % p;
        int recvfrom = (id + dist) % p;

        int send_count = std::min(dist * n, n * p - dist * n);
        int recv_count = send_count;

        MPI_Request reqs[2];
        if (recv_count > 0)
            MPI_Irecv(data + dist * n, recv_count, MPI_INT, recvfrom, 1000 + step, comm, &reqs[0]);
        if (send_count > 0)
            MPI_Isend(data, send_count, MPI_INT, sendto, 1000 + step, comm, &reqs[1]);
        if (send_count > 0 || recv_count > 0)
            MPI_Waitall(2, reqs, MPI_STATUSES_IGNORE);

        dist <<= 1;
        step++;
    }


    std::vector<int> tmp(n * p);
    for (int j = 0; j < p; ++j) {
        int src_block = (j + id) % p;
        std::memcpy(tmp.data() + j * n, data + src_block * n, n * sizeof(int));
    }
    std::memcpy(data, tmp.data(), n * p * sizeof(int));
}

void loc_bruck_allgather(MPI_Comm comm, int id, int p,
                         MPI_Comm comm_local, int id_local, int p_local,
                         int r_n, int region_id,
                         int* data, const int* init_data, int n)
{
    //  cada nodo hace Bruck dentro del comm_local
    std::vector<int> local_data(n * p_local);
    bruck_allgather(comm_local, id_local, p_local, local_data.data(), init_data, n);

    // solo el representante construye el buffer
    std::vector<int> tmp;
    if (id_local == 0) {
        tmp.assign(n * p, 0);
        int region_offset = region_id * p_local * n;
        std::memcpy(tmp.data() + region_offset, local_data.data(), n * p_local * sizeof(int));


        for (int other = 0; other < r_n; ++other) {
            if (other == region_id) continue;

            int sendto = other * p_local;   
            int recvfrom = other * p_local; 

            int send_offset = region_offset;                
            int recv_offset = other * p_local * n;          
            int count = p_local * n;                        

            MPI_Request reqs[2];
            MPI_Irecv(tmp.data() + recv_offset, count, MPI_INT, recvfrom, 3000 + region_id, comm, &reqs[0]);
            MPI_Isend(tmp.data() + send_offset, count, MPI_INT, sendto, 3000 + other, comm, &reqs[1]);
            MPI_Waitall(2, reqs, MPI_STATUSES_IGNORE);
        }

        // ahora tmp tiene los bloques de todas las regiones en las posiciones correctas
        std::memcpy(data, tmp.data(), n * p * sizeof(int));
    }

    // broadcast dentro del nodo

    MPI_Bcast(data, n * p, MPI_INT, 0, comm_local);
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    MPI_Comm comm = MPI_COMM_WORLD;
    int p, id;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &id);

    // detectar procesos por nodo y num de regiones
    char hostname[MPI_MAX_PROCESSOR_NAME];
    int len;
    MPI_Get_processor_name(hostname, &len);
    std::string host_str(hostname, len);

    int color = std::hash<std::string>{}(host_str) & 0x7FFFFFFF;
    MPI_Comm comm_local;
    MPI_Comm_split(comm, color, id, &comm_local);

    int id_local, p_local;
    MPI_Comm_rank(comm_local, &id_local);
    MPI_Comm_size(comm_local, &p_local);

    int r_n = p / p_local;
    int region_id = id / p_local;

    std::cout << "Proceso " << id << " en host " << host_str
              << " -> id_local=" << id_local << ", p_local=" << p_local
              << ", region=" << region_id << ", r_n=" << r_n << std::endl;

    std::vector<int> init(N), out1(N * p), out2(N * p), gold(N * p);
    for (int k = 0; k < N; ++k)
        init[k] = id * 100000 + k;

    // medicion bruck
    for (int w = 0; w < WARMUP; ++w) {
        bruck_allgather(comm, id, p, out1.data(), init.data(), N);
        MPI_Barrier(comm);
    }
    double t0 = MPI_Wtime();
    for (int it = 0; it < ITERS; ++it) {
        bruck_allgather(comm, id, p, out1.data(), init.data(), N);
        MPI_Barrier(comm);
    }
    double t1 = MPI_Wtime();
    double time_bruck = (t1 - t0) / ITERS;

    // medicion locality
    for (int w = 0; w < WARMUP; ++w) {
        loc_bruck_allgather(comm, id, p, comm_local, id_local, p_local, r_n, region_id, out2.data(), init.data(), N);
        MPI_Barrier(comm);
    }
    double t2 = MPI_Wtime();
    for (int it = 0; it < ITERS; ++it) {
        loc_bruck_allgather(comm, id, p, comm_local, id_local, p_local, r_n, region_id, out2.data(), init.data(), N);
        MPI_Barrier(comm);
    }
    double t3 = MPI_Wtime();
    double time_locbruck = (t3 - t2) / ITERS;

    
    if (CHECK) {
        MPI_Allgather(init.data(), N, MPI_INT, gold.data(), N, MPI_INT, comm);
        bool ok1 = (std::memcmp(out1.data(), gold.data(), sizeof(int) * N * p) == 0);
        bool ok2 = (std::memcmp(out2.data(), gold.data(), sizeof(int) * N * p) == 0);
        if (id == 0) {
            std::cout << "Check vs MPI_Allgather:\n";
            std::cout << "  Bruck: " << (ok1 ? "OK" : "MISMATCH") << "\n";
            std::cout << "  Locality-Aware Bruck: " << (ok2 ? "OK" : "MISMATCH") << "\n";
        }
    }

    // Resultados 
    double sum1 = 0.0, sum2 = 0.0;
    MPI_Reduce(&time_bruck, &sum1, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
    MPI_Reduce(&time_locbruck, &sum2, 1, MPI_DOUBLE, MPI_SUM, 0, comm);

    if (id == 0) {
        double avg_bruck = sum1 / p;
        double avg_locbruck = sum2 / p;
        std::cout << "\n===== Comparacion de rendimiento =====\n";
        std::cout << "P=" << p << "  P_local=" << p_local << "  Regiones=" << r_n << "\n";
        std::cout << "Bruck promedio: " << avg_bruck << " s\n";
        std::cout << "Locality-Aware Bruck promedio: " << avg_locbruck << " s\n";
        std::cout << "Aceleracion: " << (avg_bruck / avg_locbruck) << "x\n";
    }

    MPI_Finalize();
    return 0;
}
