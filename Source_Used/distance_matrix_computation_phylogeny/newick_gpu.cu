#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cuda.h>
#include <math.h>
#include <time.h>
#include <curand_kernel.h>
#include <sys/time.h>

#define TRUE 1

unsigned int EMPTY = UINT_MAX;

char str[200];
int i, j, k, e, p, na;
int nb; // numero de blocos a serem usados no kernel
FILE *fp;
int line=1;
int nnos, idx_ni, nfol; // numero de nos, indice de nos internos, numero de folhas
int hnnos; // tamanho da tabela hash
int ennos; // tamanho do vetor com as distancias entre as especies (matriz triangular superior)
int pos_ins, n_ins; // posicao de insercoes e numero de insercoes
int *nz; // contem indice do no; para os nos a serem inseridos, contem o indice do ponto de insercao
		 // para os nos internos a serem usados na insercao, contem -2
float *nz_br; // distancia do ramo (branch)
float *nz_dr; // distancias ate o no raiz
float *nz_de; // distancias entre especies
int *nz_qf; // altura do no
int *nz_qe; // quantidade de especies abaixo do no
int *nz_p; // pai dono
int *nz_f1; // filho da esquerda do no
int *nz_f2; // filho da direita do no
unsigned int *nz_sig; // assinatura do no - da o caminho em bits ate o raiz
unsigned int *nz_hsig; // hash da assinatura do no
unsigned int *nz_hval; // indice do no na tabela hash
long long GPU_start_time;
long long GPU_time;

// pointers to GPU memory
int *nz_d;
float *nz_br_d;
float *nz_dr_d;
float *nz_de_d;
int *nz_qf_d;
int *nz_qe_d;
int *nz_p_d;
int *nz_f1_d;
int *nz_f2_d;
unsigned int *nz_sig_d;
unsigned int *nz_hsig_d;
unsigned int *nz_hval_d;
//int pos_ins_d, idx_ni_d;
//
char *symb, **nz_sy, **n_arq;
char str_tmp[100];
char str_float[30];
int nbint, nbuint, nbhuint, nbfloat, nbefloat; // tamanho em bytes dos tipos basicos
curandState *seed_d;
float zero = 0.0; // para facilitar impressao da matriz de distancias
char arquivo[100];

// Forward function declarations
long long start_timer();
long long stop_timer(long long start_time, char *name);

// print tree in newick format
char *toNewick(int raiz);

// find next prime number greater than n
int nextprime( int n );

// kernel
__global__ void Mutate_tree_Gpu(int nnos, int hnnos, int pos_ins, int idx_ni, int *nz, float *nz_br, float *nz_dr, float *nz_de, int *nz_qf, int *nz_qe, int *nz_p, int *nz_f1, int *nz_f2, unsigned int *nz_sig, unsigned int *nz_hsig, unsigned int *nz_hval, curandState *states, unsigned long seed);

// auxiliary kernel functions
__device__ int quadratic_probing_insert(unsigned int *nz_hsig, unsigned int *nz_hval, unsigned int sig, int val, int hnnos);
__device__ int quadratic_probing_search(unsigned int *nz_hsig, unsigned int *nz_hval, unsigned int sig, int hnnos);

// Main program
int main()
{

	symb = (char *) malloc(50);

	fp = fopen("wellParser.out", "r");
	if (fp == NULL) {
		printf("\nCannot open file\n");
		exit(0);
	}
	
	fscanf(fp,"%d %d", &nb, &nnos); // numero de (arvores) blocos a serem usados no kernel
	printf("Num. arvores: %d, cada uma com %d nós\n", nb, nnos);

	nfol = nnos / 2;
	nz = (int *) malloc(nnos * sizeof(int));
	nz_sy = (char **) malloc(nnos * sizeof(char *));
	nz_dr = (float *) malloc(nnos * sizeof(float));	
	ennos = (nfol * (nfol - 1)) / 2;
	nz_de = (float *) malloc(ennos * sizeof(float));
	nz_br = (float *) malloc(nnos * sizeof(float));
	nz_qf = (int *) malloc(nnos * sizeof(int));
	nz_qe = (int *) malloc(nnos * sizeof(int));
	nz_p = (int *) malloc(nnos * sizeof(int));
	nz_f1 = (int *) malloc(nnos * sizeof(int));
	nz_f2 = (int *) malloc(nnos * sizeof(int));	
	hnnos = nextprime(2*nnos);
	nz_sig = (unsigned int *) malloc(hnnos * sizeof(unsigned int));
	nz_hsig = (unsigned int *) malloc(hnnos * sizeof(unsigned int));
	nz_hval = (unsigned int *) malloc(hnnos * sizeof(unsigned int));
	n_arq = (char **) malloc(nb * sizeof(char *));	  // guarda nome dos arquivos das arvores

	nbint = nnos * sizeof(int);
	nbuint = nnos * sizeof(unsigned int);
	nbhuint = hnnos * sizeof(unsigned int);
	nbfloat = nnos * sizeof(float);
	nbefloat = ennos * sizeof(float);
	
	cudaMalloc((void **)&nz_d, nb * nbint);
    cudaMalloc((void **)&nz_br_d, nb * nbfloat);
    cudaMalloc((void **)&nz_dr_d, nb * nbfloat);
    cudaMalloc((void **)&nz_de_d, nb * nbefloat);
    cudaMalloc((void **)&nz_qf_d, nb * nbint);
    cudaMalloc((void **)&nz_qe_d, nb * nbint);
    cudaMalloc((void **)&nz_p_d, nb * nbint);
    cudaMalloc((void **)&nz_f1_d, nb * nbint);
    cudaMalloc((void **)&nz_f2_d, nb * nbint);
    cudaMalloc((void **)&nz_sig_d, nb * nbuint);
    cudaMalloc((void **)&nz_hsig_d, nb * nbhuint);
    cudaMalloc((void **)&nz_hval_d, nb * nbhuint);
    cudaMalloc((void **)&seed_d, nb * nnos*sizeof(curandState));
    
    if( nz_d==0 || nz_br_d==0 || nz_dr_d==0 || nz_de_d==0 || nz_qf_d==0 || nz_qe_d==0 || nz_p_d==0 || nz_f1_d==0 || nz_f2_d==0 || nz_sig_d==0 || nz_hsig_d==0 || nz_hval_d==0 ) {
      printf("couldn't allocate memory\n"); 
      return 1;
	} 

////
	
	for(na=0; na < nb; na++) { 			// na = numero de arvores
	
		fscanf(fp,"%s", symb);
		if (na == 0) { printf("Arquivo: %s\n", symb); }
		n_arq[na] = (char *) malloc(50);
		strcpy(n_arq[na], symb);
	
		fscanf(fp,"%d %d", &nnos, &idx_ni);
		if (na == 0) { printf("No nos: %d, Indice no interno: %d\n", nnos, idx_ni); }
		fscanf(fp,"%d %d", &pos_ins, &n_ins);
		if (na == 0) { printf("Pos inic: %d, No insercoes: %d\n", pos_ins, n_ins); }

		if (na == 0) { printf("Arvore: "); }
		for(i=0; i<nnos; i++) {
			fscanf(fp,"%d", &nz[i]);
			if (na == 0) { printf("%d ", nz[i]); }
		}
		if (na == 0) { printf("\n"); }

		if (na == 0) { printf("Simbolos: "); }
		for(i=0; i<nnos; i++) {
			fscanf(fp,"%s", symb);
			nz_sy[i] = (char *) malloc(50);
			strcpy(nz_sy[i], symb);
			if (na == 0) { printf("%s ", nz_sy[i]); }
		}
		if (na == 0) { printf("\n"); }

		if (na == 0) { printf("Ramos: "); }
		for(i=0; i<nnos; i++) {
			fscanf(fp,"%f", &nz_br[i]);
			if (na == 0) { printf("%f ", nz_br[i]); }
		}
		if (na == 0) { printf("\n"); }

		if (na == 0) { printf("No Filhos: "); }
		for(i=0; i<nnos; i++) {
			fscanf(fp,"%d", &nz_qf[i]);
			if (na == 0) { printf("%d ", nz_qf[i]); }
		}
		if (na == 0) { printf("\n"); }
	
		if (na == 0) { printf("No Especies: "); }
		for(i=0; i<nnos; i++) {
			fscanf(fp,"%d", &nz_qe[i]);
			if (na == 0) { printf("%d ", nz_qe[i]); }
		}
		if (na == 0) { printf("\n"); }

		if (na == 0) { printf("Pais: "); }
		for(i=0; i<nnos; i++) {
			fscanf(fp,"%d", &nz_p[i]);
			if (na == 0) { printf("%d ", nz_p[i]); }
		}
		if (na == 0) { printf("\n"); }

		if (na == 0) { printf("Filhos 1: "); }
		for(i=0; i<nnos; i++) {
			fscanf(fp,"%d", &nz_f1[i]);
			if (na == 0) { printf("%d ", nz_f1[i]); }
		}
		if (na == 0) { printf("\n"); }

		if (na == 0) { printf("Filhos 2: "); }
		for(i=0; i<nnos; i++) {
			fscanf(fp,"%d", &nz_f2[i]);
			if (na == 0) { printf("%d ", nz_f2[i]); }
		}
		if (na == 0) { printf("\n"); }
	
		for(i=0; i<hnnos; i++) {
			nz_sig[i] = 0;
		}
	
		for(i=0; i<hnnos; i++) {
			nz_hsig[i] = (unsigned int) EMPTY;
			nz_hval[i] = (unsigned int) EMPTY;
		}

		if (na == 0) { toNewick(nnos-1); printf(";\n"); }

		if (na > 0) {
			for(j = 0; j < nnos; j++) {
				if (nz[j] >= 0) nz[j] = nz[j] + na*nnos;
				if (nz_p[j] >= 0) nz_p[j] = nz_p[j] + na*nnos;
				if (nz_f1[j] >= 0) nz_f1[j] = nz_f1[j] + na*nnos;
				if (nz_f2[j] >= 0) nz_f2[j] = nz_f2[j] + na*nnos;
			}
			nz[nfol] = -na;
		}
			// move data to GPU
		cudaMemcpy(nz_d + (na*nnos), nz, nbint, cudaMemcpyHostToDevice);
		cudaMemcpy(nz_br_d + (na*nnos), nz_br, nbfloat, cudaMemcpyHostToDevice);
		cudaMemcpy(nz_dr_d + (na*nnos), nz_dr, nbfloat, cudaMemcpyHostToDevice);
		cudaMemcpy(nz_de_d + (na*ennos), nz_de, nbefloat, cudaMemcpyHostToDevice);
		cudaMemcpy(nz_qf_d + (na*nnos), nz_qf, nbint, cudaMemcpyHostToDevice);
		cudaMemcpy(nz_qe_d + (na*nnos), nz_qe, nbint, cudaMemcpyHostToDevice);
		cudaMemcpy(nz_p_d + (na*nnos), nz_p, nbint, cudaMemcpyHostToDevice);
		cudaMemcpy(nz_f1_d + (na*nnos), nz_f1, nbint, cudaMemcpyHostToDevice);
		cudaMemcpy(nz_f2_d + (na*nnos), nz_f2, nbint, cudaMemcpyHostToDevice);
		cudaMemcpy(nz_sig_d + (na*nnos), nz_sig, nbuint, cudaMemcpyHostToDevice);
		cudaMemcpy(nz_hsig_d + (na*hnnos), nz_hsig, nbhuint, cudaMemcpyHostToDevice);
		cudaMemcpy(nz_hval_d + (na*hnnos), nz_hval, nbhuint, cudaMemcpyHostToDevice);
	}

	fclose(fp);

//	
	GPU_start_time = start_timer();
	
	// call kernel
	Mutate_tree_Gpu<<<nb, nfol>>>(nnos, hnnos, pos_ins, idx_ni, nz_d, nz_br_d, nz_dr_d, nz_de_d, nz_qf_d, nz_qe_d, nz_p_d, nz_f1_d, nz_f2_d, nz_sig_d, nz_hsig_d, nz_hval_d, seed_d, time(NULL));

	cudaThreadSynchronize(); // this is only needed for timing purposes	
	GPU_time = stop_timer(GPU_start_time, "\t Total");
//	
	
	p = nb-(nb/2)+(nb/4)-1;

	// copy data back to the CPU
	cudaMemcpy(nz, nz_d+p*nnos, nbint, cudaMemcpyDeviceToHost);
	cudaMemcpy(nz_br, nz_br_d+p*nnos, nbfloat, cudaMemcpyDeviceToHost);
	cudaMemcpy(nz_dr, nz_dr_d+p*nnos, nbfloat, cudaMemcpyDeviceToHost);
	cudaMemcpy(nz_de, nz_de_d+p*ennos, nbefloat, cudaMemcpyDeviceToHost);
	cudaMemcpy(nz_qf, nz_qf_d+p*nnos, nbint, cudaMemcpyDeviceToHost);
	cudaMemcpy(nz_qe, nz_qe_d+p*nnos, nbint, cudaMemcpyDeviceToHost);
	cudaMemcpy(nz_p, nz_p_d+p*nnos, nbint, cudaMemcpyDeviceToHost);
	cudaMemcpy(nz_f1, nz_f1_d+p*nnos, nbint, cudaMemcpyDeviceToHost);
	cudaMemcpy(nz_f2, nz_f2_d+p*nnos, nbint, cudaMemcpyDeviceToHost);
	cudaMemcpy(nz_sig, nz_sig_d+p*nnos, nbuint, cudaMemcpyDeviceToHost);
	cudaMemcpy(nz_hsig, nz_hsig_d+p*hnnos, nbhuint, cudaMemcpyDeviceToHost);
	cudaMemcpy(nz_hval, nz_hval_d+p*hnnos, nbhuint, cudaMemcpyDeviceToHost);
	
	if (p > 0) {
		for(j = 0; j < nnos; j++) {
			if (nz[j] >= 0) nz[j] = nz[j] - p*nnos;
			if (nz_p[j] >= 0) nz_p[j] = nz_p[j] - p*nnos;
			if (nz_f1[j] >= 0) nz_f1[j] = nz_f1[j] - p*nnos;
			if (nz_f2[j] >= 0) nz_f2[j] = nz_f2[j] - p*nnos;
		}
	}

	printf("Arquivo: %s\n", n_arq[p]);

	toNewick(nnos-1);
	printf(";\n");
//	
	printf("Arvore: ");
	for(i=0; i<nnos; i++) {
		printf("%d ", nz[i]);
	}
	printf("\n");
//
	printf("Pais: ");
	for(i=0; i<nnos; i++) {
		printf("%d ", nz_p[i]);
	}
	printf("\n");

	printf("f1: ");
	for(i=0; i<nnos; i++) {
		printf("%d ", nz_f1[i]);
	}
	printf("\n");

	printf("f2: ");
	for(i=0; i<nnos; i++) {
		printf("%d ", nz_f2[i]);
	}
	printf("\n");

	printf("Dst Raiz: ");
	for(i=0; i<nnos; i++) {
		if (i == nfol) continue; // desconta o no da posicao nfol
		printf("%.2f ", nz_dr[i]); // pois este nao e usado
	}
	printf("\n");

	printf("Assinatura: ");
	for(i=0; i<nnos; i++) {
		if (i == nfol) continue;
		printf("%u ", nz_sig[i]);
	}
	printf("\n");
	
	printf("Hash Sign: ");
	for(i=0; i<hnnos; i++) {
		if (i == nfol) continue;
		printf("%u ", nz_hsig[i]);
	}
	printf("\n");
	
	printf("Hash Val: ");
	for(i=0; i<hnnos; i++) {
		if (i == nfol) continue;
		printf("%u ", nz_hval[i]);
	}
	printf("\n");

	e = 0; // indexa a matriz triangular superior (representada num array) que contem a distancia
		   // entre as especies
	printf("Distancias: \n");
	printf("%7s ", nz_sy[0]);
	for(i=1; i<nfol; i++)
		printf("%4s\t", nz_sy[i]);
	printf("\n");
	for(i=0; i<nfol; i++) {
		for(j=0; j<=i; j++)
			printf("%.5f\t", zero);
		for(k=i+1; k<nfol; k++) {
			printf("%.5f\t", nz_de[e]);
			e++;
		}
		printf("\n");
	}
	printf("\n");
//		
	free(nz);
	free(nz_br);
	free(nz_dr);
	free(nz_de);
	free(nz_p);
	free(nz_f1);
	free(nz_f2);
	free(nz_sig);
	free(nz_hsig);
	free(nz_hval);
	free(symb);
	free(nz_sy);
	free(nz_qf);
	free(nz_qe);
//	
	cudaFree(nz_d);    
    cudaFree(nz_br_d);
    cudaFree(nz_dr_d);
    cudaFree(nz_de_d);
    cudaFree(nz_qf_d);
    cudaFree(nz_qe_d);
    cudaFree(nz_p_d);    
    cudaFree(nz_f1_d);
    cudaFree(nz_f2_d);
    cudaFree(nz_sig_d);
    cudaFree(nz_hsig_d);
    cudaFree(nz_hval_d);
//    
    return 0;
}

// Returns the current time in microseconds
long long start_timer() {
struct timeval tv;
gettimeofday(&tv, NULL);
return tv.tv_sec * 1000000 + tv.tv_usec;
}


// Prints the time elapsed since the specified time
long long stop_timer(long long start_time, char *name) {
struct timeval tv;
gettimeofday(&tv, NULL);
long long end_time = tv.tv_sec * 1000000 + tv.tv_usec;
printf("%s: %.5f sec\n", name, ((float) (end_time - start_time)) / (1000 * 1000));
return end_time - start_time;
}

char *toNewick(int idRaiz) {
  
    strcpy(str_tmp,"");
    strcpy(str_float,"");
		
	if (nz_f1[idRaiz] < 0) { // Não tem filhos 
		strcat (str_tmp, nz_sy[idRaiz]);
		strcat (str_tmp, ":");
//		sprintf(str_float,"%0.2f", nz_br[idRaiz]); 
		sprintf(str_float,"%f", nz_br[idRaiz]); 
		strcat (str_tmp, str_float);
		return str_tmp;
	} else { // Tem filhos 
		printf("(");
		printf("%s", toNewick(nz_f1[idRaiz]));
		printf(",");
		printf("%s", toNewick(nz_f2[idRaiz]));
		printf(")");
		printf("%s", nz_sy[idRaiz]);
		printf(":");
//		sprintf(str_float,"%0.2f", nz_br[idRaiz]);
		sprintf(str_float,"%f", nz_br[idRaiz]); 
		printf("%s", str_float);
		return "";
	}
}

int nextprime( int n ) {
	int Divisor, PossiblePrime;
	int FoundPrime;

	PossiblePrime = n;
	if( PossiblePrime <= 2 )
		PossiblePrime = 2;
	else
		if( PossiblePrime != 3 ) {
			if( PossiblePrime % 2 == 0 )
				PossiblePrime++;	/* Need An Odd Number */
			for( ; ; PossiblePrime += 2 ) {
				FoundPrime = !TRUE;
				for( Divisor = 3; PossiblePrime % Divisor; Divisor += 2 )
					if( Divisor * Divisor > PossiblePrime ) {
						FoundPrime = TRUE;
						break;
					}
				if( FoundPrime )
					break;
			}
		}
	return PossiblePrime;
}

__device__ int quadratic_probing_insert(unsigned int *nz_hsig, unsigned int *nz_hval, unsigned int sig, int val, int hnnos) {
    unsigned int j, hk, old;
    
    int ib = blockIdx.x; // identificador do bloco
    j = 0;
    hk = sig  % hnnos;
    while(j < hnnos) {
    	old = atomicCAS(&nz_hsig[hk+ib*hnnos], UINT_MAX, sig); // se posicao estiver vazia (UINT_MAX = EMPTY)
		if (old == UINT_MAX) {
			nz_hval[hk+ib*hnnos] = val;
			return (hk+ib*hnnos);
    	}
        j++;
        hk = (hk + j * j) % hnnos;
//        hk = (hk + j) % hnnos;
    }
    return (-1);
}

__device__ int quadratic_probing_search(unsigned int *nz_hsig, unsigned int *nz_hval, unsigned int sig, int hnnos) {
    unsigned int j, hk;
    
    int ib = blockIdx.x; // identificador do bloco
    j = 0;
    hk = sig  % hnnos;
    while(j < hnnos) {
		if (nz_hsig[hk+ib*hnnos] == sig) {
			return (nz_hval[hk+ib*hnnos]);
    	}
        j++;
        hk = (hk + j * j) % hnnos;
//        hk = (hk + j) % hnnos;
    }
    return (-1);
}

// estas duas funcoes sao usada para mapear os indices de um array para uma matriz triangular 
// superior correspondente (sem a diagonal). para uma matriz nxn, o array terá n(n-1)/2 elementos

__host__ __device__ int row_index( int i, int M ){ // retorna o indice da linha
	M--;
    float m = M;
    float row = (-2*m - 1 + sqrt( (4*m*(m+1) - 8*(float)i - 7) )) / -2;
    if( row == (float)(int) row ) row -= 1;
    return (int) row;
}

__host__ __device__ int column_index( int i, int M ){ // retorna o indice da coluna
    int row = row_index( i, M);
    M--;
    return 1 + (i - M * row + row*(row+1) / 2);
}


__global__ void Mutate_tree_Gpu(int nnos, int hnnos, int pos_ins, int idx_ni, int *nz, float *nz_br, float *nz_dr, float *nz_de, int *nz_qf,int *nz_qe, int *nz_p, int *nz_f1, int *nz_f2, unsigned int *nz_sig, unsigned int *nz_hsig, unsigned int *nz_hval, curandState *states, unsigned long seed) {
	float y; // acumula soma das arestas
	int nfol; // numero de folhas da arvore
	int j; // indice para thread ativa
	int a, b; // usados no calculo da faixa de elementos (da matriz triangular) a serem considerados
	unsigned int sig1, sig2, sig3, sig4; // assinaturas de tres nos - da o caminho em bits ate o raiz
	int bit; // contem bit sendo analizado
	int ancc; // indice do ancestral comum
	int nthreads; // numero de threads ativas
	int r, c; // linha e coluna da matriz triangular superior
	int bits; // conta quantos bits sao iguais
	int i = threadIdx.x; // identificador da thread
	int ib = blockIdx.x; // identificador do bloco
	int it; // indice de acesso global das threads
	int ennos; // tamanho da matriz de distancias
	
	nfol = nnos / 2; // folhas estao na metade inferior
	ennos = (nfol * (nfol - 1)) / 2;
	it = i + ib*nnos;
	
	if (i < nfol) { // nos folhos calculam distancia ate a raiz e armazena o caminho (assinatura
					// em bits) até a raiz
		y = 0;
		j = it;  // associa threads com nos folhas
		nz_sig[it] = 1;
		while (j != -1) {
			y = y + nz_br[j]; // acumula a distancia
			if (nz_p[j] == -1) 
				break;
			nz_sig[it] <<= 1;  // acumula o caminho
			if (nz_f1[nz_p[j]] == j) // acrescenta 0 se vier da direita (f2)
				nz_sig[it]++;         // ou 1 se vier da esquerda (f1)
			j = nz_p[j];
		}
		quadratic_probing_insert(nz_hsig, nz_hval, nz_sig[it], it, hnnos);
		nz_dr[it] = y;
	}

	__syncthreads(); // espera todas as threads chegarem até aqui

	if (i < (nfol-1)) { // nos internos calculam distancia ate a raiz e armazena o caminho
						// (assinatura em bits) até o raiz
		y = 0;
		j = it+nfol+1;  // associa threads com os nos internos
		nz_sig[j] = 1;
		if (nz_p[j] == -1) j = -1;
		while (j != -1) {
			y = y + nz_br[j]; // acumula a distancia
			if (nz_p[j] == -1) 
				break;
			nz_sig[it+nfol+1] <<= 1; // acumula o caminho
			if (nz_f1[nz_p[j]] == j)    // acrescenta 0 se vier da direita (f2)
				nz_sig[it+nfol+1]++; // ou 1 se vier da esquerda (f1)
			j = nz_p[j];
		}
		quadratic_probing_insert(nz_hsig, nz_hval, nz_sig[it+nfol+1], (it+nfol+1), hnnos);
		nz_dr[it+nfol+1] = y;
	}
	
	__syncthreads(); // espera todas as threads chegarem até aqui

	// se nfol (numero de especies) for impar, usamos nfol threads
	// se nfol (numero de especies) for par, usamos nfol-1 threads
	// isso evita termos que tratar de elementos restantes
	
	if ( (nfol % 2) == 0) {
		nthreads = nfol - 1; // nfol é par: cada thread calcula nfol/2 distancias 
		a = nfol / 2; // quantidade de elementos por thread
	} else {
		nthreads = nfol;	 // nfol é ímpar: cada thread calcula (nfol-1)/2 distancias
		a = (nfol - 1) / 2;  // quantidade de elementos por thread
	}

	if (i < nthreads) {
		for( b = i*a; b < a+(i*a); b++) {
			r = row_index(b, nfol);
			c = column_index(b, nfol);
			sig1 = nz_sig[r+ib*nnos];
			sig2 = nz_sig[c+ib*nnos];
			sig3 = 1; // inicia com 1 para diferenciar das demais assinaturas, i.e., 10, 100 etc
			bits = 0; // conta quantos bits sao iguais
			sig4 = 1; // recebe assinatura invertida
			while ( (sig1 & 1) == (sig2 & 1) && bits < 32) { // compara bit menos significativo
				bit = (sig1 & 1);
				bits++;
				sig1 >>= 1; // avanca para proximo bit
				sig2 >>= 1; // avanca para proximo bit
				sig3 <<= 1; // armazena bits coincidentes - caminho do ancestral comum
				if (bit) 
					sig3++; // soma 1 ou 0
			}
			while (bits>0) { // inverte a assinatura coincidente incluindo um 1 mais a esquerda
				sig4 <<= 1;
				if (sig3 & 1)
					sig4++;
				sig3 >>= 1;
				bits--;
			}
			ancc = quadratic_probing_search(nz_hsig, nz_hval, sig4, hnnos);
			nz_de[b+ib*ennos] = nz_dr[r+ib*nnos] + nz_dr[c+ib*nnos] - 2*nz_dr[ancc];
		}
	}
}
