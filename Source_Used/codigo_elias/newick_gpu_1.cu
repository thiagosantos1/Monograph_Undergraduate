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
int i, j, k, e;
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
int *nz_p; // pai do no
int *nz_f1; // filho da esquerda do no
int *nz_f2; // filho da direita do no
float *nz_trait; //característica a ser comparada com cada espécie 
float *nz_class_range; //Faixa para as classe de distância
float *nz_class_value; //Valores medios de I de Moran por classe de distância
float *nz_class_media; //Valores medios de I de Moran por classe de distância
float *nz_class_variance; //Variancia para cada classe de distância
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
float *nz_trait_d; 
float *nz_class_range_d; 
float *nz_class_value_d; 
unsigned int *nz_sig_d;
unsigned int *nz_hsig_d;
unsigned int *nz_hval_d;
//int pos_ins_d, idx_ni_d;


//
char *symb, **nz_sy;
char str_tmp[100];
char str_float[30];
int nbint, nbuint, nbhuint, nbfloat, nbefloat; // tamanho em bytes dos tipos basicos
curandState *seed_d;
float zero = 0.0; // para facilitar impressao da matriz de distancias


// Forward function declarations
long long start_timer();
long long stop_timer(long long start_time, char *name);

// print tree in newick format
char *toNewick(int raiz, int base);

// find next prime number greater than n
int nextprime( int n );

// kernel
__global__ void Load_memory_global_Gpu(int nnos, int *nz, float *nz_br, float *nz_dr, int *nz_qf,int *nz_qe, int *nz_p, int *nz_f1, int *nz_f2);
__global__ void My_Load_memory_global_Gpu(int nnos, int *nz, float *nz_br, float *nz_dr, int *nz_qf,int *nz_qe, int *nz_p, int *nz_f1, int *nz_f2);

__global__ void Load_memory_shared_Gpu(int nnos, int *nz, float *nz_br, float *nz_dr, int *nz_qf,int *nz_qe, int *nz_p, int *nz_f1, int *nz_f2);

__global__ void Insert_tree_Gpu(int nnos, int hnnos, int pos_ins, int idx_ni, int *nz, float *nz_br, int *nz_qf, int *nz_qe, int *nz_p, int *nz_f1, int *nz_f2, curandState *states, unsigned long seed);

__global__ void Matrix_distance_Gpu(int nnos, int hnnos, int *nz, float *nz_br, float *nz_dr, float *nz_de, int *nz_qf,int *nz_qe, int *nz_p, int *nz_f1, int *nz_f2, unsigned int *nz_sig, unsigned int *nz_hsig, unsigned int *nz_hval);


__global__ void I_moran_Gpu(int nnos, int nrClass, float *nz_de, float *nz_trait, float *nz_class_range, float *nz_class_value, float MeanY, float Variance);

// auxiliary kernel functions
__device__ int quadratic_probing_insert(unsigned int *nz_hsig, unsigned int *nz_hval, unsigned int sig, int val, int hnnos);
__device__ int quadratic_probing_search(unsigned int *nz_hsig, unsigned int *nz_hval, unsigned int sig, int hnnos);

__device__ inline void atomicFloatAdd(float *address, float val);

// Main program
int main(int argc, char *argv[])
{
	int qtdArvores = 1;
	int qtdBlock = 1;
	int qtdThreadsPerBlock = 1;
	int qualArvore = 0;	
	int tipoTransferencia = 1;
	long long *vetorTempo;
	
	printf("\nSyntax: newick <#qtdBlocos #qtdThreads <#numeroArvoreImprimir <#tipoTransferencia >>>");
	printf("\n   tipoTransferencia:");
	printf("\n       1: replicação feita na GPU, utilizando a memória GLOBAL como origem da copia");
	printf("\n       2: replicação feita na GPU, utilizando a memória COMPARTILHADA para acelerar transferência");
	printf("\n       3: replicação feita na CPU, em seguida, todos os dados são copiados para a GPU\n\n");



	if (argc >= 2)
		sscanf(argv[1], "%d", &qtdBlock);
	if (argc >= 3)
		sscanf(argv[2], "%d", &qtdThreadsPerBlock);
	if (argc >= 4)
		sscanf(argv[3], "%d", &qualArvore);
	if (argc >= 5)
		sscanf(argv[4], "%d", &tipoTransferencia);

	printf("qtdBlock: %d qtdThreadsPerBlock: %d => qtdArvores: %d) \n", qtdBlock, qtdThreadsPerBlock, qtdBlock*qtdThreadsPerBlock);
	printf("qualArvore: %d\n", qualArvore);
	
	vetorTempo = (long long *) malloc(5 * sizeof(long long));
	GPU_start_time = start_timer();
	dim3 grid(qtdBlock), block(qtdThreadsPerBlock);
	//Total de threads a serem criadas, conforme a quantidade de blocos e threads por bloco
	qtdArvores = grid.x * block.x;
 
	printf("qtdArvores  %d. \n", qtdArvores);
	fp = fopen("wellParser.out", "r");
	if (fp == NULL) {
		printf("\nCannot open file\n");
		exit(0);
	}

	fscanf(fp,"%d %d", &nnos, &idx_ni);
	printf("No nos: %d, Indice no interno: %d\n", nnos, idx_ni);

	nfol = nnos / 2;

	fscanf(fp,"%d %d", &pos_ins, &n_ins);
	printf("Inserir %d especies a partir de %d\n", n_ins, pos_ins);

	printf("Arvore: ");
	nz = (int *) malloc(nnos * sizeof(int) * qtdArvores);
	for(i=0; i<nnos; i++) {
		fscanf(fp,"%d", &nz[i]);
		printf("%d ", nz[i]);
	}
	printf("\n");

	printf("Simbolos: ");
	symb = (char *) malloc(50);
	nz_sy = (char **) malloc(nnos * sizeof(char *));
	for(i=0; i<nnos; i++) {
		fscanf(fp,"%s", symb);
		nz_sy[i] = (char *) malloc(50);
		strcpy(nz_sy[i], symb);
		printf("%s ", nz_sy[i]);
	}
	printf("\n");

	nz_dr = (float *) malloc(nnos * sizeof(float) * qtdArvores);
	
	ennos = (nfol * (nfol - 1)) / 2;
	size_t tamanho_ui = (unsigned int) (ennos * sizeof(float) * qtdArvores) ;
	nz_de = (float *) malloc( tamanho_ui );
	printf("\n\nennos * sizeof(float) * qtdArvores(%d): %u\n\n", qtdArvores, tamanho_ui );

	printf("Ramos: ");
	nz_br = (float *) malloc(nnos * sizeof(float) * qtdArvores);

	for(i=0; i<nnos; i++) {
		fscanf(fp,"%f", &nz_br[i]);
		printf("%.2f ", nz_br[i]);
	}
	printf("\n");

	printf("No Filhos: ");
	nz_qf = (int *) malloc(nnos * sizeof(int) * qtdArvores);
	for(i=0; i<nnos; i++) {
		fscanf(fp,"%d", &nz_qf[i]);
		printf("%d ", nz_qf[i]);
	}
	printf("\n");
	
	printf("No Especies: ");
	nz_qe = (int *) malloc(nnos * sizeof(int) * qtdArvores);
	for(i=0; i<nnos; i++) {
		fscanf(fp,"%d", &nz_qe[i]);
		printf("%d ", nz_qe[i]);
	}
	printf("\n");

	printf("Pais: ");
	nz_p = (int *) malloc(nnos * sizeof(int) * qtdArvores);
	for(i=0; i<nnos; i++) {
		fscanf(fp,"%d", &nz_p[i]);
		printf("%d ", nz_p[i]);
	}
	printf("\n");

	printf("Filhos 1: ");
	nz_f1 = (int *) malloc(nnos * sizeof(int) * qtdArvores);
	for(i=0; i<nnos; i++) {
		fscanf(fp,"%d", &nz_f1[i]);
		printf("%d ", nz_f1[i]);
	}
	printf("\n");

	printf("Filhos 2: ");
	nz_f2 = (int *) malloc(nnos * sizeof(int) * qtdArvores);
	for(i=0; i<nnos; i++) {
		fscanf(fp,"%d", &nz_f2[i]);
		printf("%d ", nz_f2[i]);
	}
	printf("\n");
	
	printf("Traits: ");
	nz_trait = (float *) malloc(nnos * sizeof(float) );
	for(i=0; i<nnos; i++) {
		fscanf(fp,"%f", &nz_trait[i]);
		printf("%f ", nz_trait[i]);
	}
	printf("\n");
	
	hnnos = nextprime(2*nnos);
	printf("\nPrimo/hnnos = %d, qtdArvores = %d, hnnos*qtdArvores = %d\n", hnnos, qtdArvores, hnnos * qtdArvores);
	nz_sig = (unsigned int *) malloc(hnnos * sizeof(unsigned int)  * qtdArvores);
	for(i=0; i<(hnnos * qtdArvores); i++) {
		nz_sig[i] = 0;
	}
	
	nz_hsig = (unsigned int *) malloc(hnnos * sizeof(unsigned int) * qtdArvores);
	nz_hval = (unsigned int *) malloc(hnnos * sizeof(unsigned int) * qtdArvores);
	
	for(i=0; i<(hnnos * qtdArvores); i++) {
		nz_hsig[i] = (unsigned int) EMPTY;
		nz_hval[i] = (unsigned int) EMPTY;
	}

	fclose(fp);

	toNewick(nnos-1, 0);
	printf(";\n");
	
	// move data to GPU
	nbint = nnos * sizeof(int);
	nbuint = nnos * sizeof(unsigned int);
	nbhuint = hnnos * sizeof(unsigned int);
	nbfloat = nnos * sizeof(float);
	nbefloat = ennos * sizeof(float);
	GPU_time = stop_timer(GPU_start_time, "\t Tempo: Preencher dados nas estruturas da CPU");


	GPU_start_time = start_timer();
	//cudaMalloc((void **)&pos_ins_d, sizeof(int));
	//cudaMalloc((void **)&idx_ni_d, sizeof(int));
/*	:printf("\n ________ \t tipos \t\t tipos*qtdArvores, nnos %d ", nnos);
	printf("\n nbint    \t %d \t\t %d ", nbint, nbint * qtdArvores);
	printf("\n nbfloat  \t %d \t\t %d ", nbfloat, nbfloat * qtdArvores);
	printf("\n nbefloat \t %d \t\t %d ", nbefloat, nbefloat * qtdArvores);
	printf("\n nbuint   \t %d \t\t %d ", nbuint, nbuint * qtdArvores);
	printf("\n ");
*/

	cudaDeviceReset();

	printf("\ncurandState: %d\n", sizeof(curandState));

	cudaMalloc((void **)&nz_d, nbint * qtdArvores);
    	cudaMalloc((void **)&nz_br_d, nbfloat * qtdArvores);
    	cudaMalloc((void **)&nz_dr_d, nbfloat * qtdArvores);
    	cudaMalloc((void **)&nz_qf_d, nbint * qtdArvores);
    	cudaMalloc((void **)&nz_qe_d, nbint * qtdArvores);
    	cudaMalloc((void **)&nz_p_d, nbint * qtdArvores);
    	cudaMalloc((void **)&nz_f1_d, nbint * qtdArvores);
    	cudaMalloc((void **)&nz_f2_d, nbint * qtdArvores);
    	cudaMalloc((void **)&seed_d, nnos*sizeof(curandState)*qtdArvores);
	GPU_time = stop_timer(GPU_start_time, "\t Tempo: Alocar memória na GPU");

	if( nz_d==0 ) {
		printf("couldn't allocate memory nz_d\n"); 
		return 1;
   	}
 	if( nz_br_d==0 ) {
		printf("couldn't allocate memory nz_br_d\n"); 
		return 1;
   	}
 	if( nz_dr_d==0 ) {
		printf("couldn't allocate memory nz_dr_d\n"); 
		return 1;
   	}
	if( nz_qf_d==0  ) {
		printf("couldn't allocate memory nz_qf_d\n"); 
		return 1;
   	}
 	if( nz_qe_d==0 ) {
		printf("couldn't allocate memory nz_qe_d\n"); 
		return 1;
   	} 
	if( nz_p_d==0 || nz_f1_d==0 || nz_f2_d==0 ) {
		printf("couldn't allocate memory 2\n"); 
		return 1;
   	} 
	if(seed_d ==0 ) {
		printf("couldn't allocate memory seed_d\n"); 
		return 1;
   	}


	GPU_start_time = start_timer();
	if (tipoTransferencia == 1 || tipoTransferencia == 2){
		cudaMemcpy(nz_d, nz, nbint, cudaMemcpyHostToDevice);
		cudaMemcpy(nz_br_d, nz_br, nbfloat, cudaMemcpyHostToDevice);
		cudaMemcpy(nz_qf_d, nz_qf, nbint, cudaMemcpyHostToDevice);
		cudaMemcpy(nz_qe_d, nz_qe, nbint, cudaMemcpyHostToDevice);
		cudaMemcpy(nz_p_d, nz_p, nbint, cudaMemcpyHostToDevice);
		cudaMemcpy(nz_f1_d, nz_f1, nbint, cudaMemcpyHostToDevice);
		cudaMemcpy(nz_f2_d, nz_f2, nbint, cudaMemcpyHostToDevice);

/*

		//cudaMemcpy(pos_ins_d, pos_ins, sizeof(int), cudaMemcpyHostToDevice);
		//cudaMemcpy(idx_ni_d, idx_ni, sizeof(int), cudaMemcpyHostToDevice);
*/
	}

	GPU_time = stop_timer(GPU_start_time, "\t Tempo para copiar dados (bases) para memória");
	vetorTempo[0] = GPU_time;
	int aux = sizeof(int)*nnos;

	/* OPÇÕES PARA GERAR OS DADOS NA GPU:
		1. Copiar os elementos das estruturas e replica-los na gpu, utilizando memória GLOBAL
		2. Copiar os elementos das estruturas e replica-los na gpu, utilizando memória COMPARTILHADA
		3. Replicar os elementos na CPU e copia-los para A GPU (memória global)
	*/
	if (tipoTransferencia == 1) {
		GPU_start_time = start_timer();
//		Load_memory_global_Gpu<<<grid, block, aux>>>(nnos, nz_d, nz_br_d, nz_dr_d, nz_qf_d, nz_qe_d, nz_p_d, nz_f1_d, nz_f2_d);
		My_Load_memory_global_Gpu<<<qtdArvores, nnos>>>(nnos, nz_d, nz_br_d, nz_dr_d, nz_qf_d, nz_qe_d, nz_p_d, nz_f1_d, nz_f2_d);
		cudaDeviceSynchronize();
		GPU_time = stop_timer(GPU_start_time, "\t Tempo para copiar memória GPU (transferencia via memoria global)");
		vetorTempo[0] += GPU_time;
	}else{
		if (tipoTransferencia == 2){
			GPU_start_time = start_timer();
			Load_memory_shared_Gpu<<<grid, block, aux>>>(nnos, nz_d, nz_br_d, nz_dr_d, nz_qf_d, nz_qe_d, nz_p_d, nz_f1_d, nz_f2_d);
			cudaDeviceSynchronize();
			GPU_time = stop_timer(GPU_start_time, "\t Tempo para copiar memória GPU (transferencia via memoria compartilhada)");
			vetorTempo[0] += GPU_time;
		}else{
				GPU_start_time = start_timer();
				int base = 0;
				for(i = 0; i < qtdArvores; i++) {
					if (i > 0) {
						for(j = 0; j < nnos; j++) {
							base = i * nnos;
							nz[base+j] = nz[j] + (nz[j] >= 0 ? base : 0);
							nz_br[base+j] = nz_br[j];
							nz_dr[base+j] = 0;
							//nz_de[base+j] = nz_de[j];
							nz_qf[base+j] = nz_qf[j];
							nz_qe[base+j] = nz_qe[j];
							nz_p[base+j] = nz_p[j] + (nz_p[j] >= 0 ? base : 0);
							nz_f1[base+j] = nz_f1[j] + (nz_f1[j] >= 0 ? base : 0);
							nz_f2[base+j] = nz_f2[j] + (nz_f2[j] >= 0 ? base : 0);
						}
						nz[nfol] = -i;
					}
				}
				cudaMemcpy(nz_d, nz, nbint * qtdArvores, cudaMemcpyHostToDevice);
				cudaMemcpy(nz_br_d, nz_br, nbfloat * qtdArvores, cudaMemcpyHostToDevice);
				cudaMemcpy(nz_dr_d, nz_dr, nbfloat * qtdArvores, cudaMemcpyHostToDevice);
				//cudaMemcpy(nz_de_d, nz_de, nbefloat * qtdArvores, cudaMemcpyHostToDevice);
				cudaMemcpy(nz_qf_d, nz_qf, nbint * qtdArvores, cudaMemcpyHostToDevice);
				cudaMemcpy(nz_qe_d, nz_qe, nbint * qtdArvores, cudaMemcpyHostToDevice);
				cudaMemcpy(nz_p_d, nz_p, nbint * qtdArvores, cudaMemcpyHostToDevice);
				cudaMemcpy(nz_f1_d, nz_f1, nbint * qtdArvores, cudaMemcpyHostToDevice);
				cudaMemcpy(nz_f2_d, nz_f2, nbint * qtdArvores, cudaMemcpyHostToDevice);

				GPU_time = stop_timer(GPU_start_time, "\t Tempo para copiar da CPU->GPU (carregar dados)");
				vetorTempo[0] += GPU_time;
		}
	}







	cudaDeviceSynchronize();
	/**************************************************
	*
	* I N S E R I R   E S P E C I E S   P E R D I D A S 
	*
	******************************************************/
	// call kernel

	GPU_start_time = start_timer();
	if (n_ins > 0){ //se houver nós a inserir
		Insert_tree_Gpu<<<grid, block, aux>>>(nnos, hnnos, pos_ins, idx_ni, nz_d, nz_br_d, nz_qf_d, nz_qe_d, nz_p_d, nz_f1_d, nz_f2_d, seed_d, time(NULL));
		printf("Erro (inserir): %s\n", cudaGetErrorString( cudaGetLastError() ) );
	}
	cudaDeviceSynchronize();
	GPU_time = stop_timer(GPU_start_time, "\t Tempo para incluir nós na árvore");
	vetorTempo[1] = GPU_time;


				
	//alocar memoria para outras vetores

    	cudaMalloc((void **)&nz_de_d, nbefloat * qtdArvores);
    	cudaMalloc((void **)&nz_sig_d, nbhuint * qtdArvores);
    	cudaMalloc((void **)&nz_hsig_d, nbhuint * qtdArvores);
    	cudaMalloc((void **)&nz_hval_d, nbhuint * qtdArvores);
 	if( nz_de_d==0  ) {
		printf("couldn't allocate memory nz_de_d\n"); 
		return 1;
   	} 
	if( nz_sig_d==0) {
		printf("couldn't allocate memory nz_sig_d\n"); 
		return 1;
   	} 
	if( nz_hsig_d==0 ) {
		printf("couldn't allocate memory nz_hsig_d\n"); 
		return 1;
   	} 
	if( nz_hval_d==0 ) {
		printf("couldn't allocate memory nz_hval_d\n"); 
		return 1;
   	} 

	cudaMemcpy(nz_sig_d, nz_sig, nbhuint * qtdArvores, cudaMemcpyHostToDevice);
	cudaMemcpy(nz_hsig_d, nz_hsig, nbhuint * qtdArvores, cudaMemcpyHostToDevice);
	cudaMemcpy(nz_hval_d, nz_hval, nbhuint * qtdArvores, cudaMemcpyHostToDevice);

	/**************************************************
	*
	* C A L C U L A R   A   M A T R I Z   D E   D I S T A N C I A 
	*
	******************************************************/
	cudaDeviceSynchronize();
	GPU_start_time = start_timer();
//	int nb = qtdArvores;
	Matrix_distance_Gpu<<<qtdArvores, nfol>>>(nnos, hnnos, nz_d, nz_br_d, nz_dr_d, nz_de_d, nz_qf_d, nz_qe_d, nz_p_d, nz_f1_d, nz_f2_d, nz_sig_d, nz_hsig_d, nz_hval_d);
	cudaDeviceSynchronize();
	printf("Erro (matrix distancia): %s\n", cudaGetErrorString( cudaGetLastError() ) );
	GPU_time = stop_timer(GPU_start_time, "\t Tempo total para calcular a matriz de distância");

	vetorTempo[2] = GPU_time;

	GPU_start_time = start_timer();

	// copy data back to the CPU
	//cudaMemcpy(pos_ins, pos_ins_d, sizeof(int), cudaMemcpyDeviceToHost);
	//cudaMemcpy(idx_ni, idx_ni_d, sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(nz, nz_d, nbint * qtdArvores, cudaMemcpyDeviceToHost);
	cudaMemcpy(nz_br, nz_br_d, nbfloat * qtdArvores, cudaMemcpyDeviceToHost);
	cudaMemcpy(nz_dr, nz_dr_d, nbfloat * qtdArvores, cudaMemcpyDeviceToHost);
	cudaMemcpy(nz_de, nz_de_d, nbefloat * qtdArvores, cudaMemcpyDeviceToHost);
	cudaMemcpy(nz_qf, nz_qf_d, nbint * qtdArvores, cudaMemcpyDeviceToHost);
	cudaMemcpy(nz_qe, nz_qe_d, nbint * qtdArvores, cudaMemcpyDeviceToHost);
	cudaMemcpy(nz_p, nz_p_d, nbint * qtdArvores, cudaMemcpyDeviceToHost);
	cudaMemcpy(nz_f1, nz_f1_d, nbint * qtdArvores, cudaMemcpyDeviceToHost);
	cudaMemcpy(nz_f2, nz_f2_d, nbint * qtdArvores, cudaMemcpyDeviceToHost);
	cudaMemcpy(nz_sig, nz_sig_d, nbhuint * qtdArvores, cudaMemcpyDeviceToHost);
	cudaMemcpy(nz_hsig, nz_hsig_d, nbhuint * qtdArvores, cudaMemcpyDeviceToHost);
	cudaMemcpy(nz_hval, nz_hval_d, nbhuint * qtdArvores , cudaMemcpyDeviceToHost);
	GPU_time = stop_timer(GPU_start_time, "\t Tempo copiar dados de volta (GPU -> cpu): ");
	vetorTempo[3] = GPU_time; //Copiar dados da Gpu para cpu
	
/*
	printf("\n\nImprimir uma arvore: \n");
	toNewick(nnos-1, 0);
	printf(";\n");
*/

	cudaDeviceSynchronize();
	//Desalocar memoria da GPU para utilizar no próximo kernel
	GPU_start_time = start_timer();

	cudaFree(nz_d);    
	cudaFree(nz_br_d);
	cudaFree(nz_dr_d);
	cudaFree(nz_qf_d);
	cudaFree(nz_qe_d);
	cudaFree(nz_p_d);    
	cudaFree(nz_f1_d);
	cudaFree(nz_f2_d);

	cudaFree(nz_sig_d);
	cudaFree(nz_hsig_d);
	cudaFree(nz_hval_d);
	cudaFree(seed_d);

	free(nz_qf);
	free(nz_qe);
	free(nz_f1);
	free(nz_f2);
	free(nz_p);
	free(nz_dr);

	free(nz_sig);
	free(nz_hsig);
	free(nz_hval);

	GPU_time = stop_timer(GPU_start_time, "\t Tempo: Liberar memória GPU");


	/**************************************************
	*
	* C A L C U L A R   I   D E   M O R A N
	*
	******************************************************/
	//Aloca posicoes em memoria para armazenar as classes de distância
	int nrClass = 4;
	float maiorDistancia=0, menorDistancia = nz_de[0], salto;
	nz_class_range = (float *) malloc((nrClass+1) * sizeof(float));
	nz_class_value = (float *) malloc(nrClass * sizeof(float) * qtdArvores);
	nz_class_media = (float *) malloc(nrClass * sizeof(float) );
	nz_class_variance = (float *) malloc(nrClass * sizeof(float) );
	//As classes são definidas de forma igual, entre o maior e menor valor	
	for (i=0;i<ennos;i++){
		if (maiorDistancia < nz_de[i])
			maiorDistancia = nz_de[i];	
		if (menorDistancia > nz_de[i])
			menorDistancia = nz_de[i];
	}

	//nz_class_range[0] = menorDistancia;
	salto = (maiorDistancia - menorDistancia)/nrClass;
	for(i=0;i<nrClass;i++){
		nz_class_range[i] = menorDistancia;
		nz_class_value[i] = 0.0;
		menorDistancia += salto;
	}
	nz_class_range[0] -= nz_class_range[0]/2; //para incluir distancias iguais ao menor valor
	nz_class_range[i] = maiorDistancia;

	//realiza uma cópia do vetor de características (são as mesmas para todas as especies, independente da posição na árvore)
    	cudaMalloc((void **)&nz_trait_d, nbfloat);
    	cudaMalloc((void **)&nz_class_range_d, sizeof(float) * (nrClass+1)); //+1 para guardar a faixa final da classe
    	cudaMalloc((void **)&nz_class_value_d, sizeof(float) * nrClass * qtdArvores);

	cudaMemcpy(nz_trait_d, nz_trait, nbfloat, cudaMemcpyHostToDevice);
	cudaMemcpy(nz_class_range_d, nz_class_range, sizeof(float) * (nrClass+1), cudaMemcpyHostToDevice);
	cudaMemcpy(nz_class_value_d, nz_class_value, sizeof(float) * nrClass * qtdArvores, cudaMemcpyHostToDevice);



	float Variance, MeanY, SumW;

	SumW = 0;
  	Variance = 0;

  	for (int d=0;d<nfol;d++){
    		SumW = SumW + nz_trait[d];
    		Variance = Variance + pow(nz_trait[d],2);
	}
  	MeanY = SumW / nfol;
	Variance = Variance - (pow(SumW, 2) / nfol);

	cudaDeviceSynchronize();
	GPU_start_time = start_timer();
	aux = sizeof(float)*(nrClass+1);
	I_moran_Gpu<<<qtdArvores, nfol, aux>>>(nnos, nrClass, nz_de_d, nz_trait_d, nz_class_range_d, nz_class_value_d, MeanY, Variance);
	cudaDeviceSynchronize();
	printf("Erro (I_moran_Gpu): %s\n", cudaGetErrorString( cudaGetLastError() ) );
	GPU_time = stop_timer(GPU_start_time, "\t Tempo para calcular o Indice de Moran): ");
	vetorTempo[4] = GPU_time; //Copiar dados da Gpu para cpu


	//Traz os resultados de volta (GPU para Host), as medias são armazenadas no início do vetor
	cudaMemcpy(nz_class_value, nz_class_value_d, nrClass * sizeof(float) * qtdArvores, cudaMemcpyDeviceToHost);

	//Calcula a media por classe e a variancia
	float media;
	int nrArvore;
	for(i=0;i<nrClass;i++){
		media = 0;
		for (nrArvore=i;nrArvore<(qtdArvores*nrClass);nrArvore+=nrClass){
			media += nz_class_value[nrArvore];
		}
		nz_class_media[i] = media / qtdArvores;
	}
	//calculo da variancia	
	for(i=0;i<nrClass;i++){
		media = 0;
		for (nrArvore=i;nrArvore<(qtdArvores*nrClass);nrArvore+=nrClass){
			media += pow((nz_class_value[nrArvore] -  nz_class_media[i]), 2);
		}
		nz_class_variance[i] = media / qtdArvores;
	}

	GPU_start_time = start_timer();


	/**************************************************
	*
	* E X I B I R   R E S U L T A D O S 
	*
	******************************************************/

/*

	printf("\nnz_sy, ");
	for (int jx=0;jx<qtdArvores;jx++)	
		for(i=0;i<nnos;i++){
			printf("%s,", nz_sy[i]);
		}

	printf("\nnz, ");
	for(i=0;i<(qtdArvores*nnos);i++){
		printf("%d,", nz[i]);
	}
*/
	printf("\nnz_br,");
	for(i=0;i<(qtdArvores*nnos);i++){
		printf("%f,", nz_br[i]);
		if (i == 10000) //limitar impressao para nao deixar os arquivos muito grandes
			break;
	}
/*
	printf("\nnz_dr,");
	for(i=0;i<(qtdArvores*nnos);i++){
		printf("%f,", nz_dr[i]);
	}

	printf("\nnz_qf,");
	for(i=0;i<(qtdArvores*nnos);i++){
		printf("%d,", nz_qf[i]);
	}


	printf("\nnz_qe,");
	for(i=0;i<(qtdArvores*nnos);i++){
		printf("%d,", nz_qe[i]);
	}


	printf("\nnz_p,");
	for(i=0;i<(qtdArvores*nnos);i++){
		printf("%d,", nz_p[i]);
//		if (i == 10000)
//			break;
	}

	printf("\nnz_f1,");
	for(i=0;i<(qtdArvores*nnos);i++){
		printf("%d,", nz_f1[i]);
		if (i == 10000)
			break;
	}

	printf("\nnz_f2,");
	for(i=0;i<(qtdArvores*nnos);i++){
		printf("%d,", nz_f2[i]);
		if (i == 10000)
			break;
	}
*/
	printf("\nnz_class,");
	for(i=0;i<(nrClass);i++){
		printf("\t\n [%d] %f => value: %f ; media: %f ; variance: %f ", i, nz_class_range[i], nz_class_value[i], nz_class_media[i], nz_class_variance[i]);
	}

	printf("\n");
	printf("\n");

/*
	for(i=1;i<=qtdArvores;i++){
		//toNewick((nnos)-1);
		toNewick((i*nnos)-1, (nnos*(i-1)));
		printf(";\n");
	}

	printf("Pais: ");
	for(i=0; i<(nnos*qtdArvores); i++) {
		printf("%d ", nz_p[i]);
	}
	printf("\n");

	printf("Dst Raiz: ");
	for(i=0; i<(nnos*qtdArvores); i++) {
		if ((i-((i/nnos)*nnos)) == nfol) continue; // desconta o no da posicao nfol
		printf("%.2f ", nz_dr[i]); // pois este nao e usado
	}
	printf("\n");



	printf("Assinatura: ");
	for(i=0; (i<(hnnos * qtdArvores)); i++) {
		if ((i-((i/nnos)*nnos)) == nfol) continue; // desconta o no da posicao nfol
		if (i == nfol) continue;
		printf("%u ", nz_sig[i]);
	}
	printf("\n");
	
	printf("Hash Sign: ");
	for(i=0; (i<(hnnos * qtdArvores)); i++) {
		if (i == nfol) continue;
		if ((i-((i/nnos)*nnos)) == nfol) continue; // desconta o no da posicao nfol
		printf("%u ", nz_hsig[i]);
	}
	printf("\n");
	
	printf("Hash Val: ");
	for(i=0; (i<(hnnos * qtdArvores)); i++) {
		if ((i-((i/nnos)*nnos)) == nfol) continue; // desconta o no da posicao nfol
		if (i == nfol) continue;
		printf("%u ", nz_hval[i]);
	}
	printf("\n");
*/

	e = 0; // indexa a matriz triangular superior (representada num array) que contem a distancia
		   // entre as especies
	printf("Distancias: \n");
	printf("%7s ", nz_sy[0]);


	for(i=1; i<nfol; i++)
		printf("%4s ", nz_sy[i]);
	printf("\n");
	for(i=0; i<nfol; i++) {
		printf("%3s ", nz_sy[i]);
//		if (i >= (nfol-3)){
			for(j=0; j<=i; j++)
				printf("%.2f ", zero);
			for(k=i+1; k<nfol; k++) {
					printf("%.2f ", nz_de[e+(qualArvore*ennos)]);
				e++;
			}
//		}
		printf("\n");
	}

/*
	printf("\n\nnz_de: ");
	for(k=0; k<(ennos*qtdArvores); k++) {
			printf("%.2f ", nz_de[k]);
	}	
*/
	GPU_time = stop_timer(GPU_start_time, "\t Tempo mostrar dados em tela");

	GPU_start_time = start_timer();
	printf("\n");
	free(nz);
	free(nz_br);
	free(nz_de);
	free(symb);
	free(nz_sy);
	GPU_time = stop_timer(GPU_start_time, "\t Tempo: Liberar memoria da CPU");

	GPU_start_time = start_timer();
	cudaFree(nz_de_d);
	GPU_time = stop_timer(GPU_start_time, "\t Tempo: Liberar memória GPU (matriz de distancia)");

	printf("\n\n===================   R e s u m o    d o s    T e m p o s   =======================");

	printf("\nCPU -> GPU\tIncluir Esp. \tMatriz dist.\tGPU-> CPU\tI de Moran (em sec)");
	printf("\n%.5f \t", ((float) vetorTempo[0]) / (1000 * 1000));
	printf("%.5f \t", ((float) vetorTempo[1]) / (1000 * 1000));
	printf("%.5f \t", ((float) vetorTempo[2]) / (1000 * 1000));
	printf("%.5f \t", ((float) vetorTempo[3]) / (1000 * 1000));
	printf("%.5f \t\n", ((float) vetorTempo[4]) / (1000 * 1000));
	
	free(vetorTempo);
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



char *toNewick(int idRaiz, int base) {
    strcpy(str_tmp,"");
    strcpy(str_float,"");
	if (nz_f1[idRaiz] < 0) { // Não tem filhos 
		if ((idRaiz-base) < 0 || (idRaiz-base) > (nnos-1))
			printf("ERRO %d\n", (idRaiz-base));
		else
			strcat (str_tmp, nz_sy[idRaiz-base]);
		strcat (str_tmp, ":");
		sprintf(str_float,"%0.2f", nz_br[idRaiz]); 
		strcat (str_tmp, str_float);
		return str_tmp;
	} else { // Tem filhos 
		printf("(");
		printf("%s", toNewick(nz_f1[idRaiz], base));
		printf(",");
		printf("%s", toNewick(nz_f2[idRaiz], base));
		printf(")");
		printf("%s", nz_sy[idRaiz-base]);
		printf(":");
		sprintf(str_float,"%0.2f", nz_br[idRaiz]); 
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


__global__ void Load_memory_global_Gpu(int nnos, int *nz, float *nz_br, float *nz_dr, int *nz_qf,int *nz_qe, int *nz_p, int *nz_f1, int *nz_f2) {
	int i = threadIdx.x; // identificador da thread
	int index, new_index;
	int base = (blockIdx.x * blockDim.x * nnos) + nnos * i;

	//Todas as threads copiam os dados para suas respectivas áreas
	for(index = 0; index < nnos; index++){
		new_index = base+index;		
		nz[new_index] = nz[index] + (nz[index] >= 0 ? base : 0);
		nz_br[new_index] = nz_br[index];
		nz_dr[new_index] = 0;
		nz_qf[new_index] = nz_qf[index];
		nz_qe[new_index] = nz_qe[index];
		nz_p[new_index] = nz_p[index] + (nz_p[index] >= 0 ? base : 0);
		nz_f1[new_index] = nz_f1[index] + (nz_f1[index] >= 0 ? base : 0);
		nz_f2[new_index] = nz_f2[index] + (nz_f2[index] >= 0 ? base : 0);
	}

}

__global__ void My_Load_memory_global_Gpu(int nnos, int *nz, float *nz_br, float *nz_dr, int *nz_qf,int *nz_qe, int *nz_p, int *nz_f1, int *nz_f2) {
	int index = threadIdx.x; // identificador da thread
	int new_index;
	int base = (blockIdx.x * blockDim.x);

	//Todas as threads copiam os dados para suas respectivas áreas
		new_index = base+index;		
		nz[new_index] = nz[index] + (nz[index] >= 0 ? base : 0);
		nz_br[new_index] = nz_br[index];
		nz_dr[new_index] = 0;
		nz_qf[new_index] = nz_qf[index];
		nz_qe[new_index] = nz_qe[index];
		nz_p[new_index] = nz_p[index] + (nz_p[index] >= 0 ? base : 0);
		nz_f1[new_index] = nz_f1[index] + (nz_f1[index] >= 0 ? base : 0);
		nz_f2[new_index] = nz_f2[index] + (nz_f2[index] >= 0 ? base : 0);

}

__global__ void Load_memory_shared_Gpu(int nnos, int *nz, float *nz_br, float *nz_dr, int *nz_qf, int *nz_qe, int *nz_p, int *nz_f1, int *nz_f2) {
	extern __shared__ float nzTemp[];

	int i = threadIdx.x; // identificador da thread
	int index;


	int base = (blockIdx.x * blockDim.x * nnos) + nnos * i;

	//Copiar dados do vetor NZ
		if (threadIdx.x == 0)
			for(index = 0; index < nnos; index++)
				nzTemp[index] = nz[index];
		__syncthreads(); 

		for(index = 0; index < nnos; index++)
			nz[base+index] = (int) (nzTemp[index] + (nzTemp[index] >= 0 ? base : 0));
		__syncthreads(); 

	//Copiar dados do vetor BR
		if (threadIdx.x == 0)
			for(index = 0; index < nnos; index++)
				nzTemp[index] = nz_br[index];
		__syncthreads(); 

		for(index = 0; index < nnos; index++)
			nz_br[base+index] = nzTemp[index];
		__syncthreads(); 

	//Copiar dados do vetor QF
		if (threadIdx.x == 0)
			for(index = 0; index < nnos; index++)
				nzTemp[index] = nz_qf[index];
		__syncthreads(); 

		for(index = 0; index < nnos; index++){
			nz_dr[base+index] = 0;			
			nz_qf[base+index] = nzTemp[index];
		}
		__syncthreads(); 

	//Copiar dados do vetor QE
		if (threadIdx.x == 0)
			for(index = 0; index < nnos; index++)
				nzTemp[index] = nz_qe[index];
		__syncthreads(); 

		for(index = 0; index < nnos; index++)
			nz_qe[base+index] = nzTemp[index];
		__syncthreads(); 

	//Copiar dados do vetor P
		if (threadIdx.x == 0)
			for(index = 0; index < nnos; index++)
				nzTemp[index] = nz_p[index];
		__syncthreads(); 

		for(index = 0; index < nnos; index++)
			nz_p[base+index] = (int) (nzTemp[index] + (nzTemp[index] >= 0 ? base : 0));
		__syncthreads(); 

	//Copiar dados do vetor F1
		if (threadIdx.x == 0)
			for(index = 0; index < nnos; index++)
				nzTemp[index] = nz_f1[index];
		__syncthreads(); 

		for(index = 0; index < nnos; index++)
			nz_f1[base+index] = (int) (nzTemp[index] + (nzTemp[index] >= 0 ? base : 0));
		__syncthreads(); 

	//Copiar dados do vetor F2
		if (threadIdx.x == 0)
			for(index = 0; index < nnos; index++)
				nzTemp[index] = nz_f2[index];
		__syncthreads(); 

		for(index = 0; index < nnos; index++)
			nz_f2[base+index] = (int) (nzTemp[index] + (nzTemp[index] >= 0 ? base : 0));
//		__syncthreads(); 

}




__global__ void Insert_tree_Gpu(int nnos, int hnnos, int pos_ins, int idx_ni, int *nz, float *nz_br, int *nz_qf,int *nz_qe, int *nz_p, int *nz_f1, int *nz_f2, curandState *states, unsigned long seed) {
	int i = threadIdx.x; // identificador da thread
	float x; // valor gerado aleatoriamente
	unsigned int valor2; // numero entre 1 e maximo inteiro sem sinal
	unsigned int valor1; // numero entre 1 e altura da sub-arvore
	unsigned int shift = 8*sizeof(unsigned int)-1; // bits estao na faixa 0-31, e nao em 1-32
	unsigned int mask=1<<shift; // recebe 1 deslocado 31 vezes p/ direita 
								// (10000000 00000000 00000000 00000000)
	__shared__ int nfol; // numero de folhas da arvore
	int indMdcc; // no a partir do qual sera inserido uma especie
	int indNewNode; // aponta para o no internto a ser inserido junto com a especie a ser inserida
				 // idx_ni e o indice inicial dos nos internos a serem inseridos. Este indice
				 // cresce da direita para a esquerda. Veja que pos_ins aponta para a primeira
				 // especie a ser inserida. Serao inseridas nfol-pos_ins+1 especies.
	int indSisterSpecies;	
	int index;
	int indNewSpecies;

	int base = (blockIdx.x * blockDim.x * nnos) + nnos * i;

	index = 0;
	indSisterSpecies = 0;

	nfol = nnos / 2; // folhas estao na metade inferior
	curand_init(seed+i, base, 0, &states[base]);  // 	Initialize CURAND
	for(indNewSpecies=(base+pos_ins);indNewSpecies < (base+nfol);indNewSpecies++){
		curand(&states[base]);
		x = curand_uniform (&states[base]);      // gera numero aleatorio
		indNewNode = base + (idx_ni - ((indNewSpecies-base) - pos_ins)); // recebe um no interno a ser usado na insercao das especies
       								  
		indMdcc = nz[indNewSpecies]; // a posicao species [pos_ins <= species < nfol] contem o indice do no interno que
				   // sera usado para inserir a especie, i.e., ponto inicial de insercao (MDCC-most derived consensus clade)
					   
		valor1 = (int) (1 + x*nz_qf[indMdcc]); // numero entre 1 e altura da sub-arvore
		valor2 = (unsigned int) (1 + x*UINT_MAX);    // numero entre 1 e maximo inteiro sem sinal

		// a insercao e feita a partir do ponto de insercao mas seguindo os bits de valor2
		// se o bit for 1 avanca para a esquerda (f1) e se for 0 avanca para a direita (f2)
		if (indMdcc <= (base+nfol)) //Se o ponto de inserção for uma folha, então sobe um nível
			indMdcc = nz_p[indMdcc];
		else
			while (valor1 > 0) { // faca enquando nao alcancar a altura do no em questao ou um no
						 // folha seja alcancado.
				if(valor2 & mask)
					if (nz_f1[indMdcc] <= (base+nfol)) break;
					else indMdcc = nz_f1[indMdcc]; // avanca para proximo filho
				else
					if (nz_f2[indMdcc] <= (base+nfol)) break;
					else indMdcc = nz_f2[indMdcc]; // avanca para proximo filho
				valor2 <<= 1; // avanca para proximo bit
				valor1--; // diminui altura da arvore
			}

		//
		// convencao: f1 aa esquerda e f2 aa direita
		//
		x = curand_uniform (&states[base]);      // gera numero aleatorio - reuso de x
		if(valor2 & mask) { // insere no aa direita (f2) do no folha (especie) atual (f1)
			indSisterSpecies = nz_f1[indMdcc]; //nó a partir do qual o calculo do brach para a nova especie será realizado
			nz_f1[indNewNode] = nz_f1[indMdcc];
			nz_f2[indNewNode] = indNewSpecies;
			nz_p[nz_f1[indNewNode]] = indNewNode;
			nz_f1[indMdcc] = indNewNode;
			nz_br[indNewNode] = x * nz_br[nz_f1[indNewNode]];
			nz_br[nz_f1[indNewNode]] -= nz_br[indNewNode];
			nz_qf[indNewNode] = nz_qf[nz_f1[indNewNode]]++;
			nz_qe[indNewNode] = nz_qe[nz_f1[indNewNode]]++;
		}
		else { // insere no aa esquerda (f1) do no folha (especie) atual (f2)
			indSisterSpecies = nz_f2[indMdcc]; //nó a partir do qual o calculo do brach para a nova especie será realizado
			nz_f1[indNewNode] = indNewSpecies;
			nz_f2[indNewNode] = nz_f2[indMdcc];
			nz_p[nz_f2[indNewNode]] = indNewNode;
			nz_f2[indMdcc] = indNewNode;
			//Dividir o branch do nó "quebrado", de forma proporcional para o novo nó PAI (indNewNode)
			nz_br[indNewNode] = x * nz_br[nz_f2[indNewNode]];
			nz_br[nz_f2[indNewNode]] -= nz_br[indNewNode];
			nz_qf[indNewNode] = nz_qf[nz_f2[indNewNode]]++;	
			//Atualizar informacoes de quantidade de especies
			nz_qe[indNewNode] = nz_qe[nz_f2[indNewNode]]++;
		}			
		//atualiza vetor de pais

		nz_p[indNewSpecies] = indNewNode;
		nz_p[indNewNode] = indMdcc;
		nz_qe[indNewSpecies] = 1;

		//atualizar a qtde de especies e qtd de filhos
		index = nz_p[indNewNode];
		x = nnos/2;
		while( index > -1 || x <= 0 ){
			nz_qe[index] += 1; 
//			if (nz_f1[index] == -2 || nz_f2[index] == -2) break;
			if ( (nz_f1[index] >= (base+nfol) && nz_qf[nz_f1[index]] >= nz_qf[index]) || (nz_f2[index] >= (base+nfol) && nz_qf[nz_f2[index]] >= nz_qf[index])) 
				nz_qf[index] += 1;

			index = nz_p[index]; 
			x--;

		}

		//Calcular distancia para o nó inserido
		x = curand_uniform (&states[base]);      // gera numero aleatorio - reuso de x
		if (indSisterSpecies < (base+nfol)) //se irma eh folha, então branch deve possuir tamanho igual a irma
			nz_br[indNewSpecies] = nz_br[indSisterSpecies];
		else {
			valor2 = (unsigned int) (1 + x*UINT_MAX);    // numero entre 1 e maximo inteiro sem sinal
			nz_br[indNewSpecies] = 0.0;
			index = indSisterSpecies;

			while (true){
				nz_br[indNewSpecies] += nz_br[index];
				if (valor2 & mask){
					if (nz_f1[index] == -2) break;
					index = nz_f1[index];
				}else{
					if (nz_f2[index] == -2) break;
					index = nz_f2[index];

				}
			}

		}
	
	}
}



__global__ void Matrix_distance_Gpu(int nnos, int hnnos, int *nz, float *nz_br, float *nz_dr, float *nz_de, int *nz_qf,int *nz_qe, int *nz_p, int *nz_f1, int *nz_f2, unsigned int *nz_sig, unsigned int *nz_hsig, unsigned int *nz_hval) {
	float y; // acumula soma das arestas
	__shared__ int nfol; // numero de folhas da arvore
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





/*
Calcular o I de Moran para cada classe. Sao diversas arvores, cada uma tera o I de Moran para cada classe (nz_class),
em seguida faz-se a media e calcula a variancia entre elas.

Return: I de Moran por classe e a variancia para cada classe.
*/
__global__ void I_moran_Gpu(int nnos, int nrClass, float *nz_de, float *nz_trait, float *nz_class_range, float *nz_class_value, float MeanY, float Variance){
  	int d, r, c, a, b;
	int nfol, nthreads;
 	float SumProdCross, SumW, w;
  	short int p;
	int i = threadIdx.x; // identificador da thread
	int ib = blockIdx.x; // identificador do bloco
	int ennos;
	int base;
	__shared__ float sumTotal, sumTotalProdCross;
	extern __shared__ float nzClass[];

	for(d=0;d<nrClass;d++){
		nzClass[d] = nz_class_range[d];
	}

	nfol = nnos/2;
	ennos = (nfol * (nfol - 1)) / 2;
	base = ib * ennos;

	SumW = 0;

	if ( (nfol % 2) == 0) {
		nthreads = nfol - 1; // nfol é par: cada thread calcula nfol/2 distancias 
		a = nfol / 2; // quantidade de elementos por thread
	} else {
		nthreads = nfol;	 // nfol é ímpar: cada thread calcula (nfol-1)/2 distancias
		a = (nfol - 1) / 2;  // quantidade de elementos por thread
	}

	w = 1;
      	p = 2; //Symetric

	//Inicializa variaveis compartilhadas
	sumTotalProdCross = 0;
	sumTotal = 0;
	__syncthreads();//aguarda inicializacao das variaveis para continuar execução
	if (i < nthreads) {
  		for(d=0;d<nrClass;d++){
    			SumProdCross = 0;
	    		SumW = 0;
	
			for( b = i*a; b < a+(i*a); b++) {
				if (nz_de[b+base] > nzClass[d] && nz_de[b+base] <= nzClass[d+1]){
					r = row_index(b, nfol);
					c = column_index(b, nfol);
			
					SumW += (w*p);
					SumProdCross += (((nz_trait[r] - MeanY) * (nz_trait[c] - MeanY))*p);
				}
			}

			//Utilizar operacao atomica
				atomicFloatAdd(&sumTotalProdCross, SumProdCross);
				atomicFloatAdd(&sumTotal, SumW);
			__syncthreads(); // espera todas as threads chegarem até aqui
			//apenas uma thread calcula o I de Moran
			if (threadIdx.x == 0){ 
				nz_class_value[(ib*nrClass)+d] = (nfol / sumTotal) * (sumTotalProdCross / Variance);       // I de Moran
				sumTotalProdCross = 0;
				sumTotal = 0;
			}
			__syncthreads(); // espera todas as threads chegarem até aqui
   		}
	}
}


__device__ inline void atomicFloatAdd(float *address, float val)
{
	int tmp0 = *address;
	int i_val = __float_as_int(val + __int_as_float(tmp0));
	int tmp1;
	// compare and swap v = (old == tmp0) ? i_val : old;
	// returns old 
	while( (tmp1 = atomicCAS((int *)address, tmp0, i_val)) != tmp0 )
	{
		tmp0 = tmp1;
		i_val = __float_as_int(val + __int_as_float(tmp1));
	}
}
