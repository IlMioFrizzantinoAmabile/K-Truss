#include <stdio.h>
#include <stdlib.h>
#include <math.h>
		

typedef struct nodo {
	struct lista *inizioListaAdiacenti;
	int numAdiacenti;
} nodo_t;
typedef struct edge {
	int nodo1, nodo2;
	int Delta;
	int* G;
	int posizioneNellHeap;
	int tau;
} edge_t;
typedef struct lista {
	struct lista *next;
	int indice;
} lista_t;
typedef struct heaptree {
	edge_t** pos;
	int numElementi;
} heaptree_t;

lista_t* newElementoLista (int elemento, lista_t* next) {
	lista_t* new = malloc(sizeof(lista_t*));
	new->indice = elemento;
	new->next = next;
	return new;
}
int randomSubset (int Kmax) {
	if (random() % (100*Kmax)) return 0;
	return 1;
}
int ricercaBinaria(int ricercato, int* array, int start, int end) {
	if (end-start < 1) return 0;
	if (end-start == 1) {
		if (array[start] == ricercato)	return 1;
		else							return 0;
	}
	int centro = start + (end-start)/2;
	if (array[centro] > ricercato)	return ricercaBinaria(ricercato, array, start, centro);
	if (array[centro] < ricercato)	return ricercaBinaria(ricercato, array, centro+1, end);
	return 1;
}
int comparaArchi_Delta (const void * elem1, const void * elem2) {
	//compara due archi in base al valore di Delta, utile per sortarli alla fine
	edge_t** x1 = (edge_t**) elem1;
	edge_t** x2 = (edge_t**) elem2;
	if ((*x1)->Delta > (*x2)->Delta) return -1;
	if ((*x1)->Delta < (*x2)->Delta) return 1;
	return 0;
}
int comparaArchi_Tau (const void * elem1, const void * elem2) {
	//compara due archi in base al valore di Tau, utile per sortarli alla fine
	edge_t** x1 = (edge_t**) elem1;
	edge_t** x2 = (edge_t**) elem2;
	if ((*x1)->tau > (*x2)->tau) return -1;
	if ((*x1)->tau < (*x2)->tau) return 1;
	return 0;
}
int comparaArchi_Indice (const void * elem1, const void * elem2) {
	//compara due archi in base al valore degli indici dei nodi, utile per sortarli alla fine     x1 > x2 = -1
	edge_t** x1 = (edge_t**) elem1;
	edge_t** x2 = (edge_t**) elem2;
	if ( (*x1)->nodo1 < (*x2)->nodo1 ) return -1;
	if ( (*x1)->nodo1 > (*x2)->nodo1 ) return 1;
	if ( (*x1)->nodo2 < (*x2)->nodo2 ) return -1;
	if ( (*x1)->nodo2 > (*x2)->nodo2 ) return 1;
	return 0;
}
void matrixMultiplication (int* A, int* B, int* Result, int n, int m, int r, int dim) {
	//moltiplica la matrice A [n x m] per la matrice B [m x	r] e scrive il risultato in Result [n x r]
	int i, j, k;
	for (i=0;i<n;i++) for(j=0;j<r;j++) {
		Result[dim*i+j] = 0;
		for (k=0;k<m;k++)	Result[dim*i+j] += A[dim*i+k]*B[dim*k+j];
	}
}


//operazioni heap
int padre(int i) {
	return (i-1)/2;
}
int figlioSx(int i) {
	return 2*i+1;
}
int figlioDx(int i) {
	return 2*i+2;
}
void scambia (heaptree_t* heap, int n, int m) {
	int posProvv;
	posProvv = heap->pos[n]->posizioneNellHeap;
	heap->pos[n]->posizioneNellHeap = heap->pos[m]->posizioneNellHeap;
	heap->pos[m]->posizioneNellHeap = posProvv;
	edge_t* edgeProvv;
	edgeProvv = heap->pos[n];
	heap->pos[n] = heap->pos[m];
	heap->pos[m] = edgeProvv;
}
void aggiusta (heaptree_t* heap, int daAggiustare) {
	while (heap->pos[daAggiustare]->Delta  <  heap->pos[padre(daAggiustare)]->Delta) {
		scambia(heap, daAggiustare, padre(daAggiustare));
		daAggiustare = padre(daAggiustare);
	}
	int sx,dx;
	sx = dx = 1;
	while ( dx || sx ) {
		sx = dx = 0;
		if (figlioSx(daAggiustare) < heap->numElementi)	if (heap->pos[daAggiustare]->Delta  >  heap->pos[figlioSx(daAggiustare)]->Delta)	sx=1;
		if (figlioDx(daAggiustare) < heap->numElementi)	if (heap->pos[daAggiustare]->Delta  >  heap->pos[figlioDx(daAggiustare)]->Delta)	dx=1;
		if (sx==1 && dx==1) {
			if (heap->pos[figlioDx(daAggiustare)]->Delta  >  heap->pos[figlioSx(daAggiustare)]->Delta)		{
				scambia(heap, daAggiustare, figlioSx(daAggiustare));
				daAggiustare = figlioSx(daAggiustare);
			}
			else																							{
				scambia(heap, daAggiustare, figlioDx(daAggiustare));
				daAggiustare = figlioDx(daAggiustare);
			}
		}
		if (sx==1 && dx==0)	{
			scambia(heap, daAggiustare, figlioSx(daAggiustare));
			daAggiustare = figlioSx(daAggiustare);
		}
		if (sx==0 && dx==1) {
			scambia(heap, daAggiustare, figlioDx(daAggiustare));
			daAggiustare = figlioDx(daAggiustare);
		}
	}
}
edge_t* getMin (heaptree_t* heap) {
	if (heap->numElementi > 0)	return heap->pos[0];
	else return NULL;
}
edge_t* extractMin (heaptree_t* heap) {
	edge_t* min;
	if (heap->numElementi > 0) {
		min = heap->pos[0];
		heap->numElementi--;
		heap->pos[0] = heap->pos[heap->numElementi];
		heap->pos[0]->posizioneNellHeap = 0;
		aggiusta(heap,0);
		return min;
	}
	else return NULL;
}
void aggiungiElemento (heaptree_t* heap, edge_t* new) {
	heap->pos[heap->numElementi] = new;
	new->posizioneNellHeap = heap->numElementi;
	heap->numElementi++;
	aggiusta(heap,heap->numElementi-1);
} 




int main() {
	char NOME_FILE[50];
	printf("Inserisci il nome del file: ");
	scanf("%s",NOME_FILE);
	printf("\n");
	

	int v,u,w,xl;	//queste variabili indicheranno sempre nodi
	//scorro il file la prima volta per vedere il range in cui variano gli indici dei nodi
	int IndiceNodoMin = 1000000;
	int IndiceNodoMax = 0;
	FILE *archi = fopen(NOME_FILE,"r");
	while (	fscanf(archi, "%d %d", &u, &v)!=EOF ) {
		if (u>IndiceNodoMax) IndiceNodoMax=u;
		if (v>IndiceNodoMax) IndiceNodoMax=v;
		if (u<IndiceNodoMin) IndiceNodoMin=u;
		if (v<IndiceNodoMin) IndiceNodoMin=v;
    }
	printf("I nodi in %s vanno da %d a %d\n\n",NOME_FILE, IndiceNodoMin , IndiceNodoMax);
	
	
	
	int counter;
	int i,j;
	edge_t* e;
	lista_t* trianglesInvolvingE;
	int numArchi = 0;
	int numNodi = IndiceNodoMax - IndiceNodoMin + 2;
	nodo_t* nodes = calloc(numNodi,sizeof(nodo_t));
	edge_t** edges = calloc(numNodi*numNodi,sizeof(edge_t*));
	
	int Kmax;
	printf("Inserisci K max cercato: ");
	scanf("%d",&Kmax);
	int c;					// la costante c usata nell'algoritmo per determinare L
							// maggiore è c maggiore sarà la probabilità di ottenere un risultato corretto
	printf("Inserisci la costante C: ");
	scanf("%d",&c);
	double constB;			// b deve stare nell'intervallo [0,1-a]
							// prendo b = constB * (1-a)
	printf("Inserisci ConstB per stabilire la soglia nodi Heavy/Light ( 0<ConstB<1 ): ");
	scanf("%lf",&constB);
	printf("Hint:\tSe il preprocessing Light impiega troppo tempo, scegli una costante minore\n");
	printf("\tSe il preprocessing Heavy impiega troppo tempo, scegli una costante maggiore\n\n\n");
	
	int L = (int) c*Kmax*log((double)numNodi);
	int l;
	
	//scorro il file la seconda volta per aggiungere gli archi al mio grafo
	rewind(archi);
	while (	fscanf(archi, "%d %d", &u, &v)!=EOF ) {
		u -= IndiceNodoMin; u++;		// così non ho nodi con indice=0 che
		v -= IndiceNodoMin; v++;		// potrebbero dare problemi vista la definizione di G(e,l)
		if (edges[numNodi*u+v]==NULL) {
			nodes[u].inizioListaAdiacenti = newElementoLista(v, nodes[u].inizioListaAdiacenti);
			nodes[u].numAdiacenti++;
			nodes[v].inizioListaAdiacenti = newElementoLista(u, nodes[v].inizioListaAdiacenti);
			nodes[v].numAdiacenti++;
			edges[numNodi*u+v] = edges[numNodi*v+u] = calloc(1,sizeof(edge_t));
			if (u>v) {i=u; u=v; v=i;} //cosi nodo1 < nodo2 per ogni arco
			edges[numNodi*u+v]->nodo1 = u;
			edges[numNodi*u+v]->nodo2 = v;
			edges[numNodi*u+v]->G = calloc(L,sizeof(int));
			edges[numNodi*u+v]->tau = 1000000000;
			numArchi++;
		}
    }
	printf ("numNodi = %d, numArchi = %d\n\n",numNodi-1,numArchi);
	
	printf("L = c*Kmax*log(numNodi) \t-->\t L = %d\n",L);
	
	
	
	
	
	//--------
	//STEP  1
	//--------
	int* Xprovv = malloc(numNodi*sizeof(int));
	int* Xlenght = calloc(L,sizeof(int));				//array di cardinalità degli X[l]
	int** X = calloc(L,sizeof(int*));					//array di insiemi già ordinati
	for (l=0;l<L;l++) {
		counter=0;
		for (v=0;v<numNodi;v++)		if (randomSubset(Kmax))		Xprovv[counter++] = v;
		Xlenght[l] = counter;
		X[l] = malloc(counter*sizeof(int));
		for (i=0;i<counter;i++)		(X[l])[i] = Xprovv[i];
	}
	
	
	
	
	//--------
	//STEP  2
	//--------
	//Light
	lista_t *provv1,*provv2;
	printf("Kmax = m^a \t\t\t-->\t a = %lf\n",log((double)Kmax)/log((double)numArchi));
	printf("b>0 && b<1-a && b=constB*(1-a)\t-->\t b = %lf\n",constB*(1 - log((double)Kmax)/log((double)numArchi) ) );
	int s = (int) pow(  numArchi, constB*(1 - log((double)Kmax)/log((double)numArchi) )   );
	printf ("s = m^b \t\t\t-->\t s = %d\n",s);
	
	
	counter=0;
	for (v=0;v<numNodi;v++)	if (nodes[v].numAdiacenti <= s) counter++;
	printf("\nCi sono %d nodi light e %d nodi heavy (con più o meno di %d nodi)\n\n\n",counter,numNodi-counter,s);
	
	
	printf("Inizio preprocessing Light\n");
	counter=0;
	for (v=0;v<numNodi;v++)
		if (nodes[v].numAdiacenti <= s) {
			counter++;
			for (provv1=nodes[v].inizioListaAdiacenti; provv1!=NULL; provv1=provv1->next)	for (provv2=provv1->next; provv2!=NULL; provv2=provv2->next) {
				u = provv1->indice;
				w = provv2->indice;
				if (edges[numNodi*u+w] != NULL)
					if ((nodes[u].numAdiacenti<=s && u>v) || nodes[u].numAdiacenti>s)	if ((nodes[w].numAdiacenti<=s && w>v) || nodes[w].numAdiacenti>s) {
						edges[numNodi*v+u]->Delta++;
						edges[numNodi*v+w]->Delta++;
						edges[numNodi*u+w]->Delta++;
						for (l=0;l<L;l++) {
							if (ricercaBinaria(v,X[l],0,Xlenght[l]))	edges[numNodi*u+w]->G[l] += v;
							if (ricercaBinaria(u,X[l],0,Xlenght[l]))	edges[numNodi*v+w]->G[l] += u;
							if (ricercaBinaria(w,X[l],0,Xlenght[l]))	edges[numNodi*u+v]->G[l] += w;
						}
					}
			}
		}
	printf("Fine preprocessing Light\n\n");
	
	
	
	
	
	//Heavy
	printf("Inizio preprocessing Heavy\n");
	int* bindingHeavy = malloc(numNodi*sizeof(int));
	int numNodiHeavy;
	counter=0;
	for (v=0;v<numNodi;v++)																//Bigezione tra i nodi heavy e gli interi {1, ... , numNodiHeavy}
		if (nodes[v].numAdiacenti > s) {
			bindingHeavy[counter] = v;
			counter++;
		}
	numNodiHeavy = counter;
	int* Atrasp = malloc(numNodiHeavy*numNodiHeavy*sizeof(int));
	int* Aprimo = malloc(numNodiHeavy*numNodiHeavy*sizeof(int));
	int* Result = malloc(numNodiHeavy*numNodiHeavy*sizeof(int));
	int* bindingHeavyL = malloc(numNodiHeavy*sizeof(int));
	int numNodiHeavyL;
	
	for (l=0;l<L;l++) {		//Calcolo G(e,l)
		counter=0;
		for (i=0;i<Xlenght[l];i++) {													//Bigezione tra i nodi heavy appartenenti a X(l) e gli interi {1, ... , numNodiHeavyL}
			v = (X[l])[i];
			if (nodes[v].numAdiacenti > s) {
				bindingHeavyL[counter] = v;
				counter++;
			}
		}
		numNodiHeavyL = counter;
		for (i=0;i<numNodiHeavy;i++)	for(j=0;j<numNodiHeavyL;j++) {					
			v = bindingHeavy[i];
			w = bindingHeavyL[j];
			if (edges[numNodi*v+w] != NULL) {											//Assegno i giusti valori alle matrici A trasposta e A primo
				Atrasp[numNodiHeavy*j+i] = 1;
				Aprimo[numNodiHeavy*i+j] = w;
			}
			else {
				Atrasp[numNodiHeavy*j+i] = 0;
				Aprimo[numNodiHeavy*i+j] = 0;
			}
		}
		matrixMultiplication(Aprimo, Atrasp, Result, numNodiHeavy, numNodiHeavyL, numNodiHeavy, numNodiHeavy);
		for (i=0;i<numNodiHeavy;i++)	for(j=i+1;j<numNodiHeavy;j++) {					//Aggiorno i valori di G(e,l)
			v = bindingHeavy[i];
			u = bindingHeavy[j];
			if (edges[numNodi*v+u] != NULL) 	edges[numNodi*v+u]->G[l] += Result[numNodiHeavy*i+j];
		}
	}
	
	for (i=0;i<numNodiHeavy;i++)	for(j=0;j<numNodiHeavy;j++) { //Calcolo Delta(e)
		v = bindingHeavy[i];
		w = bindingHeavy[j];
		if (edges[numNodi*v+w] != NULL) {												//Assegno i giusti valori alle matrici A trasposta e A primo
			Atrasp[numNodiHeavy*j+i] = 1;
			Aprimo[numNodiHeavy*i+j] = 1;
		}
		else {
			Atrasp[numNodiHeavy*j+i] = 0;
			Aprimo[numNodiHeavy*i+j] = 0;
		}
	}
	matrixMultiplication(Aprimo, Atrasp, Result, numNodiHeavy, numNodiHeavy, numNodiHeavy, numNodiHeavy);
	for (i=0;i<numNodiHeavy;i++)	for(j=i+1;j<numNodiHeavy;j++) {						//Aggiorno i valori di Delta(e)
		v = bindingHeavy[i];
		u = bindingHeavy[j];
		if (edges[numNodi*v+u] != NULL) 	edges[numNodi*v+u]->Delta += Result[numNodiHeavy*i+j];
	}
	printf("Fine preprocessing Heavy\n");

	
	
	
	
	

	
	
	
	
	//metto gli archi nell'heap 	(ordinato per Delta)
	heaptree_t* heap = malloc(sizeof(heaptree_t));
	heap->pos = calloc(numArchi,sizeof(edge_t*));
	heap->numElementi = 0;
	for (v=0;v<numNodi;v++) for(u=v+1;u<numNodi;u++) if(edges[numNodi*v+u]!=NULL) {
		aggiungiElemento(heap, edges[numNodi*v+u]);
	}
	//metto gli archi in un array 	(per non perderli dopo averli rimossi dal grafo)
	edge_t** edgesArray = malloc(numArchi*sizeof(edge_t*));
	counter=0;
	for (v=0;v<numNodi;v++) for(u=v+1;u<numNodi;u++) if(edges[numNodi*v+u]!=NULL)	edgesArray[counter++] = edges[numNodi*v+u];
	if (counter!=numArchi) printf ("\n\nNON VA BENE %d != %d\n\n",counter,numArchi);
	
	//stampo i valori di Delta su file
	FILE *OutputFILEDelta = fopen("Risultati_Delta.txt","w");
	fprintf(OutputFILEDelta, "Delta\tEdge\n");
	for (i=0;i<numArchi;i++) 	fprintf(OutputFILEDelta, "%d\t:\t%d - %d\n", edgesArray[i]->Delta, edgesArray[i]->nodo1 + IndiceNodoMin-1, edgesArray[i]->nodo2 + IndiceNodoMin-1);
	
	qsort(edgesArray, numArchi, sizeof(edge_t*), comparaArchi_Delta);
	FILE *OutputFILEDeltaSorted = fopen("Risultati_Delta_Sorted.txt","w");
	fprintf(OutputFILEDeltaSorted, "Delta\tEdge\n");
	for (i=0;i<numArchi;i++)		fprintf(OutputFILEDeltaSorted, "%d\t:\t%d - %d\n", edgesArray[i]->Delta, edgesArray[i]->nodo1 + IndiceNodoMin-1, edgesArray[i]->nodo2 + IndiceNodoMin-1);



	
	//------------
	//STEP  3 e 4
	//------------
	lista_t* provv;
	int k;
	int giaInLista;
	for (k=1; k<=Kmax; k++) {
		while (getMin(heap)!=NULL && getMin(heap)->Delta < k) {
		
			//-------
			//STEP 5
			//-------
			e = extractMin(heap);
			e->tau = k-1;
		
			//-------
			//STEP 6
			//-------
			v = e->nodo1;
			u = e->nodo2;
			trianglesInvolvingE = NULL;
			for (l=0;l<L;l++) {
				xl = e->G[l];
				if (xl>=0 && xl < numNodi)		if (edges[numNodi*xl+v]!=NULL && edges[numNodi*xl+u]!=NULL) {
					giaInLista=0;
					for (provv=trianglesInvolvingE; provv!=NULL; provv=provv->next) 	if (xl == provv->indice) giaInLista=1;
					if (!giaInLista)	trianglesInvolvingE = newElementoLista(xl, trianglesInvolvingE);
				}
			}
		
			//-------
			//STEP 7
			//-------
			edges[numNodi*u+v] = edges[numNodi*v+u] = NULL;
		
			//-------
			//STEP 8
			//-------
			for (provv=trianglesInvolvingE; provv!=NULL; provv=provv->next) {
				w = provv->indice;
				edges[numNodi*v+w]->Delta--;
				aggiusta(heap, edges[numNodi*v+w]->posizioneNellHeap);
				edges[numNodi*u+w]->Delta--;
				aggiusta(heap, edges[numNodi*u+w]->posizioneNellHeap);
			}
			for (l=0;l<L;l++) {
				if (ricercaBinaria(v,X[l],0,Xlenght[l]))	for(provv=trianglesInvolvingE; provv!=NULL; provv=provv->next) {
					w = provv->indice;
					edges[numNodi*v+w]->G[l] -= u;
				}
				if (ricercaBinaria(v,X[l],0,Xlenght[l]))	for(provv=trianglesInvolvingE; provv!=NULL; provv=provv->next) {
					w = provv->indice;
					edges[numNodi*u+w]->G[l] -= v;
				}
			}
		}
	}
	//-------------
	//STEP 11 e 12
	//-------------
	/*FILE *controllo = fopen("Sottografo_di_controllo.txt","w");
	while (getMin(heap)!=NULL)		{
		e = extractMin(heap);
		fprintf(controllo,"%d %d\n",e->nodo1 + IndiceNodoMin-1,e->nodo2 + IndiceNodoMin-1);
	}*/
	
	counter=0;
	while (getMin(heap)!=NULL)		{
		extractMin(heap)->tau = k-1;
		counter++;
	}
	printf("\n\nCi sono %d archi con tau >= %d",counter,k-1);

	
	
	qsort(edgesArray, numArchi, sizeof(edge_t*), comparaArchi_Indice);
	FILE *OutputFILETau = fopen("Risultati_Tau.txt","w");
	fprintf(OutputFILETau, "Tau\t\tEdge\n");
	for (i=0;i<numArchi;i++) 	fprintf(OutputFILETau, "%d\t:\t%d - %d\n", edgesArray[i]->tau, edgesArray[i]->nodo1 + IndiceNodoMin-1, edgesArray[i]->nodo2 + IndiceNodoMin-1);
	
	qsort(edgesArray, numArchi, sizeof(edge_t*), comparaArchi_Tau);
	FILE *OutputFILETauSorted = fopen("Risultati_Tau_Sorted.txt","w");
	fprintf(OutputFILETauSorted, "Tau\t\tEdge\n");
	for (i=0;i<numArchi;i++) 	fprintf(OutputFILETauSorted, "%d\t:\t%d - %d\n", edgesArray[i]->tau, edgesArray[i]->nodo1 + IndiceNodoMin-1, edgesArray[i]->nodo2 + IndiceNodoMin-1);
	
	return 0;
}
