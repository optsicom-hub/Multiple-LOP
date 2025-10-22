#include <fstream>
#include <string>
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <math.h> 
#include <vector>
#include <algorithm>
#include <random>



using namespace std; 	
int N=0;	
int P=15;		
		
long** matriz;




//pos_i Posicion donde voy
//val_j Valor que quiero cambiar

void CargarMatriz(ifstream& f,int fila, int columna,long**& matriz);
int  eva_move(vector<int> sol,int val_sol, long**& matriz,int pos_i,int pos_j);
void make_move(vector<int>& sol,int pos_i, int pos_j);

int greedy(int pos_i,long**& matriz );
void make_sol(vector<int>& sol,vector<int> val, long**& matriz, float alpha);
int random_in(int a, int b);
void local_search(vector<int>& sol,int& val_act,long**& matriz );
int eva_sol(vector<int> sol, long**& matriz );
void normalize(long**& matriz);
int hash(vector<int> sol);
bool is_equal(vector<int> sol1, vector<int> sol2);
int kendall_tau(vector<int> sol1, vector<int> sol2);
vector<int>  shake(vector<int> inputVector, int k, int& val_act, long**& matriz);
int eva_dispersion(vector<vector<int>> max_dispersion, vector<vector<int>> kendall_matrix );
bool inclusion_gain(vector<int> sol,vector<vector<int>>& max_dispersion, vector<vector<int>>& kendall_matrix, int& best_MDP);
bool are_equal(vector<int> sol, vector<vector<int>> max_dispersion);
float linear_degree( long**& matriz, int best_value);
vector<int> max_range(vector<vector<int>> max_dispersion);
int main()
{
	int num_arch=21;
 	double time_elapsed=0;
 	double time_elapsed_wide=0;
 	double wide_time_limit=10;
	ofstream final_results;
	final_results.open("Results_C++_VNS_k15_largeOpt");
	string name[num_arch]={"N-be75eec", "N-be75oi", "N-stabu70", "N-t59b11xx", "N-t59d11xx", "N-t59f11xx", "N-t59i11xx", "N-t59n11xx", "N-t65b11xx", "N-t65d11xx", "N-t65l11xx", "N-t65n11xx", "N-t69r11xx", "N-t70b11xx", "N-t70k11xx", "N-t70l11xx", "N-t70n11xx", "N-t70u11xx", "N-t75k11xx", "N-t75n11xx", "N-t75u11xx"};
	for(int iter_arch=0; iter_arch<num_arch; iter_arch++)
	{
	time_elapsed=0;
	time_elapsed_wide=0;
	clock_t start=clock();
	srand(time(NULL));
	
	int num_solutions=0;
	int val_act=0;
	int best_val=0;
	int MaxIterOpt=5000;
	int MaxIterMDP=10000000;
	
	
	float threshold=0;
	//Leer Matriz
	ifstream f;
	f.open(name[iter_arch].c_str());

	if (!f){
		cout << "El fichero  no se ha podido abrir";
	}
	if (f)
	{
		f>>N;	
	}
	matriz= new long*[N];
	vector<int> sol;
	sol.resize(N);
	vector<int> best_sol;
	best_sol.resize(N);
	vector<int> sol_shaken;
	sol_shaken.resize(N);
	CargarMatriz(f,N,N,matriz);
	normalize(matriz);
	int kmax=N/2;
	f.close();	
	

	
	//Calculo de la funcion greedy, maximos y minimos
	float alpha=0.2;
	vector<int> greedy_val;
	greedy_val.resize(N);
	for(int i=0; i<N ; i++)
	{
		
		greedy_val[i]=greedy(i,matriz);
	}
	 //Busqueda mejor valor
	cout<<"mejor valor"<<endl;
	int iter=0;
	int k=0;
	make_sol( sol, greedy_val, matriz,  alpha);
	val_act=eva_sol(sol,matriz);
	while( iter<MaxIterOpt )
	{

		
		k=0;
		
		while(k<kmax)		
		{
			
			
			sol_shaken=shake(sol,k,val_act,matriz);
			val_act=eva_sol(sol_shaken,matriz);
			local_search(sol_shaken,val_act,matriz);
			
			
			if(val_act>best_val)
			{
				best_val=val_act;
				best_sol=sol_shaken;
				sol=sol_shaken;
				iter=0;
				k=0;
				//cout<< " Con un valor de "<<best_val<< " Calculado: "<<eva_sol(sol,matriz)<<endl ;
				/*
				for(int i=0; i<N; i++)
				{
					cout<<sol[i]<<" ";
				}
				cout <<endl;
				*/
			}
			else
			{
				k++;//val_act=val_backup;
				iter++;
			}
			
		}

	}
	cout<< best_val<<endl;
	//Busqueda de soluciones mas alejadas
	
	cout<<"Mas diverso"<<endl;
	
	int best_MDP=0;
	
	
	vector<vector<int>> max_dispersion(P, vector<int>(N));
	vector<vector<int>> max_dispersion_shaken(P, vector<int>(N));
	vector<vector<int>> kendall_matrix(P, vector<int>(P+1));
	
	//Initialize solutions.
	for(int i=0; i <P ;i++)
	{
		
		max_dispersion[i]=best_sol;
		
		for(int j=0; j<P ; j++)
		{
			kendall_matrix[i][j]=0;
		}
		
		kendall_matrix[i][P]=0;
	}

	iter=0;

	make_sol( sol, greedy_val, matriz,  alpha);
	
	val_act=eva_sol(sol,matriz);
	
	int best_val_2=0;
	vector<int> best_sol_2;
	bool flag_next=false;
	bool flag_far=false;
	int rnd;
	int par=1;
	while(iter<MaxIterMDP && time_elapsed_wide<wide_time_limit)
	{
		
		clock_t start_wide=clock();
		k=0;
		//par=par*(-1);
		//make_sol( sol, greedy_val, matriz,  alpha);
		//(par>0)
		//{
			rnd=random_in(0,P-1);
			sol=max_dispersion[rnd];	
		//}

		//best_val_2=0;
	
		
		while(k<kmax)
		{
			//cout<<k<<endl;
			sol_shaken=shake(sol,k,val_act,matriz);
			val_act=eva_sol(sol_shaken,matriz);
			local_search(sol_shaken,val_act,matriz);
			//cout<<"val act "<<val_act<<endl;
			if(abs(val_act-best_val)<=threshold*best_val)
			{
				//cout<<"optimo"<<endl;
				if(!are_equal(sol_shaken,max_dispersion))
				{
					//cout<<"diff"<<endl;
					flag_next=true;
					num_solutions++;
					if(inclusion_gain(sol_shaken,max_dispersion, kendall_matrix,best_MDP))
					{
						flag_far=true;

						//cout<<"far: "<< best_MDP<<endl;

						
					}
				}
					
			}
			if(val_act>best_val_2  )
			{
				
				best_val_2=val_act;
				best_sol_2=sol_shaken;
				sol=sol_shaken;
				iter=0;
				k=0;

				
				//cout<< " Con un valor de "<<best_val_2<< " Calculado: "<<eva_sol(sol,matriz)<<endl ;
				
				
				
			}
			/*
			else if(flag_next && !flag_far)
			{	

					
					iter++;
					k++;
					flag_next=false;
					//make_sol( sol, greedy_val, matriz,  alpha);
					//val_act=eva_sol(sol,matriz);

					//best_val_2=0;
					
					
			}
			*/
			else if(flag_far)
			{
				//make_sol( sol, greedy_val, matriz,  alpha);
				//val_act=eva_sol(sol,matriz);
				flag_far=false;
				flag_next=false;
			
				
				k=0;
				iter++;
				sol=sol_shaken;
				//best_val_2=0;
	
			}
			
			else
			{
				k++;
				iter++;	
			}	
		}
		
		clock_t end_wide = clock();
		time_elapsed_wide +=  double(end_wide - start_wide)/CLOCKS_PER_SEC;
		
	}
	
	

		clock_t end = clock();
		cout<<"far: "<<best_MDP<<endl;
		
		for(int i=0; i<P ; i++)
		{
			for(int j=0; j<N ;j++)
			{
				cout<<max_dispersion[i][j]+1<<" ";
			}
			cout<<endl;
		}
		cout<<endl;
					
		time_elapsed +=  double(end - start)/CLOCKS_PER_SEC;
		//cout<< "LD: "<<linear_degree(matriz,best_val)<<endl;
		//final_results<< name[iter_arch].c_str() << " , " <<best_val<<" , " <<num_solutions<< " , "<< best_MDP<<", " <<linear_degree(matriz,best_val) <<" , "<<time_elapsed<<endl;
		final_results<< "[ ";
		vector<int> res=max_range(max_dispersion);
		for(int i=0; i<N ; i++)
		{
			final_results<<res[i]<<". ";
			cout<<res[i]<<" ";
		}
		final_results<< " ]"<<endl;
	}
	
	final_results.close();
	

	
	

}

int  eva_move(vector<int> sol,int val_sol, long**& matriz,int pos_i,int pos_j)
{
	int new_val;
	int diff=0;
	if(pos_i<pos_j)
	{
		for(int i=pos_i; i<pos_j; i++)
		{
			diff=diff+(*(*(matriz+sol[pos_j])+sol[i])-*(*(matriz+sol[i])+sol[pos_j]));
		}
	}
	else
	{
		for(int i=pos_j+1; i<=pos_i; i++)
		{
			diff=diff-(*(*(matriz+sol[pos_j])+sol[i])-*(*(matriz+sol[i])+sol[pos_j]));
		}
	}
	new_val=val_sol+diff;
	return new_val;
}
void make_move(vector<int>& sol,int pos_i, int pos_j)
{
	
	if(pos_i<pos_j)
	{
		sol.insert(sol.begin()+pos_i,sol[pos_j]);
		sol.erase(sol.begin()+pos_j+1);
	}
	else
	{
		sol.insert(sol.begin()+pos_i+1,sol[pos_j]);
		sol.erase(sol.begin()+pos_j);
	}

	
	return ;		
}
int greedy(int pos_i,long**& matriz )
{
	int suma=1;
	for(int i=0; i<N ; i++)
	{
		suma=suma+ *(*(matriz+pos_i)+i);
	}
	return suma;
}		
void make_sol(vector<int>& sol,vector<int> val,long**& matriz, float alpha)
{
	//cout<<"hola"<<endl;
	vector<int> RCL;

	int max;
	float threshold=0;
	int entra;
	
	for(int i=0 ; i<N; i++)
	{	
		max=0;
		for(int j=0; j<N; j++)
		{
			if(val[j]>max)
			{
				max=val[j];	
				//cout<<"Max: "<<max<<endl; 
			}
		}
		
		threshold=alpha*max;
		
		for(int j=0; j<N ; j++)
		{
			if(val[j]>threshold)
			{
				
				RCL.insert(RCL.begin(),j);
			}
		}
		entra=random_in(0,RCL.size()-1);
		sol[i]=RCL[entra];
		val[sol[i]]=0;	
		RCL.clear();	
	}
		
}
int random_in(int a, int b)
{

	int random=1;
	if(b>=0)
	{
		int range=b-a+1;
		random=  a + (rand() % range);	
	}

	return random;
}
void CargarMatriz(ifstream& f,int fila, int columna,long**& matriz)
{
	
	long valor;
	double basura;
	
	matriz = new long*[fila];
	for(int c=0; c<fila; c++ )
	{
		matriz[c]= new long[columna];
	}
	
	
	
	
	for (int i=0; i<N ; i++ )
	{
		for(int j=0; j<N; j++)
		{
			f>>valor;
			*(*(matriz+i)+j)=valor;
		}
	}
}
void local_search(vector<int>& sol,int& val_act,long**& matriz )
{
	int iter_1=0;
	int iter_2=0;
	int val_candidate=0;
	int best_row=0;
	int best_iter_1=0;
	int best_iter_2=0;

	bool flag_2=true;
	bool flag_1=true;
	while(flag_1)
	{
		flag_2=true;
		iter_2=0;
		while(iter_2<N )
		{			
			
			val_candidate=eva_move(sol,val_act,matriz,iter_2,iter_1);
			if(val_candidate>best_row)
			{
				best_row=val_candidate;
				best_iter_1=iter_1;
				best_iter_2=iter_2;
			}
			iter_2++;
		}
		if(best_row>val_act)
		{
			val_act=best_row;
				//cout<<"cmabia valor: "<<val_act<<endl;
				//cout<< "Iter1: "<<iter_1 <<" iter_2: "<<iter_2<<endl;
			make_move(sol,best_iter_2,best_iter_1);
			flag_2=false;
		}
		best_row=0;
		val_candidate=0;
		iter_1++;
		
		if(iter_1==N && flag_2==true)
		{
			flag_1=false;
		}
		else if(iter_1==N)
		{
			iter_1=0;
		}
	}
	
}
int eva_sol(vector<int> sol, long**& matriz )
{
	int suma=0;
	for(int i=0; i<N; i++)
	{
		for(int j=i+1; j<N; j++)
		{
			suma=suma+*(*(matriz+sol[i])+sol[j]);
		}
	}
	return suma;
}
void normalize(long**& matriz)
{
	for(int i=0 ; i<N; i++)
	{
		for(int j=i+1; j<N; j++)
		{
			if(*(*(matriz+i)+j)> *(*(matriz+j)+i))
			{
				*(*(matriz+i)+j)=*(*(matriz+i)+j)-*(*(matriz+j)+i);
				*(*(matriz+j)+i)=0;
			}
			else
			{
				*(*(matriz+j)+i)=*(*(matriz+j)+i)-*(*(matriz+i)+j);
				*(*(matriz+i)+j)=0;
			}
		}
		*(*(matriz+i)+i)=0;
	}	
}
int hash(vector<int> sol)
{
	int suma=0;
	for(int i=0; i<N; i++)
	{
		suma=suma+ sol[i]*i;
	}
	return suma;
}
bool is_equal(vector<int> sol1, vector<int> sol2)
{
	bool flag=true;
	int iter=0;
	int size=sol1.size();
	while(flag && iter<size)
	{
		
		if(sol1[iter]!=sol2[iter])
		{
			flag=false;
		}
		iter++;
	}
	return flag;
}
int kendall_tau(vector<int> arr1, vector<int> arr2)
{
	int distancia = 0;
    int n = arr1.size();
	vector<int> index;
	index.resize(n);
	int element=0;
	for(int i=0; i<n ; i++)
	{
		element=arr2[i];
		index[element]=i;
	}
    for (int i = 0; i < n - 1; ++i) 
	{
    	for (int j = i + 1; j < n; ++j)
		{
            if(index[arr1[i]]>index[arr1[j]])
            {
            	distancia++;
			}
        }
    }


    return distancia;
}
vector<int> shake(vector<int> inputVector, int k, int& val_act, long**& matriz) 
{
	int pos_i;
	int pos_j;
	int flip;
	if(k==1)
	{
		pos_i=random_in(1,N-2);
		flip=random_in(0,1);
		if(flip==0)
		{
			pos_j=pos_i+1;	
		}
		else
		{
			pos_j=pos_i-1;
		}
		//val_act=eva_move(inputVector,val_act,matriz,pos_i,pos_j);
		make_move(inputVector,pos_i,pos_j);
		
	}
	else
	{
		for(int i=1; i<k ; i++)
		{
			pos_i=random_in(0,N-1);
			pos_j=random_in(0,N-1);
			//val_act=eva_move(inputVector,val_act,matriz,pos_i,pos_j);
			make_move(inputVector,pos_i,pos_j);
		}	
	}
	
	

    return inputVector;
}
int eva_dispersion(vector<vector<int>> max_dispersion, vector<vector<int>> kendall_matrix )
{
	int suma=0;
	int dist=0;
	for(int i=0; i<P ; i++)
	{
		for(int j=0; j<P ; j++)
		{
			dist=kendall_tau(max_dispersion[i],max_dispersion[j]);
			suma=suma+dist;
			kendall_matrix[i][j]=dist;
			
		}
	}	
	return suma;	
}
bool inclusion_gain(vector<int> sol,vector<vector<int>>& max_dispersion, vector<vector<int>>& kendall_matrix, int& best_MDP)
{
	int best=0;
	int best_where=0;
	int act=0; 
	int dist=0;
	bool flag=false;
	
	vector<int> new_dist(P);
	
	for(int i=0 ; i<P ;i++)
	{
		new_dist[i]=kendall_tau(sol, max_dispersion[i]);
	}

	for(int i=0; i<P ; i++)
	{
		act=0;
		dist=0;
		for(int j=0; j<P ; j++)
		{
			if(j!=i)
			{
				dist=kendall_tau(sol,max_dispersion[j]);
				
				act=act+ dist;	
			}	
		}
		
		
		if((act-kendall_matrix[i][P])>best)
		{
			
			best=act-kendall_matrix[i][P];
			best_where=i;
		}	
		
		
	}
	if( best>0)
	{
		flag=true;
		best_MDP=best_MDP+best;
		max_dispersion[best_where]=sol;
		for(int i=0 ; i< P; i++)
		{
			if(i!=best_where)
			{
				kendall_matrix[best_where][i]=new_dist[i];
				kendall_matrix[i][best_where]=new_dist[i];
			}

		}
		for(int i=0; i<P; i++)
		{
			kendall_matrix[i][P]=0;
			for(int j=0; j<P; j++)
			{
				kendall_matrix[i][P]=kendall_matrix[i][P]+kendall_matrix[i][j];
			}
		}
		
	}
	return flag;
}
bool are_equal(vector<int> sol, vector<vector<int>> max_dispersion)
{
	bool flag=false;
	int i=0;
	while(!flag && i<P)
	{
		if(is_equal(sol,max_dispersion[i]))
		{
			flag=true;
			
		}
		i++;	
	}
	return flag;
}
float linear_degree( long**& matriz, int best_value)
{
	float suma=0;
	for(int i=0; i<N; i++)
	{
		for(int j=0; j<N ;j++)
		{
			suma= suma + *(*(matriz+j)+i);
		}
	}
	float res=best_value / suma;

	return res;
}
vector<int> max_range(vector<vector<int>> max_dispersion)
{
	vector<int> max_range(N, 0);
	vector<int> min_range(N, N);
	vector<int> res(N, 0);
	
	for(int i=0; i<P ; i++)
	{
		for(int j=0 ; j< N ; j++ )
		{
			if(j>max_range[max_dispersion[i][j]])
			{
				max_range[max_dispersion[i][j]]=j;
			}
			if(j<min_range[max_dispersion[i][j]])
			{
				min_range[max_dispersion[i][j]]=j;
			}
		}
	}
	for(int i=0 ; i<N; i++)
	{
		res[i]=max_range[i]-min_range[i];
	}
	return res;
}
