package assignment4exercise7;

import java.util.Arrays;

/**
 * @author lucio.yz E-mail:lucio.yang@qq.com 
 * @date 2016年11月14日 下午6:55:46
 * @version 1.0
*/
public class Simplex {
	/**
	 * simplex算法执行类
	 * @param A 原问题松弛形式的系数矩阵
	 * @param b 松弛形式Ax=b
	 * @param c 目标函数的系数
	 * @return
	 */
	private Solution simplex(double[][] A,double[] b,double[] c,boolean needInitialize,boolean[] BI_input){
		//A是m*n的矩阵
		int m=A.length;
		int n=A[0].length;
		double[] x=new double[n];//最优解
		double z=0;//最优目标函数值
		double[] delta=new double[m];
		//集合BI，代表相应位置的x是否是基
		boolean[] BI;
		if( needInitialize ){
			BI=new boolean[n];
			Arrays.fill(BI, false);
			z=initializeSimplex(A, b, c, BI);
		}
		else
			BI=BI_input;
		if( z==Double.MAX_VALUE )//无解
			return new Solution();
		while( true ){
			//找到一个<0的c
			int e=-1;
			for( int i=0;i<n;i++ )
				if( c[i]<0 ){
					e=i;//break;
				}					
			if( e==-1 ){//没有找到一个<0的c
				calculateX(A, b, c, BI, x);//计算此时的最优解
				return new Solution(-z, x);
			}
			//选择e所在的向量为入基向量
			//li是出基的行,lj是出基的列
			Arrays.fill(delta, Double.MAX_VALUE);
			int li=0;
			for( int i=0;i<m;i++ ){//找基向量变换的倍数
				if( A[i][e]>0 ){
					delta[i]=b[i]/A[i][e];
					li=(delta[li]>delta[i])?i:li;
				}
			}
			if( delta[li]==Double.MAX_VALUE )//没有找到最小倍数，意味着解空间无界
				return new Solution();
			//在基中，找到最小的delta，对应的基向量的列lj
			int lj=0;
			for( int j=0;j<n;j++ ){
				if( BI[j]==true && A[li][j]==1 ){
					lj=j;
				}
			}
			z=pivot(A,b,c,BI,e,li,lj,z);//第l列是出基，第e列是入基，进行pivot操作
		}
	}
	/**
	 * 测试一个线性规划是否可行，如果可行，找到一个解。
	 * @param A 原问题松弛形式的系数矩阵
	 * @param b 松弛形式Ax=b
	 * @param c 目标函数的系数
	 * @param BI 相应位置的向量是否为基
	 * @return
	 */
	private double initializeSimplex(double[][] A,double[] b,double[] c,boolean[] BI){
		//简单方法：假定没有负值，让松弛变量=b，非松弛变量=0，就可以构造初始解
		for( int i=BI.length-1;i>=BI.length-b.length;i-- )
			BI[i]=true;
		//找到b中最小的，bi保存最小值位置
		int bi=0;
		for( int i=1;i<b.length;i++ )
			bi=(b[bi]>b[i])?i:bi;
		if( b[bi]>=0 )
			return 0;//因为所有的非松弛变量均取0，所以目标函数解为0
		//存在b<0，则不能通过以上的简单方法得到一个基本可行解
		//构造辅助线性规划，求解初始解和目标;新增的松弛变量x0放在第一个元素。
		int m=A.length;
		int n=A[0].length;
		double[] c_aux=new double[n+1];
		Arrays.fill(c_aux, 0);
		c_aux[0]=1;//目标函数为x0
		double[] b_aux=new double[m];
		double[][] A_aux=new double[m][n+1];
		for( int i=0;i<m;i++ ){
			for( int j=0;j<n+1;j++ ){
				if( j==0 ){//x0的系数全为-1
					A_aux[i][j]=-1;
				}else
					A_aux[i][j]=A[i][j-1];
			}
			b_aux[i]=b[i];//b的值没有变化
		}				
		boolean[] BI_aux=new boolean[n+1];
		BI_aux[0]=false;//x0不是基向量
		for( int j=1;j<n+1;j++ )
			BI_aux[j]=BI[j-1];
		//找到最小的b对应的基，作为出基
		int bj=0;
		for( int j=0;j<BI_aux.length;j++ ){
			if( BI_aux[j]==true && A_aux[bi][j]==1 ){
				bj=j;
			}
		}
		double z=0;
		//使得所有的b为正
		z=pivot(A_aux,b_aux,c_aux,BI_aux,0,bi,bj,z);//b为负的基向量作为出基，x0的向量作为入基，进行pivot操作
		Solution solution=simplex(A_aux, b_aux, c_aux,false,BI_aux);
		//之前的一次pivot有z;simplex算法对z在pivot中每次减去一个量
		solution.z=z-solution.z;
		if( solution.z==0 ){//x0=0,找到了一个可行解
			//移除x0，用aux的simplex表格作为初始化表格,c不变
			for( int j=1;j<n+1;j++ ){
				BI[j-1]=BI_aux[j];
			}
			for( int i=0;i<m;i++ ){
				 for( int j=1;j<n+1;j++ ){
					 A[i][j-1]=A_aux[i][j];
				 }
				 b[i]=b_aux[i];
			}
			z=solution.z;
			//确保所有基向量的c为0，使用的是pivot中调整c的方法
			for( int e=0;e<n;e++ ){
				if( BI[e]==true && c[e]!=0 ){
					int li=0;
					//基向量对应的行
					for( int i=0;i<m;i++ ){
						if( A[i][e]==1 )
							li=i;
					}
					//c减去第l行
					z=z-b[li]*c[e];
					double ce=c[e];//这个值是变换前c[e]的值，需要保留
					for( int j=0;j<n;j++ )
						c[j]=c[j]-ce*A[li][j];
					
				}
			}
			return z;
		}else{//无解
			return Double.MAX_VALUE; 
		}
	}
	/**
	 * 计算当前的解
	 * @param A
	 * @param b
	 * @param c
	 * @param BI
	 * @param x
	 */
	private void calculateX(double[][] A,double[] b,double[] c,boolean[] BI,double[] x){
		//基所在的变量指定为相应的b，非基的变量指定为0
		for( int j=0;j<BI.length;j++ ){
			if( BI[j]==false ){
				x[j]=0;
			}else{
				for( int i=0;i<A.length;i++ ){
					if( A[i][j]==1 )
						x[j]=b[i];
				}
			}
		}
	}
	/**
	 * pivot操作
	 * @param A 松弛形式的系数
	 * @param b 松弛形式Ax=b
	 * @param c 目标函数系数
	 * @param BI 基集合
	 * @param e 入基列数
	 * @param li 出基行数
	 * @param lj 出基列数
	 * @param z 目标函数值
	 * @return
	 */
	private double pivot(double[][] A,double[] b,double[] c,boolean[] BI,int e,int li,int lj,double z){
		//A是m*n的矩阵
		int m=A.length;
		int n=A[0].length;
		b[li]=b[li]/A[li][e];
		//对第l行按Ale缩小
		double Alie=A[li][e];
		for( int j=0;j<n;j++ )
			A[li][j]=A[li][j]/Alie;
		//初等行变换，使得第e行的值为第l行的值
		for( int i=0;i<m;i++ ){
			if( i!=li ){
				double Aie=A[i][e];//这个值是变换前A[i][e]的值，需要保留
				b[i]=b[i]-Aie*b[li];
				for( int j=0;j<n;j++ )
					A[i][j]=A[i][j]-Aie*A[li][j];
			}				
		}
		//c减去第l行
		z=z-b[li]*c[e];
		double ce=c[e];//这个值是变换前c[e]的值，需要保留
		for( int j=0;j<n;j++ )
			c[j]=c[j]-ce*A[li][j];
		//重计算BI
		BI[lj]=false;
		BI[e]=true;
		return z;
	}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		Simplex s=new Simplex();
		Solution solution;
		//非松弛变量在前，松弛变量在后
		double[][] A={{1,1,1,1,0,0,0},{1,0,0,0,1,0,0},{0,0,1,0,0,1,0},{0,3,1,0,0,0,1}};
		double[] b={4,2,3,6};
		double[] c={-1,-14,-6,0,0,0,0};
		solution=s.simplex(A,b,c,true,null);
		solution.print();
		
		//没有松弛解
//		double[][] A={{-1,-1,1,0},{1,1,0,1}};
//		double[] b={-2,1};
//		double[] c={1,2,0,0};
//		solution=s.simplex(A,b,c,true,null);
//		solution.print();
		//有松弛解
//		double[][] A={{-1,-1,1,0},{1,1,0,1}};
//		double[] b={-1,2};
//		double[] c={1,2,0,0};
//		solution=s.simplex(A,b,c,true,null);
//		solution.print();
	}
	/*
	 * 辅助类，用来表示存储最优解和最优目标函数值
	 */
	private class Solution{
		double z;
		double[] x;
		boolean solved;
		public Solution(){
			this.solved=false;
		}
		public Solution(double z,double[] x){
			this.z=z;
			this.x=x;
			this.solved=true;
		}
		public void print(){
			if( solved ){
				System.out.println("找到最优解！");
				System.out.print("最优解为：");
				for( int i=0;i<x.length;i++ ){
					System.out.print(x[i]+" ");
				}
				System.out.println();
				System.out.println("最优目标函数值为："+z);
			}else{
				System.out.println("没有最优解！");
			}
		}
	}
}
