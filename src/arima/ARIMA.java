package arima;


import arima.ARMAMath;

import java.util.*;


public class ARIMA {

	double[] originalData={};
	ARMAMath armamath=new ARMAMath();
	double stderrDara=0;
	double avgsumData=0;
	Vector<double[]> armaARMAcoe=new Vector<double[]>();
	Vector<double[]> bestarmaARMAcoe=new Vector<double[]>();
	
/**
 * 构造函数
 * @param originalData 原始时间序列数据
 */
	public ARIMA(double [] originalData)
	{
		this.originalData=originalData;
	}
/**
 * 原始数据标准化处理：一阶季节性差分
 * @return 差分过后的数据
 */ 
	public double[] preDealDif()
	{
		
		//seasonal Difference:Peroid=7
		double []tempData=new double[originalData.length-7];
		for(int i=0;i<originalData.length-7;i++)
		{
			tempData[i]=originalData[i+7]-originalData[i];
		}

		return tempData;
	}
/**
 * 原始数据标准化处理：Z-Score归一化
 * @param 待处理数据
 * @return 归一化过后的数据
 */
	public double[] preDealNor(double[] tempData)
	{
		//Z-Score
		avgsumData=armamath.avgData(tempData);
		stderrDara=armamath.stderrData(tempData);
		
		for(int i=0;i<tempData.length;i++)
		{
			tempData[i]=(tempData[i]-avgsumData)/stderrDara;
		}
		
		return tempData;
	}
/**
* 得到ARMA模型=[p,q]
 * @return ARMA模型的阶数信息
 */
	public int[] getARIMAmodel()
	{
		double[] stdoriginalData=this.preDealDif();//原始数据差分处理
		
		int paraType=0;
		double minAIC=9999999;
		int bestModelindex=0;
		int[][] model=new int[][]{{0,1},{1,0},{1,1},{0,2},{2,0},{2,2},{1,2},{2,1}};//,{3,0},{0,3},{3,1},{1,3},{3,2},{2,3},{3,3}};//,{4,0},{0,4},{4,1},{1,4},{4,2},{2,4},{4,3},{3,4},{4,4}};
		//对8种模型进行迭代，选出AIC值最小的模型作为我们的模型
		for(int i=0;i<model.length;i++)
		{
			if(model[i][0]==0)
			{
				MA ma=new MA(stdoriginalData, model[i][1]);
				armaARMAcoe=ma.MAmodel(); //拿到ma模型的参数
				paraType=1;
			}
			else if(model[i][1]==0)
			{
				AR ar=new AR(stdoriginalData, model[i][0]);
				armaARMAcoe=ar.ARmodel(); //拿到ar模型的参数
				paraType=2;
			}
			else
			{
				ARMA arma=new ARMA(stdoriginalData, model[i][0], model[i][1]);
				armaARMAcoe=arma.ARMAmodel();//拿到arma模型的参数
				paraType=3;
			}
			

			double temp=getmodelAIC(armaARMAcoe,stdoriginalData,paraType);
			System.out.println("AIC of these model="+temp);
			if (temp<minAIC)
			{
				bestModelindex=i;
				minAIC=temp;
				bestarmaARMAcoe=armaARMAcoe;
			}
		}
		
		return model[bestModelindex];
 	}
/**
 * 计算ARMA模型的AIC
 * @param para 装载模型的参数信息
 * @param stdoriginalData   预处理过后的原始数据
 * @param type 1：MA；2：AR；3：ARMA
 * @return 模型的AIC值
 */
	public double getmodelAIC(Vector<double[]> para,double[] stdoriginalData,int type)
	{
		double temp=0;
		double temp2=0;
		double sumerr=0;
		int p=0;//ar1,ar2,...,sig2
		int q=0;//sig2,ma1,ma2...
		int n=stdoriginalData.length;
		Random random=new Random();
		
		if(type==1)
		{
			double[] maPara=para.get(0);
			q=maPara.length;
			double[] err=new double[q];  //error(t),error(t-1),error(t-2)...
			for(int k=q-1;k<n;k++)
			{
				temp=0;
				
				for(int i=1;i<q;i++)
				{
					temp+=maPara[i]*err[i];
				}
			
				//产生各个时刻的噪声
				for(int j=q-1;j>0;j--)
				{
					err[j]=err[j-1];
				}
				err[0]=random.nextGaussian()*Math.sqrt(maPara[0]);
				
				//估计的方差之和
				sumerr+=(stdoriginalData[k]-(temp))*(stdoriginalData[k]-(temp));
				
			}
			//return  (n-(q-1))*Math.log(sumerr/(n-(q-1)))+(q)*Math.log(n-(q-1));//AIC 最小二乘估计
			return (n-(q-1))*Math.log(sumerr/(n-(q-1)))+(q+1)*2;
		}
		else if(type==2)
		{
			double[] arPara=para.get(0);
			p=arPara.length;
			for(int k=p-1;k<n;k++)
			{
				temp=0;
				for(int i=0;i<p-1;i++)
				{
					temp+=arPara[i]*stdoriginalData[k-i-1];
				}
				//估计的方差之和
				sumerr+=(stdoriginalData[k]-temp)*(stdoriginalData[k]-temp);
			}
			return (n-(q-1))*Math.log(sumerr/(n-(q-1)))+(p+1)*2;
			//return (n-(p-1))*Math.log(sumerr/(n-(p-1)))+(p)*Math.log(n-(p-1));//AIC 最小二乘估计
		}
		else
		{
			double[] arPara=para.get(0);
			double[] maPara=para.get(1);
			p=arPara.length;
			q=maPara.length;
			double[] err=new double[q];  //error(t),error(t-1),error(t-2)...
			
			for(int k=p-1;k<n;k++)
			{
				temp=0;
				temp2=0;
				for(int i=0;i<p-1;i++)
				{
					temp+=arPara[i]*stdoriginalData[k-i-1];
				}
			
				for(int i=1;i<q;i++)
				{
					temp2+=maPara[i]*err[i];
				}
			
				//产生各个时刻的噪声
				for(int j=q-1;j>0;j--)
				{
					err[j]=err[j-1];
				}
				//System.out.println("predictBeforeDiff="+1);
				err[0]=random.nextGaussian()*Math.sqrt(maPara[0]);
				//估计的方差之和
				sumerr+=(stdoriginalData[k]-(temp2+temp))*(stdoriginalData[k]-(temp2+temp));
			}
			return (n-(q-1))*Math.log(sumerr/(n-(q-1)))+(p+q)*2;
			//return (n-(p-1))*Math.log(sumerr/(n-(p-1)))+(p+q-1)*Math.log(n-(p-1));//AIC 最小二乘估计
		}
	}
/**
 * 对预测值进行反差分处理
 * @param predictValue 预测的值
 * @return 反差分过后的预测值
 */
	public int aftDeal(int predictValue)
	{
		//System.out.println("predictBeforeDiff="+predictValue);
		return (int)(predictValue+originalData[originalData.length-7]);
	}
/**
 * 进行一步预测
 * @param p ARMA模型的AR的阶数
 * @param q ARMA模型的MA的阶数
 * @return 预测值
 */
	public int predictValue(int p,int q)
	{
		int predict=0;
		double[] stdoriginalData=this.preDealDif();
		int n=stdoriginalData.length;
		double temp=0,temp2=0;
		double[] err=new double[q+1];
	
		Random random=new Random();
		if(p==0)
		{
			double[] maPara=bestarmaARMAcoe.get(0);
			for(int k=q;k<n;k++)
			{
				temp=0;
				for(int i=1;i<=q;i++)
				{
					temp+=maPara[i]*err[i];
				}
				//产生各个时刻的噪声
				for(int j=q;j>0;j--)
				{
					err[j]=err[j-1];
				}
				err[0]=random.nextGaussian()*Math.sqrt(maPara[0]);
			}
			predict=(int)(temp); //产生预测
		}
		else if(q==0)
		{
			double[] arPara=bestarmaARMAcoe.get(0);
			for(int k=p;k<n;k++)
			{
				temp=0;
				for(int i=0;i<p;i++)
				{
					temp+=arPara[i]*stdoriginalData[k-i-1];
				}
			}
			predict=(int)(temp);
		}
		else
		{

			double[] arPara=bestarmaARMAcoe.get(0);
			double[] maPara=bestarmaARMAcoe.get(1);
			err=new double[q+1];  //error(t),error(t-1),error(t-2)...
			for(int k=p;k<n;k++)
			{
				temp=0;
				temp2=0;
				for(int i=0;i<p;i++)
				{
					temp+=arPara[i]*stdoriginalData[k-i-1];
				}
			
				for(int i=1;i<=q;i++)
				{
					temp2+=maPara[i]*err[i];
				}
			
				//产生各个时刻的噪声
				for(int j=q;j>0;j--)
				{
					err[j]=err[j-1];
				}
				
				err[0]=random.nextGaussian()*Math.sqrt(maPara[0]);
			}
			
			predict=(int)(temp2+temp);
			
		}
	
		
		return predict;
	}
/**
 * 计算MA模型的参数
 * @param autocorData 自相关系数Grma
 * @param q MA模型的阶数
 * @return 返回MA模型的参数
 */
	public double[] getMApara(double[] autocorData,int q)
	{
		double[] maPara=new double[q+1];//第一个存放噪声参数，后面q个存放ma参数sigma2,ma1,ma2...
		double[] tempmaPara=maPara;
		double temp=0;
		boolean iterationFlag=true;
		//解方程组
		//迭代法解方程组
		System.out.println("autocorData[0]"+autocorData[0]);
		while(iterationFlag)
		{
			for(int i=1;i<maPara.length;i++)
			{
				temp+=maPara[i]*maPara[i];
			}
			tempmaPara[0]=autocorData[0]/(1+temp);
		
			for(int i=1;i<maPara.length;i++)
			{
				temp=0;
				for(int j=1;j<maPara.length-i;j++)
				{
					temp+=maPara[j]*maPara[j+i];
				}
				tempmaPara[i]=-(autocorData[i]/tempmaPara[0]-temp);
			}
			iterationFlag=false;
			for(int i=0;i<maPara.length;i++)
			{
				if(maPara[i]!=tempmaPara[i])
				{
					iterationFlag=true;
					break;
				}
			}
			
			maPara=tempmaPara;
		}
		
		return maPara;
	}

}
