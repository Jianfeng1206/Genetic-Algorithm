
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace GA
{
    class GeneticAlgorithm
    {
        // to create the Genetic algorithms
        private int popsize = 20;  //population size 种群数目

		
		
        public int Popsize { get => popsize; set => popsize = value; }
       // 这些都是启发式的算法
        private double crossOver = 0.8;
        private double mutationRate = 0.03; // mutation probability 

        //  private int generations = 20;
        // suppose that the variable range is 0-5 
        public int maxbound = 5;
        private int DNAsize = 10;

        // 其实种群的数目也不必过多
        double min = 0;


        public individual[] InitialPopulation(individual[] Population)
        {
            Population = new individual[popsize];
            //但是每个种群的 DNA 编码还是没有进行编码
            for (int i = 0; i < Popsize; i++)
            {
                // uq create the object 
                Population[i] = new individual();
                //构造一个数据去保存随机数
                int[] num = new int[DNAsize];
                for (int j = 0; j < DNAsize; j++)
                {
                  
                    // 解决方案
                   byte[] buffer = Guid.NewGuid().ToByteArray(); int iSeed = BitConverter.ToInt32(buffer, 0);
                    Random random = new Random(iSeed);
                  int mm= random.Next(0,2);
                    num[j] = mm;

                }

                // 去赋值,这是个引用并非赋值
                Array.Copy(num, Population[i].num,DNAsize);
              //  Population[i].num = num;
                Population[i].fitness = 0;
                Population[i].selectP = 0;
            }
         return Population;
        }
       public double Function(double x)
        {
            double result = 0;
            //np.sin(10*x)*x + np.cos(2*x)*x 
            result = Math.Sin(10 * x) * x + Math.Cos(2 * x) * x;
          //  result = x * x;
            return result;
         }
        //传入种群数去获取相对应的具体的数值
         public double TranslateDNA(individual ind)
        {
            int result =ind.num [1] * 2 + ind.num[2] * 4 + ind.num[3] * 8 + ind.num[4] * 16 + ind.num[5] * 32 +
               ind.num[6] * 64 + ind.num[7] * 128 + ind.num[8] * 256 + ind.num[9] * 512;
            //正数;
            return result / 1023.0 * 5;
         
        }
        // 去输入种群去得到其适应度的大小
        public void GetFitness(individual[] Population)
          { 
          for (int i=0;i<popsize;i++)
           {
                // to translate information of DNA to the number 
          
                // 适应度返回的是个函数的大小
               double variableX = TranslateDNA(Population[i]);
                Population[i].fitness = 0;
                //we define the fitness is the value of the function 
               double fitness = Function(variableX);
                // 看到一个fitness 

                //选择的fitness 有可能有新的问题
                if (fitness < 0 || fitness == 0)
                {
                    Population[i].fitness = 0;
                }
                else
                {
                    Population[i].fitness = fitness;
                }
             //   Console.WriteLine("{0},{1},{2}", variableX, fitness, Population[i].fitness);

            }
            //去找到一个最好的

            // 重置fitness 这样才是比较好的了

        
            //for (int i = 0; i < popsize; i++)
            //{
            //    if (min > Population[i].fitness) min = Population[i].fitness;
            //}
            //// 对fitness做出一个改变 去得到概率不为0 的数字,否则概率为负的一个数值
            //double mm= 0;
            //for (int i = 0; i < popsize; i++)
            //{
               

            //    mm = Population[i].fitness;
            // // fitness 选择的不太好的了,这样计算出来的数据可能效果不太好
            //    Population[i].fitness = mm - min + 0.00005;

            //    //  Console.WriteLine("{0},{1}", mm, Population[i].fitness);

            //}
            //// 根据计算出来的





        }
        //  根据种群和fitness，优胜劣汰
        //
        public void CalculateProbability(individual[] Population)
        {
            //个体概率比较大的容易被选择上
            // 为了保证每个种群的fitness 必须大于0，概率不会出现0 的状态，可以减去一个最小的数字
            // find a non-zero fitness to choose

            double sum = 0;
            
            for (int i = 0; i < popsize; i++)
            {
               sum+= Population[i].fitness;

            }
            // 为了不让最小的选择的概率为00000
            for (int i = 0; i < popsize; i++)
            {
             //计算每个个体被选择的概率 整数     
             Population[i].selectP= Population[i].fitness/sum;

            }
            // 去计算累积的概率
         }
        //去根据 累计概率，返回被选择种群的索引
          private int SelectIndex(double[] sp)
        {
            byte[] buffer = Guid.NewGuid().ToByteArray();
            int iSeed = BitConverter.ToInt32(buffer, 0);
            Random rd = new Random(iSeed);
            double mm = rd.NextDouble();
          //挑选的效果不太明显   double mm = rd.Next(0,101)/100;
            //记住一定要类型转化，否则概率不准确
            // 看一般如何去选择了，选择的多样性
            //  Console.WriteLine("{0}", mm);
            //相当于轮盘法
            // 说明这一块的选择是不怎么对的了
            for (int i=0;i<popsize;i++)
            {
                if (mm < sp[i]) return i;
               
          }

           return 0;
         }
        // 去产生优秀的种群 猴子
      
         public individual[] PopProduce(individual[] Population)
        {



            //  选择种群
            CalculateProbability(Population);


            // 定义了一个新的对象
            individual[] pp = new individual[popsize];

           for(int i=0;i<popsize;i++)
           {
            pp[i] = new individual();
            }


            // to calcualte the 累计的概率
            // 这样的写法会增加算法的运算时间的

            // to get the minimun number 




            // 去计算累计的概率数值
            double[] sp = new double[popsize];
            double r = 0.0;
            for (int i = 0; i < popsize; i++)
            {
                //去计算累加的概率
                r += Population[i].selectP;
                // 选择一个最接近的
                sp[i] = r;
           //     Console.WriteLine(r);
            }
            // 这一块的内容如何去写了

            for (int i=0;i<popsize;i++)
            {
                int k = SelectIndex(sp);
               // Console.WriteLine(k);
                //去进行一个复制
                pp[i].fitness = Population[k].fitness;
                pp[i].selectP = Population[k].selectP;
                for(int j=0;j<DNAsize;j++)
                {
                    pp[i].num[j] = Population[k].num[j];
                }
                //  Console.WriteLine(k);
            }
            // 为什么要赋给其他的一个数组了 
          return pp;
        }





        //去返回来一个种群的对象
        public individual CrossOver(individual parent,individual[] pop)
        {
            individual child = new individual();
            child = parent;
            //如何去写遗传算法的交叉概率

            byte[] buffer = Guid.NewGuid().ToByteArray(); int iSeed = BitConverter.ToInt32(buffer, 0);
            Random random = new Random(iSeed);

            int r1 = random.Next(0,101);
            double r = r1 / 100.0;
          //  Console.WriteLine(r);
           //去判断是否进行交叉
           //这样如何去交叉配对了
           if(r1<=crossOver)
            {
                //去进行DNA片段信息的交叉  去决定选择与哪个去进行交叉
                byte[] buffer1 = Guid.NewGuid().ToByteArray(); int iSeed1 = BitConverter.ToInt32(buffer1, 0);
                Random random1 = new Random(iSeed1);
                int index = random1.Next(0,popsize);
                //to get the information of this population
                individual ig = null;
                ig = pop[index];
                // to decide which position to cross over.
                for (int j=1;j<DNAsize-1; j++)
                {
                    byte[] buffer2 = Guid.NewGuid().ToByteArray(); int iSeed2 = BitConverter.ToInt32(buffer2, 0);
                    Random random2 = new Random(iSeed2);
                    int mm = random2.Next(0,2);
                    if(mm==1)
                    {
                        //to record the index
                       // 
                     child.num[j] = pop[index].num[j];
                      
                    }
       
              //to change the position  
              
               }
        // 可以直接返回parent
            return child;
            }
           return child;
        }

        public individual Mutation(individual child)
        {
            individual newchild = null;
            // 如何去更好的变异
            // 去控制概率的
            newchild = new individual();
            newchild = child;
            for (int j=0;j<DNAsize;j++)

           //变异后计算fitness 这样会怎么样的了。。。。。
            {
                byte[] buffer3 = Guid.NewGuid().ToByteArray(); int iSeed3 = BitConverter.ToInt32(buffer3, 0);
                Random random = new Random(iSeed3);
                //这个函数可以更好一点的了
                // 说明变异选择的是比较好的了
                double p = random.Next(0,101)/100.0;
           //  每一个DNA变异的概率  Console.Write("pp{0}\n",p);
                if(p<mutationRate)
                {
            // 这个就是个条件运算符号
               newchild.num[j]=(child.num[j]) == 0 ? 1 :0;

            //去得到其 DMA 的信息
                }

            }
          return newchild;
      
        }

    } 























    }
