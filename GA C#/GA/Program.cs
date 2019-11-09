using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GA
{
    class Program
    {   

          static int generations = 200;
          static int DNAsize = 10;
          static int popsize = 20;
          // 去测试基因遗传算法的优越性能


           

       
         //去写一个copy 的参数
          static void Copy(individual[] source, individual[] destination)
         {
            // 定义种群的大小
            for (int i = 0; i < popsize; i++)
            {
           
                for (int j = 0; j < DNAsize; j++)
                {
                     destination[i].num[j]=source[i].num[j];
                }
                 destination[i].selectP=source[i].selectP ;
                 destination[i].fitness=source[i].fitness ;
            }
          }
          static void Copychild(individual source, individual destination)
          {

            // 去新建一个对象 按照常规失去了作用域
         // 这个是错误的，无法计算出应有的结果数值
         // destination = new individual();
            for (int j = 0; j < DNAsize; j++)
            {
             destination.num[j]  = source.num[j] ;
            }
            destination.selectP= source.selectP;
            destination.fitness =source.fitness ;
        }
     //
        static void Main(string[] args)
        {

         Console.WriteLine("遗传算法");
         //传入种群数，DNA 数 迭代次数
         int NumberGenerations = 200;
         individual[] Population=null;
         GeneticAlgorithm ga = new GeneticAlgorithm();
         //生成初始化的种群。。。
         //生成一系列的DNA序列数,返回种群
            Console.WriteLine("初始化: ");
            individual[] pop= ga.InitialPopulation(Population);
         // translate DNA

            //根据迭代次数去迭代计算
            int count = 0;
            for (int iiii=0;iiii< NumberGenerations; iiii++)
            {


                ga.GetFitness(pop);
                 //这里得到的也是一个引用的类型
                // 这个会交叉引用会变数值的，有可能
                pop = ga.PopProduce(pop);
                //去得到每次的函数的数值

            
                Console.WriteLine("第{0}次迭代", count);
                count++;
                for (int m = 0; m < ga.Popsize; m++)
                {
                    //去得到函数值
                    // most fitted DNA is 
                    double variableX = ga.TranslateDNA(pop[m]);
                    //we define the fitness is the value of the function 
                    double functionValue = ga.Function(variableX);
                    Console.WriteLine("{0}, {1}, {2}, {3}", variableX, functionValue, pop[m].fitness, pop[m].selectP);
                }
                
               //给每个群体一个适应度
               // 根据 fitness 和轮盘转法确定出优秀的群体信息

               //  去copy 一份
                individual[] PopCopy = new individual[popsize];
                for (int iii=0;iii<popsize;iii++)
                {
                    PopCopy[iii] = new individual();

                }
                //这个依旧是引用，不是复制了    
                Copy(pop, PopCopy);


                //父代的信息已经保存到了 parentcopy了


                individual child = new individual();
                // 选择交叉和基因变异

                // 新建一个child的对象的类

                individual[] Child= new individual[popsize];
                // 进行初始化
                for (int ii =0;ii<popsize;ii++)
                {
                Child[ii] = new individual();
                }


                for (int j=0;j<popsize;j++)
                {

                child=ga.CrossOver(pop[j],PopCopy);
                    child = ga.Mutation(child); 
                 //对象引用了而不是复制了
                 //进行复制 child 
                 // parentcopy[j] = child;
                  // 会丢失一段DNA的信息
                Copychild(child, Child[j]);
                }
                pop = Child;

            }
            // 去打印出最大的数值，这样才是比较好的了。。。

            double maxValue = 0;
            for(int i=0;i<popsize;i++)
            {
             //如何去处理一下
             if (maxValue < pop[i].fitness) maxValue = pop[i].fitness;


            }
            Console.WriteLine("最大的数值{0}",maxValue);


        }
    }
}
