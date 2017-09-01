using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using libPSO;
using libOptimTests;

namespace PSOTester
{
    class Program
    {
        static void Main(string[] args)
        {
            OptimTests myGoal = new StyblinskiTang(2);// (4, 100);

            double[] lb = new double[] { -100,  -100,-100,-100,-100};
            double[] ub = new double[] { 100 , 100,100,100,100};
            PSO myPSO = new PSO(myGoal.numDims, 10, myGoal.goalFunc, myGoal.lbounds, myGoal.ubounds, 500, true);//, .5, 2, 2);
            double best = myPSO.Optimize();

            Console.WriteLine("Best: " + best);
            Console.ReadLine();
            
            for(int i = 0; i < myPSO.gBestHistory.Count; i++)
            {
                Console.WriteLine(i + " " + myPSO.gBestHistory[i]);

            }
            Console.ReadLine();
        }

        static double g(double[] x)
        {
            double val=0;
            for (int i = 0; i < x.Length; i++)
                val += (x[i] * x[i]);
            return val;
        }
    }
}
