using System;

namespace libOptimTests
{
    public class OptimTests
    {
        public double[] lbounds;
        public double[] ubounds;
        public int numDims;
        public double exactSolution;
        public double[] minPoint;

        public virtual double goalFunc(double[] x)
        {
            return 0;
        }
    }

    public class Spherical : OptimTests
    {
        public Spherical(int dims, double absBound)
        {
            numDims = dims;
            lbounds = new double[dims];
            ubounds = new double[dims];
            minPoint = new double[numDims];
            for (int i = 0; i < numDims; i++)
            {
                lbounds[i] = -absBound;
                ubounds[i] = absBound;
                minPoint[i] = 0.0;
            }

            exactSolution = 0.0;
        }

        public override double goalFunc(double[] x)
        {
            double val = 0;
            for (int i = 0; i < numDims; i++) val += x[i] * x[i];
            return val;
        }
    }

    public class Rastrigin : OptimTests
    {
        public Rastrigin()
        {
            numDims = 2;
            lbounds = new double[numDims];
            ubounds = new double[numDims];
            minPoint = new double[numDims];
            for (int i = 0; i < numDims; i++)
            {
                lbounds[i] = -5.12;
                ubounds[i] = 5.12;
                minPoint[i] = 0.0;
            }

            exactSolution = 0.0;
        }

        public override double goalFunc(double[] x)
        {
            double val = 0;
            for (int i = 0; i < numDims; i++) val += (x[i] * x[i] - 10 * Math.Cos(2 * Math.PI * x[i]));
            return val + 20;
        }
    }

    public class Ackley : OptimTests
    {
        public Ackley()
        {
            numDims = 2;
            lbounds = new double[numDims];
            ubounds = new double[numDims];
            minPoint = new double[numDims];
            for (int i = 0; i < numDims; i++)
            {
                lbounds[i] = -5;
                ubounds[i] = 5;
                minPoint[i] = 0.0;
            }

            exactSolution = 0.0;
        }

        public override double goalFunc(double[] x)
        {
            return -20 * Math.Exp(-.2 * Math.Sqrt(.5 * (x[0] * x[0] - x[1] * x[1]))) - Math.Exp(.5 * (Math.Cos(2 * Math.PI * x[0]) + Math.Cos(2 * Math.PI * x[1]))) + Math.E + 20;
        }
    }

    public class Rosenbrock : OptimTests
    {
        public Rosenbrock(int dims, double absBound)
        {
            numDims = dims;
            lbounds = new double[dims];
            ubounds = new double[dims];
            minPoint = new double[numDims];
            for (int i = 0; i < numDims; i++)
            {
                lbounds[i] = -absBound;
                ubounds[i] = absBound;
                minPoint[i] = 1.0;
            }

            exactSolution = 0.0;
        }

        public override double goalFunc(double[] x)
        {
            double val = 0;
            for (int i = 0; i < numDims - 1; i++)
                val += (100 * (x[i + 1] - x[i] * x[i]) * (x[i + 1] - x[i] * x[i]) + (x[i] * x[i] - 1) * (x[i] * x[i] - 1));
            return val;
        }
    }

    public class Beale : OptimTests
    {
        public Beale()
        {
            numDims = 2;
            lbounds = new double[numDims];
            ubounds = new double[numDims];
            minPoint = new double[] { 3, .5 };
            for (int i = 0; i < numDims; i++)
            {
                lbounds[i] = -4.5;
                ubounds[i] = 4.5;
            }

            exactSolution = 0.0;
        }

        public override double goalFunc(double[] x)
        {
            double val = (1.5-x[0]+x[0]*x[1])* (1.5 - x[0] + x[0] * x[1]);
            val += (2.25 - x[0] + x[0] * x[1] * x[1]) * (2.25 - x[0] + x[0] * x[1] * x[1]);
            val += (2.625 - x[0] + x[0] * x[1] * x[1]*x[1]) * (2.625 - x[0] + x[0] * x[1] * x[1]*x[1]);
            return val;
        }
    }


}
