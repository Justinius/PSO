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
            return -20 * Math.Exp(-.2 * Math.Sqrt(.5 * (x[0] * x[0] - x[0] * x[0]))) - Math.Exp(.5 * (Math.Cos(2 * Math.PI * x[0]) + Math.Cos(2 * Math.PI * x[0]))) + Math.E + 20;
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
            double val = (1.5-x[0]+x[0]*x[0])* (1.5 - x[0] + x[0] * x[0]);
            val += (2.25 - x[0] + x[0] * x[0] * x[0]) * (2.25 - x[0] + x[0] * x[0] * x[0]);
            val += (2.625 - x[0] + x[0] * x[0] * x[0]*x[0]) * (2.625 - x[0] + x[0] * x[0] * x[0]*x[0]);
            return val;
        }
    }

    public class GoldmanPrice : OptimTests
    {
        public GoldmanPrice()
        {
            numDims = 2;
            lbounds = new double[numDims];
            ubounds = new double[numDims];
            minPoint = new double[] { 0, -1 };
            for (int i = 0; i < numDims; i++)
            {
                lbounds[i] = -2;
                ubounds[i] = 2;
            }

            exactSolution = 3.0;
        }

        public override double goalFunc(double[] x)
        {
            double fact1a = (x[0] + x[1] + 1) * (x[0] + x[1] + 1);
            double fact1b = 19 - 14 * x[0] + 3 * x[0] * x[0] - 14 * x[1] + 6 * x[0] * x[1] + 3 * x[1] * x[1];
            double fact1 = 1 + fact1a * fact1b;

            double fact2a = (2 * x[0] - 3 * x[1])* (2 * x[0] - 3 * x[1]);
            double fact2b = 18 - 32 * x[0] + 12 * x[0] * x[0] + 48 * x[1] - 36 * x[0] * x[1] + 27 * x[1] * x[1];
            double fact2 = 30 + fact2a * fact2b;

            return fact1 * fact2;
        }
    }

    public class Booth : OptimTests
    {
        public Booth()
        {
            numDims = 2;
            lbounds = new double[numDims];
            ubounds = new double[numDims];
            minPoint = new double[] { 1, 3 };
            for (int i = 0; i < numDims; i++)
            {
                lbounds[i] = -10;
                ubounds[i] = 10;
            }

            exactSolution = 0.0;
        }

        public override double goalFunc(double[] x)
        {
            double term1 = (x[0] + 2 * x[1] - 7);
            double term2 = (2 * x[0] + x[1] - 5);

            return term1*term1 + term2*term2;

        }
    }

    public class Bukin6 : OptimTests
    {
        public Bukin6()
        {
            numDims = 2;
            lbounds = new double[] { -15, -3 };
            ubounds = new double[] { -5, 3 };
            minPoint = new double[] { -10, 1 };
            exactSolution = 0.0;
        }

        public override double goalFunc(double[] x)
        {
            double term1 = 100 * Math.Sqrt(Math.Abs(x[1] - 0.01 * x[0]*x[0]));
            double term2 = 0.01 * Math.Abs(x[0] + 10);

            return term1 + term2;
            
        }
    }

    public class Matyas : OptimTests
    {
        public Matyas()
        {
            numDims = 2;
            lbounds = new double[numDims];
            ubounds = new double[numDims];
            minPoint = new double[] { 0, 0 };
            for (int i = 0; i < numDims; i++)
            {
                lbounds[i] = -10;
                ubounds[i] = 10;
            }

            exactSolution = 0.0;
        }

        public override double goalFunc(double[] x)
        {
            return .26 * (x[0] * x[0] + x[1] * x[1]) - .48 * x[0] * x[1];
        }
    }

    public class Levi : OptimTests
    {
        public Levi()
        {
            numDims = 2;
            lbounds = new double[numDims];
            ubounds = new double[numDims];
            minPoint = new double[] { 1, 1 };
            for (int i = 0; i < numDims; i++)
            {
                lbounds[i] = -10;
                ubounds[i] = 10;
            }

            exactSolution = 0.0;
        }

        public override double goalFunc(double[] x)
        {
            double term1 = Math.Pow(Math.Sin(3 * Math.PI * x[0]),2);
            double term2 = Math.Pow(x[0] - 1, 2);
            double term3 = 1 + Math.Pow(Math.Sin(3 * Math.PI * x[1]), 2);
            double term4 = Math.Pow(x[1] - 1, 2);
            double term5 = 1 + Math.Pow(Math.Sin(2 * Math.PI * x[1]), 2);

            return term1 + term2 * term3 + term4 * term5;
        }
    }

    public class Camel3 : OptimTests
    {
        public Camel3()
        {
            numDims = 2;
            lbounds = new double[numDims];
            ubounds = new double[numDims];
            minPoint = new double[] { 0, 0 };
            for (int i = 0; i < numDims; i++)
            {
                lbounds[i] = -5;
                ubounds[i] = 5;
            }

            exactSolution = 0.0;
        }

        public override double goalFunc(double[] x)
        {
            double term1 = 2 * x[0] * x[0];
            double term2 = -1.05 * Math.Pow(x[0],4);
            double term3 = Math.Pow(x[0],6)/6;
            double term4 = x[0] * x[1];
            double term5 = x[1] * x[1];

            return term1 + term2 + term3 + term4 + term5;
            
        }
    }

    public class Easom : OptimTests
    {
        public Easom()
        {
            numDims = 2;
            lbounds = new double[numDims];
            ubounds = new double[numDims];
            minPoint = new double[] { Math.PI, Math.PI };
            for (int i = 0; i < numDims; i++)
            {
                lbounds[i] = -100;
                ubounds[i] = 100;
            }

            exactSolution = -1.0;
        }

        public override double goalFunc(double[] x)
        {
            double fact1 = -Math.Cos(x[0]) * Math.Cos(x[1]);
            double fact2 = Math.Exp(-Math.Pow(x[0] - Math.PI,2) - Math.Pow(x[1] - Math.PI,2));

            return fact1 * fact2;
        }
    }

    public class Eggholder : OptimTests
    {
        public Eggholder()
        {
            numDims = 2;
            lbounds = new double[numDims];
            ubounds = new double[numDims];
            minPoint = new double[] { 512, 404.2319};
            for (int i = 0; i < numDims; i++)
            {
                lbounds[i] = -512;
                ubounds[i] = 512;
            }

            exactSolution = -959.6407;
        }

        public override double goalFunc(double[] x)
        {
            double term1 = -(x[1] + 47) * Math.Sin(Math.Sqrt(Math.Abs(x[1] + x[0] / 2 + 47)));
            double term2 = -x[0] * Math.Sin(Math.Sqrt(Math.Abs(x[0] - (x[1] + 47))));

            return term1 + term2;
        }
    }

    public class StyblinskiTang : OptimTests
    {
        public StyblinskiTang(int dims)
        {
            numDims = dims;
            lbounds = new double[numDims];
            ubounds = new double[numDims];
            minPoint = new double[numDims];
            for (int i = 0; i < numDims; i++)
            {
                minPoint[i] = -2.903534;
                lbounds[i] = -5;
                ubounds[i] = 5;
            }

            exactSolution = -39.16599*numDims;
        }

        public override double goalFunc(double[] x)
        {
            double sum = 0;
            for(int i = 0; i < numDims; i++)
            {
                double curr_x = x[i];
                double term = Math.Pow(curr_x, 4) - 16 * curr_x * curr_x + 5 * curr_x;
                sum += term;
            }
            return sum / 2;
        }
    }

}
