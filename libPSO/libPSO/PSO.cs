using System;
using System.Collections.Generic;
using System.Text;

namespace libPSO
{
    public class PSO
    {
        int numDims;
        int numParticles;
        double[] lbounds;
        double[] ubounds;
        double[] gBestPos;
        double gBest;

        double[][] particles;
        double[][] velocities;
        double[][] lbestPos;
        double[] scores;
        double[] localBestScores;

        bool keepHistory;
        public List<double> gBestHistory;
        public List<double[][]> particleHistory;

        Func<double[], double> goalFunc;
        Random rnd = new Random();

        enum EdgeEffect { Wrap, Clamp, Reflect };

        //Ioan Cristian Trelea. 
        //The particle swarm optimization algorithm: convergence analysis and parameter selection.
        //Inf.Process.Lett., 85(6):317-325, 2003. 
        double w;
        double lWeight;
        double gWeight;

        int maxIters;
        EdgeEffect edgeEffect = EdgeEffect.Reflect;
        double velEdgeEffectMul = .8;
        
        public PSO(int numDims, int numParticles, 
                   Func<double[], double> goal, 
                   double[] lbounds, double[] ubounds, 
                   int iters,
                   bool keepHistory = false,
                   double velocityWeight = .7968, double localBestWeight = 1.4962, double globalBestWeight = 1.4962)
        {
            
            if (numDims <= 0)
                throw new ArgumentException("Number of dimensions should be greater than 0.");

            if(numParticles <= 0)
                throw new ArgumentException("Number of particles should be greather than 0.");

            if (iters <= 0)
                throw new ArgumentException("Number of iterations should be greather than 0.");

            if (lbounds.Length != numDims || ubounds.Length != numDims)
                throw new ArgumentException("Bounds should have a length of NUMDIMS.");

            w = velocityWeight;
            lWeight = localBestWeight;
            gWeight = globalBestWeight;
            maxIters = iters;
            
            for(int i = 0; i < numDims; i++)
            {
                if (lbounds[i] >= ubounds[i])
                    throw new ArgumentException("Bounda are degenerate.");
            }

            this.numDims = numDims;
            this.numParticles = numParticles;
            this.lbounds = new double[numDims];
            this.ubounds = new double[numDims];
            for(int i = 0; i < numDims; i++)
            {
                this.lbounds[i] = lbounds[i];
                this.ubounds[i] = ubounds[i];
            }

            goalFunc = goal;

            particles = new double[numParticles][];
            velocities = new double[numParticles][];
            lbestPos = new double[numParticles][];
            scores = new double[numParticles];
            localBestScores = new double[numParticles];
            for(int i = 0; i < numParticles; i++)
            {
                particles[i] = new double[numDims];
                velocities[i] = new double[numDims];
                lbestPos[i] = new double[numDims];
            }
                        
            this.keepHistory = keepHistory;
            gBestPos = new double[numDims];
            gBestHistory = new List<double>();
            particleHistory = new List<double[][]>();
            
        }

        private int GetMinScoreIdx()
        {
            int minIdx = 0;
            double minVal = scores[minIdx];
            for(int i = 1; i < numParticles; i++)
            {
                if(scores[i]< minVal)
                {
                    minIdx = i;
                    minVal = scores[i];
                }
            }
            return minIdx;
        }

        private double RndRange(double min, double max)
        {
            return rnd.NextDouble() * (max - min) + min;
        }

        private void Init()
        {
            for (int i = 0; i < numParticles; i++)
            {
                for (int j = 0; j < numDims; j++)
                {
                    double range = Math.Abs(ubounds[j] - lbounds[j]);

                    particles[i][j] = RndRange(lbounds[j],ubounds[j]);
                    velocities[i][j] = .5*range*(rnd.NextDouble()-rnd.NextDouble());
                    lbestPos[i][j] = particles[i][j];
                }
                scores[i] = goalFunc(particles[i]);
                localBestScores[i] = scores[i];
            }

            int minIdx = GetMinScoreIdx();
            gBest = scores[minIdx];
            particles[minIdx].CopyTo(gBestPos,0);

            if(keepHistory)
            {
                gBestHistory.Add(gBest);
                particleHistory.Add(DuplicateParticles());
            }
        }

        private void CalcEdgeEffects(int partIdx, int dimIdx)
        {
            if (particles[partIdx][dimIdx] < lbounds[dimIdx])
            {
                if(edgeEffect == EdgeEffect.Clamp)
                {
                    particles[partIdx][dimIdx] = lbounds[dimIdx];
                    velocities[partIdx][dimIdx] = 0;
                }
                else if(edgeEffect == EdgeEffect.Wrap)
                {
                    particles[partIdx][dimIdx] = ubounds[dimIdx] - ((lbounds[dimIdx] - particles[partIdx][dimIdx]) % (ubounds[dimIdx]- lbounds[dimIdx]));
                    velocities[partIdx][dimIdx] *= velEdgeEffectMul;
                }
                else
                {
                    particles[partIdx][dimIdx] = lbounds[dimIdx] + ((lbounds[dimIdx] - particles[partIdx][dimIdx]) % (ubounds[dimIdx] - lbounds[dimIdx]));
                    velocities[partIdx][dimIdx] *= velEdgeEffectMul;
                }                                
            }

            if (particles[partIdx][dimIdx] > ubounds[dimIdx])
            {
                if (edgeEffect == EdgeEffect.Clamp)
                {
                    particles[partIdx][dimIdx] = ubounds[dimIdx];
                    velocities[partIdx][dimIdx] = 0;
                }
                else if (edgeEffect == EdgeEffect.Wrap)
                {
                    particles[partIdx][dimIdx] = lbounds[dimIdx] + ((particles[partIdx][dimIdx] - ubounds[dimIdx]) % (ubounds[dimIdx] - lbounds[dimIdx]));
                    velocities[partIdx][dimIdx] *= velEdgeEffectMul;
                }
                else
                {
                    particles[partIdx][dimIdx] = ubounds[dimIdx] - ((particles[partIdx][dimIdx] - ubounds[dimIdx]) % (ubounds[dimIdx] - lbounds[dimIdx]));
                    velocities[partIdx][dimIdx] *= velEdgeEffectMul;
                }                                
            }
        }

        private void Update()
        {
            for(int i = 0; i < numParticles; i++)
            {
                //update velocity
                double localR, globalR;
                for(int j = 0; j < numDims; j++)
                {
                    localR = rnd.NextDouble();
                    globalR = rnd.NextDouble();
                    velocities[i][j] = w * velocities[i][j] +
                                       lWeight * localR * (lbestPos[i][j] - particles[i][j]) +
                                       gWeight * globalR * (gBestPos[j] - particles[i][j]);
                    particles[i][j] = particles[i][j] + velocities[i][j];

                    CalcEdgeEffects(i, j);
                }

                scores[i] = goalFunc(particles[i]);
                if(scores[i] < localBestScores[i])
                {
                    particles[i].CopyTo(lbestPos[i], 0);
                    localBestScores[i] = scores[i];
             
                }
            }

            int minIdx = GetMinScoreIdx();
            if(scores[minIdx] < gBest)
            {
                gBest = scores[minIdx];
                particles[minIdx].CopyTo(gBestPos, 0);
            }

            if(keepHistory)
            {
                gBestHistory.Add(gBest);
                particleHistory.Add(DuplicateParticles());
            }
        }
        
        private double[][] DuplicateParticles()
        {
            double[][] newHistory = new double[numParticles][];
            for(int i = 0; i < numParticles; i++)
            {
                newHistory[i] = new double[numDims];
                for(int j = 0; j < numDims; j++)
                {
                    newHistory[i][j] = particles[i][j];
                }
            }
            return newHistory;

        }

        public double Optimize()
        {
            Init();
            for(int i = 0; i < maxIters; i++)
            {
                Update();
            }
            return gBest;
        }
    }
}
