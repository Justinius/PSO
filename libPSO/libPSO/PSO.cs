using System;
using System.Collections.Generic;
using System.Text;

namespace libPSO
{
    class PSO
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
        List<double> gBestHistory;
        List<double[][]> particleHistory;

        Func<double[], double> goalFunc;
        Random rnd = new Random();


        //Ioan Cristian Trelea. 
        //The particle swarm optimization algorithm: convergence analysis and parameter selection.
        //Inf.Process.Lett., 85(6):317-325, 2003. 
        double w;
        double lWeight;
        double gWeight;
        
        public PSO(int numDims, int numParticles, 
                   Func<double[], double> goal, 
                   double[] lbounds, double[] ubounds, 
                   double velocityWeight = .7968, double localBestWeight = 1.4962, double globalBestWeight = 1.4962,
                   bool keepHistory = false)
        {
            
            if (numDims <= 0)
                throw new ArgumentException("Number of dimensions should be greater than 0.");

            if(numParticles <= 0)
                throw new ArgumentException("Number of particles should be greather than 0.");

            if (lbounds.Length != numDims || ubounds.Length != numDims)
                throw new ArgumentException("Bounds should have a length of NUMDIMS.");

            w = velocityWeight;
            lWeight = localBestWeight;
            gWeight = globalBestWeight;
            
            for(int i = 0; i < numDims; i++)
            {
                if (lbounds[i] >= ubounds[i])
                    throw new ArgumentException("Bounda are degenerate.");
            }

            this.numDims = numDims;
            this.numParticles = numParticles;
            this.lbounds = new double[numDims];
            this.lbounds = new double[numDims];
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
                    velocities[i][j] = RndRange(-range, range);
                    lbestPos[i][j] = particles[i][j];
                }
                scores[i] = goalFunc(particles[i]);
                localBestScores[i] = scores[i];
            }

            int minIdx = GetMinScoreIdx();
            gBest = scores[minIdx];
            gBestPos = particles[minIdx];

            if(keepHistory)
            {
                gBestHistory.Add(gBest);
                particleHistory.Add(particles);
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

                    if (particles[i][j] < lbounds[j])
                    {
                        velocities[i][j] = -velocities[i][j];
                        double updatedPos = lbounds[j] + (lbounds[j] - particles[i][j]);
                        if(updatedPos > ubounds[j])
                        {
                            particles[i][j] = RndRange(lbounds[j], ubounds[j]);
                        }
                        else
                        {
                            particles[i][j] = updatedPos;
                        }
                     
                    }

                    if (particles[i][j] > ubounds[j])
                    {
                        velocities[i][j] = -velocities[i][j];
                        double updatedPos = ubounds[j] - (particles[i][j] - ubounds[j]);
                        if (updatedPos < lbounds[j])
                        {
                            particles[i][j] = RndRange(lbounds[j], ubounds[j]);
                        }
                        else
                        {
                            particles[i][j] = updatedPos;
                        }
                    }
                }

                scores[i] = goalFunc(particles[i]);
                if(scores[i] < localBestScores[i])
                {
                    localBestScores[i] = scores[i];
                    lbestPos[i] = particles[i];
                }
            }

            int minIdx = GetMinScoreIdx();
            if(scores[minIdx] < gBest)
            {
                gBest = scores[minIdx];
                gBestPos = particles[minIdx];
            }

            if(keepHistory)
            {
                gBestHistory.Add(gBest);
                particleHistory.Add(particles);
            }
        }

    }
}
