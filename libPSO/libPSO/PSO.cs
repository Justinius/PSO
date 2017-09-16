using System;
using System.Collections.Generic;
using System.Text;

namespace libPSO
{

    public enum EdgeEffect { Wrap, Clamp, Reflect };
    public enum Topology { Global, Ring, Random };
    public enum WeightScaling { Constant, Linear }

    public class PSOSettings
    {
        int numParticles;
        bool keepHistory;
        int historySpan;
        
        //Ioan Cristian Trelea. 
        //The particle swarm optimization algorithm: convergence analysis and parameter selection.
        //Inf.Process.Lett., 85(6):317-325, 2003. 
        double w;
        double lWeight;
        double gWeight;

        int maxIters;
        int numItersNoImprovement;
        double noImprovementThreshold;

        EdgeEffect edgeEffect = EdgeEffect.Reflect;
        double velEdgeEffectMul = .8;

        Topology topology;
        int randomTopologyNeighborhoodSize;

        WeightScaling velScalingType;
        double startVelocityWeight;
        double endVelocityWeight;
        int startVelocityIter;
        int endVelocityIter;

        WeightScaling neighborScalingType;
        double startNeighborhoodWeight;
        double endNeighborhoodWeight;
        int startNeighborhoodIter;
        int endNeighborhoodIter;

        WeightScaling localScalingType;
        double startLocalWeight;
        double endLocalWeight;
        int startLocalIter;
        int endLocalIter;

        public int NumParticles { get => numParticles; set => numParticles = value; }

        public bool KeepHistory { get => keepHistory; set => keepHistory = value; }
        public int HistorySpan { get => historySpan; set => historySpan = value; }

        public int MaxIters { get => maxIters; set => maxIters = value; }
        public int NumItersNoImprovement { get => numItersNoImprovement; set => numItersNoImprovement = value; }
        public double NoImprovementThreshold { get => noImprovementThreshold; set => noImprovementThreshold = value; }

        public EdgeEffect EdgeEffect { get => edgeEffect; set => edgeEffect = value; }
        public double VelEdgeEffectMul { get => velEdgeEffectMul; set => velEdgeEffectMul = value; }

        public Topology SelTopology { get => topology; set => topology = value; }
        public int RandomTopologyNeighborhoodSize { get => randomTopologyNeighborhoodSize; set => randomTopologyNeighborhoodSize = value; }

        public WeightScaling VelScalingType { get => velScalingType; set => velScalingType = value; }
        public double StartVelocityWeight { get => startVelocityWeight; set => startVelocityWeight = value; }
        public double EndVelocityWeight { get => endVelocityWeight; set => endVelocityWeight = value; }
        public int StartVelocityIter { get => startVelocityIter; set => startVelocityIter = value; }
        public int EndVelocityIter { get => endVelocityIter; set => endVelocityIter = value; }

        public WeightScaling NeighborScalingType { get => neighborScalingType; set => neighborScalingType = value; }
        public double StartNeighborhoodWeight { get => startNeighborhoodWeight; set => startNeighborhoodWeight = value; }
        public double EndNeighborhoodWeight { get => endNeighborhoodWeight; set => endNeighborhoodWeight = value; }
        public int StartNeighborhoodIter { get => startNeighborhoodIter; set => startNeighborhoodIter = value; }
        public int EndNeighborhoodIter { get => endNeighborhoodIter; set => endNeighborhoodIter = value; }

        public WeightScaling LocalScalingType { get => localScalingType; set => localScalingType = value; }
        public double StartLocalWeight { get => startLocalWeight; set => startLocalWeight = value; }
        public double EndLocalWeight { get => endLocalWeight; set => endLocalWeight = value; }
        public int StartLocalIter { get => startLocalIter; set => startLocalIter = value; }
        public int EndLocalIter { get => endLocalIter; set => endLocalIter = value; }

        public double W { get => w; set => w = value; }
        public double LWeight { get => lWeight; set => lWeight = value; }
        public double GWeight { get => gWeight; set => gWeight = value; }

        public PSOSettings()
        {
            NumParticles = 20;
            KeepHistory = false;
            W = .7968;
            LWeight = 1.4962;
            GWeight = 1.4962;
            MaxIters = 500;
            NumItersNoImprovement = (int)(MaxIters * .1);
            NoImprovementThreshold = .05;
            EdgeEffect = EdgeEffect.Reflect;
            VelEdgeEffectMul = .8;
            SelTopology = Topology.Global;
        }

        public void EnableHistory(int span)
        {
            KeepHistory = true;
            HistorySpan = span;
        }

        public void SetTopology(Topology T, int param = 0)
        {
            SelTopology = T;
            if(SelTopology == Topology.Random)
            {
                RandomTopologyNeighborhoodSize = param;
                if(RandomTopologyNeighborhoodSize == 0)
                {
                    throw new ArgumentException("Random Topology needs neighborhood size.");
                }
            }
        }

        public void SetNumParticles(int num)
        {
            if (num <= 0) throw new ArgumentException("Number of Particles must be greater than 0.");
            NumParticles = num;
        }

        public void SetEdgeEffect(EdgeEffect E, double param = .8)
        {
            EdgeEffect = E;
            if (EdgeEffect != EdgeEffect.Clamp)
            {
                VelEdgeEffectMul = param;
            }
        }

        public void SetVelocityWeight(double velWeight)
        {
            StartVelocityWeight = velWeight;
            VelScalingType = WeightScaling.Constant;
        }

        public void SetVelocityWeight(int startIter, double startWeight, int endIter, double endWeight)
        {
            VelScalingType = WeightScaling.Linear;
            StartVelocityWeight = startWeight;
            EndVelocityWeight = endWeight;
            StartVelocityIter = startIter;
            EndVelocityIter = endIter;
        }

        public void SetNeighborhoodWeight(double neighborWeight)
        {
            StartNeighborhoodWeight = neighborWeight;
            NeighborScalingType = WeightScaling.Constant;
        }

        public void SetNeightborhoodWeight(int startIter, double startWeight, int endIter, double endWeight)
        {
            NeighborScalingType = WeightScaling.Linear;
            StartNeighborhoodWeight = startWeight;
            EndNeighborhoodWeight = endWeight;
            StartNeighborhoodIter = startIter;
            EndNeighborhoodIter = endIter;
        }

        public void SetLocalWeight(double GlobalWeight)
        {
            StartLocalWeight = GlobalWeight;
            LocalScalingType = WeightScaling.Constant;
        }

        public void SetLocalWeight(int startIter, double startWeight, int endIter, double endWeight)
        {
            LocalScalingType = WeightScaling.Linear;
            StartLocalWeight = startWeight;
            EndLocalWeight = endWeight;
            StartLocalIter = startIter;
            EndLocalIter = endIter;
        }
    }
    
    public class PSO
    {
        PSOSettings currSettings;

        int numDims;
        int numParticles;
        double[] lbounds;
        double[] ubounds;
        double[] gBestPos;
        double gBest;

        double[][] particles;
        double[][] velocities;
        double[][] lBestPos;
        double[][] nBestPos;
        bool[][] topologyGraph;
        double[] scores;
        double[] localBestScores;

        public List<double> gBestHistory;
        public List<double[][]> particleHistory;

        Func<double[], double> goalFunc;
        Random rnd = new Random();

        double w_b, w_m;
        double lW_b, lW_m;
        double gW_b, gW_m;

        Func<int, double> updateVelocityWeight;
        Func<int, double> updateLocalWeight;
        Func<int, double> updateGlobalWeight;

        Action GraphUpdate;

        public PSO(PSOSettings settings, int inNumDims, 
                   Func<double[], double> goal, 
                   double[] lbounds, double[] ubounds 
                  )
        {

            currSettings = settings;
            numDims = inNumDims;
            numParticles = currSettings.NumParticles;

            if (numDims <= 0)
                throw new ArgumentException("Number of dimensions should be greater than 0.");

            if(currSettings.NumParticles <= 0)
                throw new ArgumentException("Number of particles should be greather than 0.");

            if (currSettings.MaxIters <= 0)
                throw new ArgumentException("Number of iterations should be greather than 0.");

            if (lbounds.Length != numDims || ubounds.Length != numDims)
                throw new ArgumentException("Bounds should have a length of NUMDIMS.");
                                   
            for(int i = 0; i < numDims; i++)
            {
                if (lbounds[i] >= ubounds[i])
                    throw new ArgumentException("Bounds are degenerate.");
            }

            this.lbounds = new double[numDims];
            this.ubounds = new double[numDims];
            for (int i = 0; i < numDims; i++)
            {
                this.lbounds[i] = lbounds[i];
                this.ubounds[i] = ubounds[i];
            }

            goalFunc = goal;

            Init();
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
            //create all the particles and associated data
            particles = new double[numParticles][];
            velocities = new double[numParticles][];
            lBestPos = new double[numParticles][];
            nBestPos = new double[numParticles][];
            scores = new double[numParticles];
            localBestScores = new double[numParticles];
            for (int i = 0; i < numParticles; i++)
            {
                particles[i] = new double[numDims];
                velocities[i] = new double[numDims];
                lBestPos[i] = new double[numDims];
                nBestPos[i] = new double[numDims];
            }

            gBestPos = new double[numDims];
            gBestHistory = new List<double>();
            particleHistory = new List<double[][]>();
            
            //algorithm starts with random init of pos and vel
            //and the scores that go along with it
            for (int i = 0; i < numParticles; i++)
            {
                for (int j = 0; j < numDims; j++)
                {
                    double range = Math.Abs(ubounds[j] - lbounds[j]);

                    particles[i][j] = RndRange(lbounds[j],ubounds[j]);
                    velocities[i][j] = .5*range*(rnd.NextDouble()-rnd.NextDouble());
                    lBestPos[i][j] = particles[i][j];
                }
                scores[i] = goalFunc(particles[i]);
                localBestScores[i] = scores[i];
            }

            int minIdx = GetMinScoreIdx();
            gBest = scores[minIdx];
            particles[minIdx].CopyTo(gBestPos,0);

            if(currSettings.KeepHistory)
            {
                gBestHistory.Add(gBest);
                particleHistory.Add(DuplicateParticles());
            }

            //somewhat advanced variants of algorithm have various information exchange
            //global, random, ring
            if (currSettings.SelTopology != Topology.Global)
            {
                topologyGraph = new bool[numParticles][];
                for (int i = 0; i < numParticles; i++)
                {
                    topologyGraph[i] = new bool[numParticles];
                    for (int j = 0; j < numParticles; j++)
                    {
                        topologyGraph[i][j] = false;
                    }
                }

                if (currSettings.SelTopology == Topology.Random)
                {
                    InitRandomTopology();
                }
                else if (currSettings.SelTopology == Topology.Ring)
                {
                    InitRingTopology();
                }

                GraphUpdate = UpdateNeighborhood;
            }
            else
            {
                GraphUpdate = UpdateGlobalNeighborhood;
            }
            
            if (currSettings.VelScalingType == WeightScaling.Constant)
            {
                updateVelocityWeight = UpdateVelWeight;
            }
            else
            {
                w_m = (currSettings.StartVelocityWeight - currSettings.EndVelocityWeight) / (currSettings.StartVelocityIter - currSettings.EndVelocityIter);
                w_b = currSettings.StartVelocityWeight - w_m * currSettings.StartVelocityIter;
                updateVelocityWeight = UpdateVelWeightLinear;
            }

            if (currSettings.LocalScalingType == WeightScaling.Constant)
            {
                updateLocalWeight = UpdateLWeight;
            }
            else
            {
                lW_m = (currSettings.StartLocalWeight - currSettings.EndLocalWeight) / (currSettings.StartLocalIter - currSettings.EndLocalIter);
                lW_b = currSettings.StartLocalWeight - lW_m * currSettings.StartLocalIter;
                updateLocalWeight = UpdateLWeightLinear;
            }

            if (currSettings.NeighborScalingType == WeightScaling.Constant)
            {
                updateGlobalWeight = UpdateGWeight;
            }
            else
            {
                gW_m = (currSettings.StartNeighborhoodWeight - currSettings.EndNeighborhoodWeight) / (currSettings.StartNeighborhoodIter - currSettings.EndNeighborhoodIter);
                gW_b = currSettings.StartNeighborhoodWeight - gW_m * currSettings.StartNeighborhoodIter;
                updateGlobalWeight = UpdateGWeightLinear;
            }
            
        }

        public void UpdateGlobalNeighborhood()
        {
            //decided to do this with logic in main loop than possibly waste all these memory copies
            return;
            //for(int i = 0; i < numParticles; i++)
            //{
            //    gBestPos.CopyTo(nBestPos[i], 0);
            //}
        }

        public void UpdateNeighborhood()
        {
            for(int i = 0; i < numParticles; i++)
            {
                int bestParticleInNeighborhood = i;
                for(int j = 0; j < numParticles; j++)
                {
                    if(topologyGraph[i][j] && (localBestScores[j] < localBestScores[bestParticleInNeighborhood]))
                    {
                        bestParticleInNeighborhood = j;
                    }
                }
                lBestPos[bestParticleInNeighborhood].CopyTo(nBestPos[i], 0);
            }
        }

        public double UpdateVelWeight(int iter)
        {
            return currSettings.W;
        }

        public double UpdateVelWeightLinear(int iter)
        {
            if(iter < currSettings.StartVelocityIter)
            {
                return currSettings.StartVelocityWeight;
            }
            if(iter > currSettings.EndVelocityIter)
            {
                return currSettings.EndVelocityWeight;
            }

            //wStart = m*wStartIter + b
            //wEnd   = m*wEndIter + b
            //wStart - wEnd = m(wStartIter-wEndIter)
            //double m = (wStart - wEnd) / (wStartIter - wEndIter);
            //double b = wStart - m*wStartIter

            return w_m * iter + w_b;
        }

        public double UpdateLWeight(int iter)
        {
            return currSettings.LWeight;
        }

        public double UpdateLWeightLinear(int iter)
        {
            if (iter < currSettings.StartLocalIter)
            {
                return currSettings.StartLocalWeight;
            }
            if (iter > currSettings.EndLocalIter)
            {
                return currSettings.EndLocalWeight;
            }

            //wStart = m*wStartIter + b
            //wEnd   = m*wEndIter + b
            //wStart - wEnd = m(wStartIter-wEndIter)
            //double m = (wStart - wEnd) / (wStartIter - wEndIter);
            //double b = wStart - m*wStartIter

            return lW_m * iter + lW_b;
        }

        public double UpdateGWeight(int iter)
        {
            return currSettings.GWeight;
        }

        public double UpdateGWeightLinear(int iter)
        {
            if (iter < currSettings.StartVelocityIter)
            {
                return currSettings.StartNeighborhoodWeight;
            }
            if (iter > currSettings.EndNeighborhoodIter)
            {
                return currSettings.EndNeighborhoodWeight;
            }

            //wStart = m*wStartIter + b
            //wEnd   = m*wEndIter + b
            //wStart - wEnd = m(wStartIter-wEndIter)
            //double m = (wStart - wEnd) / (wStartIter - wEndIter);
            //double b = wStart - m*wStartIter

            return gW_m * iter + gW_b;
        }

        public void InitRingTopology()
        {
            for(int i = 0; i < numParticles; i++)
            {
                topologyGraph[i][i] = true;
                if(i == 0)
                {
                    topologyGraph[i][i + 1] = true;
                    topologyGraph[i][numParticles - 1] = true;
                }
                else if(i == (numParticles - 1))
                {
                    topologyGraph[i][i - 1] = true;
                    topologyGraph[i][0] = true;
                }
                else
                {
                    topologyGraph[i][i + 1] = true;
                    topologyGraph[i][i - 1] = true;
                }
            }
        }

        public void InitRandomTopology()
        {
            for (int i = 0; i < numParticles; i++)
            {
                topologyGraph[i][i] = true;
                for(int j = 0; j < currSettings.RandomTopologyNeighborhoodSize; j++)
                {
                    int randIdx = rnd.Next(0, numParticles);
                    topologyGraph[i][randIdx] = true;
                }
            }
        }

        private void CalcEdgeEffects(int partIdx, int dimIdx)
        {
            if (particles[partIdx][dimIdx] < lbounds[dimIdx])
            {
                if(currSettings.EdgeEffect == EdgeEffect.Clamp)
                {
                    particles[partIdx][dimIdx] = lbounds[dimIdx];
                    velocities[partIdx][dimIdx] = 0;
                }
                else if(currSettings.EdgeEffect == EdgeEffect.Wrap)
                {
                    particles[partIdx][dimIdx] = ubounds[dimIdx] - ((lbounds[dimIdx] - particles[partIdx][dimIdx]) % (ubounds[dimIdx]- lbounds[dimIdx]));
                    velocities[partIdx][dimIdx] *= currSettings.VelEdgeEffectMul;
                }
                else
                {
                    particles[partIdx][dimIdx] = lbounds[dimIdx] + ((lbounds[dimIdx] - particles[partIdx][dimIdx]) % (ubounds[dimIdx] - lbounds[dimIdx]));
                    velocities[partIdx][dimIdx] *= currSettings.VelEdgeEffectMul;
                }                                
            }

            if (particles[partIdx][dimIdx] > ubounds[dimIdx])
            {
                if (currSettings.EdgeEffect == EdgeEffect.Clamp)
                {
                    particles[partIdx][dimIdx] = ubounds[dimIdx];
                    velocities[partIdx][dimIdx] = 0;
                }
                else if (currSettings.EdgeEffect == EdgeEffect.Wrap)
                {
                    particles[partIdx][dimIdx] = lbounds[dimIdx] + ((particles[partIdx][dimIdx] - ubounds[dimIdx]) % (ubounds[dimIdx] - lbounds[dimIdx]));
                    velocities[partIdx][dimIdx] *= currSettings.VelEdgeEffectMul;
                }
                else
                {
                    particles[partIdx][dimIdx] = ubounds[dimIdx] - ((particles[partIdx][dimIdx] - ubounds[dimIdx]) % (ubounds[dimIdx] - lbounds[dimIdx]));
                    velocities[partIdx][dimIdx] *= currSettings.VelEdgeEffectMul;
                }                                
            }
        }

        public double Optimize()
        {
            for (int iter = 1; iter < currSettings.MaxIters; iter++)
            {
                GraphUpdate();
                for (int i = 0; i < numParticles; i++)
                {
                    //update velocity
                    double localR, globalR;
                    for (int j = 0; j < numDims; j++)
                    {
                        localR = rnd.NextDouble();
                        globalR = rnd.NextDouble();

                        velocities[i][j] = updateVelocityWeight(iter) * velocities[i][j] +
                                           updateLocalWeight(iter) * localR * (lBestPos[i][j] - particles[i][j]);

                        //if global topology then use gBest, if not then has to use nBest
                        //if global, which is probably most use cases, save numParticles*numIters array copy
                        if (currSettings.SelTopology == Topology.Global)
                        {
                            velocities[i][j] += updateGlobalWeight(iter) * globalR * (gBestPos[j] - particles[i][j]);
                        }
                        else
                        {
                            velocities[i][j] += updateGlobalWeight(iter) * globalR * (nBestPos[i][j] - particles[i][j]);
                        }
                                                
                        particles[i][j] = particles[i][j] + velocities[i][j];

                        CalcEdgeEffects(i, j);
                    }

                    scores[i] = goalFunc(particles[i]);
                    if (scores[i] < localBestScores[i])
                    {
                        particles[i].CopyTo(lBestPos[i], 0);
                        localBestScores[i] = scores[i];
                    }
                }

                int minIdx = GetMinScoreIdx();
                if (scores[minIdx] < gBest)
                {
                    gBest = scores[minIdx];
                    particles[minIdx].CopyTo(gBestPos, 0);
                }

                if (currSettings.KeepHistory && (iter % currSettings.HistorySpan) == 0)
                {
                    gBestHistory.Add(gBest);
                    particleHistory.Add(DuplicateParticles());
                }
            }
            return gBest;
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
                
    }
}
