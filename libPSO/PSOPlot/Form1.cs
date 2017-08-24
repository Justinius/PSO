using libPSO;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace PSOPlot
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        private static int ScaleDim(double val, double lb, double ub, int imageDim)
        {
            //y = mx + b
            //0 = m(lb) + b
            //imageDim = m(ub) + b
            //b = -m(lb)
            //imageDim = m(ub) - m(lb)
            //imageDim = m(ub - lb)
            //m = imageDim/(ub-lb)

            double slope = imageDim / (ub - lb);
            double intercept = -1 * slope * lb;
            return (int)(slope * val + intercept);

        }

        private void button1_Click(object sender, EventArgs e)
        {
            double[] lb = new double[] { -100, -100 };
            double[] ub = new double[] { 100, 100 };
            PSO myPSO = new PSO(2, 20, g, lb, ub, 1000, true, .5, 2, 2);
            double best = myPSO.Optimize();

            Bitmap bmp = new Bitmap(750,750);

            

            for (int numIters = 0; numIters < 1000; numIters++)
            {
                for (int i = 0; i < 750; i++)
                {
                    for (int j = 0; j < 750; j++)
                    {
                        bmp.SetPixel(i, j, Color.White);
                    }
                }


                for (int j = 0; j < 20; j++)
                {
                    int x = ScaleDim(myPSO.particleHistory[numIters][j][0], lb[0], ub[0], 750);
                    int y = ScaleDim(myPSO.particleHistory[numIters][j][1], lb[1], ub[1], 750);
                         
                    if(x > 2 && x < 748 && y > 2 && y < 748)
                    {
                        for(int q = x-2; q < x+2; q++)
                        {
                            for(int r = y-2; r < y+2; r++)
                            {
                                bmp.SetPixel(q, r, Color.Red);
                            }
                        }
                    }
                    else
                        bmp.SetPixel(x, y, Color.Red);
                }

                //bmp.Save("C:\\Users\\fista\\Documents\\pso\\" + numIters.ToString() + ".png", System.Drawing.Imaging.ImageFormat.Png);

                pictureBox1.Image = bmp;
                pictureBox1.Invalidate();
                pictureBox1.Refresh();
                textBox1.Text = numIters.ToString();
            }
            
        }

        static double g(double[] x)
        {
            double val = 0;
            for (int i = 0; i < x.Length; i++)
                val += (x[i] * x[i]);
            return val;
        }
    }
}
