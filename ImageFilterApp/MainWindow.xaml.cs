using System;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Media.Imaging;
using Microsoft.Win32;
using System.Drawing;
using System.IO;
using System.Drawing.Imaging;
using System.Windows.Input;
using System.Text.RegularExpressions;
using System.Collections.Generic;
using System.Linq;

namespace ImageFilterApp
{
    public partial class MainWindow : Window
    {

        private BitmapImage original_Image;
        private BitmapImage previous_Image;

        public MainWindow()
        {
            InitializeComponent();
        }

        // The following 5 functions were taken from docs.microsoft or stackoverflow to help with some basic functionalities

        private void Load_Image(object sender, RoutedEventArgs e)
        {
            OpenFileDialog op = new OpenFileDialog();
            op.Title = "Select a picture";
            op.Filter = "All supported graphics|    *.jpg;*.jpeg;*.png|" +
              "JPEG (*.jpg;*.jpeg)|*.jpg;*.jpeg|" +
              "Portable Network Graphic (*.png)|*.png";
            if (op.ShowDialog() == true)
            {
                input_Image.Source = new BitmapImage(new Uri(op.FileName));
                original_Image = new BitmapImage(new Uri(op.FileName));
            }

        }

        private void Save_Image(object sender, RoutedEventArgs e)
        {
            SaveFileDialog dlg = new SaveFileDialog();
            dlg.FileName = "Untitled";
            dlg.Filter = "JPeg Image|*.jpg|Bitmap Image|*.bmp";
            if (dlg.ShowDialog() == true)
            {
                var encoder = new JpegBitmapEncoder(); // Or PngBitmapEncoder, or whichever encoder you want
                encoder.Frames.Add(BitmapFrame.Create(input_Image.Source as BitmapImage));
                using (var stream = dlg.OpenFile())
                {
                    encoder.Save(stream);
                }
            }
        }

        private Bitmap BitmapImage2Bitmap(BitmapImage bitmapImage)
        {
            // BitmapImage bitmapImage = new BitmapImage(new Uri("../Images/test.png", UriKind.Relative));

            using (MemoryStream outStream = new MemoryStream())
            {
                BitmapEncoder enc = new BmpBitmapEncoder();
                enc.Frames.Add(BitmapFrame.Create(bitmapImage));
                enc.Save(outStream);
                System.Drawing.Bitmap bitmap = new System.Drawing.Bitmap(outStream);

                return new Bitmap(bitmap);
            }
        }

        private BitmapImage ToBitmapImage(Bitmap bitmap)
        {
            using (var memory = new MemoryStream())
            {
                bitmap.Save(memory, ImageFormat.Png);
                memory.Position = 0;

                var bitmapImage = new BitmapImage();
                bitmapImage.BeginInit();
                bitmapImage.StreamSource = memory;
                bitmapImage.CacheOption = BitmapCacheOption.OnLoad;
                bitmapImage.EndInit();
                bitmapImage.Freeze();

                return bitmapImage;
            }
        }

        private static T Clamp<T>(T val, T min, T max) where T : IComparable<T>
        {
            if (val.CompareTo(min) < 0) return min;
            else if (val.CompareTo(max) > 0) return max;
            else return val;
        }

        //===============================================================================
        // Get pixels into a Color array
        Color[,] TranslatePicture(Bitmap pic)
        {
            Color[,] ColorPic = new Color[pic.Width, pic.Height];
            for (int y = 0; (y <= (pic.Height - 1)); y++)
            {
                for (int x = 0; (x <= (pic.Width - 1)); x++)
                {
                    ColorPic[x, y] = pic.GetPixel(x, y);
                }
            }
            return ColorPic;
        }

        // Compare 3 numbers and find the greatest. 1 - first number 2 - second number 3 - third number (R,G,B)
        int CompareNumbers(double first, double second, double third)
        {
            if (first >= second && first >= third)
                return 1;
            else if (second >= first && second >= third)
                return 2;
            else
                return 3;
        }
        // Calculate Min - Max differences in a side of a box
        Dictionary<double, double> GetDiff(Dictionary<double, List<int>> side)
        {
            Dictionary<double, double> diff = new Dictionary<double, double>();

            foreach (var keyValuePair in side)
            {
                double Max = 0;
                double Min = 255;

                foreach (var listItem in keyValuePair.Value)
                {
                    if (Max < listItem)
                        Max = listItem;
                    if (Min > listItem)
                        Min = listItem;
                }
                diff.Add(keyValuePair.Key, Max - Min);
            }
                return diff;
        }

        Dictionary<double, double> GetMedians(Dictionary<double, List<int>> side)
        {
            Dictionary<double, double> medians = new Dictionary<double, double>();
            double median;
            foreach (var keyValuePair in side)
            {
                keyValuePair.Value.Sort();
                if (keyValuePair.Value.Count() % 2 == 0)
                    median = (keyValuePair.Value[(keyValuePair.Value.Count() / 2) - 1] + keyValuePair.Value[keyValuePair.Value.Count() / 2]) / 2;
                else
                    median = keyValuePair.Value[(keyValuePair.Value.Count() - 1) / 2];
                medians.Add(keyValuePair.Key, median);
            }

            return medians;
        }



        //===============================================================================

        // Function filters
        private void Inversionfilter()
        {
            Bitmap pic = BitmapImage2Bitmap(input_Image.Source as BitmapImage);

            for (int y = 0; (y <= (pic.Height - 1)); y++)
            {
                for (int x = 0; (x <= (pic.Width - 1)); x++)
                {
                    Color pixelColor = pic.GetPixel(x, y);
                    pic.SetPixel(x, y, Color.FromArgb(pixelColor.A, 255 - pixelColor.R, 255 - pixelColor.G, 255 - pixelColor.B));
                }
            }
            output_Image.Source = ToBitmapImage(pic);
        }

        private void BrightnesCorrection()
        {
            Bitmap pic = BitmapImage2Bitmap(input_Image.Source as BitmapImage);

            int correction = 30;
            int red, green, blue;

            for (int y = 0; (y <= (pic.Height - 1)); y++)
            {
                for (int x = 0; (x <= (pic.Width - 1)); x++)
                {
                    Color pixelColor = pic.GetPixel(x, y);

                    red = Clamp(pixelColor.R + correction, 0, 255);
                    green = Clamp(pixelColor.G + correction, 0, 255);
                    blue = Clamp(pixelColor.B + correction, 0, 255);

                    pic.SetPixel(x, y, Color.FromArgb(pixelColor.A, red, green, blue));
                }
            }
            output_Image.Source = ToBitmapImage(pic);
        }

        private void ContrastEnhancement()
        {
            Bitmap pic = BitmapImage2Bitmap(input_Image.Source as BitmapImage);
            double aCoeff = 1.65, bCoeff = -82.46;
            int red, green, blue;

            for (int y = 0; (y <= (pic.Height - 1)); y++)
            {
                for (int x = 0; (x <= (pic.Width - 1)); x++)
                {
                    Color pixelColor = pic.GetPixel(x, y);

                    red = Clamp(Convert.ToInt32(aCoeff * pixelColor.R + bCoeff), 0, 255);
                    green = Clamp(Convert.ToInt32(aCoeff * pixelColor.G + bCoeff), 0, 255);
                    blue = Clamp(Convert.ToInt32(aCoeff * pixelColor.B + bCoeff), 0, 255);

                    pic.SetPixel(x, y, Color.FromArgb(pixelColor.A, red, green, blue));
                }
            }
            output_Image.Source = ToBitmapImage(pic);
        }

        private void GammaCorrection()
        {
            Bitmap pic = BitmapImage2Bitmap(input_Image.Source as BitmapImage);
            double gamma = 1.5;

            for (int y = 0; (y <= (pic.Height - 1)); y++)
            {
                for (int x = 0; (x <= (pic.Width - 1)); x++)
                {
                    Color pixelColor = pic.GetPixel(x, y);
                    pic.SetPixel(x, y, Color.FromArgb(pixelColor.A, 
                        GammaCorrectionCalculation(gamma,pixelColor.R), 
                        GammaCorrectionCalculation(gamma, pixelColor.G),
                        GammaCorrectionCalculation(gamma, pixelColor.B)));
                }
            }
            output_Image.Source = ToBitmapImage(pic);
        }

        private int GammaCorrectionCalculation(double gamma, double color)
        {
            return Convert.ToInt32(Math.Pow(color / 255, gamma) * 255);
        }



        // ===============================================================================

        // Convolution filters

        private void Blur3x3()
        {
            Bitmap pic = BitmapImage2Bitmap(input_Image.Source as BitmapImage);
            Bitmap output = new Bitmap(pic.Width, pic.Height);
            int kernelSize = 3;
            int matrixLimits = Convert.ToInt32((kernelSize - 1) / 2);

            for (int y = 0; (y <= (pic.Height - 1)); y++)
            {
                for (int x = 0; (x <= (pic.Width - 1)); x++)
                {
                    int avgR = 0, avgG = 0, avgB = 0;
                    int blurPixelCount = 0;

                    for (int ky = -matrixLimits; ky <= matrixLimits; ky++)
                    {
                        for(int kx = -matrixLimits; kx <= matrixLimits; kx++)
                        {
                            if (y + ky >= 0 && y + ky < pic.Height && x + kx >= 0 && x + kx < pic.Width) 
                            {
                                Color pixel = pic.GetPixel(x + kx, y + ky);

                                avgR += pixel.R;
                                avgG += pixel.G;
                                avgB += pixel.B;

                                blurPixelCount++;
                            }
                        }
                    }

                    avgR = avgR / blurPixelCount;
                    avgG = avgG / blurPixelCount;
                    avgB = avgB / blurPixelCount;

                    output.SetPixel(x, y, Color.FromArgb(pic.GetPixel(x, y).A, avgR, avgG, avgB));
                }
            }

            output_Image.Source = ToBitmapImage(output);
        }

        private void GaussianSmoothing3x3()
        {
            Bitmap pic = BitmapImage2Bitmap(input_Image.Source as BitmapImage);
            Bitmap output = new Bitmap(pic.Width, pic.Height);
            int kernelSize = 3;
            int matrixLimits = Convert.ToInt32((kernelSize - 1) / 2);
            int[,] kernelWeighs = new int[3, 3] { { 0, 1, 0 }, 
                                                  { 1, 4, 1 },
                                                  { 0, 1, 0 } };

            for (int y = 0; (y <= (pic.Height - 1)); y++)
            {
                for (int x = 0; (x <= (pic.Width - 1)); x++)
                {
                    int avgR = 0, avgG = 0, avgB = 0;
                    int blurPixelWeighs = 0;

                    int kernelCoeffX = 0, kernelCoeffY = 0;

                    for (int ky = -matrixLimits; ky <= matrixLimits; ky++)
                    {
                        for (int kx = -matrixLimits; kx <= matrixLimits; kx++)
                        {
                            if (y + ky >= 0 && y + ky < pic.Height && x + kx >= 0 && x + kx < pic.Width)
                            {
                                Color pixel = pic.GetPixel(x + kx, y + ky);

                                avgR += kernelWeighs[kernelCoeffX, kernelCoeffY] * pixel.R;
                                avgG += kernelWeighs[kernelCoeffX, kernelCoeffY] * pixel.G;
                                avgB += kernelWeighs[kernelCoeffX, kernelCoeffY] * pixel.B;

                                blurPixelWeighs += kernelWeighs[kernelCoeffX, kernelCoeffY];
                                kernelCoeffX++;
                            }
                        }
                        kernelCoeffX = 0;
                        kernelCoeffY++;
                    }

                    avgR = avgR / blurPixelWeighs;
                    avgG = avgG / blurPixelWeighs;
                    avgB = avgB / blurPixelWeighs;

                    output.SetPixel(x, y, Color.FromArgb(pic.GetPixel(x, y).A, avgR, avgG, avgB));
                }
            }

            output_Image.Source = ToBitmapImage(output);
        }

        private void Sharpen3x3()
        {
            Bitmap pic = BitmapImage2Bitmap(input_Image.Source as BitmapImage);
            Bitmap output = new Bitmap(pic.Width, pic.Height);

            int b = 5, a = 1;
            int s = b - 4*a;

            int kernelSize = 3;
            int matrixLimits = Convert.ToInt32((kernelSize - 1) / 2);
            double[,] kernelWeighs = new double[3, 3] { {   0, -a/s,   0 },
                                                        { -a/s, b/s, -a/s },
                                                        {   0, -a/s,   0 } };

            for (int y = 0; (y <= (pic.Height - 1)); y++)
            {
                for (int x = 0; (x <= (pic.Width - 1)); x++)
                {
                    int avgR = 0, avgG = 0, avgB = 0;
                    double blurPixelWeighs = 0;

                    int kernelCoeffX = 0, kernelCoeffY = 0;

                    for (int ky = -matrixLimits; ky <= matrixLimits; ky++)
                    {
                        for (int kx = -matrixLimits; kx <= matrixLimits; kx++)
                        {
                            if (y + ky >= 0 && y + ky < pic.Height && x + kx >= 0 && x + kx < pic.Width)
                            {
                                Color pixel = pic.GetPixel(x + kx, y + ky);

                                Clamp(avgR += Convert.ToInt32(kernelWeighs[kernelCoeffX, kernelCoeffY] * pixel.R), 0, 255);
                                Clamp(avgG += Convert.ToInt32(kernelWeighs[kernelCoeffX, kernelCoeffY] * pixel.G), 0, 255);
                                Clamp(avgB += Convert.ToInt32(kernelWeighs[kernelCoeffX, kernelCoeffY] * pixel.B), 0, 255);

                                blurPixelWeighs += kernelWeighs[kernelCoeffX, kernelCoeffY];
                                kernelCoeffX++;
                        }
                            }
                        kernelCoeffX = 0;
                        kernelCoeffY++;
                    }

                    avgR = Clamp(Convert.ToInt32(avgR / blurPixelWeighs),0,255);
                    avgG = Clamp(Convert.ToInt32(avgG / blurPixelWeighs),0,255);
                    avgB = Clamp(Convert.ToInt32(avgB / blurPixelWeighs),0,255);

                    output.SetPixel(x, y, Color.FromArgb(pic.GetPixel(x, y).A, avgR, avgG, avgB));
                }
            }

            output_Image.Source = ToBitmapImage(output);
        }

        private void EdgeDetect3x3()
        {
            Bitmap pic = BitmapImage2Bitmap(input_Image.Source as BitmapImage);
            Bitmap output = new Bitmap(pic.Width, pic.Height);
            int kernelSize = 3;
            int matrixLimits = Convert.ToInt32((kernelSize - 1) / 2);
            int[,] kernelWeighs = new int[3, 3] { { -1,  0,  0 },
                                                  { 0,  1,  0 },
                                                  { 0,  0,  0 } };

            for (int y = 0; (y <= (pic.Height - 1)); y++)
            {
                for (int x = 0; (x <= (pic.Width - 1)); x++)
                {
                    int avgR = 0, avgG = 0, avgB = 0;
                    int blurPixelWeighs = 0;

                    int kernelCoeffX = 0, kernelCoeffY = 0;

                    for (int ky = -matrixLimits; ky <= matrixLimits; ky++)
                    {
                        for (int kx = -matrixLimits; kx <= matrixLimits; kx++)
                        {
                            if (y + ky >= 0 && y + ky < pic.Height && x + kx >= 0 && x + kx < pic.Width)
                            {
                                Color pixel = pic.GetPixel(x + kx, y + ky);

                                avgR += kernelWeighs[kernelCoeffX, kernelCoeffY] * pixel.R;
                                avgG += kernelWeighs[kernelCoeffX, kernelCoeffY] * pixel.G;
                                avgB += kernelWeighs[kernelCoeffX, kernelCoeffY] * pixel.B;

                                blurPixelWeighs += kernelWeighs[kernelCoeffX, kernelCoeffY];
                                kernelCoeffX++;
                            }
                        }
                        kernelCoeffX = 0;
                        kernelCoeffY++;
                    }

                    avgR = Clamp(Convert.ToInt32(avgR / 1), 0, 255);
                    avgG = Clamp(Convert.ToInt32(avgG / 1), 0, 255);
                    avgB = Clamp(Convert.ToInt32(avgB / 1), 0, 255);

                    output.SetPixel(x, y, Color.FromArgb(pic.GetPixel(x, y).A, avgR, avgG, avgB));
                }
            }

            output_Image.Source = ToBitmapImage(output);
        }

        private void Emboss3x3()
        {
            Bitmap pic = BitmapImage2Bitmap(input_Image.Source as BitmapImage);
            Bitmap output = new Bitmap(pic.Width, pic.Height);
            int kernelSize = 3;
            int matrixLimits = Convert.ToInt32((kernelSize - 1) / 2);
            int[,] kernelWeighs = new int[3, 3] { { -1, 0, 1 },
                                                  { -1, 1, 1 },
                                                  { -1, 0, 1 } };

            for (int y = 0; (y <= (pic.Height - 1)); y++)
            {
                for (int x = 0; (x <= (pic.Width - 1)); x++)
                {
                    int avgR = 0, avgG = 0, avgB = 0;
                    int blurPixelWeighs = 0;

                    int kernelCoeffX = 0, kernelCoeffY = 0;

                    for (int ky = -matrixLimits; ky <= matrixLimits; ky++)
                    {
                        for (int kx = -matrixLimits; kx <= matrixLimits; kx++)
                        {
                            if (y + ky >= 0 && y + ky < pic.Height && x + kx >= 0 && x + kx < pic.Width)
                            {
                                Color pixel = pic.GetPixel(x + kx, y + ky);

                                avgR += kernelWeighs[kernelCoeffX, kernelCoeffY] * pixel.R;
                                avgG += kernelWeighs[kernelCoeffX, kernelCoeffY] * pixel.G;
                                avgB += kernelWeighs[kernelCoeffX, kernelCoeffY] * pixel.B;

                                blurPixelWeighs += kernelWeighs[kernelCoeffX, kernelCoeffY];
                                kernelCoeffX++;
                            }
                        }
                        kernelCoeffX = 0;
                        kernelCoeffY++;
                    }

                    avgR = Clamp(Convert.ToInt32(avgR / 1), 0, 255);
                    avgG = Clamp(Convert.ToInt32(avgG / 1), 0, 255);
                    avgB = Clamp(Convert.ToInt32(avgB / 1), 0, 255);

                    output.SetPixel(x, y, Color.FromArgb(pic.GetPixel(x, y).A, avgR, avgG, avgB));
                }
            }

            output_Image.Source = ToBitmapImage(output);
        }

        // ===============================================================================

        // Lab 2 functions

        private void GreyScaleConvertion()
        {
            Bitmap pic = BitmapImage2Bitmap(input_Image.Source as BitmapImage);

            for (int y = 0; (y <= (pic.Height - 1)); y++)
            {
                for (int x = 0; (x <= (pic.Width - 1)); x++)
                {
                    int sum = 0;
                    Color pixelColor = pic.GetPixel(x, y);
                    sum = Convert.ToInt32((pixelColor.R + pixelColor.G + pixelColor.B) / 3);
                    pic.SetPixel(x, y, Color.FromArgb(pixelColor.A, sum, sum, sum));
                }
            }
            output_Image.Source = ToBitmapImage(pic);
        }

        // Fixed average dithering
        //private void Average_Dithering()
        //{
        //    Bitmap pic = BitmapImage2Bitmap(input_Image.Source as BitmapImage);
        //    int average = 0;

        //    int number = Convert.ToInt32(slValue_Dither.Value);

        //    int[] intervals = new int[number];
        //    int step = (255 / number);

        //    for (int i = 0; i < number - 1; i++)
        //        intervals[i] = step * (i + 1);

        //    intervals[number - 1]=255;

        //    for (int y = 0; (y <= (pic.Height - 1)); y++)
        //    {
        //        for (int x = 0; (x <= (pic.Width - 1)); x++)
        //        {
        //            Color pixelColor = pic.GetPixel(x, y);
        //            average += Convert.ToInt32((pixelColor.R + pixelColor.G + pixelColor.B) / 3);
        //        }
        //    }

        //    average = Convert.ToInt32(average / (pic.Height * pic.Width));

        //    for (int y = 0; (y <= (pic.Height - 1)); y++)
        //    {
        //        for (int x = 0; (x <= (pic.Width - 1)); x++)
        //        {
        //            Color pixelColor = pic.GetPixel(x, y);

        //            for (int i = 0; i < number; i++)
        //            {
        //                if (Convert.ToInt32((pixelColor.R + pixelColor.G + pixelColor.B) / 3) < intervals[i])
        //                {
        //                    pic.SetPixel(x, y, Color.FromArgb(pixelColor.A, intervals[i], intervals[i], intervals[i]));
        //                    break;
        //                }
        //            }
        //        }
        //    }

        //    output_Image.Source = ToBitmapImage(pic);
        //}

        private void RandomDithering()
        {
            GreyScaleConvertion();
            Bitmap pic = BitmapImage2Bitmap(output_Image.Source as BitmapImage);
            Color[,] ColorPic = TranslatePicture(pic);

            int k = Convert.ToInt32(slValue_Random.Value) - 1;
            Random rnd = new Random();

            List<double> levels = new List<double>();
            double step = 256 / k;
            for (int i = 1; i <= k; i++)
                levels.Add(step * i);

            for (int y = 0; (y <= (pic.Height - 1)); y++)
            {
                for (int x = 0; (x <= (pic.Width - 1)); x++)
                {
                    for(int i = 0; i < levels.Count; i++)
                    {
                        if(ColorPic[x, y].R < levels[i])
                        {
                            int threshold = rnd.Next(Convert.ToInt32(levels[i] - step), Convert.ToInt32(levels[i]));
                            pic.SetPixel(x, y, Color.FromArgb(ColorPic[x, y].A,
                                Clamp(Convert.ToInt32((ColorPic[x, y].R > threshold ? levels[i] : levels[i] - step)),0,255),
                                Clamp(Convert.ToInt32((ColorPic[x, y].G > threshold ? levels[i] : levels[i] - step)), 0, 255),
                                Clamp(Convert.ToInt32((ColorPic[x, y].B > threshold ? levels[i] : levels[i] - step)), 0, 255)));
                            break;
                        }
                    }
                }
            }
            output_Image.Source = ToBitmapImage(pic);
        }
        
        private void MedianCut()
        {
            Bitmap pic = BitmapImage2Bitmap(input_Image.Source as BitmapImage);
            Color[,] ColorPic = TranslatePicture(pic);

            int numberOfCuts = Convert.ToInt32(slValue.Value);

            //Dictionary (Upper value bound of the box, List of values in the box)
            Dictionary<double, List<int>> RSide = new Dictionary<double, List<int>>();
            Dictionary<double, List<int>> GSide = new Dictionary<double, List<int>>();
            Dictionary<double, List<int>> BSide = new Dictionary<double, List<int>>();

            List<int> Rtmp = new List<int>();
            List<int> Gtmp = new List<int>();
            List<int> Btmp = new List<int>();
            // Create one big box with all Colors
            for (int y = 0; (y <= (pic.Height - 1)); y++)
            {
                for (int x = 0; (x <= (pic.Width - 1)); x++)
                {
                    Rtmp.Add(ColorPic[x, y].R);
                    Gtmp.Add(ColorPic[x, y].G);
                    Btmp.Add(ColorPic[x, y].B);
                }
            }
            RSide.Add(256, Rtmp);
            GSide.Add(256, Gtmp);
            BSide.Add(256, Btmp);

            for (int i = 0; i < numberOfCuts; i++)
            { 
                // Determine the range of colors in each box
                // Dictionary (Upper value bound of the box, Min - Max difference in the box)
                Dictionary<double, double> Rdiff = GetDiff(RSide);
                Dictionary<double, double> Gdiff = GetDiff(GSide);
                Dictionary<double, double> Bdiff = GetDiff(BSide);

                // Longest side of the box, and split the box along that axis into two smaller boxes
                double median;
                double MaxValueKey;
                // Lists of values that will split the box with the biggest difference
                List<int> box1 = new List<int>();
                List<int> box2 = new List<int>();
                switch (CompareNumbers(Rdiff.Values.Max(), Gdiff.Values.Max(), Bdiff.Values.Max()))
                {
                    case 1: // Case 1: R channel has the biggest value spread
                        MaxValueKey = Rdiff.Aggregate((x, y) => x.Value > y.Value ? x : y).Key;
                        RSide[MaxValueKey].Sort();
                        if (RSide[MaxValueKey].Count() % 2 == 0)
                            median = (RSide[MaxValueKey][(RSide[MaxValueKey].Count() / 2) - 1] + RSide[MaxValueKey][RSide[MaxValueKey].Count() / 2]) / 2;
                        else
                            median = RSide[MaxValueKey][(RSide[MaxValueKey].Count() - 1) / 2];

                        foreach (var value in RSide[MaxValueKey])
                        {
                            if (value < median)
                                box1.Add(value);
                            else
                                box2.Add(value);
                        }
                        RSide.Remove(MaxValueKey);
                        RSide.Add(median, box1);
                        RSide.Add(MaxValueKey, box2);

                        break;
                    case 2:// Case 2: G channel has the biggest value spread
                        MaxValueKey = Gdiff.Aggregate((x, y) => x.Value > y.Value ? x : y).Key;
                        GSide[MaxValueKey].Sort();
                        if (GSide[MaxValueKey].Count() % 2 == 0)
                            median = (GSide[MaxValueKey][(GSide[MaxValueKey].Count() / 2) - 1] + GSide[MaxValueKey][GSide[MaxValueKey].Count() / 2]) / 2;
                        else
                            median = GSide[MaxValueKey][(GSide[MaxValueKey].Count() - 1) / 2];

                        foreach (var value in GSide[MaxValueKey])
                        {
                            if (value < median)
                                box1.Add(value);
                            else
                                box2.Add(value);
                        }
                        GSide.Remove(MaxValueKey);
                        GSide.Add(median, box1);
                        GSide.Add(MaxValueKey, box2);
                        break;
                    case 3: // Case 3: B channel has the biggest value spread
                        MaxValueKey = Bdiff.Aggregate((x, y) => x.Value > y.Value ? x : y).Key;
                        BSide[MaxValueKey].Sort();
                        if (BSide[MaxValueKey].Count() % 2 == 0)
                            median = (BSide[MaxValueKey][(BSide[MaxValueKey].Count() / 2) - 1] + BSide[MaxValueKey][BSide[MaxValueKey].Count() / 2]) / 2;
                        else
                            median = BSide[MaxValueKey][(BSide[MaxValueKey].Count() - 1) / 2];

                        foreach (var value in BSide[MaxValueKey])
                        {
                            if (value < median)
                                box1.Add(value);
                            else
                                box2.Add(value);
                        }
                        BSide.Remove(MaxValueKey);
                        BSide.Add(median, box1);
                        BSide.Add(MaxValueKey, box2);
                        break;
                }

            }
            // After the box of all colors has been split find median of each one
            // Dictionary (Upper value bound of the box, Median of the box)
            Dictionary<double, double> Rmedians = GetMedians(RSide);
            Dictionary<double, double> Gmedians = GetMedians(GSide);
            Dictionary<double, double> Bmedians = GetMedians(BSide);

            // Iterate over the full picture applying color medians with respect to the box any color is in
            for (int y = 0; (y <= (pic.Height - 1)); y++)
            {
                for (int x = 0; (x <= (pic.Width - 1)); x++)
                {
                    int R = 0, G = 0, B = 0;
                    foreach(var keyValuePair in Rmedians)
                    {
                        if(ColorPic[x,y].R<=keyValuePair.Key)
                        {
                            R = Convert.ToInt32(keyValuePair.Value);
                            break;
                        }
                    }
                    foreach (var keyValuePair in Gmedians)
                    {
                        if (ColorPic[x, y].G <= keyValuePair.Key)
                        {
                            G = Convert.ToInt32(keyValuePair.Value);
                            break;
                        }
                    }
                    foreach (var keyValuePair in Bmedians)
                    {
                        if (ColorPic[x, y].B <= keyValuePair.Key)
                        {
                            B = Convert.ToInt32(keyValuePair.Value);
                            break;
                        }
                    }
                    pic.SetPixel(x, y, Color.FromArgb(ColorPic[x, y].A, R, G, B));
                }
            }
            
            output_Image.Source = ToBitmapImage(pic);

        }

        private void Lab1()
        {
            Bitmap pic = BitmapImage2Bitmap(input_Image.Source as BitmapImage);
            Bitmap output = new Bitmap(pic.Width, pic.Height);
            int kernelSize = 5;
            int matrixLimits = Convert.ToInt32((kernelSize - 1) / 2);
            int M = 200;

            for (int y = 0; (y <= (pic.Height - 1)); y++)
            {
                for (int x = 0; (x <= (pic.Width - 1)); x++)
                {
                    int avgR = 0, avgG = 0, avgB = 0;
                    int blurPixelCount = 0;
                    Color centrePixel = pic.GetPixel(x, y);

                    for (int ky = -matrixLimits; ky <= matrixLimits; ky++)
                    {
                        for (int kx = -matrixLimits; kx <= matrixLimits; kx++)
                        {
                            if (y + ky >= 0 && y + ky < pic.Height && x + kx >= 0 && x + kx < pic.Width)
                            {
                                Color pixel = pic.GetPixel(x + kx, y + ky);

                                if(Math.Abs(centrePixel.R-pixel.R) + Math.Abs(centrePixel.G-pixel.G)+Math.Abs(centrePixel.B-pixel.B) < M)
                                {
                                    avgR += pixel.R;
                                    avgG += pixel.G;
                                    avgB += pixel.B;

                                    blurPixelCount++;
                                }
                            }
                        }
                    }

                    avgR = avgR / blurPixelCount;
                    avgG = avgG / blurPixelCount;
                    avgB = avgB / blurPixelCount;

                    output.SetPixel(x, y, Color.FromArgb(pic.GetPixel(x, y).A, avgR, avgG, avgB));
                }
            }

            output_Image.Source = ToBitmapImage(output);
        }

        private void Hist_Stretch()
        {
            Bitmap pic = BitmapImage2Bitmap(input_Image.Source as BitmapImage);
            Color[,] ColorPic = TranslatePicture(pic);

            int[] Rchannel = new int[256];
            int[] Gchannel = new int[256];
            int[] Bchannel = new int[256];

            for (int i = 0; i < 256; i++)
            {
                Rchannel[i] = 0;
                Gchannel[i] = 0;
                Bchannel[i] = 0;
            }

            for (int y = 0; (y <= (pic.Height - 1)); y++)
            {
                for (int x = 0; (x <= (pic.Width - 1)); x++)
                {
                    Rchannel[ColorPic[x, y].R] += 1;
                    Gchannel[ColorPic[x, y].G] += 1;
                    Bchannel[ColorPic[x, y].B] += 1;
                }
            }

            int Rmax = 0, Rmin = 0, Gmax = 0, Gmin = 0, Bmax = 0, Bmin = 0;

            for (int i = 0; i < 256; i++)
            {
                if (Rchannel[i] != 0 && Rmin == 0)
                    Rmin = i;
                if (Gchannel[i] != 0 && Gmin == 0)
                    Gmin = i;
                if (Bchannel[i] != 0 && Bmin == 0)
                    Bmin = i;
                if (Rmin != 0 && Gmin != 0 && Bmin != 0)
                    break;
            }
            for (int i = 255; i >= 0; i--)
            {
                if (Rchannel[i] != 0 && Rmax == 0)
                    Rmax = i;
                if (Gchannel[i] != 0 && Gmax == 0)
                    Gmax = i;
                if (Bchannel[i] != 0 && Bmax == 0)
                    Bmax = i;
                if (Rmax != 0 && Gmax != 0 && Bmax != 0)
                    break;
            }

            for (int y = 0; (y <= (pic.Height - 1)); y++)
            {
                for (int x = 0; (x <= (pic.Width - 1)); x++)
                {
                    
                    pic.SetPixel(x, y, Color.FromArgb(ColorPic[x,y].A,
                        Clamp(Convert.ToInt32((255 / (Rmax - Rmin)) * (ColorPic[x, y].R - Rmin)),0,255),
                        Clamp(Convert.ToInt32((255 / (Gmax - Gmin)) * (ColorPic[x, y].G - Gmin)),0,255),
                        Clamp(Convert.ToInt32((255 / (Bmax - Bmin)) * (ColorPic[x, y].B - Rmin)),0,255)));
                }
            }
            output_Image.Source = ToBitmapImage(pic);
        }

        private void Hist_Eq()
        {
            Bitmap pic = BitmapImage2Bitmap(input_Image.Source as BitmapImage);
            Color[,] ColorPic = TranslatePicture(pic);

            int[] Rchannel = new int[256];
            int[] Gchannel = new int[256];
            int[] Bchannel = new int[256];

            for (int i = 0; i < 256; i++)
            {
                Rchannel[i] = 0;
                Gchannel[i] = 0;
                Bchannel[i] = 0;
            }

            for (int y = 0; (y <= (pic.Height - 1)); y++)
            {
                for (int x = 0; (x <= (pic.Width - 1)); x++)
                {
                    Rchannel[ColorPic[x, y].R] += 1;
                    Gchannel[ColorPic[x, y].G] += 1;
                    Bchannel[ColorPic[x, y].B] += 1;
                }
            }

            double[] DR = new double[256];
            double[] DG = new double[256];
            double[] DB = new double[256];

            int Npixels = pic.Width * pic.Height;
            for (int i = 0; i < 256; i++)
            {
                double hR = 0, hG = 0, hB = 0;
                for (int j = 0; j <= i; j++)
                {
                    hR += Rchannel[j];
                    hG += Gchannel[j];
                    hB += Bchannel[j];
                }
                DR[i] = hR / Npixels;
                DG[i] = hG / Npixels;
                DB[i] = hB / Npixels;
            }

            double  DRmin = 0, DGmin = 0, DBmin = 0;

            for (int i = 0; i < 256; i++)
            {
                if (DR[i] != 0 && DRmin == 0)
                    DRmin = DR[i];
                if (DG[i] != 0 && DGmin == 0)
                    DGmin = DG[i];
                if (DB[i] != 0 && DBmin == 0)
                    DBmin = DB[i];
                if (DRmin != 0 && DGmin != 0 && DBmin != 0)
                    break;
            }

            for (int y = 0; (y <= (pic.Height - 1)); y++)
            {
                for (int x = 0; (x <= (pic.Width - 1)); x++)
                {

                    pic.SetPixel(x, y, Color.FromArgb(ColorPic[x, y].A,
                        Clamp(Convert.ToInt32(((DR[ColorPic[x, y].R] - DRmin) / 1 - DRmin) * (256 - 1)), 0, 255),
                        Clamp(Convert.ToInt32(((DG[ColorPic[x, y].R] - DGmin) / 1 - DGmin) * (256 - 1)), 0, 255),
                        Clamp(Convert.ToInt32(((DB[ColorPic[x, y].R] - DBmin) / 1 - DBmin) * (256 - 1)), 0, 255)));
                }
            }
            output_Image.Source = ToBitmapImage(pic);
        }


        // ===============================================================================

        // Button events 

        private void Apply_Function(object sender, RoutedEventArgs e)
        {
            previous_Image = input_Image.Source as BitmapImage;
            input_Image.Source = output_Image.Source;
        }

        private void Undo_Click(object sender, RoutedEventArgs e)
        {
            input_Image.Source = previous_Image;
        }

        private void UndoAll_Click(object sender, RoutedEventArgs e)
        {
            input_Image.Source = original_Image;
        }

        private void Inversion_Click(object sender, RoutedEventArgs e)
        {
            Inversionfilter();
        }

        private void Brightnes_Click(object sender, RoutedEventArgs e)
        {
            BrightnesCorrection();
        }

        private void Contrast_Click(object sender, RoutedEventArgs e)
        {
            ContrastEnhancement();
        }

        private void Gamma_Click(object sender, RoutedEventArgs e)
        {
            GammaCorrection();
        }

        private void Blur_Click(object sender, RoutedEventArgs e)
        {
            Blur3x3();
        }

        private void GaussianBlur_Click(object sender, RoutedEventArgs e)
        {
            GaussianSmoothing3x3();
        }

        private void Sharpen_Click(object sender, RoutedEventArgs e)
        {
            Sharpen3x3();
        }

        private void Edge_Click(object sender, RoutedEventArgs e)
        {
            EdgeDetect3x3();
        }

        private void Emboss_Click(object sender, RoutedEventArgs e)
        {
            Emboss3x3();
        }

        private void Grey_Click(object sender, RoutedEventArgs e)
        {
            GreyScaleConvertion();
        }

        private void Dithering_Click(object sender, RoutedEventArgs e)
        {
            //Average_Dithering();
        }

        private void MedCut_Click(object sender, RoutedEventArgs e)
        {
            MedianCut();
        }

        private void Lab1_Click(object sender, RoutedEventArgs e)
        {
            Lab1();
        }

        private void Hist_Stretch_Click(object sender, RoutedEventArgs e)
        {
            Hist_Stretch();
        }

        private void Hist_Eq_Click(object sender, RoutedEventArgs e)
        {
            Hist_Eq();
        }

        private void Random_dithering_button_Click(object sender, RoutedEventArgs e)
        {
            RandomDithering();
        }



        // =============================================================

    }
}
